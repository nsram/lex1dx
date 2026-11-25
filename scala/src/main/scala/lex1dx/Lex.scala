package lex1dx

import java.util.Arrays

enum Mode:
  case Tuned, Sig5

final case class Options(
    mode: Mode = Mode.Tuned,
    tolerance: Double = 1e-11,
    coarseGroups: Int = 0 // 0 means defaultM(k)
)

object Options:
  def defaultM(k: Int): Int =
    if k < 1 then 4096
    else math.max(4096, 32 * k)

final case class Result(
    bounds: IArray[Int],
    counts: IArray[Int],
    means: IArray[Double],
    upper: IArray[Double],
    sse: Double,
    D: Int,
    Dprime: Int,
    sweeps: Int,
    moves: Int,
    presorted: Boolean
):
  def k: Int = counts.length
  def assign(x: Double): Int =
    var (lo, hi) = (0, k - 1)
    while lo < hi do
      val mid = (lo + hi) >>> 1
      if upper(mid) >= x then hi = mid else lo = mid + 1
    lo

object Lex:
  def cluster(xs: Array[Double], k: Int, options: Options = Options()): Result =
    if xs.isEmpty then throw new IllegalArgumentException("input is empty")
    if k < 1 then throw new IllegalArgumentException(s"k must be at least 1, got $k")
    if !(options.tolerance > 0 && options.tolerance < 1) then
      throw new IllegalArgumentException(s"tolerance must be in (0, 1), got ${options.tolerance}")
    var i = 0
    while i < xs.length do
      val v = xs(i)
      if v.isNaN then throw new IllegalArgumentException(s"input value at index $i is NaN; only finite values are allowed")
      if v.isInfinite then
        val sp = if v > 0 then "Infinity" else "-Infinity"
        throw new IllegalArgumentException(s"input value at index $i is $sp; only finite values are allowed")
      i += 1

    val sorted = xs.map(_ + 0.0)
    val presorted = isNonDecreasing(sorted)
    if !presorted then Arrays.sort(sorted)
    val D = countDistinct(sorted)
    val sig5 = options.mode == Mode.Sig5

    val (atoms, elemStart) =
      if sig5 then
        val g = SigGroup.group(sorted, 5)
        (Atoms(g.mean, g.w, 0.0), g.start)
      else
        (Atoms(sorted, 0.0), null: Array[Int])

    val Dprime = atoms.d
    var kk = k
    if kk > atoms.d then kk = atoms.d
    val Mreq = if options.coarseGroups == 0 then Options.defaultM(k) else options.coarseGroups
    val m = math.max(kk, math.min(Mreq, atoms.d))

    val groups =
      if m >= atoms.d then atoms.start.clone()
      else Collapse(atoms, m)

    val b0: Array[Int] =
      if kk <= 1 then Array(0, atoms.n)
      else if kk >= groups.length - 1 then groups
      else
        val agg = aggAtoms(atoms, groups)
        val gb = Wilber.cluster(agg._1, agg._2, agg._3, kk)
        // gb is in group-index space [0, G]; map to atom/element indices.
        Array.tabulate(gb.length)(i => groups(gb(i)))

    val boundsWork = b0.clone()
    val refine = Refine(atoms, boundsWork, options.tolerance)
    refine.run()
    val raw = mapRaw(boundsWork, elemStart, sig5)
    summarize(sorted, raw, D, Dprime, refine.sweeps, refine.moves, presorted)

  /** Paper-faithful exact DP, used only as a quality reference in tests. */
  def exact(xs: Array[Double], k: Int): Result =
    val sorted = xs.map(_ + 0.0)
    Arrays.sort(sorted)
    val atoms = Atoms(sorted, 0.0)
    var kk = k
    if kk > atoms.d then kk = atoms.d
    val b = Exact(atoms, kk)
    val refine = Refine(atoms, b, 1e-11)
    refine.run()
    refine.sweeps = 0 // keep lex1d exact_largest_cuts convention? tests only need SSE
    summarize(sorted, b, atoms.d, atoms.d, 0, refine.moves, false)

  private def isNonDecreasing(a: Array[Double]): Boolean =
    var i = 1
    while i < a.length do
      if a(i) < a(i - 1) then return false
      i += 1
    true

  private def countDistinct(sorted: Array[Double]): Int =
    if sorted.isEmpty then 0
    else
      var d = 1
      var i = 1
      while i < sorted.length do
        if sorted(i) != sorted(i - 1) then d += 1
        i += 1
      d

  private def aggAtoms(a: Atoms, groups: Array[Int]): (Array[Double], Array[Double], Array[Double]) =
    val G = groups.length - 1
    val w = new Array[Double](G)
    val wx = new Array[Double](G)
    val wx2 = new Array[Double](G)
    var i = 0
    while i < G do
      val lo = groups(i)
      val hi = groups(i + 1)
      // point mass at the group mean (wx2 = W mu^2); see agg_atoms in lex1dx.cpp
      if a.weighted then
        val W = a.wsum(lo, hi)
        val mu = a.mean(lo, hi)
        w(i) = W; wx(i) = W * mu; wx2(i) = W * mu * mu
      else
        val c = (hi - lo).toDouble
        val mu = a.mean(lo, hi)
        w(i) = c; wx(i) = c * mu; wx2(i) = c * mu * mu
      i += 1
    (w, wx, wx2)

  private def mapRaw(b: Array[Int], elemStart: Array[Int], sig5: Boolean): Array[Int] =
    if !sig5 then b
    else Array.tabulate(b.length)(c => elemStart(b(c)))

  private def summarize(
      sorted: Array[Double],
      b: Array[Int],
      D: Int,
      Dprime: Int,
      sweeps: Int,
      moves: Int,
      presorted: Boolean
  ): Result =
    val k = b.length - 1
    val n = sorted.length
    var center = 0.0
    var i = 0
    while i < n do { center += sorted(i); i += 1 }
    center /= n
    val counts = new Array[Int](k)
    val means = new Array[Double](k)
    val upper = new Array[Double](k)
    var sse = 0.0
    var bin = 0
    while bin < k do
      val lo = b(bin)
      val hi = b(bin + 1)
      var s = 0.0
      var e = lo
      while e < hi do { s += sorted(e) - center; e += 1 }
      val m = s / (hi - lo).toDouble
      e = lo
      while e < hi do
        val dv = sorted(e) - center - m
        sse += dv * dv
        e += 1
      counts(bin) = hi - lo
      means(bin) = center + m
      upper(bin) = sorted(hi - 1)
      bin += 1
    Result(IArray.unsafeFromArray(b.clone()), IArray.unsafeFromArray(counts),
      IArray.unsafeFromArray(means), IArray.unsafeFromArray(upper),
      sse, D, Dprime, sweeps, moves, presorted)
end Lex
