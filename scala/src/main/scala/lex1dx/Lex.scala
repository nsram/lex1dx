package lex1dx

import java.util.Arrays

enum Mode:
  case Tuned, Sig5, Exact

enum Solver:
  case Wilber, Ckmeans

final case class Options(
    mode: Mode = Mode.Tuned,
    tolerance: Double = 1e-11,
    coarseGroups: Int = 0, // 0 means defaultM(k)
    sig: Int = -1, // -1 follows mode; 0 = full precision; 1..9 = significant digits
    solver: Solver = Solver.Wilber
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
    sig: Int, // significant digits used for grouping; 0 = full precision
    sweeps: Int,
    moves: Int,
    presorted: Boolean,
    solverUsed: String = "wilber" // "wilber" or "ckmeans"; mirrors the request if the DP is skipped
):
  def k: Int = counts.length
  def assign(x: Double): Int =
    var (lo, hi) = (0, k - 1)
    while lo < hi do
      val mid = (lo + hi) >>> 1
      if upper(mid) >= x then hi = mid else lo = mid + 1
    lo

object Lex:
  def defaultM(k: Int): Int = Options.defaultM(k)

  def cluster(xs: Array[Double], k: Int, options: Options = Options()): Result =
    clusterImpl(xs, k, null, options)

  def cluster(xs: Array[Double], k: Int, weights: Array[Double]): Result =
    cluster(xs, k, weights, Options())

  def cluster(xs: Array[Double], k: Int, weights: Array[Double], options: Options): Result =
    if weights == null then throw new IllegalArgumentException("weights must not be null")
    clusterImpl(xs, k, weights, options)

  private def clusterImpl(xs: Array[Double], k: Int, weights: Array[Double], options: Options): Result =
    if xs.isEmpty then throw new IllegalArgumentException("input is empty")
    if k < 1 then throw new IllegalArgumentException(s"k must be at least 1, got $k")
    if !(options.tolerance > 0 && options.tolerance < 1) then
      throw new IllegalArgumentException(s"tolerance must be in (0, 1), got ${options.tolerance}")
    if options.coarseGroups != 0 && options.coarseGroups < 1 then
      throw new IllegalArgumentException(s"M must be at least 1, got ${options.coarseGroups}")
    var i = 0
    while i < xs.length do
      val v = xs(i)
      if v.isNaN then throw new IllegalArgumentException(s"input value at index $i is NaN; only finite values are allowed")
      if v.isInfinite then
        val sp = if v > 0 then "Infinity" else "-Infinity"
        throw new IllegalArgumentException(s"input value at index $i is $sp; only finite values are allowed")
      i += 1

    val weighted = weights != null
    if weighted then checkWeights(weights, xs.length)

    val sorted = xs.map(_ + 0.0)
    val wts = if weighted then weights.clone() else null
    val didSort =
      if weighted then sortPairs(sorted, wts)
      else
        val already = isNonDecreasing(sorted)
        if !already then Arrays.sort(sorted)
        !already
    val D = countDistinct(sorted)
    var usedSig = resolveSig(options)

    var mapGroups = false
    var elemStart: Array[Int] = null
    val atoms =
      if usedSig > 0 then
        // Raise precision until there are at least k groups, or we reach raw atoms.
        // This is the v0.1.0 k-preserving fallback, starting from the requested sig.
        var sig = usedSig
        var g = SigGroup.group(sorted, sig)
        if weighted then applyUserWeights(g, sorted, wts)
        while g.w.length < k && g.w.length < D && sig < 15 do
          sig += 1
          g = SigGroup.group(sorted, sig)
          if weighted then applyUserWeights(g, sorted, wts)
        if g.w.length < k && g.w.length < D then
          usedSig = 0
          if weighted then Atoms(sorted, wts, 0.0) else Atoms(sorted, 0.0)
        else
          usedSig = sig
          mapGroups = true
          elemStart = g.start
          sortPairs(g.mean, g.w)
          Atoms(g.mean, g.w, 0.0)
      else if weighted then Atoms(sorted, wts, 0.0)
      else Atoms(sorted, 0.0)

    val Dprime = atoms.d
    var kk = k
    if kk > atoms.d then kk = atoms.d

    val groups =
      if options.mode == Mode.Exact then atoms.start.clone()
      else
        val Mreq = if options.coarseGroups == 0 then Options.defaultM(k) else options.coarseGroups
        val m = math.max(kk, math.min(Mreq, atoms.d))
        if m >= atoms.d then atoms.start.clone()
        else Collapse(atoms, m)

    val requested = if options.solver == Solver.Ckmeans then "ckmeans" else "wilber"
    var solverUsed = requested
    val b0: Array[Int] =
      if kk <= 1 then
        solverUsed = requested
        Array(0, atoms.n)
      else if kk >= groups.length - 1 then
        solverUsed = requested
        groups
      else if options.solver == Solver.Ckmeans then
        solverUsed = "ckmeans"
        ckmeansOnGroups(atoms, groups, kk)
      else
        wilberOnGroupsOrFallback(atoms, groups, kk) match
          case (b, true)  => solverUsed = "ckmeans"; b
          case (b, false) => solverUsed = "wilber"; b

    // Refine on the atoms actually used (grouped at the chosen sig, or raw
    // distinct-value atoms after fallback). Then map grouped bounds to raw
    // element indices. Tuned already Refines on raw atoms.
    val boundsWork = b0.clone()
    val refine = Refine(atoms, boundsWork, options.tolerance)
    refine.run()
    val raw = if mapGroups then mapRaw(boundsWork, elemStart, true) else boundsWork
    if weighted then
      summarizeWeighted(sorted, wts, raw, D, Dprime, usedSig, refine.sweeps, refine.moves, !didSort, solverUsed)
    else
      summarize(sorted, raw, D, Dprime, usedSig, refine.sweeps, refine.moves, !didSort, solverUsed)

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
    summarize(sorted, b, atoms.d, atoms.d, 0, 0, refine.moves, false, "ckmeans")

  private def resolveSig(options: Options): Int =
    if options.sig == -1 then
      if options.mode == Mode.Sig5 then 5 else 0
    else if options.sig < 0 || options.sig > 9 then
      throw new IllegalArgumentException(
        s"sig must be an integer in 0..9, got ${options.sig}. Use sig=0 for full precision"
      )
    else options.sig

  private def checkWeights(weights: Array[Double], n: Int): Unit =
    if weights.length != n then
      throw new IllegalArgumentException(
        s"weights size must match the input length, got ${weights.length} and $n"
      )
    var W = 0.0
    var i = 0
    while i < n do
      val w = weights(i)
      if w.isNaN || w.isInfinite || w < 0 then
        val shown = if w.isNaN then "NaN" else if w.isInfinite then (if w > 0 then "Infinity" else "-Infinity") else w.toString
        throw new IllegalArgumentException(
          s"weight at index $i is $shown; weights must be finite and >= 0"
        )
      W += w
      i += 1
    if !(W > 0) || W.isNaN || W.isInfinite then
      throw new IllegalArgumentException(s"weights must sum to a positive finite total, got $W")

  /** Returns true if a sort was performed. Ties break on weight. */
  private def sortPairs(x: Array[Double], w: Array[Double]): Boolean =
    val n = x.length
    var i = 0
    while i < n do { x(i) = x(i) + 0.0; i += 1 }
    def less(xa: Double, wa: Double, xb: Double, wb: Double): Boolean =
      xa < xb || (xa == xb && wa < wb)
    var sorted = true
    i = 1
    while i < n && sorted do
      if less(x(i), w(i), x(i - 1), w(i - 1)) then sorted = false
      i += 1
    if sorted then false
    else
      val z = Array.tabulate(n)(j => (x(j), w(j)))
      Arrays.sort(z, (a: (Double, Double), b: (Double, Double)) =>
        if less(a._1, a._2, b._1, b._2) then -1
        else if less(b._1, b._2, a._1, a._2) then 1
        else 0
      )
      i = 0
      while i < n do { x(i) = z(i)._1; w(i) = z(i)._2; i += 1 }
      true

  private def applyUserWeights(g: SigGroup.Groups, sorted: Array[Double], wts: Array[Double]): Unit =
    val G = g.start.length - 1
    var i = 0
    while i < G do
      val lo = g.start(i)
      val hi = g.start(i + 1)
      var W = 0.0
      var S = 0.0
      var e = lo
      while e < hi do
        W += wts(e)
        S += wts(e) * sorted(e)
        e += 1
      g.w(i) = W
      if W > 0 then g.mean(i) = S / W
      i += 1


  private def validGroupBounds(gb: Array[Int], kk: Int, G: Int): Boolean =
    if gb.length != kk + 1 then false
    else
      var ok = gb(0) >= 0 && gb(kk) <= G
      var i = 0
      while ok && i < kk do
        if !(gb(i) < gb(i + 1)) then ok = false
        if gb(i) < 0 || gb(i) > G then ok = false
        i += 1
      ok && gb(kk) >= 0 && gb(kk) <= G

  // C++ layout: (k-1)*(G+1)*4. Scala unweighted Exact uses (k+1)*(G+1)*4.
  private def ckmeansTableBytes(k: Int, G: Int, unweighted: Boolean): Long =
    if unweighted then
      if k < 0 || G < 0 then 0L
      else
        val rows = k.toLong + 1L
        val cols = G.toLong + 1L
        if cols != 0 && rows > Long.MaxValue / cols then Long.MaxValue
        else
          val cells = rows * cols
          if cells > Long.MaxValue / 4L then Long.MaxValue else cells * 4L
    else if k <= 1 || G < 0 then 0L
    else
      val rows = k.toLong - 1L
      val cols = G.toLong + 1L
      if cols != 0 && rows > Long.MaxValue / cols then Long.MaxValue
      else
        val cells = rows * cols
        if cells > Long.MaxValue / 4L then Long.MaxValue else cells * 4L

  private val ckmeansFallbackCapBytes: Long = 2L * 1024L * 1024L * 1024L

  private def checkCkmeansFallbackBudgetImpl(k: Int, G: Int, unweighted: Boolean): Unit =
    val bytes = ckmeansTableBytes(k, G, unweighted)
    if bytes > ckmeansFallbackCapBytes then
      val gib = bytes.toDouble / (1024.0 * 1024.0 * 1024.0)
      // RuntimeException, not IllegalStateException, so the Wilber catch does not swallow it.
      throw new RuntimeException(
        s"Wilber interpolation search failed; the Ckmeans fallback would allocate a (k-1)x(G+1) back-pointer table of $gib GiB, which exceeds the 2 GiB cap. Use mode=\"tuned\", or solver=\"ckmeans\" with a smaller k. (k=$k, G=$G)"
      )

  private def wilberOnGroupsOrFallback(a: Atoms, groups: Array[Int], kk: Int): (Array[Int], Boolean) =
    val G = groups.length - 1
    def fallback(): (Array[Int], Boolean) =
      checkCkmeansFallbackBudgetImpl(kk, G, !a.weighted)
      (ckmeansOnGroups(a, groups, kk), true)
    try
      val agg = aggAtoms(a, groups)
      val gb = Wilber.cluster(agg._1, agg._2, agg._3, kk)
      // Validate in group-index space before groups(gb(i)).
      if !validGroupBounds(gb, kk, G) then fallback()
      else (Array.tabulate(gb.length)(i => groups(gb(i))), false)
    catch
      case _: IllegalStateException | _: ArrayIndexOutOfBoundsException =>
        fallback()

  /** Public so tests can assert the 2 GiB cap without allocating the table. */
  def checkCkmeansFallbackBudget(k: Int, G: Int, unweighted: Boolean = true): Unit =
    checkCkmeansFallbackBudgetImpl(k, G, unweighted)

  private def ckmeansOnGroups(a: Atoms, groups: Array[Int], kk: Int): Array[Int] =
    if a.weighted then Exact.weighted(a, groups, kk)
    else if groups.length - 1 == a.d && Arrays.equals(groups, a.start) then Exact(a, kk)
    else
      val G = groups.length - 1
      val means = Array.tabulate(G)(g => a.mean(groups(g), groups(g + 1)))
      Exact(groups, means, kk)

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

  private def mapRaw(b: Array[Int], elemStart: Array[Int], sigGrouped: Boolean): Array[Int] =
    if !sigGrouped then b
    else Array.tabulate(b.length)(c => elemStart(b(c)))

  private def summarize(
      sorted: Array[Double],
      b: Array[Int],
      D: Int,
      Dprime: Int,
      sig: Int,
      sweeps: Int,
      moves: Int,
      presorted: Boolean,
      solverUsed: String
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
      sse, D, Dprime, sig, sweeps, moves, presorted, solverUsed)

  private def summarizeWeighted(
      sorted: Array[Double],
      w: Array[Double],
      b: Array[Int],
      D: Int,
      Dprime: Int,
      sig: Int,
      sweeps: Int,
      moves: Int,
      presorted: Boolean,
      solverUsed: String
  ): Result =
    val k = b.length - 1
    val n = sorted.length
    var Wtot = 0.0
    var Stot = 0.0
    var e0 = 0
    while e0 < n do
      Wtot += w(e0)
      Stot += w(e0) * sorted(e0)
      e0 += 1
    val center = Stot / Wtot
    val counts = new Array[Int](k)
    val means = new Array[Double](k)
    val upper = new Array[Double](k)
    var sse = 0.0
    var bin = 0
    while bin < k do
      val lo = b(bin)
      val hi = b(bin + 1)
      var W = 0.0
      var S = 0.0
      var e = lo
      while e < hi do
        W += w(e)
        S += w(e) * (sorted(e) - center)
        e += 1
      val m = if W == 0 then 0.0 else S / W
      if W > 0 then
        e = lo
        while e < hi do
          val dv = sorted(e) - center - m
          sse += w(e) * dv * dv
          e += 1
      counts(bin) = hi - lo
      means(bin) = center + m
      upper(bin) = sorted(hi - 1)
      bin += 1
    Result(IArray.unsafeFromArray(b.clone()), IArray.unsafeFromArray(counts),
      IArray.unsafeFromArray(means), IArray.unsafeFromArray(upper),
      sse, D, Dprime, sig, sweeps, moves, presorted, solverUsed)
end Lex
