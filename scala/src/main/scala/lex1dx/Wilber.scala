package lex1dx

import scala.collection.mutable.ArrayBuffer

/** Weighted wilber-interp DP, ported from Grønlund kmeans1d commit 585386d
  * with the weighted interval cost cost = max(0, S2 - S1^2/W).
  */
private[lex1dx] object Wilber:

  final class IntervalSum(w: Array[Double], wx: Array[Double], wx2: Array[Double]):
    private val n = w.length
    private val prefixW = new Array[Double](n + 1)
    private val prefixS = new Array[Double](n + 1)
    private val prefixQ = new Array[Double](n + 1)
    locally {
      var W, S, Q = 0.0
      var i = 1
      while i <= n do
        W += w(i - 1); S += wx(i - 1); Q += wx2(i - 1)
        prefixW(i) = W; prefixS(i) = S; prefixQ(i) = Q
        i += 1
    }
    def queryW(i: Int, j: Int): Double = prefixW(j) - prefixW(i)
    def query(i: Int, j: Int): Double = prefixS(j) - prefixS(i)
    def querySq(i: Int, j: Int): Double = prefixQ(j) - prefixQ(i)
    /** cost of points [start, i] inclusive, as in C++ cost_interval_l2. */
    def costIntervalL2(start: Int, i: Int): Double =
      if i <= start then 0.0
      else
        val S = query(start, i + 1)
        val S2 = querySq(start, i + 1)
        val W = queryW(start, i + 1)
        if W == 0 then 0.0
        else
          val c = S2 - (S * S) / W
          if c < 0 then 0.0 else c
    def mean(i: Int, j: Int): Double =
      val W = queryW(i, j + 1)
      if W == 0 then 0.0 else query(i, j + 1) / W
    def size: Int = n

  /** @return (bounds in group index space of length k+1, cost) */
  def cluster(w: Array[Double], wx: Array[Double], wx2: Array[Double], k: Int): Array[Int] =
    val n = w.length
    if k >= n then
      return Array.tabulate(n + 1)(i => i)
    if k <= 1 then
      return Array(0, n)
    val is = new IntervalSum(w, wx, wx2)
    val solver = new Solver(is)
    val path = solver.compute(k)
    pathToBounds(path, n)

  private def pathToBounds(path: Array[Int], n: Int): Array[Int] =
    val p = path.sorted.distinct
    val b = ArrayBuffer.empty[Int]
    for gi <- p if gi >= 0 && gi <= n do b += gi
    if b.isEmpty || b(0) != 0 then b.insert(0, 0)
    if b.last != n then b += n
    b.toArray

  private final class Solver(is: IntervalSum):
    private val n = is.size
    private var lambda = 0.0
    private var f = new Array[Double](n + 1)
    private var bestleft = new Array[Int](n + 1)

    def compute(k: Int): Array[Int] =
      val (ok, _) = interpolationSearch(k)
      if !ok then
        throw new IllegalStateException("wilber interpolation search failed (empty lambda interval)")
      val path = ArrayBuffer.empty[Int]
      var m = n
      path += m
      while m != 0 do
        val prev = bestleft(m)
        path += prev
        m = prev
      path.toArray

    private def weight(i: Int, j: Int): Double =
      if i >= j then Double.MaxValue
      else is.costIntervalL2(i, j - 1) + lambda

    private def g(i: Int, j: Int): Double = f(i) + weight(i, j)

    private def interpolationSearch(k: Int): (Boolean, Int) =
      var loIntercept = 0.0
      var loK = n
      var hi = is.costIntervalL2(0, n - 1)
      var hiIntercept = hi
      var hiK = 1
      var cnt = 0
      while true do
        cnt += 1
        lambda = (hiIntercept - loIntercept) / (loK - hiK).toDouble
        val (valFound, kFound) = wilber()
        if kFound <= hiK || kFound >= loK then
          return (false, kFound)
        if kFound > k then
          loK = kFound
          loIntercept = valFound - loK * lambda
        else if kFound < k then
          hi = lambda
          hiK = kFound
          hiIntercept = valFound - hiK * lambda
        else
          return (true, kFound)
        if cnt > 1000 then
          return (false, kFound)
      (false, 0)

    private def wilber(): (Double, Int) =
      f = Array.fill(n + 1)(0.0)
      bestleft = Array.fill(n + 1)(0)
      f(0) = 0
      var c = 0
      var r = 0
      while c < n do
        val p = math.min(2 * c - r + 1, n)
        val (colMin, bl) = smawk(r, c, c + 1, p)
        var j = c + 1
        while j <= p do
          f(j) = colMin(j - (c + 1))
          bestleft(j) = bl(j - (c + 1))
          j += 1
        if c + 1 <= p - 1 then
          val (hCol, bl2) = smawk(c + 1, p - 1, c + 2, p)
          var j0 = p + 1
          j = p
          while j >= c + 2 do
            if hCol(j - (c + 2)) < f(j) then j0 = j
            j -= 1
          if j0 == p + 1 then c = p
          else
            f(j0) = hCol(j0 - (c + 2))
            bestleft(j0) = bl2(j0 - (c + 2))
            r = c + 1
            c = j0
        else c = p
      var m = n
      var length = 0
      while m > 0 do
        m = bestleft(m)
        length += 1
      (f(n), length)

    private def smawk(i0: Int, i1: Int, j0: Int, j1: Int): (Array[Double], Array[Int]) =
      val rows = Array.tabulate(i1 - i0 + 1)(t => i0 + t)
      val cols = Array.tabulate(j1 - j0 + 1)(t => j0 + t)
      val columnMinimaIdx = smawkInner(cols, 1, rows)
      val res = new Array[Double](columnMinimaIdx.length)
      val idxes = new Array[Int](columnMinimaIdx.length)
      var i = 0
      while i < res.length do
        idxes(i) = rows(columnMinimaIdx(i))
        res(i) = g(idxes(i), cols(i))
        i += 1
      (res, idxes)

    private def smawkInner(columns: Array[Int], e: Int, rows: Array[Int]): Array[Int] =
      val nn = columns.length
      val resultSize = (nn + e - 1) / e
      if rows.length == 1 then
        return Array.fill(resultSize)(0)
      val newRows = ArrayBuffer.empty[Int]
      val translate = ArrayBuffer.empty[Int]
      if resultSize < rows.length then
        var i = 0
        while i < rows.length do
          val r = rows(i)
          while newRows.nonEmpty &&
              g(r, columns(e * (newRows.length - 1))) <= g(newRows.last, columns(e * (newRows.length - 1))) do
            newRows.dropRightInPlace(1)
            translate.dropRightInPlace(1)
          if e * newRows.length < nn then
            newRows += r
            translate += i
          i += 1
      else
        var i = 0
        while i < rows.length do
          newRows += rows(i)
          translate += i
          i += 1
      if resultSize == 1 then
        return Array(translate(0))
      val rec = smawkInner(columns, 2 * e, newRows.toArray)
      val columnMinima = ArrayBuffer.empty[Int]
      columnMinima += translate(rec(0))
      var i = 1
      while i < rec.length do
        val from = rec(i - 1)
        val to = rec(i)
        val newColumn = 2 * i - 1
        columnMinima += from
        var r = from
        while r <= to do
          if g(newRows(r), columns(newColumn * e)) <= g(rows(columnMinima(newColumn)), columns(newColumn * e)) then
            columnMinima(newColumn) = translate(r)
          r += 1
        columnMinima += translate(to)
        i += 1
      if columnMinima.length < resultSize then
        val from = rec.last
        val newColumn = columnMinima.length
        columnMinima += translate(from)
        var r = from
        while r < newRows.length do
          if g(newRows(r), columns(newColumn * e)) <= g(rows(columnMinima(newColumn)), columns(newColumn * e)) then
            columnMinima(newColumn) = translate(r)
          r += 1
      columnMinima.toArray
  end Solver
end Wilber
