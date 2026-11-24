package lex1dx

import scala.collection.mutable.ArrayBuffer

/** Tuned batch Collapse: 0.75-quantile threshold per pass, heap finish once
  * (bins - target) <= 16384. Mirrors cpp/src/tuned_collapse.hpp.
  */
private[lex1dx] object Collapse:
  private val Cutover = 16384
  private val Q = 0.75

  def apply(a: Atoms, target: Int): Array[Int] =
    val d = a.d
    if target >= d then return a.start.clone()
    val prev = Array.tabulate(d)(j => j - 1)
    val next = Array.tabulate(d)(j => j + 1)
    def cost(x: Int): Double = a.mergeCost(a.start(prev(x)), a.start(x), a.start(next(x)))

    var bins = d

    def heapFinish(): Unit =
      val heap = MinHeapDir(d)
      var x = next(0)
      while x < d do
        heap.append(x, cost(x))
        x = next(x)
      heap.heapify()
      while bins > target && heap.size > 0 do
        val y = heap.pop()
        if prev(y) >= 0 && next(prev(y)) == y then
          val p = prev(y)
          val q = next(y)
          next(p) = q
          if q < d then prev(q) = p
          bins -= 1
          if p > 0 && heap.has(p) then heap.update(p, cost(p))
          if q < d && heap.has(q) then heap.update(q, cost(q))

    var done = false
    while bins > target && !done do
      val need = bins - target
      if need <= Cutover then
        heapFinish()
        done = true
      else
        val ids = ArrayBuffer.empty[Int]
        val costs = ArrayBuffer.empty[Double]
        var x = next(0)
        while x < d do
          ids += x
          costs += cost(x)
          x = next(x)
        if ids.isEmpty then done = true
        else
          val tmp = costs.toArray
          val qi = math.min(tmp.length - 1, (Q * tmp.length.toDouble).toInt)
          java.util.Arrays.sort(tmp)
          val thr = tmp(qi)
          var merged = 0
          var skip = false
          var i = 0
          while i < ids.length && bins > target do
            if skip then skip = false
            else
              val y = ids(i)
              if prev(y) >= 0 && next(prev(y)) == y && costs(i) <= thr then
                val p = prev(y)
                val qq = next(y)
                next(p) = qq
                if qq < d then prev(qq) = p
                bins -= 1
                merged += 1
                skip = true
            i += 1
          if merged == 0 then
            heapFinish()
            done = true

    val b = new Array[Int](target + 1)
    var j = 0
    var bin = 0
    while bin < d && j < target do
      b(j) = a.start(bin)
      j += 1
      bin = next(bin)
    b(target) = a.n
    b

  /** Binary min-heap with decrease/increase key, same tie-break as the C++ MinHeapDir. */
  private final class MinHeapDir(capacity: Int):
    private val key = new Array[Double](capacity)
    private val heap = new Array[Int](capacity)
    private val pos = Array.fill(capacity)(-1)
    private var sz = 0
    def size: Int = sz
    def has(id: Int): Boolean = id >= 0 && id < pos.length && pos(id) >= 0
    def append(id: Int, k: Double): Unit = { put(id, k, sz); sz += 1 }
    def heapify(): Unit =
      var i = sz / 2 - 1
      while i >= 0 do { down(i); i -= 1 }
    def pop(): Int =
      val top = heap(0)
      sz -= 1
      if sz > 0 then { put(heap(sz), key(sz), 0); down(0) }
      pos(top) = -1
      top
    def update(id: Int, k: Double): Unit =
      val i = pos(id)
      val old = key(i)
      key(i) = k
      if k < old then up(i)
      else if k > old then down(i)
    private def put(id: Int, k: Double, i: Int): Unit =
      heap(i) = id; key(i) = k; pos(id) = i
    private def less(ka: Double, a: Int, kb: Double, b: Int): Boolean =
      ka < kb || (ka == kb && a < b)
    private def up(i0: Int): Unit =
      var i = i0
      val id = heap(i)
      val k = key(i)
      while i > 0 && less(k, id, key((i - 1) / 2), heap((i - 1) / 2)) do
        val p = (i - 1) / 2
        put(heap(p), key(p), i)
        i = p
      put(id, k, i)
    private def down(i0: Int): Unit =
      var i = i0
      val id = heap(i)
      val k = key(i)
      var cont = true
      while cont do
        val l = 2 * i + 1
        val c =
          if l + 1 < sz && less(key(l + 1), heap(l + 1), key(l), heap(l)) then l + 1 else l
        if c < sz && less(key(c), heap(c), k, id) then { put(heap(c), key(c), i); i = c }
        else cont = false
      put(id, k, i)
end Collapse
