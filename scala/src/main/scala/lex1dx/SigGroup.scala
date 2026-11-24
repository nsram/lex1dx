package lex1dx

/** 5-significant-digit grouping of an already-sorted array.
  * The key is nearbyint(x / 10^(e-s+1)) * 10^(e-s+1), the same rule as the C++ `siggroup`.
  * Group means are accumulated in IEEE-754 double. The C++ pack uses `long double` for the
  * running sum, which may differ on some platforms; see the README.
  */
private[lex1dx] object SigGroup:
  def sigKey(x: Double, s: Int): Double =
    if x == 0.0 then 0.0
    else
      val e = math.floor(math.log10(math.abs(x))).toInt
      val q = math.pow(10.0, (e - s + 1).toDouble)
      math.rint(x / q) * q

  final case class Groups(start: Array[Int], mean: Array[Double], w: Array[Double])

  def group(sorted: Array[Double], sig: Int): Groups =
    val n = sorted.length
    val startB = Array.newBuilder[Int]
    val meanB = Array.newBuilder[Double]
    val wB = Array.newBuilder[Double]
    var i = 0
    while i < n do
      val kk = sigKey(sorted(i), sig)
      var j = i
      var s = 0.0
      while j < n && sigKey(sorted(j), sig) == kk do
        s += sorted(j)
        j += 1
      startB += i
      meanB += s / (j - i).toDouble
      wB += (j - i).toDouble
      i = j
    startB += n
    Groups(startB.result(), meanB.result(), wB.result())
end SigGroup
