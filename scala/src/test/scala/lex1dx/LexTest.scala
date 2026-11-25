package lex1dx

class LexTest extends munit.FunSuite:

  test("empty input is rejected"):
    interceptMessage[IllegalArgumentException]("input is empty"):
      Lex.cluster(Array.empty[Double], 1)

  test("NaN and inf are rejected"):
    intercept[IllegalArgumentException]:
      Lex.cluster(Array(1.0, Double.NaN), 1)
    intercept[IllegalArgumentException]:
      Lex.cluster(Array(Double.PositiveInfinity), 1)
    intercept[IllegalArgumentException]:
      Lex.cluster(Array(0.0, Double.NegativeInfinity), 1)

  test("k < 1 is rejected"):
    intercept[IllegalArgumentException]:
      Lex.cluster(Array(1.0, 2.0), 0)

  test("single value"):
    for mode <- List(Mode.Tuned, Mode.Sig5) do
      val r = Lex.cluster(Array(3.5), 1, Options(mode = mode))
      assertEquals(r.k, 1)
      assertEquals(r.bounds.toList, List(0, 1))
      assertEquals(r.means(0), 3.5)
      assertEquals(r.sse, 0.0)
      assertEquals(r.assign(3.5), 0)

  test("all ties: k is clamped to 1"):
    val xs = Array.fill(20)(2.0)
    for mode <- List(Mode.Tuned, Mode.Sig5) do
      val r = Lex.cluster(xs, 8, Options(mode = mode))
      assertEquals(r.k, 1)
      assertEquals(r.bounds.toList, List(0, 20))
      assertEquals(r.means(0), 2.0)
      assertEquals(r.sse, 0.0)
      assertEquals(r.D, 1)
      assertEquals(r.Dprime, 1)

  test("k >= D is clamped to D"):
    val r = Lex.cluster(Array(1.0, 2.0, 3.0, 1.0, 2.0), 10)
    assertEquals(r.k, 3)
    assertEquals(r.D, 3)
    assertEquals(r.sse, 0.0)

  test("round-trip: shuffled vs sorted; assign covers the range"):
    val xs = Array.tabulate(200)(i => math.sin(0.17 * i) * 10.0 + 0.01 * i)
    val r = Lex.cluster(xs, 8)
    val sorted = xs.sorted
    val q = Lex.cluster(sorted, 8)
    assertEquals(r.bounds.toList, q.bounds.toList)
    assert(q.presorted)
    assert(!r.presorted)
    for x <- xs do
      val b = r.assign(x)
      assert(b >= 0 && b < r.k)
    assertEquals(r.assign(r.upper(0)), 0)
    assertEquals(r.assign(r.upper(r.k - 1)), r.k - 1)

  test("REE vs exact is under 0.01% on small data"):
    val xs = Array.tabulate(400)(i => math.sin(0.31 * i) + 0.05 * (i % 17) + ((i % 5) - 2) * 0.2)
    for k <- List(4, 8, 16) do
      val ex = Lex.exact(xs, k)
      assert(ex.sse > 0)
      for mode <- List(Mode.Tuned, Mode.Sig5) do
        val r = Lex.cluster(xs, k, Options(mode = mode))
        val ree = (r.sse - ex.sse) / ex.sse * 100.0
        assert(ree < 0.01, s"k=$k mode=$mode ree=$ree")
        assert(ree >= -1e-9, s"SSE better than exact? ree=$ree")

  test("default M is max(4096, 32k)"):
    assertEquals(Options.defaultM(64), 4096)
    assertEquals(Options.defaultM(128), 4096)
    assertEquals(Options.defaultM(256), 8192)
    assertEquals(Options.defaultM(512), 16384)

  test("coarse DP ignores a group's internal SSE consistently (sig5, M=4)"):
    // see the same test in cpp/tests/test_quality.cpp; optimum {0,3,5,6}, SSE 1.1746
    val xs = Array(0.08, 0.63, 1.57, 3.27, 3.55, 4.66)
    for mode <- List(Mode.Tuned, Mode.Sig5) do
      val r = Lex.cluster(xs, 3, Options(mode = mode, coarseGroups = 4))
      assertEquals(r.bounds.toList, List(0, 3, 5, 6), s"mode=$mode")
      assert(math.abs(r.sse - 1.1746) < 1e-12, s"mode=$mode sse=${r.sse}")
end LexTest
