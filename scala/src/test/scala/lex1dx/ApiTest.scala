package lex1dx

class ApiTest extends munit.FunSuite:

  test("sig 10 is rejected with a hint to use sig=0"):
    val ex = intercept[IllegalArgumentException] {
      Lex.cluster(Array(1.0, 2.0, 3.0), 2, Options(sig = 10))
    }
    assert(ex.getMessage.contains("0..9"), ex.getMessage)
    assert(ex.getMessage.contains("sig=0"), ex.getMessage)

  test("negative sig other than the default follow-mode value is rejected"):
    intercept[IllegalArgumentException] {
      Lex.cluster(Array(1.0, 2.0, 3.0), 2, Options(sig = -2))
    }

  test("sig=5 matches mode Sig5 and raises precision on the tight-cluster case"):
    val xs = Array.fill(50)(0.0) ++
      Array.tabulate(50)(i => 1.000012 + 1e-10 * (i - 25)) ++
      Array.tabulate(50)(i => 1.000017 + 1e-10 * (i - 25)) ++
      Array.fill(50)(2.0)
    val r5 = Lex.cluster(xs, 4, Options(sig = 5))
    val alias = Lex.cluster(xs, 4, Options(mode = Mode.Sig5))
    assertEquals(r5.k, 4)
    assertEquals(r5.sig, 6)
    assertEquals(r5.Dprime, 4)
    assertEquals(r5.bounds.toList, List(0, 50, 100, 150, 200))
    assertEquals(alias.sig, 6)
    assertEquals(alias.bounds.toList, r5.bounds.toList)
    assertEquals(alias.sse, r5.sse)

  test("weights match expanded counts, and invalid weights are rejected"):
    val vals = Array(1.0, 2.0, 3.0)
    val wts = Array(1.0, 3.0, 1.0)
    val expanded = Array(1.0, 2.0, 2.0, 2.0, 3.0)
    val rw = Lex.cluster(vals, 2, wts)
    val re = Lex.cluster(expanded, 2)
    assertEquals(rw.sse, re.sse)
    assert(rw.means.sameElements(re.means), s"means ${rw.means.mkString(",")} vs ${re.means.mkString(",")}")
    val mass = Array.tabulate(rw.bounds.length - 1)(i => wts.slice(rw.bounds(i), rw.bounds(i + 1)).sum)
    val sizes = Array.tabulate(re.counts.length)(i => re.counts(i).toDouble)
    assert(mass.sameElements(sizes), s"mass ${mass.mkString(",")} vs sizes ${sizes.mkString(",")}")

    val neg = intercept[IllegalArgumentException] {
      Lex.cluster(vals, 2, Array(1.0, -1.0, 1.0))
    }
    assert(neg.getMessage.contains("weights must be finite and >= 0"), neg.getMessage)
    val nan = intercept[IllegalArgumentException] {
      Lex.cluster(vals, 2, Array(1.0, Double.NaN, 1.0))
    }
    assert(nan.getMessage.contains("weights must be finite and >= 0"), nan.getMessage)
    val len = intercept[IllegalArgumentException] {
      Lex.cluster(vals, 2, Array(1.0, 1.0))
    }
    assert(len.getMessage.contains("size must match"), len.getMessage)
    val zero = intercept[IllegalArgumentException] {
      Lex.cluster(vals, 2, Array(0.0, 0.0, 0.0))
    }
    assert(zero.getMessage.contains("positive finite total"), zero.getMessage)

  test("exact wilber and ckmeans solvers agree on SSE"):
    val rng = new scala.util.Random(1)
    val x = Array.fill(2000)(math.exp(rng.nextGaussian()))
    val rw = Lex.cluster(x, 8, Options(mode = Mode.Exact, solver = Solver.Wilber))
    val rc = Lex.cluster(x, 8, Options(mode = Mode.Exact, solver = Solver.Ckmeans))
    val denom = math.max(rw.sse, rc.sse)
    val rel = if denom == 0 then 0.0 else math.abs(rw.sse - rc.sse) / denom
    assert(rel < 1e-12, s"rel=$rel wilber=${rw.sse} ckmeans=${rc.sse}")

  test("sig Dprime is nondecreasing on a synthetic sample, and sig=0 matches tuned"):
    val rng = new scala.util.Random(1)
    val x = Array.fill(10000)(math.exp(rng.nextGaussian()))
    var prev = 0
    for d <- 1 to 9 do
      val r = Lex.cluster(x, 8, Options(sig = d))
      assertEquals(r.sig, d, s"sig=$d used ${r.sig}")
      assert(r.Dprime >= prev, s"sig=$d Dprime=${r.Dprime} prev=$prev")
      prev = r.Dprime
    val r0 = Lex.cluster(x, 8, Options(sig = 0))
    val rt = Lex.cluster(x, 8)
    assertEquals(r0.sig, 0)
    assertEquals(r0.Dprime, rt.Dprime)
    assertEquals(r0.bounds.toList, rt.bounds.toList)
    assertEquals(r0.sse, rt.sse)

  test("defaultM matches Options.defaultM"):
    assertEquals(Lex.defaultM(4), 4096)
    assertEquals(Lex.defaultM(64), Options.defaultM(64))
    assertEquals(Lex.defaultM(256), 8192)

  test("athletics20, k=5: wilber and ckmeans give the same partition"):
    val x = Datasets.load("athletics20")
    for (mode, label) <- List((Mode.Exact, "exact"), (Mode.Tuned, "tuned")) do
      val rw = Lex.cluster(x, 5, Options(mode = mode, solver = Solver.Wilber))
      val rc = Lex.cluster(x, 5, Options(mode = mode, solver = Solver.Ckmeans))
      assertEquals(rw.bounds.toList, rc.bounds.toList, label)
      val denom = math.max(rw.sse, rc.sse)
      val rel = if denom == 0 then 0.0 else math.abs(rw.sse - rc.sse) / denom
      assert(rel < 1e-12, s"$label rel=$rel wilber=${rw.sse} ckmeans=${rc.sse}")
end ApiTest
