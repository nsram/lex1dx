# lex1dx (Scala)

Scala 3.3 and sbt 1.10, from this directory:

```
sbt test
```

```scala
import lex1dx.{Lex, Options, Mode, Solver, Datasets}

val x = Datasets.load("athletics20")
val r = Lex.cluster(x, 4)                                 // tuned, default M
val s = Lex.cluster(x, 4, Options(sig = 5))               // group to 5 significant digits
val e = Lex.cluster(x, 4, Options(mode = Mode.Exact, solver = Solver.Ckmeans))
val w = Lex.cluster(vals, 4, counts)                      // observation weights
```

`Lex.cluster(xs, k)` and `Lex.cluster(xs, k, options)` are the unweighted entry points.

`Lex.cluster(xs, k, weights)` and `Lex.cluster(xs, k, weights, options)` take a non-negative finite weight per value.

`Options` has `mode` (`Tuned`, `Sig5`, `Exact`), `sig` (integer 0..9, or -1 to follow `mode`), `solver` (`Wilber`, `Ckmeans`), `coarseGroups` (`M`; 0 uses `defaultM(k)`), and `tolerance` (default `1e-11`).

Both solvers are exact and give the same optimum.

Wilber is the default and is faster as k grows, because wilber-interp is O(k N) in style while ckmeans is a divide-and-conquer O(k N log N) DP.

ckmeans is the original LEX exact DP and serves as a cross-check.

`Mode.Sig5` is an alias for `sig = 5`.

If grouping at the requested `sig` yields fewer than `k` groups, precision is raised so `k` is kept.

`Result.sig` is the sig actually used, and `Result.Dprime` is the group count after that grouping.

`Result.counts` and `Result.means` are the bin sizes and centers.

`Lex.defaultM(k)` is `max(4096, 32 * k)`.

`Datasets.load` reads the bundled athletics time lists (`athletics20` and `worldathletics-mens-100m-2379`).

Source: World Athletics 100m all-time lists, retrieved January 2024.

Sample Data Attribution: The sample 100m timing datasets (top 20 performances and 2379 all time bests) contain historical athletic performance records compiled from publicly available statistics published by World Athletics (https://worldathletics.org/). These datasets are provided strictly as static benchmark/example inputs for testing and demonstrating the algorithms in this package. This project is not affiliated with, endorsed by, or sponsored by World Athletics.

The Python and R bindings can draw a horizontal bin table (`histogram()`).

This Scala package has no plot dependency.

A text table with the same rows can be built from a `Result` as follows.

```scala
val x = Datasets.load("athletics20")
val r = Lex.cluster(x, 5)
val sorted = x.sorted
val sst = { val m = x.sum / x.length; x.map(v => { val d = v - m; d * d }).sum }
val r2 = if sst == 0 then 1.0 else 1.0 - r.sse / sst
val mse = r.sse / x.length
println(s"Squared correlation r² of ${r.k}-bin model with data = $r2")
println(s"MSE of ${r.k}-bin model with data = $mse")
for i <- 0 until r.k do
  val lo = sorted(r.bounds(i))
  val hi = r.upper(i)
  val n = r.counts(i)
  val pct = 100.0 * n / x.length
  val mean = r.means(i)
  println(f"$lo%.5g - $hi%.5g  [n=$n, $pct%%, $mean%.5g]")
```
