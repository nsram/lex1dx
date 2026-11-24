# lex1dx (R)

R bindings for lex1dx. From this directory, in a full repository checkout:

```
sh tools/sync_core.sh
R CMD build .
R CMD check --no-manual --no-vignettes lex1dx_0.1.0.tar.gz
```

`configure` copies the C++ core into `src/` so the package is self-contained. The core is compiled with `-ffp-contract=off`, so results are bit-identical to the C++ library.

```r
library(lex1dx)
r <- cluster(x, 64L)
s <- cluster(x, 64L, mode = "sig5")
```
