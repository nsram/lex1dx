# lex1dx 0.1.4

* When the default Wilber solver's λ search fails (an empty interval, typically
  on collinear SSE*(k), for example a regular grid or `cluster(c(0,1,2,3), 3)`),
  lex1dx falls back to the Fisher / Ckmeans DP instead of throwing. If that
  fallback's back-pointer table would exceed 2 GiB, the table is not allocated
  and a clear error is raised. The result now includes `solver_used`, which is
  `"wilber"` or `"ckmeans"` and records which DP ran (when the DP is skipped,
  it mirrors the requested solver).
* The R `tolerance` argument is now checked to lie in (0, 1), matching the C++
  core (it previously accepted `tolerance >= 0` and then failed in C++).
* The maintainer ORCID is recorded in DESCRIPTION. Example output headers that
  still said "from lex1dx 0.1.1" now say 0.1.4.
