## Test environments

* Debian GNU/Linux 13 (trixie), x86_64, R 4.5.0 (2025-04-11),
  R CMD check --as-cran (including the PDF manual)

GitHub Actions (ubuntu-latest, macos-latest, windows-latest) and win-builder
were used for 0.1.3 and have not been re-run on this 0.1.4 tree; they should
be run before the CRAN upload. A local macOS --as-cran --no-manual check
was also used for 0.1.3 and has not been repeated here.

## R CMD check results

0 errors | 0 warnings | 2 notes

* CRAN incoming feasibility reports "New submission". This is an update
  from 0.1.3, which was still in incoming/newbies at the time of this check.
  It fixes a hard error in the default Wilber solver on some evenly spaced
  inputs (for example `cluster(c(0, 1, 2, 3), 3)`), falling back to the
  Fisher / Ckmeans DP, and aligns the R `tolerance` check with the C++ core.

* "unable to verify current time" came from this check environment (the
  clock could not be verified) and is not a package issue.

A spelling note, if CRAN reports one, would list possibly misspelled words
in DESCRIPTION. They are correct: 'Sriram' is the author's surname in the
citation, 'interp' is part of the solver name 'wilber-interp', and 'sig' is
the name of an argument (significant digits).
