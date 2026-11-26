## Test environments

* local macOS, R CMD check --as-cran --no-manual
* local Debian GNU/Linux 13 (trixie), x86_64, R 4.5.0, R CMD check --as-cran
  (including the PDF manual)
* GitHub Actions: ubuntu-latest, macos-latest, windows-latest
* win-builder: R-devel (2026-11-17 r90926 ucrt)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

The note also lists possibly misspelled words in DESCRIPTION. They are
correct: 'Sriram' is the author's surname in the citation, 'interp' is
part of the solver name 'wilber-interp', and 'sig' is the name of an
argument (significant digits).
