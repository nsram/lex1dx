## Test environments

* local Debian GNU/Linux 13 (trixie), x86_64, R 4.5.0 (2025-04-11)
  R CMD check --as-cran, including the PDF manual
* GitHub Actions is already configured for ubuntu-latest, macos-latest,
  and windows-latest (R CMD check --as-cran --no-manual). Confirm a green
  run on the commit you submit.
* win-builder (R-devel) and R-hub have not been run yet; please run them
  before submitting (see the submit notes).

## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new submission.
* checking for future file timestamps: unable to verify current time

The second note is from the local machine, which could not reach a time
server. It is not expected on CRAN.
