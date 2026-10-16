## Submission

This is a new submission.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Zhigang Li <zhigang.li@ufl.edu>'
  New submission

## Test environments

* local Windows 11 x64, R 4.4.2
* win-builder: R-devel (2026-10-14 r90673 ucrt)
* R-hub: ubuntu-release (R 4.6.1)
* R-hub: donttest (Linux, R-devel, examples run with `--run-donttest`)

## Notes for reviewers

* Examples that fit the models take roughly 5-10 seconds each (they run
  delete-group jackknife replicates), so they are wrapped in `\donttest{}`.
