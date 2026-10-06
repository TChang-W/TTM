## Submission

This is a new submission.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Zhigang Li <zhigang.li@ufl.edu>'
  New submission

## Test environments

* local Ubuntu 24.04, R 4.6.1: `R CMD check --as-cran`

## Notes for reviewers

* Examples that fit the models take roughly 5-10 seconds each (they run
  delete-group jackknife replicates), so they are wrapped in `\donttest{}`.
