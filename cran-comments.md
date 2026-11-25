## Submission

islasso 1.7.0. Main changes (see NEWS.md):

* the Fortran numerical core has been replaced by a new core written in C (BLAS/LAPACK of R);
* any GLM family accepted by glm() can be used;
* several bug fixes, new vignettes and unit tests;
* glmnet moved from Depends to Imports; dependencies on cli and gridExtra removed.

## Test environments

* local macOS (aarch64), R 4.6.1
* win-builder (devel and release)
* GitHub Actions: macOS, Windows, Ubuntu (devel, release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 0 notes

(On the local macOS machine the only warning is the linker message
"search path '/opt/gfortran/lib/...' not found", due to the absence of the CRAN gfortran
toolchain on that machine; it does not appear on the CRAN-like check services.)

## Reverse dependencies

All reverse dependencies were checked with no new problems.
