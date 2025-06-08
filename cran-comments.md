## Comments

The chemodiv 0.3.1 release mainly fixes the cause of an an error that was
temporarily present for the package, due to an unavailable internet resource. 
The fix ensures that functions "fail gracefully with an informative message",
instead of causing an error in such cases.

## Test environments

* local Windows 10 install, R 4.5.0
* win-builder (R-devel)
* ubuntu 24.04.2 (on Github Actions), R 4.4.3, R 4.5.0, R-devel
* windows-latest (on Github Actions), R 4.5.0
* macOS-latest (on Github Actions), R 4.5.0

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

There are no downstream dependencies for this package.
