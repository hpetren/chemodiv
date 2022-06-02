## **Second re-submission**

This is a re-submission. The following changes have been made:

* The use of options(warn=-1) has now been completely removed.

The test environments and check results are identical to the initial submission.


## **First re-submission**

This is a re-submission. The following changes have been made:

> I think there are too many blank spaces between some words in your
description text. Please check that.

Thank you for spotting this. The accidental extra spaces between words in the
description text in DESCRIPTION have been removed.

> If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file. [...]

There are currently no further references describing the package
(but this will be added once a paper describing the package is published).

> You are setting options(warn=-1) in your function. This is not allowed.
Please make sure that you do not change the user's options, par or
working directory. If you really have to do so within functions, please
ensure with an *immediate* call of on.exit() that the settings are reset
when the function is exited. [...]

The options(warn=-1) was included to manage a peculiar warning message
printed when using the function ChemmineR::pubchemCidToSDF() in our compDis()
function. We wanted to suppress that warning, as it is not relevant in our
case, and could be confusing for the user. The only way we managed to solve
this was to use options(warn=-1) before calling the function. As instructed,
we have now added an immediate call of on.exit() where applicable, and hope
that this is an acceptable solution.

The test environments and check results are identical to the initial submission.


## **Initial submission**

## Test environments

* local Windows 10 install, R 4.2.0
* Ubuntu 22.04 LTS, R 4.2.0
* win-builder (R-devel)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Reverse dependencies

None.
