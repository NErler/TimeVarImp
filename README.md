# Bayesian imputation of time-varying covariates in linear mixed models
This repository contains example syntax accompanying the paper and the syntax used for the simulation study

## Example trimester-specific weight gain
The file `trimester-spec-wg.R` contains the JAGS syntax used for the example in Section 5.2

## Syntax for the simulation study
The files ending with ...-batch.R can be run with the help of an .bat file in which first the working directory is specified
for example:

`cd C:\Documents\Simulation`

and then the commands

`"<path of R.exe>" CMD BATCH --vanilla "--args seedset=1" "<name of the ...-batch.R file>"`

`"<path of R.exe>" CMD BATCH --vanilla "--args seedset=2" "<name of the ...-batch.R file>"`

...

`"<path of R.exe>" CMD BATCH --vanilla "--args seedset=20" "<name of the ...-batch.R file>"`
