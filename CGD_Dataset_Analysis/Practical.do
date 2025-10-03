* -------------------------------------------------------------------
* DO File for CGD Dataset Practical Session
* -------------------------------------------------------------------

* Load the dataset
use "C:\Users\AMCD\Downloads\fwddataforpracticalsession\CGD dataset.dta", clear

* Browse the data
br

* Describe the dataset structure
describe

* Tabulate the 'event' variable
tab event

* Execute temporary scripts (if needed)
do "C:\Users\AMCD\AppData\Local\Temp\STD44a8_000000.tmp"
do "C:\Users\AMCD\AppData\Local\Temp\STD44a8_000000.tmp"
do "C:\Users\AMCD\AppData\Local\Temp\STD44a8_000000.tmp"

* Reload dataset to ensure clean start
clear
use "C:\Users\AMCD\Downloads\fwddataforpracticalsession\CGD dataset.dta", clear

* Tabulate 'sequence' variable and cross-tabulate with 'event'
tab sequence
tab sequence event

* Keep only the first sequence
keep if sequence == 1

* Inspect dataset
describe
codebook
describe

* Set survival-time data
stset T1

* Try setting failure variable (some incorrect syntax attempts removed)
stset T1, failure(event)

* Reload dataset again for clean start
clear
use "C:\Users\AMCD\Downloads\fwddataforpracticalsession\CGD dataset.dta", clear
keep if sequence == 1

* Proper survival-time setup
stset T1, failure(event)

* Open help for stset if needed
help stset

* Cox proportional hazards models
stcox
stcox trtmt
stcox trtmt age

* Tabulate age for inspection
tab age

* Reload dataset again for extended stset example
clear
use "C:\Users\AMCD\Downloads\fwddataforpracticalsession\CGD dataset.dta", clear

* Setting survival data with entry and exit times
stset T1, failure(event) id(ID) enter(T0) exit(Time.)

* Cox regression with robust standard errors and strata
stcox trtmt, robust
stcox trtmt, robust strata(sequence)
stcox trtmt age, robust strata(sequence)
stcox trtmt age weight, robust strata(sequence)
stcox trtmt##sequence, robust strata(sequence)

* Cox regression with shared frailty
stcox trtmt, shared(ID)

* Attempt Weibull model with shared frailty
streg trtmt, dist(weibull) frailty(gamma) shared(ID)
