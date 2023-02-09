/****************************************************************
ECON 237: Problem Set 2 Question 1
Estimating an Idiosyncratic Income Labor Process- Getting residual income
****************************************************************/

clear all
global in "/Users/reiko/Documents/heterogeneity/PSet2/data"
global out "/Users/reiko/Documents/heterogeneity/PSet2/data"

********* PSID data *********
use "$in/PSID_sampleC.dta", clear

* Keep relevant variables
keep year id68 hpersno hrace hyrsed hage hlabinc wlabinc hwage

* Get log of wage of head of household
gen lhwage = ln(hwage)
gen lhlabinc = ln(hlabinc)

*EDU OF THE HEAD;
gen hedu =.
replace hedu = 1 if hyrsed<12
replace hedu = 2 if hyrsed==12
replace hedu = 3 if hyrsed>=12 & hyrsed<16
replace hedu = 4 if hyrsed>=16

* Get a residual income. Estimate with a fourth order polynomial in age 
egen panelid = group(id68 hpersno)
reg lhlabinc i.year i.hrace i.hedu hyrsed hage c.hage#c.hage c.hage#c.hage#c.hage c.hage#c.hage#c.hage#c.hage
predict reslhlabinc, residuals

* Additional manipulations
drop if reslhlabinc== .
gen eff_hage = hage - 25

keep panelid reslhlabinc eff_hage
gsort panelid eff_hage
reshape wide reslhlabinc, i(panelid) j(eff_hage)
drop panelid
export delimited "$out/residual_logincome.csv", replace

clear
set obs 36
gen age = .
local i 1
foreach a of numlist 25/60 {
	replace age = `a' in `i'
	local ++i
}
gen a = exp(_b[_cons] + _b[hage] * age + _b[c.hage#c.hage] * age^2 + _b[c.hage#c.hage#c.hage] * age^3 + _b[c.hage#c.hage#c.hage#c.hage] * age^4 + _b[hyrsed] * 12 + _b[3.hedu])
drop age
export delimited "$out/age_profile.csv", replace


********* CSF data ***********
/* 
Pull x7001 (number of HH members) from the Full Public Data Set "p19i6.dta",
downloaded from: https://www.federalreserve.gov/econres/scfindex.htm

We did not include this file in the replication because it caused our zip file 
to exceed email attachment file limits.
*/
// use x7001 using "$in/p19i6.dta", clear
// rename x7001 PEOPLE
//
// export delimited "$out/people.csv", replace
