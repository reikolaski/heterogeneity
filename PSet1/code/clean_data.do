/**********************
ECON 237: Problem Set 1 Question 1
Estimating an Idiosyncratic Income Labor Process- Getting residual income
********************/

clear all
global in "/Users/reiko/Documents/heterogeneity/PSet1/data"
global out "/Users/reiko/Documents/heterogeneity/PSet1/data"

use "$in/PSID_sampleC.dta", clear

* Keep relevant variables
keep year id68 hpersno hrace hyrsed hage hlabinc wlabinc hwage

* Get log of wage of head of household
gen lhwage = ln(hwage)

*EDU OF THE HEAD;
gen hedu =.
replace hedu = 1 if hyrsed<12
replace hedu = 2 if hyrsed==12
replace hedu = 3 if hyrsed>=12 & hyrsed<16
replace hedu = 4 if hyrsed>=16

* Get a residual income. Estimate with a fourth order polynomial in age 
egen panelid = group(id68 hpersno)
reg lhwage i.year i.hrace i.hedu hyrsed hage c.hage#c.hage c.hage#c.hage#c.hage c.hage#c.hage#c.hage#c.hage
predict reslhwage, residuals

* Additional manipulations
drop if reslhwage== .
gen eff_hage = hage - 25

keep panelid reslhwage eff_hage
gsort panelid eff_hage
reshape wide reslhwage, i(panelid) j(eff_hage)
drop panelid
export delimited "$out/residual_logwage.csv", replace
