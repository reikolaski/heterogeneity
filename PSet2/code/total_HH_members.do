// use x7001 x101 using "/Users/reiko/Documents/heterogeneity/PSet2/data/p19i6.dta", clear
// rename x101 PEOPLE_HHL
// rename x7001 PEOPLE

use x7001 using "/Users/reiko/Documents/heterogeneity/PSet2/data/p19i6.dta", clear
rename x7001 PEOPLE

export delimited "/Users/reiko/Documents/heterogeneity/PSet2/data/people.csv", replace
