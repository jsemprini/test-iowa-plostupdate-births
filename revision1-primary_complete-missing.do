
clear all

import delimited "C:\Users\jsemprini\OneDrive - Des Moines University\4-DMU-Research\a-Water\final-plosupdate\analysis_ready_T1_allvars_v3.csv"

			      
****final restrictions****
drop if t1_max_finished_water>=10 & t1_max_finished_water!=.
drop if t1_mean_complete>=10 & t1_mean_complet!=.
drop if t1_mean_observed>=10 & t1_mean_observed!=.

drop if t1_mean_complete==.
drop if t1_mean_observed==.

drop if birthweight_g==.
keep if maternal_race_broad=="White"
keep if prenatal_by5==1

gen pnc5=0
replace pnc5=1 if prenatal_by5==1


* Primary categorical or binary covariates for all models
global controls  ///
				parity2 parity3plus ///
              meduc_hs meduc_gt_hs ///
               infant_male married 

****outcomes
gen lbw=0
replace lbw=1 if birthweight_g<2500

gen ptb=0
replace ptb=1 if gest_age_weeks<37

*******extreme outcomes*******
gen vlbw=0
replace vlbw=1 if birthweight_g<1500

gen vptb=0
replace vptb=1 if gest_age_weeks<32


****create cat*****
gen n_cat=0 if  t1_mean_complete!=.
replace n_cat=1 if t1_mean_complete>.1  & t1_mean_complete!=.
replace n_cat=2 if t1_mean_complete>=5  & t1_mean_complete!=.

reghdfe ptb i.n_cat i.($controls) c.maternal_age , vce(cluster county_fips) absorb(county_fips#birth_year  birth_year#conception_quarter)


cd "C:\Users\jsemprini\OneDrive - Des Moines University\4-DMU-Research\a-Water\final-plosupdate\stata"
****summary****


egen race_eth=group(maternal_race_broad)


log using "missing-mylog_tests-allyr-primary.log", replace text
*****test exposure on outcomes*****
* Start AFTER all other primary eligibility restrictions,
* but BEFORE excluding births with missing control variables.
* Unknown source codes must already be converted to missing.

local binary parity2 parity3plus meduc_hs meduc_gt_hs ///
             infant_male married

* Identify complete and incomplete covariate records.
egen nmiss_cov = rowmiss(`binary' maternal_age)
gen byte complete_cov = (nmiss_cov == 0)

tab complete_cov
misstable summarize `binary' maternal_age

* Describe whether missingness differs across exposure groups.
tab n_cat complete_cov, column
tabstat lbw vlbw ptb vptb, by(complete_cov) statistics(n mean)

* Construct filled controls and variable-specific missing indicators.
local cc_controls
local missing_controls

foreach v of local binary {
    assert inlist(`v', 0, 1) | missing(`v')

    gen byte miss_`v' = missing(`v')
    clonevar fill_`v' = `v'
    replace fill_`v' = 0 if missing(`v')

    local cc_controls `cc_controls' i.`v'
    local missing_controls `missing_controls' ///
        i.fill_`v' i.miss_`v'
}

* Continuous maternal age: placeholder plus missingness indicator.
gen byte miss_age = missing(maternal_age)
gen double fill_age = maternal_age
replace fill_age = 0 if missing(maternal_age)

* Compare the original categorical exposure models.
foreach y in ptb vptb lbw vlbw {

    display as text "`y': COMPLETE-CASE ANALYSIS"

    reghdfe `y' ib0.n_cat `cc_controls' c.maternal_age ///
        if complete_cov == 1, ///
        absorb(county_fips#birth_year birth_year#conception_quarter) ///
        vce(cluster county_fips)

    margins, dydx(n_cat)
    margins, dydx(n_cat) mcompare(sidak)

    display as text "`y': MISSING-INDICATOR ANALYSIS"

    reghdfe `y' ib0.n_cat `missing_controls' ///
        c.fill_age i.miss_age, ///
        absorb(county_fips#birth_year birth_year#conception_quarter) ///
        vce(cluster county_fips)

    margins, dydx(n_cat)
    margins, dydx(n_cat) mcompare(sidak)
}

log close