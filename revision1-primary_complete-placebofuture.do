
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

gen n_five_c=0
replace n_five_c=1 if t1_mean_complete>=5

reghdfe ptb i.n_cat i.($controls) c.maternal_age , vce(cluster county_fips) absorb(county_fips#birth_year  birth_year#conception_quarter)

keep if e(sample)


*==============================================================================
* FUTURE-EXPOSURE PLACEBO (CORRECTED)
*
* Each focal birth is assigned the first-trimester exposure category of the
* earliest eligible pregnancy that STARTED in the same county in the quarter
* that falls future_offset quarters after the focal birth's estimated
* delivery (4 = one year after delivery).
*
* Correction to the previous version: the lookup was merged on future_q,
* which in the focal data was the quarter of the focal birth's OWN
* gestational start. Every birth therefore received the exposure of a
* pregnancy that began in its own start quarter (a contemporaneous exposure),
* and the target quarter that had been computed (future_qtr) was never used.
* The merge key is now target_q, defined below for both files, and an assert
* stops the do-file if the target quarter is not after the focal pregnancy.
*==============================================================================

local future_offset 4

* For text dates such as 1983-01-15.
generate double gest_start = daily(t1_start, "YMD")
format gest_start %td
assert !missing(gest_start) if strtrim(t1_start) != ""

* 1. Build the lookup from ALL eligible births (including 1988) BEFORE
*    restricting the focal birth years.
tempfile future_lookup
preserve
    keep county_fips gest_start n_cat n_five_c
    drop if missing(gest_start)

    * Quarter in which each pregnancy's first-trimester window begins.
    generate target_q = qofd(gest_start)
    format target_q %tq

    * Keep the earliest gestational start in each county-quarter.
    bysort county_fips target_q (gest_start): ///
        keep if gest_start == gest_start[1]

    * Remove identical records sharing that start date.
    drop gest_start
    duplicates drop
    isid county_fips target_q

    rename n_cat    n_cat_future
    rename n_five_c n_five_c_future
    save `future_lookup'
restore

* 2. Target quarter for each focal birth.
generate double birth_date_est = gest_start + 7 * gest_age_weeks
format birth_date_est %td

generate target_q = qofd(birth_date_est) + `future_offset'
format target_q %tq

* Guard against the original error: the target must come after the focal
* pregnancy's own start quarter.
generate own_start_q = qofd(gest_start)
format own_start_q %tq
assert target_q > own_start_q

generate birth_year_est = year(birth_date_est)

* Retain 1988 records in the lookup, but analyze births through 1987.
keep if inrange(birth_year, 1983, 1987)

* 3. Pull the future exposure category from the same county.
merge m:1 county_fips target_q using `future_lookup', keep(master match) generate(future_match)


cd "C:\Users\jsemprini\OneDrive - Des Moines University\4-DMU-Research\a-Water\final-plosupdate\stata"

log using "future-corrected-mylog_tests-allyr-primary.log", replace text

* Diagnostics for the supplement notes.
tabulate birth_year birth_year_est, missing
tabulate birth_year future_match, row

* 4. Analyze only matched births so the current- and future-exposure models
*    use the identical sample.
keep if future_match == 3
count

* How often current and future categories agree (should be well below 100%).
tabulate n_cat n_cat_future, row
correlate n_five_c n_five_c_future

*****test exposure on outcomes*****
qui: reghdfe lbw i.n_cat i.($controls) c.maternal_age , vce(cluster county_fips) absorb(county_fips#birth_year  birth_year#conception_quarter)
margins , dydx(n_cat)
margins , dydx(n_cat) mcompare(sidak)

qui: reghdfe lbw i.n_cat_future i.($controls) c.maternal_age , vce(cluster county_fips) absorb(county_fips#birth_year  birth_year#conception_quarter)
margins , dydx(n_cat_future)
margins , dydx(n_cat_future) mcompare(sidak)


qui: reghdfe ptb i.n_cat i.($controls) c.maternal_age , vce(cluster county_fips) absorb(county_fips#birth_year  birth_year#conception_quarter)
margins , dydx(n_cat)
margins , dydx(n_cat) mcompare(sidak)

qui: reghdfe ptb i.n_cat_future i.($controls) c.maternal_age , vce(cluster county_fips) absorb(county_fips#birth_year  birth_year#conception_quarter)
margins , dydx(n_cat_future)
margins , dydx(n_cat_future) mcompare(sidak)


qui: reghdfe vlbw i.n_cat i.($controls) c.maternal_age , vce(cluster county_fips) absorb(county_fips#birth_year  birth_year#conception_quarter)
margins , dydx(n_cat)
margins , dydx(n_cat) mcompare(sidak)

qui: reghdfe vlbw i.n_cat_future i.($controls) c.maternal_age , vce(cluster county_fips) absorb(county_fips#birth_year  birth_year#conception_quarter)
margins , dydx(n_cat_future)
margins , dydx(n_cat_future) mcompare(sidak)


qui: reghdfe vptb i.n_cat i.($controls) c.maternal_age , vce(cluster county_fips) absorb(county_fips#birth_year  birth_year#conception_quarter)
margins , dydx(n_cat)
margins , dydx(n_cat) mcompare(sidak)

qui: reghdfe vptb i.n_cat_future i.($controls) c.maternal_age , vce(cluster county_fips) absorb(county_fips#birth_year  birth_year#conception_quarter)
margins , dydx(n_cat_future)
margins , dydx(n_cat_future) mcompare(sidak)

log close


log using "future-corrected-emylog_tests-allyr-five-sens.log", replace text
*****test exposure on outcomes*****
estimates clear
***sens***
foreach y in ptb  lbw   {
	foreach x in n_five_c n_five_c_future {
qui: reghdfe `y' i.`x' i.($controls) c.maternal_age, vce(cluster county_fips) absorb(county_fips#birth_year  birth_year#conception_quarter)
eststo: margins , dydx(`x') post
qui: reghdfe `y' i.`x' i.($controls) c.maternal_age , vce(cluster county_fips) absorb(county_fips#birth_year  birth_year#conception_quarter county_fips#conception_quarter)
eststo: margins , dydx(`x') post
}
}
esttab , b(4) ci(4)
log close
