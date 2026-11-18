
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

reghdfe ptb i.n_cat i.($controls) c.maternal_age , vce(cluster county_fips) absorb(county_fips#birth_year  birth_year#conception_quarter)

keep if e(sample)

* For text dates such as 1983-01-15.
generate double gest_start = daily(t1_start, "YMD")
format gest_start %td

* Check that nonblank dates converted successfully.
assert !missing(gest_start) if strtrim(t1_start) != ""

generate future_q = qofd(gest_start)
format future_q %tq
gen n_five_c=0
replace n_five_c=1 if t1_mean_complete>=5


* Start with your current birth dataframe.
tempfile future_lookup

* Build the lookup BEFORE restricting birth years.
preserve

    keep county_fips gest_start n_cat n_five_c
    drop if missing(gest_start)

    * Quarter when the existing exposure window begins.
    generate future_q = qofd(gest_start)

    * Select the earliest gestational start in each county-quarter.
    bysort county_fips future_q (gest_start): ///
        keep if gest_start == gest_start[1]

    * Remove identical records sharing that start date.
    drop gest_start
    duplicates drop
    isid county_fips future_q

    rename n_cat    n_cat_future
    rename n_five_c n_five_c_future

    save `future_lookup'

restore

* Convert text gestational-start date.
generate double gest_start_date = daily(t1_start, "YMD")
format gest_start_date %td

* Estimate delivery date from gestational age.
generate double birth_date_est = ///
    gest_start_date + 7 * gest_age_weeks
format birth_date_est %td

* Target the quarter after estimated delivery.
generate future_qtr = qofd(birth_date_est) + 1
format future_qtr %tq

* Check agreement with recorded birth year.
generate birth_year_est = year(birth_date_est)
tabulate birth_year birth_year_est, missing

* Retain 1988 records in the lookup, but analyze births through 1987.
keep if inrange(birth_year, 1983, 1987)

* Pull existing future exposure categories from the same county.
merge m:1 county_fips future_q using `future_lookup', keep(master match) generate(future_match)

* Check availability; unmatched future exposures remain missing.
tabulate birth_year future_match, row
tabulate n_cat_future, missing
tabulate n_five_c_future, missing

cd "C:\Users\jsemprini\OneDrive - Des Moines University\4-DMU-Research\a-Water\final-plosupdate\stata"


log using "future-mylog_tests-allyr-primary.log", replace text
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








log using "future-emylog_tests-allyr-five-sens.log", replace text
*****test exposure on outcomes*****

estimates clear
***sens***

foreach y in ptb  lbw   {
	foreach x in n_five_c n_five_c_future {
qui: reghdfe `y' i.`x' i.($controls) c.maternal_age, vce(cluster county_fips) absorb(county_fips#birth_year  birth_year#conception_quarter )

eststo: margins , dydx(`x') post

qui: reghdfe `y' i.`x' i.($controls ) c.maternal_age , vce(cluster county_fips) absorb(county_fips#birth_year  birth_year#conception_quarter birth_year#conception_quarter county_fips#conception_quarter)

eststo: margins , dydx(`x') post

}

} 

esttab , b(4) ci(4)



log close