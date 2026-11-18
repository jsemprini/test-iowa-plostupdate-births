
clear all

import delimited "C:\Users\jsemprini\OneDrive - Des Moines University\4-DMU-Research\a-Water\final-plosupdate\analysis_ready_T1_allvars_v3.csv"

			   
			   
****final restrictions****
drop if birthweight_g==.
keep if maternal_race_broad=="White"
*keep if prenatal_by5==1

drop if t1_mean_complete>=10 & t1_mean_complet!=.
drop if t1_mean_observed>=10 & t1_mean_observed!=.
drop if t1_mean_complete==.
drop if t1_mean_observed==.
drop if t1_max_finished_water>=10

gen pnc5=0
replace pnc5=1 if prenatal_by5==1

* Primary categorical or binary covariates for all models
global controls  ///
				parity2 parity3plus ///
              meduc_hs meduc_gt_hs ///
               infant_male married pnc5

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

gen n_five_o=0
replace n_five_o=1 if t1_mean_observed>=5

gen n_five_c=0
replace n_five_c=1 if t1_mean_complete>=5

cd "C:\Users\jsemprini\OneDrive - Des Moines University\4-DMU-Research\a-Water\final-plosupdate"
****summary****

egen id=group(county_fips birth_year conception_quarter)
gen white=0
replace white=1 if maternal_race_broad=="White"

log using "pnc5-emylog_tests-allyr-five-sens.log", replace text
*****test exposure on outcomes*****

estimates clear
***sens***

foreach y in ptb  lbw   {
	foreach x in n_five_c {
qui: reghdfe `y' i.`x' i.($controls) c.maternal_age if  pnc5==1, vce(cluster county_fips) absorb(county_fips#birth_year  birth_year#conception_quarter)

eststo: margins , dydx(`x') post

qui: reghdfe `y' i.`x' i.($controls ) c.maternal_age , vce(cluster county_fips) absorb(county_fips#birth_year  birth_year#conception_quarter)

eststo: margins , dydx(`x') post

qui: reghdfe `y' i.`x' i.($controls pnc5) c.maternal_age , vce(cluster county_fips) absorb(county_fips#birth_year  birth_year#conception_quarter)

eststo: margins , dydx(`x') post
}

} 

esttab , b(4) ci(4)



log close