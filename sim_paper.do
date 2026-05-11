********************************************************************************
* sim_paper.do
* Causal Estimators Performance Simulation for the eltmle paper.
* Run with the updated eltmle.ado v4.0.3.  Outputs:
*   - sim_results.dta            : raw per-rep ATE/SE estimates
*   - sim_summary.csv            : per-method bias/relbias/coverage table
*   - ../Boxplot.pdf             : ATE distribution by method
*   - ../Relbias.pdf             : per-rep relative bias scatter + means
* Sample size 1000, 100 repetitions, true ATE computed analytically from a
* 10-million-observation panel (same DGP).  Set scheme sj for paper format.
********************************************************************************

clear all
set more off

local outdir "sim"
local simdir "`outdir'/sim"
cd "`simdir'"

* ----- True ATE from 10M-obs panel (deterministic with seed) -----
clear
set obs 10000000
set seed 777
gen w1 = round(runiform(1, 5))
gen w2 = rbinomial(1, 0.45)
gen w3 = round(runiform(0, 1) + 0.75*(w2) + 0.8*(w1))
recode w3 (5/6=1)
gen w4 = round(runiform(0, 1) + 1.2*(w2) + 0.2*(w1))
gen Y1 = (invlogit(-3 + 1 + 0.25*(w4) + 0.75*(w3) + 0.8*(w2)*(w4) + 0.05*(w1)))
gen Y0 = (invlogit(-3 + 0 + 0.25*(w4) + 0.75*(w3) + 0.8*(w2)*(w4) + 0.05*(w1)))
gen psi = Y1 - Y0
quietly mean psi
scalar true_ate = r(table)[1,1]
display "True Population ATE: " true_ate

local obs  1000
local reps 100

capture postclose ests
postfile ests int(repno) float(ATE SE Method) using sim_results, replace

scalar t1 = c(current_time)
qui {
    noi _dots 0, title("Running Monte Carlo Simulation...")
    forval rep = 1/`reps' {
        drop _all
        set seed `rep'
        set obs `obs'
        gen w1 = round(runiform(1, 5))
        gen w2 = rbinomial(1, 0.45)
        gen w3 = round(runiform(0, 1) + 0.75*(w2) + 0.8*(w1))
        recode w3 (5/6=1)
        gen w4 = round(runiform(0, 1) + 1.2*(w2) + 0.2*(w1))
        gen A  = (rbinomial(1,invlogit(-3 - 0.5*(w4) + 1.5*(w2) + 0.75*(w3) + 0.05*(w1) + 0.8*(w2)*(w4))))
        gen Y1 = (invlogit(-3 + 1 + 0.25*(w4) + 0.75*(w3) + 0.8*(w2)*(w4) + 0.05*(w1)))
        gen Y0 = (invlogit(-3 + 0 + 0.25*(w4) + 0.75*(w3) + 0.8*(w2)*(w4) + 0.05*(w1)))
        gen Y  = A*(Y1) + (1 - A)*Y0

        capture teffects ra (Y i.w1 i.w2 i.w3 i.w4) (A)
        if _rc == 0 post ests (`rep') (r(table)[1,1]) (r(table)[2,1]) (1)

        capture teffects ipw (Y) (A i.w1 i.w2 i.w3 i.w4)
        if _rc == 0 post ests (`rep') (r(table)[1,1]) (r(table)[2,1]) (2)

        capture teffects ipwra (Y i.w1 i.w2 i.w3 i.w4) (A i.w1 i.w2 i.w3 i.w4)
        if _rc == 0 post ests (`rep') (r(table)[1,1]) (r(table)[2,1]) (3)

        capture teffects aipw (Y i.w1 i.w2 i.w3 i.w4) (A i.w1 i.w2 i.w3 i.w4)
        if _rc == 0 post ests (`rep') (r(table)[1,1]) (r(table)[2,1]) (4)

        capture eltmle Y A w1 w2 w3 w4, tmle
        if _rc == 0 post ests (`rep') (r(ATEtmle)) (r(ATE_SE_tmle)) (5)

        capture eltmle Y A w1 w2 w3 w4, cvtmle cvfolds(10)
        if _rc == 0 post ests (`rep') (r(ATEtmle)) (r(ATE_SE_tmle)) (6)

        noi _dots `rep' 0
    }
}
postclose ests
scalar t2 = c(current_time)
display "Total processing time: " (clock(t2, "hms") - clock(t1, "hms")) / 1000 " seconds"

*===============================================================================
* Post-simulation analysis
*===============================================================================
use sim_results, clear
label define meth 1 "RA" 2 "IPW" 3 "IPWRA" 4 "AIPTW" 5 "TMLE" 6 "CVTMLE"
label values Method meth

scalar TRUEATE = `=true_ate'
gen double bias    = ATE - TRUEATE
gen double sq_err  = bias^2
gen double relbias = (abs(TRUEATE - ATE) / TRUEATE) * 100
gen double ci_l    = ATE - 1.96*SE
gen double ci_u    = ATE + 1.96*SE
gen byte   covered = (TRUEATE >= ci_l & TRUEATE <= ci_u)

* Summary table (printed and exported)
display ""
display "True ATE = " %9.4f TRUEATE
tabstat bias relbias covered, by(Method) stats(mean) format(%9.4f)

preserve
collapse (mean) bias relbias coverage=covered, by(Method)
list, sep(0)
export delimited using sim_summary.csv, replace
restore

*===============================================================================
* Boxplot.pdf -- ATE distribution by method
*===============================================================================
set scheme sj
graph box ATE, over(Method) ///
    yline(`=true_ate', lcolor(red) lpattern(dash)) ///
    title("ATE Estimator Distribution") ///
    subtitle("Red dashed line: true ATE") ///
    note("Sample size 1,000; `reps' repetitions") ///
    graphregion(color(white)) bgcolor(white) plotregion(fcolor(white))
graph export "`outdir'/Boxplot.pdf", as(pdf) replace

*===============================================================================
* Relbias.pdf -- per-rep relative bias scatter with method means
*===============================================================================
by Method, sort: egen meanrelbias = mean(relbias)
twoway (scatter relbias Method, msymbol(oh) mcolor(black)) ///
       (scatter meanrelbias Method, msymbol(D) mcolor(red) msize(medium)), ///
       xlabel(1 "RA" 2 "IPTW" 3 "IPTW-RA" 4 "AIPTW" 5 "TMLE" 6 "CVTMLE") ///
       legend(off) graphregion(color(white)) bgcolor(white) plotregion(fcolor(white)) ///
       ytitle("Relative bias of ATE (%)") xtitle("Method")
graph export "`outdir'/Relbias.pdf", as(pdf) replace

display "Done."
exit, clear
