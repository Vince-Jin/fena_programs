// special variables
/*
_st
_d
_t
_t0
*/

capture program drop survplot_fena
program define survplot_fena

qui {
	
	syntax [, np sp fp var(string) title(string) name(string)]
	
	// detect if stset
	// detect variable list
	if 1 {
		local set = 1
		local test _st _d _t _t0
		foreach i in `test' {
			capture confirm variable `i'
			if _rc {
				local set = 0
			}
		}
		
		// if set is 0, terminate and warn user to set up stset
		if (`set' == 0) {
			noi di as error "No survival dataset being established" ", Please use stset before plotting"
			noi di as error "If stset was called, please make sure key variables of _t, _t0, _d, _st exist"
			noi di in g "Program terminated"
			exit
		} 
		
		tempfile pls
		save `pls', replace
		global plot_set `pls'
		
		// detect discrete or continuous
		preserve
		capture drop _dif
		gen _dif = _t - _t0
		levelsof _dif, local(dif_levels)
		if (strlen("`dif_levels'") != 1) {
			global dis = 0
		}
		else if (strlen("`dif_levels'") == 1) {
			global dis = 1
		}
		
	}
	
	// detect name
	if 3 {
		if ("`name'" == "") {
			local name : di "Cumulative Incidence Graph"
		}
	}
	
	// detect variable type
	if 6 {
		qui var_type, var(`var')
			
		while ("${con}" != "") {
			noi di "This following variables were potentially continuous variables:"
			noi di "${con}"
			noi di "To generate survival plots, variables need to be categorical."
			noi di "Please confirm using following codes to proceed"
			noi di "1 - proceed with default grouping for continuous variables (quartiles)"
			noi di "2 - enforce certain variables as categorical for wrong detection results"
			noi di "3 - exit the program to manually group to categories", _request(categorical)

			if ("${categorical}" == "1") {
				global con_c
				foreach i in ${con} {
					
					local var_helper : di "`i'"
					local var_helper : di "`var_helper'_cat"
					capture program drop `var_helper'
					gen `var_helper' = .
					
					sum `i', detail
					local p25 = r(p25)
					local p50 = r(p50)
					local p75 = r(p75)
					
					replace `var_helper' = 1 if `i' < `p25'
					replace `var_helper' = 2 if inrange(`i', `p25', `p50')
					replace `var_helper' = 3 if inrange(`i', `p50', `p75')
					replace `var_helper' = 4 if `i' > `p75' & !missing(`i')
					
					global con_c : di "`var_helper'" " " "${con_c}"
				}
				global con_c : di strtrim(stritrim("${con_c}"))
			}
			else if ("${categorical}" == "2") {
				noi di "Please indicate variables to be forced as categorical variables"
				noi di "Use space to separate variables"
				noi di "e.g.: to force variable a, b, and, c -> a b c", _request(cat_force)
				
				foreach i in ${cat_force} {
					global con : di strtrim(subinstr("${con}", "`i'", "", .))
					global cat : di "${cat}" " " "`i'"
				}
				
				global cat : di strtrim(stritrim("$cat"))
			}
			else if ("${categorical}" == "3") {
				noi di as error "User Requested To Abort Plotting"
				noi di in g " "
				exit
			}
		}
		noi di "variable check cleared"
		local vars ${bin} ${cat} ${con_c}
	}
	
	// get title
	if 4 {
		local title
		if ("`title'" != "") {
			local title : di "title(" "`title'" ")"
		} 
		else if ("`title'" == "") {
			local title : di "Cumulative Incidence"
		}
	}
	
	// need to first determine if by is a categorical variable
	
	// non-parametric (Kaplan-Meier)
	if ("`np'" == "np") {
		noi di "Generating non-parametric plots"
		
		foreach i in `vars' {
			use ${plot_set}, clear
			noi di "Variable: `i'"
			sts graph, title("`title'") by(`i')
			graph export "`name'_`i'_np.png", as(png) replace
		}
		noi di "Function completed"
		noi di " "
	}
	
	// semi-parametric
	// use the largest cat as default baseline
	if ("`sp'" == "sp") {
		noi di "Generating semi-parametric plots"
		
		// detect coefplot pack
		capture coefplot
		if _rc {
			noi di "Missing required Stata package: coefplot for semi-parametric method"
			noi di "Function terminated"
			exit
		}		
		
		// discrete assumption warning
		if (${dis} == 1) {
			noi di as error "WARNING: The semi-parametric tool of this program uses Cox Proportional Hazard model for estimation."
			noi di as error "Current survival dataset was detected to have discrete timeline. Assumption checks may need to be completed."
			noi di in g " "
		}
		foreach i in `vars' {
			use ${plot_set}, clear
			noi di "Variable: `i'"
			
			levelsof `i', local(i_val)
			global i_v `i_val'
			// count unique values of i to create sector variables
			local i_n = strlen("`i_val'") + 1 - strlen(subinstr("`i_val'", " ", "", .))
			noi freq_detect, var(`i')
			
			// check if sv_pos is the first level
			if (${sv_pos} == 1) {
				noi di "The most frequent value for variable `i' was detected to be the same as reference level of model"
				noi di "The second frequent value will be used as default for prediction"
				noi di "Please confirm to proceed"
				noi di "(Enter 0 or N to force most frequent value as default)", _request(confirm)
				
				if !(("${confirm}" == "0") | (strupper(substr("${confirm}", 1, 1)) == "N")) {
					replace `i' = . if `i' == ${sv_val}
					noi freq_detect, var(`i')
					use ${plot_set}, clear
				}
			}
			
			/*
			capture drop temp
			egen temp = mode(`i')
			levelsof temp, local(sv_val)
			// identify on which place sv = 1
			local sv_pos = 1
			local sv_helper = 0
			local counter = 1
			tokenize `i_val'
			while (`sv_helper' == 0) {
				if (`sv_val' != ``sv_pos'') {
					local sv_pos = `sv_pos' + 1
				}
				else if (`sv_val' == ``sv_pos'') {
					local sv_helper = 1
				}
			} 
			*/
			capture drop s0
			stcox i.`i', basesurv(s0)
			matrix define cox=r(table)
			noi di "Semiparametric Cox Regression Result Matrix: "
			noi matrix list cox
			// make the corret place 1 and other 0 for sv
			local sv_string
			forvalues xx = 1/`i_n' {
				if (`xx' != ${sv_pos}) {
					local temp = 0
				}
				else if (`xx' == ${sv_pos}) {
					local temp = 1
				}
				local sv_string : di "`sv_string'" " " "`temp'"
			}
			local sv_string : di subinstr(strtrim("`sv_string'"), " ", ", ", .)
			matrix sv = (`sv_string')
			noi matrix list sv
			matrix cox2 = cox'
			noi matrix list cox2
			matrix risk_score = sv * cox'
			noi di "Risk Score Matrix: "
			noi matrix list risk_score
			
			/* turned off - not sure if we are presenting this for base level of regression
			// hazard ratio
			noi di "Hazard Ratio: ", _continue
			noi di risk_score[1, 1]
			*/
			capture drop f0
			gen f0 = (1 - s0) * 100
			capture drop f1
			gen f1 = f0 * risk_score[1, 1]
			
			// plot out the actual graph
			// absolute risk_score
			line f1 _t, sort connect(step step)
			graph export "`name'_`i'_absolute.png", as(png) replace
			// coefplot
			coefplot, xline(1) xlabel(.5 1 2 3) eform
			graph export "`name'_`i'_coef.png", as(png) replace
		}
		noi di "Function completed"
		noi di " "
	}
	
	// parametric
	if ("`fp'" == "fp") {
		noi di "Generating fully parametric plots: "

		foreach i in `vars' {
			use ${plot_set}, clear
			noi di "Variable: `i'"
			
			levelsof `i', local(i_val)
			global i_v `i_val'
			// count unique values of i to create sector variables
			local i_n = strlen("`i_val'") + 1 - strlen(subinstr("`i_val'", " ", "", .))
			
			freq_detect, var(`i')
			
			// check if sv_pos is the first level
			if (${sv_pos} == 1) {
				noi di "The most frequent value for variable `i' was detected to be the same as reference level of model"
				noi di "The second frequent value will be used as default for prediction"
				noi di "Please confirm to proceed"
				noi di "(Enter 0 or N to force most frequent value as default)", _request(confirm)
				
				if !(("${confirm}" == "0") | (strupper(substr("${confirm}", 1, 1)) == "N")) {
					replace `i' = . if `i' == ${sv_val}
					noi freq_detect, var(`i')
					use ${plot_set}, clear
				}
			}
			
			/*
			egen temp = mode(`i')
			levelsof temp, local(sv_val)
			capture drop temp
			
			levelsof `i', local(i_val)
			// count unique values of i to create sector variables
			local i_n = strlen("`i_val'") + 1 - strlen(subinstr("`i_val'", " ", "", .))
			// identify on which place sv = 1
			local sv_pos = 1
			local sv_helper = 0
			local counter = 1
			tokenize `i_val'
			while (`sv_helper' == 0) {
				if (`sv_val' != ``sv_pos'') {
					local sv_pos = `sv_pos' + 1
				}
				else if (`sv_val' == ``sv_pos'') {
					local sv_helper = 1
				}
			}
			*/
			
			streg i.`i', d(weibull)
			matrix define w = r(table)
			noi di "Weibull regression result matrix: "
			noi matrix list w
			stcurve, failure at(`i' == ${sv_val})
			graph export "`name'_`i'_fp.png", as(png) replace
			gen lambda = w[1, `i_n' + 2]
			gen p = w[1, `i_n' + 3]
			
			replace _t = round(_t)
			predict s if `i' == ${sv_val}, surv
			preserve
			collapse (max) s, by(_t)
			gen f = (1 - s) * 100
			
			noi list
			restore
		}
		noi di "Function completed"
		noi di " "
	}
}

end

capture program drop var_type
program define var_type

	qui {
		
		// determine the type of each variable being taken in
		// should take list of variables
		// default should be everything
		// end product:
		// three global macro to for binary, categorical, continuous variables
		syntax [, var(string) catt(int 15)]
		
		// sets up proper variable list for action
		if ("`var'" != "") {
			local vlist `var'
		} 
		else if ("`var'" == "") {
			ds
			local vlist `r(varlist)'
		}
		
		// determine if any of the variable specified is not in the dataset
		noi missing_detect, v(`vlist')

		if (strupper("${terminator}") == "EXIT") {
			noi di as error "WARNING: User Requested To Terminate The Program"
			exit
		}
		else {
			local vlist ${terminator}
		}

		// create returning global macros
		global bin
		global cat
		global con
		global tran
		
		// start to detect
		foreach i in `vlist' {
			
			levelsof `i', local(l1) s("!*!")
			local l2 : di subinstr("`l1'", "!*!", "", .)
			local ln1 = strlen("`l1'")
			local ln2 = strlen("`l2'")
			local ln_diff = (`ln1' - `ln2') / 3

			if (`ln_diff' == 1) {
				global bin : di "${bin}" " " "`i'"
			}
			else if (`ln_diff' > 1 & `ln_diff' < (`catt')) {
				global cat : di "${cat}" " " "`i'"
			}
			else if (inrange(`ln_diff', (`catt' - 1), (`catt' - 1 + 15))) {
				global tran : di "${tran}" " " "`i'"
				global cat : di "${cat}" " " "`i'"
			}
			else {
				global con : di "${con}" " " "`i'"
			}
			
		}
		
		// display results
		noi di "Detected binary variables: "
		noi di "${bin}"
		noi di "Detected categorical variables: "
		noi di "${cat}"
		noi di "Detected continuous variables: "
		noi di "${con}"
		
		if ("${tran}" != "") {
			noi di as error "WARNING: Variables Require User Attention"
			noi di in g "${tran}"
		}
		
		
	}

end

capture program drop missing_detect
program define missing_detect
	qui {
		syntax, v(string)
		ds
		local testing `r(varlist)'
		local missing
		foreach i in `v'{	
			if (strpos("`testing'", "`i'") == 0) {
				local missing : di "`missing'" " " "`i'"
				
			}
		}
		
		if ("`missing'" == "") {
			noi di in g "Variable check cleared"
			global terminator `v'
		}
		else if ("`missing'" != "") {				
			noi di ""
			noi di in g "Variables below not found in current dataset: "
			noi di "`missing'"
			noi di "Please re-enter variable list for the program"
			noi di "(Enter exit to terminate the program)", _request(terminator)
			
			if (strupper("${terminator}") == "EXIT") {
				exit
			}				
			if (strupper("${terminator}") != "EXIT") {
				noi missing_detect, v(${terminator})
			}
			
		}
	}
end

capture program drop freq_detect
program define freq_detect

	qui {
		
		syntax, var(varname)

		capture drop temp
		egen temp = mode(`var')
		levelsof temp, local(sv_v)
		if ("`sv_v'" == "") {
			noi di "Variable `var' had multiple most frequent levels"
			noi di "Please indicate if the max or min mode should be used"
			noi di "(Enter 1 for max mode and 2 for min mode)", _request(mode)
			
			if ("${mode}" == "1") {
				capture drop temp
				egen temp = mode(`var'), maxmode
				levelsof temp, local(sv_v)
			}
			else if ("${mode}" == "2") {
				capture drop temp
				egen temp = mode(`var'), minmode
				levelsof temp, local(sv_v)
			}
		}
		global sv_val `sv_v'
		capture drop temp

		use ${plot_set}, clear
		
		// identify on which place sv = 1
		global sv_pos = 1
		local sv_helper = 0
		local counter = 1
		tokenize ${i_v}

		while (`sv_helper' == 0) {
			if (${sv_val} != `${sv_pos}') {
				global sv_pos = ${sv_pos} + 1
			}
			else if (${sv_val} == `${sv_pos}') {
				local sv_helper = 1
			}
			
		}

	}

end
