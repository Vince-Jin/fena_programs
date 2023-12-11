*! NHANES Dataset Fast Creation
*! By Zhenghao(Vincent) Jin and Abimereki Muzaale

qui {
	
	capture program drop nhanes_fena
	program define nhanes_fena
	
		syntax , [ys(int 1999) ye(int 999999) ex(numlist) ds(int 0) format(int 0) timeout(int 1000) help mort inf(int 99999) pro s2017]
		
		qui {
			
			global display : di "noi di"
			
			// check if in program mode
			if ("`pro'" == "pro") {
				
				// modify all checkers to 1 to proceed automatically
				// substitute the macro to turn off confirmation
				global display : di "capture"
				global confirmation : di "1"
				
			}
			
			
			if 1 { // basic checks
				
				local op_list : di substr("`0'", 2, .)
				
				// check ys
				if ((`ys' < 1988) | (inrange(`ys', 1995, 1998)) | (`ys' > 2020)) {
					noi di as error "************************ERROR***************************"
					noi di as error "*******No Data Available For Requested Start Year*******"
					noi di as error "********************************************************"
					exit
				} 
				
				//check if ds is in range
				if !(inrange(`ds', 0, 4)) {
					
					noi di as error "************************ERROR***************************"
					noi di as error "*********No Data Available For Requested Dataset********"
					noi di as error "********************************************************"
					exit
					
				}
				
				// replace ye to default value
				if (`ye' == 999999) {
					local ye = `ys' + 1
				}
				
				// check if ye is greater than ys
				if !(`ye' > `ys') {
					noi di as error "************************ERROR***************************"
					noi di as error "******Requested End Year is Smaller Than Start Year*****"
					noi di as error "********************************************************"
					exit
				}
				
				// check ye
				if ((`ye' < 1988) | (inrange(`ye', 1996, 1999)) | (`ye' > 2021)) {
					noi di as error "************************ERROR***************************"
					noi di as error "*******No Data Available For Requested End Year*********"
					noi di as error "********************************************************"
					exit
				} 
				
				// modify ys for ys that starts in second year
				if (`ys' > 1998 & mod(`ys', 2) == 0) {
					local ys = `ys' - 1

					noi di as error "WARNING: Start Year Automatically Adjusted To `ys', Since Start Year Should Be An Uneven Year"

				}
				
				// modify ys to later years if requested diet or questionnaire dataset_counter
				// and ye is later than NHANES III years
				if (inrange(`ys', 1988, 1994) & (inlist(`ds', 3, 4))) {
					
					if (`ye' < 1999) {
						noi di as error "***********************ERROR*****************************"
						noi di as error "**NHANES III Data Does Not Contain Requested Information*"
						noi di as error "******************For Requested Start Year***************"
						noi di as error "*********************************************************"
						exit
					} 
					else {

						noi di as error "WARNING: NHANES III Data Does Not Contain Requested Information For Requested Start Year, Start Year Automatically Changed To 1999"

						local ys = 1999
						if (`ye' == 1999) {
							local ye = 2000
						}
					} 
					
				}
				
				// check if output format is in list
				/*
				0 - dta
				1 - excel
				2 - csv
				3 - sasxport5
				4 - sasxport8
				5 - text
				*/
				if !(inlist(`format', 0, 1, 2, 3, 4, 5)) {
					noi di as error "************************ERROR***************************"
					noi di as error "*****************Unknown Output Format******************"
					noi di as error "********************************************************"
					exit
				}
				
				// check income inflation year
				if (`inf' != 99999 & !inrange(`inf', 1913, 2023)) {
					noi di as error "************************ERROR***************************"
					noi di as error "***********No CPI Information For Year `inf'************"
					noi di as error "********************************************************"
					exit
				}
				
				if (`inf' != 99999 & `ds' != 1) {
					noi di as error "WARNING: Inflation Adjustment Will NOT Be Executed As Requested Dataset Does Not Contain Income Information"
				}
				
				// check if exclusion year is in NHANES III
				foreach i in `ex' {
					if ((inrange(`i', 1988, 1994)) & (`i' > `ys')) {

						noi di as error "WARNING: Exclusion For Year `i' Will NOT Be Executed Due To The Request Of NHANES III DATA"

					}
				}
				
				// warning for NHANES III then modify ys to 1988
				if (inrange(`ys', 1988, 1994)) {

					noi di as error "WARNING: Data For Year Between 1988 To 1994 Requested, NHANES III Data For That Period Will Be Loaded"

					local ys = 1988
				}
				
				// warning for time gap
				if (`ys' < 1995 & `ye' > 1999) {

					noi di as error "WARNING: No Data Between 1995 And 1999 Available"

				}
				
				// warning for pre-pandemic
				if ((`ys' < 2017 & `ye' > 2017) | inrange(`ys', 2017, 2020)) {
					
					if ("`s2017'" == "s2017") {
						noi di "WARNING: User Requested Data for Survey Cycle 2017-2018, Pre-Pandemic Data Will Not Be Loaded. Request End Year Will Be Modified To 2018"
						local ye = 2018
					}
					noi di as error "WARNING: Data For Year Between 2017 To 2020 Requested, NHANES Pre-Pandemic Data For That Period Will Be Load"

				}
				
				// warning for 2020
				if (inrange(2020, `ys', `ye')) {

					noi di as error "WARNING: Only Data Until March Was Included For 2020"

				}
				
				// warning for mortality
				if ((inrange(2019, `ys', `ye') | (inrange(2020, `ys', `ye')))) {
					
					noi di as error "WARNING: No mortality information available for 2019 to 2020."
					
				}
				
				// warning for harmonization of NHANES III and later datasets
				local har = 0
				if ((inrange(`ys', 1988, 1994)) & (`ye' > 1999)) {
					
					noi di as error "WARNING: " ustrtitle("Automatic harmonization will be performed to aggregate ") "NHANES III " ustrtitle("variables and later ") "NHANES " ustrtitle("dataset variables")
					local har = 1
					
				}
				
				// current ds:
				/* 
				base          1
				demographic   2
				exam          3
				diet          4
				questionnaire 5
				
				*/
				local ds_helper base demographic exam diet questionnaire
				tokenize `ds_helper'
				local ds = `ds' + 1
				macro drop ds_helper2
				global ds_helper2 : di "``ds''"
				if ("`mort'" == "mort") {
					global ds_helper3 : di "${ds_helper2}" " with mortality"
				} 
				else {
					global ds_helper3 : di "${ds_helper2}"
				}
				
				// check if help option is called
				if ("`help'" == "help") {
					
					noi nhanes1_helper
					
					noi di "Continue the inquiry? (Y/N)", _request(help_ter)
					if !((strupper(substr("${help_ter}", 1, 1)) == "Y") | (substr("${help_ter}", 1, 1) == "1")) {
						noi di as error "***********************WARNING**************************"
						noi di as error "***********User Requested To Abort Creation*************"
						noi di as error "********************************************************"
						exit
					}
					
				}
				
			}
			
			if 2 { // print out summary for dataset creation and ask for confirmation
				
				// also prepare a list for exclusion to be tested
				// screen if ex_helper contains even years
				macro drop ex_helper
				macro drop ex_helper2
				foreach i in `ex' {
					local disp_helper : di "`disp_helper'" " " `i'
					if (`i' > 1998 & mod(`i', 2) == 0) {
						local ex_helper = `i' - 1
					}
					else {
						local ex_helper = `i'
					}
					global ex_helper2 : di "`ex_helper'" ", " "${ex_helper2}"
				}
				global ex_helper2 : di substr("${ex_helper2}", 1, strlen("${ex_helper2}") - 2)

				noi di as text " "
				noi di "NHANES Dataset Creation Summary"
				noi di "Requested Start Year: `ys'"
				noi di "Requested End Year: `ye'"
				if ("`disp_helper'" != "") {
					noi di "Requested Year To Be Excluded: `disp_helper'"
				}
				else {
					noi di "Requested Year To Be Excluded: None"
				}
				noi di "Requested Dataset Information: ${ds_helper3}"
				noi di "Current Saving Directory: `c(pwd)'"
				noi di "Timeout Setting: `timeout'"
				${display} "Please Confirm To Proceed (Y/N)", _request(confirmation)
				if !((strupper(substr("${confirmation}", 1, 1)) == "Y") | (substr("${confirmation}", 1, 1) == "1")) {
					noi di as error "***********************WARNING**************************"
					noi di as error "***********User Requested To Abort Creation*************"
					noi di as error "********************************************************"
					exit
				}	
				set timeout1 `timeout'
			}
			
			if 3 { // pulling dataset
				
				local ye_helper = `ye' - 1
				forvalues i = `ys'(1)`ye_helper' {
						if ("`ex'" != "") {
							if !(inlist(`i', ${ex_helper2})) {
								local year_helper `i' `year_helper'
							}
						}
						else {
							local year_helper `i' `year_helper'
						}
				}

				foreach i in `year_helper' {
					if (`i' > 1998 & mod(`i', 2) != 0) {
						local ys_2 `i' `ys_2'
					}
					else if (`i' < 1995 & `i' > 1987) {
						local ys_1 `i' `ys_1'
					}
				}

				local dataset_counter = 0
				
				// loading the NHANES III Data
				if ("`ys_1'" != "") {
					
					local dataset_counter = `dataset_counter' + 1
					
					noi di "Loading NHANES III Data For 1988-1994...................", _continue
					
					macro drop nh3
					global nh3 https://wwwn.cdc.gov/nchs/data/nhanes3/1a/
					
					// base dataset
					*** PERHAPS CONDENSE DS CONDITIONS TO ONE MACRO ****
					*** WORK ON THIS LATER ***
					if (`ds' == 1) {
						noi infix seqn 1-5 ///
									dmarethn 12 ///
									dmaracer 13 ///
									dmaethnr 14 ///
									hssex 15 ///
									hsageir 18-19 ///
									using ${nh3}adult.DAT, clear
						tempfile demo1
						save `demo1'
						
					}
					
					// demographic dataset
					if (`ds' == 2) {
						
						noi infix seqn 1-5 ///
									dmpstat 11 /// 
									dmarethn 12 ///
									dmaethnr 14 ///
									hssex 15 ///
									hsageir	18-19 ///
									hfa8r 1256-1257 ///
									hff19r 1409-1410 ///
									sdpphase 42 ///
									using ${nh3}adult.DAT, clear
													
					}
					
					// exam dataset
					if (`ds' == 3) {
						
						infix seqn 1-5 pepmnk1r 1423-1425 pepmnk5r 1428-1430 bmpbmi 1524-1527 /// 
						 using ${nh3}exam.dat, clear
						tempfile exam
						save `exam'

						infix seqn 1-5 ghp 1861-1864 ubp 1965-1970 urp 1956-1960 sgpsi 1761-1765 /// 
						 cep 1784-1787 using ${nh3}lab.dat, clear
						tempfile lab
						save `lab'
						 
						infix seqn 1-5 hae4d2 1605 had1 1561 had6 1568 had10 1578 hak1 1855 /// 
						 hae2 1598 hae5a 1610 hat1met 2393-2395 har1 2281 ///
						 using ${nh3}adult.DAT, clear
							 
						tempfile qa 
						save `qa'  
						 
						use `exam', clear 
						merge 1:1 seqn using `lab', nogen
						merge 1:1 seqn using `qa', nogen 
						
						noi di "(`c(N)' observations read)"
						
					}					
					
					// save to a tempfile for future appending
					tempfile ds`dataset_counter'
					capture drop year_start
					gen year_start = 1988
					save `ds`dataset_counter'', replace

				}
				
				// loading NHANES Data in other years
				// having prepan_counter to ensure 2017-2020 only loaded once
				if ("`ys_2'" != "") {
					
					local prepan_counter = 0
					
					macro drop nhanes
					global nhanes "https://wwwn.cdc.gov/Nchs/Nhanes"
					
					foreach i of numlist `ys_2' {
						
						local ye_helper2 = `i' + 1
						
						local name_helper = `i' - 1999
						if (`name_helper' > 0) {
							local name_helper2 = 2 + ((`name_helper' - 2) / 2)
							tokenize "`c(ALPHA)'"
							local name_helper3 : di "_" "``name_helper2''"
							if ("`s2017'" != "s2017") {
								if (("`name_helper3'" == "_J") | ("`name_helper3'" == "_K")) {
									local name_helper3 : di "P_"
								}
							}
						}

						if (!(inrange(`i', 2017, 2020)) | (inrange(`i', 2017, 2020) & `prepan_counter' == 0)) {
								
							local dataset_counter = `dataset_counter' + 1
							if (inrange(`i', 2017, 2020)) {
								local prepan_counter = `prepan_counter' + 1
							}
							
							if (`ds' == 1) {
								
								local ds_name : di "DEMO`name_helper3'"
								local period : di "`i'" "-" "`ye_helper2'"
								if ("`s2017'" != "s2017") {
									if (inrange(`i', 2017, 2020)) {
										local ds_name : di "`name_helper3'DEMO"
										local period : di "2017-2020"
									}
								}
								noi di "Loading NHANES Data For " "`period'...................", _continue
								import sasxport5 "$nhanes/`i'-`ye_helper2'/`ds_name'.XPT", clear
								tempfile ds`dataset_counter'
								keep seqn ///
										riagendr ///
										ridageyr ///
										ridreth1
								capture drop year_start
								gen year_start = `i'
								save `ds`dataset_counter'', replace
								noi di "(" "`c(N)'" " " "observations read)"
								
							}
							
							if (`ds' == 2) {
								
								local ds_name : di "DEMO`name_helper3'"
								local period : di "`i'" "-" "`ye_helper2'"
								if ("`s2017'" != "s2017") {
									if (inrange(`i', 2017, 2020)) {
										local ds_name : di "`name_helper3'DEMO"
										local period : di "2017-2020"
									}
								}
								noi di "Loading NHANES Data For " "`period'...................", _continue
								import sasxport5 "$nhanes/`i'-`ye_helper2'/`ds_name'.XPT", clear
								tempfile ds`dataset_counter'
								capture drop year_start
								gen year_start = `i'
								drop wtm* wti*
								save `ds`dataset_counter'', replace
								noi di "(" "`c(N)'" " " "observations read)"
								
							}
							
							if (`ds' == 3) {
								
								// 1999
								local DATA1 BPX BMX LAB10 LAB16  LAB18	
								// 2001 2003
								local DATA2 BPX BMX L10   L16    L40  
								// 2005+
								local DATA3 BPX BMX GHB   ALB_CR BIOPRO 
								// 2017-2020.3
								local DATA4 BPXO BMX GHB   ALB_CR BIOPRO 
								
								if (`i' == 1999) {
									local ex_name `DATA1'
								} 
								else if (inlist(`i', 2001, 2003)) {
									local ex_name `DATA2'
								}
								else if ((inrange(`i', 2004,9999)) & !(inrange(`i', 2017, 2020))) { // leave as 9999 for future datasets
									local ex_name `DATA3'
								}
								else if (inrange(`i', 2017, 2020)) {
									local ex_name `DATA4'
								}
								
								local period : di "`i'" "-" "`ye_helper2'"
								if ("`s2017'" != "s2017") {
									if (inrange(`i', 2017, 2020)) {
										local period : di "2017-2020"
									}
								}
								noi di "Loading NHANES Data For " "`period'...................", _continue
								
								// combine to get dataset name
								local j_helper = 0
								foreach j in `ex_name' {
									
									local j_helper = `j_helper' + 1

									local ds_name : di "`j'" "`name_helper3'"

									if ("`s2017'" != "s2017") {
										if (inrange(`i', 2017, 2020)) {
											local ds_name : di "`name_helper3'" "`j'"
										}
									}
									
									import sasxport5 "$nhanes/`i'-`ye_helper2'/`ds_name'.XPT", clear
									tempfile ex_`j_helper'
									save `ex_`j_helper'', replace
									
								}
								
								// merge exam and lab data together
								use `ex_1', clear

								forvalues k = 2/`j_helper' {
									
									merge 1:1 seqn using `ex_`k'', gen(merge`k')
									
								}
								
								// save each survey year dataset
								tempfile ds`dataset_counter'
								capture drop year_start
								gen year_start = `i'
								save `ds`dataset_counter'', replace
								noi di "(" "`c(N)'" " " "observations read)"
								
							}

							if (`ds' == 4) {
								
								// 1999 & 2001 DRXTOT
								// 2003+ DR1TOT DR2TOT
								

								local period : di "`i'" "-" "`ye_helper2'"

								if ("`s2017'" != "s2017") {
									if (inrange(`i', 2017, 2020)) {
										local period : di "2017-2020"
									}
								}
								noi di "Loading NHANES Data For " "`period'...................", _continue
								
								if ((`i' == 1999) | (`i' == 2001)) {
									local ds_name : di "DRXTOT`name_helper3'"
									import sasxport5 "$nhanes/`i'-`ye_helper2'/`ds_name'.XPT", clear
								} 
								else if (inrange(`i', 2017, 2020)) {
									forvalues j = 1/2 {
										if ("`s2017'" != "s2017") {										
											local ds_name : di "`name_helper3'" "DR" "`j'" "TOT"
										}
										else if ("`s2017'" == "s2017") {
											local ds_name : di "DR" "`j'" "TOT" "`name_helper3'"
										}
										import sasxport5 "$nhanes/`i'-`ye_helper2'/`ds_name'", clear
										
										tempfile diet_`j'
										save `diet_`j'', replace
									}
									use `diet_1', clear
									merge 1:1 seqn using `diet_2'
								}
								else if (`i' > 2002) {
									forvalues j = 1/2 {
										local ds_name : di "DR" "`j'" "TOT" "`name_helper3'"
										import sasxport5 "$nhanes/`i'-`ye_helper2'/`ds_name'", clear
										
										tempfile diet_`j'
										save `diet_`j'', replace
									}
									use `diet_1', clear
									merge 1:1 seqn using `diet_2'
								}
								
								tempfile ds`dataset_counter'
								capture drop year_start
								gen year_start = `i'
								save `ds`dataset_counter'', replace
								noi di "(" "`c(N)'" " " "observations read)"
							}
							
							if (`ds' == 5) {
								
								// 1999
								local DATA1 ALQ DIQ KIQ   BPQ PAQ SMQ MCQ BPQ
								// 2001+
								local DATA2 ALQ DIQ KIQ_U BPQ PAQ SMQ MCQ BPQ
								
								if (`i' == 1999) {
									local qs_name `DATA1'
								}
								else if ((inrange(`i', 2001,9999))) { // leave as 9999 for future datasets
									local qs_name `DATA2'
								}
								
								local period : di "`i'" "-" "`ye_helper2'"
								if ("`s2017'" != "s2017") {
									if (inrange(`i', 2017, 2020)) {
										local period : di "2017-2020"
									}
								}
							
								noi di "Loading NHANES Data For " "`period'...................", _continue
								
								// combine to get dataset name
								local j_helper = 0
								foreach j in `qs_name' {
									
									local j_helper = `j_helper' + 1
									
									local ds_name : di "`j'" "`name_helper3'"
									
									if ("`s2017'" != "s2017") {
										if (inrange(`i', 2017, 2020)) {
											local ds_name : di "`name_helper3'" "`j'"
										}
									}
									
									import sasxport5 "$nhanes/`i'-`ye_helper2'/`ds_name'.XPT", clear
									tempfile qs_`j_helper'
									save `qs_`j_helper'', replace
									
								}
								
								// merge exam and lab data together
								use `qs_1', clear

								forvalues k = 2/`j_helper' {
									
									merge 1:1 seqn using `qs_`k'', gen(merge`k')
									
								}
								
								// save each survey year dataset
								tempfile ds`dataset_counter'
								capture drop year_start
								gen year_start = `i'
								save `ds`dataset_counter'', replace
								noi di "(" "`c(N)'" " " "observations read)"
								
							}
						
						}
						
					}
					
				}
				
				local ds_max = `dataset_counter' - 1
				
				noi di "Appending All Datasets..............................", _continue
				forvalues i = 1/`ds_max' {
					
					append using `ds`i''
					
				}
				
				compress
				noi di "Done"
				tempfile master_dataset
				save `master_dataset', replace
		
				if ("`mort'" == "mort") {
					
					local mort_counter = 0
					
					global nchs https://ftp.cdc.gov/pub/Health_Statistics/NCHS/
					global linkage datalinkage/linked_mortality/
					capture macro drop nh3
					global nh3 NHANES_III
					global mort_end _MORT_2019_PUBLIC.dat
					
					global mort_var ///
							 seqn 1-6 ///
							 eligstat 15 ///
							 mortstat 16 ///
							 ucod_leading 17-19 ///
							 diabetes 20 ///
							 hyperten 21 ///
							 permth_int 43-45 ///
							 permth_exm 46-48
					
					noi di "Loading Mortality Information.......................", _continue
					
					if ("`ys_1'" != "") {
						
						local mort_counter = `mort_counter' + 1
						
						infix ${mort_var} using ${nchs}${linkage}${nh3}${mort_end}, clear
						// replace seqn=-1*seqn
						tempfile mort`mort_counter'
						save `mort`mort_counter''
						
					}
					
					if ("`ys_2'" != "") {
						
						foreach i of numlist `ys_2' {
							if !(inrange(`i', 2019, 2020)) {
								local mort_counter = `mort_counter' + 1
								local ye_helper2 = `i' + 1
								local nhanes NHANES_`i'_`ye_helper2'
								infix ${mort_var} using ${nchs}${linkage}`nhanes'${mort_end}, clear
								tempfile mort`mort_counter'
								save `mort`mort_counter''
							}
						}
						
					}
					
					local mort_max = `mort_counter' - 1
					forvalues i = 1/`mort_max' {
						append using `mort`i''
					}
					
					capture lab def premiss .z "Missing"
					capture lab def eligfmt 1 "Eligible" 2 "Under age 18, not available for public release" 3 "Ineligible" 
					capture lab def mortfmt 0 "Assumed alive" 1 "Assumed deceased" .z "Ineligible or under age 18"
					capture lab def flagfmt 0 "No - Condition not listed as a multiple cause of death" ///
									  1 "Yes - Condition listed as a multiple cause of death"	///
									  .z "Assumed alive, under age 18, ineligible for mortality follow-up, or MCOD not available"
					capture lab def qtrfmt 1 "January-March" ///
									 2 "April-June" ///
									 3 "July-September" ///
									 4 "October-December" ///
									 .z "Ineligible, under age 18, or assumed alive"
					capture lab def dodyfmt .z "Ineligible, under age 18, or assumed alive"
					capture lab def ucodfmt 1 "Diseases of heart (I00-I09, I11, I13, I20-I51)"                           
					capture lab def ucodfmt 2 "Malignant neoplasms (C00-C97)"                                            , add
					capture lab def ucodfmt 3 "Chronic lower respiratory diseases (J40-J47)"                             , add
					capture lab def ucodfmt 4 "Accidents (unintentional injuries) (V01-X59, Y85-Y86)"                    , add
					capture lab def ucodfmt 5 "Cerebrovascular diseases (I60-I69)"                                       , add
					capture lab def ucodfmt 6 "Alzheimer's disease (G30)"                                                , add
					capture lab def ucodfmt 7 "Diabetes mellitus (E10-E14)"                                              , add
					capture lab def ucodfmt 8 "Influenza and pneumonia (J09-J18)"                                        , add
					capture lab def ucodfmt 9 "Nephritis, nephrotic syndrome and nephrosis (N00-N07, N17-N19, N25-N27)"  , add
					capture lab def ucodfmt 10 "All other causes (residual)"                                             , add
					capture lab def ucodfmt .z "Ineligible, under age 18, assumed alive, or no cause of death data"      , add
					
					
					replace mortstat = .z if mortstat >=.
					replace ucod_leading = .z if ucod_leading >=.
					replace diabetes = .z if diabetes >=.
					replace hyperten = .z if hyperten >=.
					replace permth_int = .z if permth_int >=.
					replace permth_exm = .z if permth_exm >=.
					
					label var seqn "NHANES Respondent Sequence Number"
					label var eligstat "Eligibility Status for Mortality Follow-up"
					label var mortstat "Final Mortality Status"
					label var ucod_leading "Underlying Cause of Death: Recode"
					label var diabetes "Diabetes flag from Multiple Cause of Death"
					label var hyperten "Hypertension flag from Multiple Cause of Death"
					label var permth_int "Person-Months of Follow-up from NHANES Interview date"
					label var permth_exm "Person-Months of Follow-up from NHANES Mobile Examination Center (MEC) Date"

					label values eligstat eligfmt
					label values mortstat mortfmt
					label values ucod_leading ucodfmt
					label values diabetes flagfmt
					label values hyperten flagfmt
					label value permth_int premiss
					label value permth_exm premiss
					
					//eligibility
					drop if inlist(eligstat,2,3)
					duplicates drop 

					keep seqn mortstat permth_int permth_exm
					tempfile mort_dataset
					save `mort_dataset', replace
					
					use `master_dataset', clear
					merge 1:1 seqn using `mort_dataset', nogen
					
					//variable labels derived from wwwn.cdc.gov/nchs/nhanes/
					capture lab var seqn "Sequence number (seqn x -1 for NHANES III)"
					capture lab var dmpstat "Examinmation/interview status"
					capture lab var dmarethn "Race/ethnicity"
					capture lab var hssex "Sex"
					capture lab var hsageir "Age at interview (screener) -qty"
					capture lab var hfa8r "Education level"
					capture lab var hff19r "Total family 12 month income group"
					capture lab var pepmnk1r "Overall average K1, systolic, BP(age 5+)"
					capture lab var pepmnk5r "Overall average K5, diastolic, BP(age5+)"
					capture lab var bmpbmi "Body mass index"
					capture lab var ghp "Glycated hemoglobin: %"
					capture lab var ubp "Urinary albumin (ug/mL)"
					capture lab var urp "Urinary creatinine (mg/dL)"
					capture lab var sgpsi "serum glucose: SI (mmol/L)"
					capture lab var cep "Serum creatine (mg/dL)"
					capture lab var had1 "Ever been told you have sugar/diabetes"
					capture lab var hae2 "Doctor ever told had hypertension/HBP"
					capture lab var har1 "Have tiy snijed 100_ cigarettes in life"
					
					noi di "Done"
					
				}
				
				// harmonization process for NHANES III and other NHANES datasets
				if (`har' == 1) {

					// base - only race needs recode
					if (`ds' == 1) {
						replace riagendr = hssex if year_start == 1988
						replace ridageyr = hsageir if year_start == 1988
						replace ridreth1 = 1 if dmarethn == 3 & year_start == 1988
						replace ridreth1 = 3 if dmarethn == 1 & year_start == 1988
						replace ridreth1 = 4 if dmarethn == 2 & year_start == 1988
						replace ridreth1 = 5 if dmarethn == 4 & year_start == 1988
						rename (riagendr ridageyr ridreth1) (gender age ethn)
						
						keep seqn age gender ethn dmaracer dmaethnr
						
						label variable dmaracer "Race information for NHANES III"
						label variable dmaethnr "Ethnicity information for NHANES III"
					}
					
					// demographic 
					/* 
					seqn
					dmpstat Examination/interview status:
					1 - interviewed, not examed
					2 - interviewed, MEC-examed
					3 - interviewed, home-examed
					dmarethn - race-ethnicity
					dmaethnr - race
					hssex Sex
					hsageir age at interview
					hfa8r Highest grade
					00 never attend
					01 - 17 grades
					01 - 05 elementary
					06 - 08 middle high
					09 - 12 high school
					13 - 16 bachelor
					17 - college graduate or above
					88 blank but applicable
					99 - unknown
					hff19r income groups:
					0 - no income
					1 - less than 1000
					2 - 1999
					3 - 2999
					4 - 3999
					...
					21 - 24999
					22 - 29999
					23 - 34999
					24 - 39999
					25 - 44999
					26 - 49999
					27 - 50000 and over
					88 blank but applicable
					99 - unknown
					sdpphase phase of nhanes iii survey
					1 - phase 1
					2 - phase 2
					*/
					if (`ds' == 2) {
						
						// dmpstat
						replace ridstatr = 1 if dmpstat == 1 & year_start == 1988
						replace ridstatr = 2 if ((dmpstat == 2 | dmpstat == 3) & year_start == 1988)
						
						// hssex
						replace riagendr = hssex if year_start == 1988
						
						// hsageir
						replace ridageyr = hsageir if year_start == 1988
						
						// dmarethn
						replace ridreth1 = 1 if dmarethn == 3 & year_start == 1988
						replace ridreth1 = 3 if dmarethn == 1 & year_start == 1988
						replace ridreth1 = 4 if dmarethn == 2 & year_start == 1988
						replace ridreth1 = 5 if dmarethn == 4 & year_start == 1988
						
						// hfa8r
						/* code to dmdeduc3, dmdeduc2, dmdeduc
						dmdeduc3 - education level - children/youth 6-19:
						0 - never attened / kindergarten only
						1 - 1th grade
						...
						12 - 12th grade, no diploma
						13 - high school graduate
						14 - GED or equivalent
						15 - more than high school
						55 - less than 5th grade
						66 - less than 9th grade
						77 - refused
						99 - don't know
						
						dmdeduc2 - education - adults 20+:
						1 - less than 9th grade
						2 - 9-11th grade
						3 - high school grad / ged or equivalent
						4 - some college or aa degree
						5 - college graduate or above
						7 - refused
						9 - unknown
						
						dmdeduc - education recode:
						1 - less than high school
						2 - high school diploma
						3 - more than high school
						7 - refused
						9 - unkown
						*/
						
						replace dmdeduc3 = 0 if ((hfa8r == 0 & inrange(ridageyr, 6, 19) & year_start == 1988))
						replace dmdeduc3 = hfa8r if (inrange(hfa8r, 1, 11) & inrange(ridageyr, 6, 19) & year_start == 1988)
						replace dmdeduc3 = 13 if (hfa8r == 12 & inrange(ridageyr, 6, 19) & year_start == 1988)
						replace dmdeduc3 = 15 if (hfa8r > 12 & inrange(ridageyr, 6, 19) & year_start == 1988)
						replace dmdeduc3 = 77 if (hfa8r == 88 & inrange(ridageyr, 6, 19) & year_start == 1988)
						replace dmdeduc3 = 99 if (hfa8r == 99 & inrange(ridageyr, 6, 19) & year_start == 1988)
						
						replace dmdeduc2 = 1 if (hfa8r < 9 & ridageyr > 19 & year_start == 1988)
						replace dmdeduc2 = 2 if (inrange(hfa8r, 9 , 11) & ridageyr > 19 & year_start == 1988)
						replace dmdeduc2 = 3 if (hfa8r == 12 & ridageyr > 19 & year_start == 1988)
						replace dmdeduc2 = 4 if (inrange(hfa8r, 13, 15) & ridageyr > 19 & year_start == 1988)
						replace dmdeduc2 = 5 if (hfa8r > 15 & ridageyr > 19 & year_start == 1988)
						replace dmdeduc2 = 7 if (hfa8r == 88 & ridageyr > 19 & year_start == 1988)
						replace dmdeduc2 = 9 if (hfa8r == 99 & ridageyr > 19 & year_start == 1988)
						
						replace dmdeduc = 1 if (hfa8r < 12 & year_start == 1988)
						replace dmdeduc = 2 if (hfa8r == 12 & year_start == 1988)
						replace dmdeduc = 3 if (hfa8r > 12 & year_start == 1988)
						replace dmdeduc = 7 if (hfa8r == 88 & year_start == 1988)
						replace dmdeduc = 9 if (hfa8r == 99 & year_start == 1988)
						
						// hff19r - recode to indfminc
						/* temporary turned off due to vague definition of 50000+
						replace indfminc = 1 if (hff19r == 0 & year_start == 1988)
						forvalues i = 1/4 {
							/*
							999 - 4999 1-5
							5999 - 9999 6-10
							10999 - 14999 11-15
							15999 - 19999 16-20
							*/
							local min = 1 + 5 * (`i' - 1)
							local max = 0 + 5 * `i'
							replace indfminc = `i' if inrange(hff19r, `min', `max') 
						}
						
						replace indfminc = 5 if hff19r == 21
						replace indfminc = 6 if inrange(hff19r, 22, 23)
						replace indfminc = 7 if inrange(hff19r, 24, 25)
						*/
						
					}
					
				}
				
				// adjusting inflation for income
				if (`inf' != 99999 & `ds' == 2) {
				
					noi di "Adjusting income for inflation......................", _continue
					global inf_yr = `inf'
					tempfile income
					save `income', replace
		
					// use "https://github.com/jhurepos/projectbeta/blob/main/nh_projectbeta_5/yearly_cpi", clear
					use yearly_cpi, clear
					capture drop cpi_helper
					gen cpi_helper = cpi if year_start == ${inf_yr}
					capture drop cpi_t
					egen cpi_t = min(cpi_helper)
					drop cpi_helper
					capture drop ifr
					gen ifr = cpi_t / cpi
					keep  year_start ifr
					tempfile cpi
					save `cpi', replace
					
					use `income', replace
					joinby year_start using `cpi', unmatched(none)

					
					// switch categories to numeric incomes
					capture drop income_n
					gen income_n = .
					/*
					indfminc:
					1 - 0-4999
					2 - 5000-9999
					3 - 10000-14999
					4 - 15000-19999
					5 - 20000-24999
					6 - 25000-34999
					7 - 35000-44999
					8 - 45000-54999
					9 - 55000-64999
					10 - 65000-74999
					11 - 75000
					12 - 20000
					13 - 20000
					treat > 20000 as 20001
					use 10000 for < 20000 assuming normal distribution
					*/
					local inc_vars indhhinc indhhin2
					foreach r in `inc_vars' {
						capture forvalues i = 1/6 {
							local inc = 5000 * (`i' - 1)
							capture replace income_n = `inc' if  `r' == `i'
						}
						capture forvalues i = 7/11 {
							local inc = 35000 + (`i' - 7) * 10000
							capture replace income_n = `inc' if  `r' == `i'
						}
						capture replace income_n = 20000 if  `r' == 12
						capture replace income_n = 10000 if  `r' == 13
						capture replace income_n = 7777777 if  `r' == 77
						capture replace income_n = 9999999 if  `r' == 99
						capture replace income_n = 75000 if `r' == 14
						capture replace income_n = 100000 if `r' == 15
					}
					/*
					hff19r income groups:
					0 - no income
					1 - less than 1000
					2 - 1999
					3 - 2999
					4 - 3999
					...
					21 - 24999
					22 - 29999
					23 - 34999
					24 - 39999
					25 - 44999
					26 - 49999
					27 - 50000 and over
					88 blank but applicable
					99 - unknown
					using 500 for 0 - 1000 assuming a normal distribution
					*/
					capture replace income_n = 0 if hff19r == 0
					capture replace income_n = 500 if hff19r == 1
					capture forvalues i = 2/21 {
						capture local inc = 1000 * (`i' - 1)
						capture replace income_n = `inc' if hff19r == `i'
					}
					capture forvalues i = 22/27 {
						capture local inc = 25000 + 5000 * (`i' - 22)
						capture replace income_n = `inc' if hff19r == `i'
					}
					capture replace income_n = 7777777 if hff19r == 88
					capture replace income_n = 9999999 if hff19r == 99
					
					capture drop inf_n
					gen inf_n = ifr * income_n
					
					
					// recode back - using new codes
					/*
					indfminc:
					1 - 0-4999
					2 - 5000-9999
					3 - 10000-14999
					4 - 15000-19999
					5 - 20000-24999
					6 - 25000-34999
					7 - 35000-44999
					8 - 45000-54999
					9 - 55000-64999
					10 - 65000-74999
					11 - 75000
					12 - 20000
					13 - 20000
					*/
					capture drop income_adjusted
					gen income_adjusted = .

					forvalues i = 1/5 {
						local min = 5000 * (`i' - 1)
						local max = 5000 * `i' - 1
						replace income_adjusted = `i' if inrange(inf_n, `min', `max')
					}
					forvalues i = 6/10 {
						local min = 25000 + 10000 * (`i' - 6)
						local max = 35000 + 10000 * (`i' - 6) - 1
						replace income_adjusted = `i' if inrange(inf_n, `min', `max')
					}
					
					replace income_adjusted = 14 if (inrange(inf_n, 75000, 99999))
					replace income_adjusted = 15 if (inf_n > 100000 & inf_n != .)
					replace income_adjusted = 12 if (income_n == 20001 & inf_n > 19999) | (income_n == 10000 & inf_n > 19999)
					replace income_adjusted = 13 if (income_n == 20001 & inf_n < 20001) | (income_n == 10000 & inf_n < 20001)
					replace income_adjusted = 77 if (income_n == 7777777)
					replace income_adjusted = 99 if income_n == 9999999
					
					// capture drop inf_n 
					// capture drop income_n
					// capture drop ifr
					
					noi di "Done"
					
				}
				
				// labeling Variables
				if inrange(`ds', 1, 2) {
					
					noi di "Labeling Demographic Variables......................", _continue
					
					if (`ds' == 1) {
						
						// gender
						capture label drop gender
						label define gender 1 "Male" 2 "Female"
						capture label values gender gender
						
						// ethn
						capture label drop ethn
						label define ethn 1 "Mexican American" 2 "Other Hispanic" ///
							3 "Non-Hispanic White" 4 "Non-Hispanic Black" 5 "Other Race - Including Multiracial"
						capture label values ethn ethn
						
					}
					
					if (`ds' == 2) {
						
						// gender
						capture label drop gender
						label define gender 1 "Male" 2 "Female"
						capture label values gender gender
						
						// ethn
						capture label drop ethn
						label define ethn 1 "Mexican American" 2 "Other Hispanic" ///
							3 "Non-Hispanic White" 4 "Non-Hispanic Black" 5 "Other Race - Including Multiracial"
						capture label values ethn ethn
						
						// veteran
						capture label drop veteran
						label define veteran 1 "Yes" 2 "No" 7 "Refused" 9 "Don't Know"
						capture label values dmqilit veteran
						
						// country of birth
						capture label drop birth
						label define birth 1 "50 US States/Washington D.C." ///
							2 "Mexico" 3 "Elsewhere" 7 "Refused" 9 "Don't Know"
						capture label values dmdborn birth
						
						// citizenship
						capture label drop citi
						label define citi 1 "Citizen By Birht or Naturalization" ///
							2 "Not a Citizen of the US" 7 "Refused" 9 "Don't Know"
						capture label values dmdcitzn citi
						
						// length of time in US
						capture label drop stay
						label define stay ///
							1 "<1 year" ///
							2 "1-5 years" ///
							3 "5-10 years" ///
							4 "10-15 years" ///
							5 "15-20 years" ///
							6 "20-30 years" ///
							7 "30-40 years" ///
							8 "40-50 years" ///
							9 ">50 years" ///
							77 "Refused" ///
							88 "Could Not Determine" ///
							99 "Don't Know"
						capture label values dmdyrsus stay
						
						// education
						capture label drop educ
						label define educ ///
							1 "Less Than High School" ///
							2 "High School Diploma (Including GED)" ///
							3 "More Than High School" ///
							7 "Refused" ///
							9 "Don't Know"
						capture label values dmdeduc educ
						
						// marital status
						capture label drop marital
						label define marital ///
							1 "Married" ///
							2 "Widowed" ///
							3 "Divorced" ///
							4 "Separated" ///
							5 "Never Married" ///
							6 "Living with Partner" ///
							77 "Refused" ///
							99 "Don't Know"
							
						// income
						capture label drop income
						label define income ///
							1 "0-4999" ///
							2 "5000-9999" ///
							3 "10000-14999" ///
							4 "15000-19999" ///
							5 "20000-24999" ///
							6 "25000-34999" ///
							7 "35000-44999" ///
							8 "45000-54999" ///
							9 "55000-64999" ///
							10 "65000-74999" ///
							11 ">75000" ///
							12 ">20000" ///
							13 "<20000" ///
							77 "Refused" ///
							99 "Don't Know"
						capture label values income_adjusted income
					}
					
					noi di "Done"
				}
				
				noi di "Saving Dataset......................................", _continue
				/*
				0 - dta
				1 - excel
				2 - csv
				3 - sasxport5
				4 - sasxport8
				5 - text
				*/
				if (`format' == 0) {
					save NHANES_${ds_helper2}_`ys'_`ye', replace
				}
				else if (`format' == 1) {
					export excel using NHANES_${ds_helper2}_`ys'_`ye', replace
				}
				else if (`format' == 2) {
					export delimited using NHANES_${ds_helper2}_`ys'_`ye'.csv, replace
				}
				else if (`format' == 3) {
					export sasxport5 NHANES_${ds_helper2}_`ys'_`ye', replace
				}
				else if (`format' == 4) {
					export sasxport8 NHANES_${ds_helper2}_`ys'_`ye', replace
				}
				else if (`format' == 5) {
					outfile using NHANES_${ds_helper2}_`ys'_`ye', replace
					outfile using NHANES_${ds_helper2}_`ys'_`ye', replace dictionary
				}
				noi di "Dataset Saved!"
				
			}
			
		}
	end
	
	capture program drop nhanes1_helper
	program define nhanes1_helper
	
		qui {
			
			cls
			noi di in g "You have entered the help option of nhanes1 program"
			noi di "Please specify the help information you would like to access: "
			noi di "0 - For general information"
			noi di "1 - For ys option"
			noi di "2 - For ye option"
			noi di "3 - For ex option"
			noi di "4 - For ds option"
			noi di "5 - For format option"
			noi di "6 - For timeout option"
			noi di "7 - For mort option"
			noi di "To request multiple information at the same time"
			noi di "Please using space to separate each number"
			noi di "(e.g. 0 1 2 3 - Accessing information for general, ys, ye, ex at the same time)"
			noi di "If you would like to terminate the help option, please enter 99", _request(help_list)
			if ("${help_list}" == "") {
				
				local option
				
			}
			else if (strpos("${help_list}", "99") != 0){
				
				noi di "The help option is terminated."
				/*
				noi di "Please confirm to proceed (Y/N)"
				noi di "Current start year: " "`ys'"
				noi di "Current end year: " "`ye'", _request(help_terminator)
				if !((strupper(substr("${help_terminator}", 1, 1)) == "Y") | (substr("${help_terminator}", 1, 1) == "1")) {
					noi di as error "***********************WARNING**************************"
					noi di as error "***********User Requested To Abort Creation*************"
					noi di as error "********************************************************"
					exit
				}
				*/
				exit
				
			}
			else {
				local help_num
				foreach i in ${help_list} {
					if (inrange(`i', 0, 7)) {
						local help_num : di "`help_num'" " " "`i'"
					}
				}
				local option `help_num'
				
				/* temporary closure for structural change
				noi di "The help option is terminated."
				noi di "Please confirm to proceed (Y/N)"
				noi di "Current start year: " "`ys'"
				noi di "Current end year: " "`ye'", _request(help_terminator)
				if !((strupper(substr("${help_terminator}", 1, 1)) == "Y") | (substr("${help_terminator}", 1, 1) == "1")) {
					noi di as error "***********************WARNING**************************"
					noi di as error "***********User Requested To Abort Creation*************"
					noi di as error "********************************************************"
					exit
				}
				*/
				
			}
			
			if 1 {
				// display the title of help page
				
				noi di ""
				noi di "NHANES Dataset Creation Program Helper"
				noi di ""
			}
			
			if ("`option'" == "") {
			
				noi di as error "General Information:"
				
				// display the general help information
				noi di in g "This is a program written by Zhenghao (Vincent) Jin and Abimereki Muzaale."
				noi di "The purpose of this program is to produce a dataset based on US NHANES Data."
				noi di "This program downloads and combines data from NHANES website."
				noi di "For more information, please visit NHANES official website at: "
				noi di "https://www.cdc.gov/nchs/nhanes/index.htm"
				noi di ""
				noi di "If you would like to know more about us, please visit our website FENA at: "
				noi di "https://jhutrc.github.io/beta/intro.html" // waiting for the FENA website
				noi di ""
				
				// basic information
				noi di "Basic Program Information:"
				noi di "Current program version: v2.1"
				noi di "This version is developed under Stata v17"
				noi di ""
				
				// option information
				noi di as error "Basic Instructions of Operation:"
				noi di in g " "
				noi di in g "To call the program in Stata:"
				noi di "nhanes1"
				noi di ""
				noi di "To ask for instructions in Stata for this program:"
				noi di "Call for help option as below:"
				noi di "nhanes1, help"
				noi di ""
				noi di "To regulate the start year of data inquiry:"
				noi di "Call for ys option (year start) as below:"
				noi di "nhanes1, ys(year_you_want)"
				noi di ""
				noi di "To regulate the end year of data inquiry:"
				noi di "Call for ye option (year end) as below:" 
				noi di "nhanes1, ye(year_you_want)"
				noi di ""
				noi di "To regulate years to be excluded from inquiry:"
				noi di "Call for ex option (exclude) as below:"
				noi di "nhanes1, ex(list_of_years_you_want)"
				noi di "***each year should be separated by space***"
				noi di ""
				noi di "To regulate the information that is being inquired:"
				noi di "Call for ds option (dataset) as below:"
				noi di "nhanes1, ds(number_of_information_set_you_want)"
				noi di ""
				noi di "To regulate type of dataset produced by the inquiry:"
				noi di "Call for format option:"
				noi di "nhanes1, format(number_of_dataset_type_you_want)"
				noi di ""
				noi di "To regulate timeout time for the inquiry:"
				noi di "Call for timeout option as below:"
				noi di "nhanes1, timeout(time_you_want_to_set)"
				noi di ""
				noi di "***For Detail Information Of Each Option***"
				noi di "Please choose corresponding option name when calling help option."
				noi di ""
			
			}
			else {
			
				/*
					noi di "1 - For ys option"
					noi di "2 - For ye option"
					noi di "3 - For ex option"
					noi di "4 - For ds option"
					noi di "5 - For format option"
					noi di "6 - For timeout option"
				*/
				
				foreach i in `option' {
					
					if (`i' == 0) {
						noi di as error "General Information:"
				
						// display the general help information
						noi di in g "This is a program written by Zhenghao (Vincent) Jin and Abimereki Muzaale."
						noi di "The purpose of this program is to produce a dataset based on US NHANES Data."
						noi di "This program downloads and combines data from NHANES website."
						noi di "For more information, please visit NHANES official website at: "
						noi di "https://www.cdc.gov/nchs/nhanes/index.htm"
						noi di ""
						noi di "If you would like to know more about us, please visit our website FENA at: "
						noi di "https://jhutrc.github.io/beta/intro.html" // waiting for the FENA website
						noi di ""
						
						// basic information
						noi di "Basic Program Information:"
						noi di "Current program version: v2.1"
						noi di "This version is developed under Stata v17"
						noi di ""
						
						// option information
						noi di as error "Basic Instructions of Operation:"
						noi di in g " "
						noi di in g "To call the program in Stata:"
						noi di "nhanes1"
						noi di ""
						noi di "To ask for instructions in Stata for this program:"
						noi di "Call for help option as below:"
						noi di "nhanes1, help"
						noi di ""
						noi di "To regulate the start year of data inquiry:"
						noi di "Call for ys option (year start) as below:"
						noi di "nhanes1, ys(year_you_want)"
						noi di ""
						noi di "To regulate the end year of data inquiry:"
						noi di "Call for ye option (year end) as below:" 
						noi di "nhanes1, ye(year_you_want)"
						noi di ""
						noi di "To regulate years to be excluded from inquiry:"
						noi di "Call for ex option (exclude) as below:"
						noi di "nhanes1, ex(list_of_years_you_want)"
						noi di "***each year should be separated by space***"
						noi di ""
						noi di "To regulate the information that is being inquired:"
						noi di "Call for ds option (dataset) as below:"
						noi di "nhanes1, ds(number_of_information_set_you_want)"
						noi di ""
						noi di "To regulate type of dataset produced by the inquiry:"
						noi di "Call for format option:"
						noi di "nhanes1, format(number_of_dataset_type_you_want)"
						noi di ""
						noi di "To regulate timeout time for the inquiry:"
						noi di "Call for timeout option as below:"
						noi di "nhanes1, timeout(time_you_want_to_set)"
						noi di ""
						noi di "***For Detail Information Of Each Option***"
						noi di "Please choose corresponding option name when calling help option."
						noi di ""
					}
					
					if (`i' == 1) {
						// help for ys option
						noi di ""
						noi di as error "Help Page: ys option"
						noi di in g "The ys option, stands for year start, is the option "
						noi di "regulates the start year of the inquiry."
						noi di ""
						noi di "The range of value for ys option for this version is: "
						noi di "1988 to 2019"
						noi di "The default value for ys option if not specified by user is: "
						noi di "1999"
						noi di ""
						noi di "There are several years are special in NHANES study and ys option: "
						noi di "1." _col(10) "For any year start that is between 1988 - 1994, "
						noi di _col(10) "the whole NHANES III dataset (1988-1994) will be included"
						noi di "2." _col(10) "There is no information in NHANES datasets for"
						noi di "3." _col(10) "Data for pandemic period (NHANES 2019 - 2020) has been"
						noi di _col(10) "aggregated with the prior dataset (NHANES 2017 - 2019)"
						noi di ""
						noi di as error "Example of use: "
						noi di in g "To request information that starts at 1988: "
						noi di "nhanes1, ys(1988)"
						noi di "To request information that starts at 2002: "
						noi di "nhanes1, ys(2002)"
						noi di ""
					
					}
					
					if (`i' == 2) {
						// help for ye option
						noi di ""
						noi di as error "Help Page: ye option"
						noi di in g "The ye option, stands for year end, is the option "
						noi di "regulates the end year of inquiry. The end year should"
						noi di "always be larger than the start year (ys). "
						noi di "In addition, since NHANES datasets are aggregated into"
						noi di "two year periods, all end year in odds number will be"
						noi di "modified to 1 year later (e.g.: 2001 -> 2002)"
						noi di ""
						noi di "The range of value for ye option for this version is: "
						noi di "1989 to 2020"
						noi di "The default value for ye option is not specified by user is: "
						noi di "One year after year start"
						noi di "(E.g.: If ys is 1999, then default for ye will be 2000)"
						noi di "(E.g.: If ys is 2011, then default for ye will be 2012)"
						noi di ""
						noi di "There are several years are special in NHANES study and ye option: "
						noi di "1." _col(10) "For any end year that is between 1989-1994, "
						noi di _col(10) "ye option will be modified to 1995 to include NHANES III."
						noi di "2." _col(10) "For any end year that is between 2017-2020, "
						noi di _col(10) "ye option will be modified to 2020 to include pandemic dataset."
						noi di ""
						noi di as error "Example of use: "
						noi di in g "To request information that ends at 2012"
						noi di "nhanes1, ye(2012)"
						noi di "To request information that ends at 2016"
						noi di "nhanes1, ye(2016)"
						noi di ""
						
					}
					
					if (`i' == 3) {
						// help for ex option
						noi di ""
						noi di as error "Help Page: ex option"
						noi di in g "The ex option, stands for exclusion, is the option "
						noi di "speficies specific years to be excluded from the inquiry."
						noi di "Exclusion of multiple years in one inquiry is ok."
						noi di "To exclude mutilple years at the same time, use space (" ")"
						noi di "to separate years in ex option."
						noi di "Since NHANES datasets are aggregated into two year periods,"
						noi di "excluding one single year will result in excluding the whole"
						noi di "dataset for the two year period."
						noi di "E.g.: excluding year of 1999 will results in exluding all"
						noi di "information for 1999 - 2000"
						noi di ""
						noi di "The ex option takes any positive integers as years"
						noi di "No default value was set for ex option"
						noi di "No actions will be taken for years specified in ex option that "
						noi di "are not between the inquiry start (ys) and end (ye) period"
						noi di ""
						noi di "There are several years are special in NHANES study and ex option: "
						noi di "1." _col(10) "For any year between 1988 - 1994 specified in ex option, "
						noi di _col(10) "the whole NHANES III dataset will be excluded"
						noi di "2." _col(10) "For any year between 2017 - 2020 specified in ex option, "
						noi di _col(10) "The whole dataset for pandemic period (2017 - 2020) will be excluded"
						noi di ""
						noi di as error "Example of use: "
						noi di in g "To exclude year of 2014 in an inqury between 2009 - 2016"
						noi di "nhanes1, ys(2009) ye(2016) ex(2014)"
						noi di "To exclude year of 2004, 2007, 2010 in an inqury between 1999 - 2016"
						noi di "nhanes1, ys(1999) ye(2016) ex(2004 2007 2010)"
						
					}
					
					if (`i' == 4) {
						// help for ds option
						noi di ""
						noi di as error "Help Page: ds option"
						noi di in g "The ds option, stands for dataset, is th option "
						noi di "specifies the type of information requested by the inqury."
						noi di ""
						noi di "The range of value for ds option is:"
						noi di "1 to 5"
						noi di "The default value for ds option if not specified by user is:"
						noi di "1 (base dataset)"
						noi di ""
						noi di "Coding for ds option: "
						noi di "0" _col(20) "base dataset"
						noi di "1" _col(20) "demographic dataset"
						noi di "2" _col(20) "exam dataset"
						noi di "3" _col(20) "diet dataset"
						noi di "4" _col(20) "questionnaire dataset"
						noi di ""
						noi di "Information included for each dataset:"
						noi di "base dataset:"
						noi di _col(20) "age"
						noi di _col(20) "sex"
						noi di _col(20) "race"
						noi di "demographic dataset:"
						noi di _col(20) "Please refer to NHANES website for specific list at" 
						noi di _col(20) "https://wwwn.cdc.gov/nchs/nhanes/default.aspx"
						noi di "exam dataset:"
						noi di _col(20) "blood pressure"
						noi di _col(20) "body measures"
						noi di _col(20) "glycohemoglobin"
						noi di _col(20) "albumin & creatinine - urine"
						noi di _col(20) "standard biochemistry profile"
						noi di "diet dataset:"
						noi di _col(20) "dietary intake - total nutrients intake"
						noi di "questionnaire dataset"
						// ALQ DIQ KIQ_U BPQ PAQ SMQ 
						noi di _col(20) "alcohol use"
						noi di _col(20) "diabetes"
						noi di _col(20) "kidney condition - urology"
						noi di _col(20) "blood pressure & cholesterol"
						noi di _col(20) "physical activity"
						noi di _col(20) "smoking - cigarette use"
						noi di ""
						noi di as error "Example of use: "
						noi di in g "To request the base dataset for age/sex/race:"
						noi di "nhanes1, ds(1)"
						noi di "To request the exam dataset for examination and lab information:"
						noi di "nhanes1, ds(3)"
						noi di ""
					}
					
					if (`i' == 5) {
						// help for format option
						noi di ""
						noi di as error "Help Page: format option"
						noi di in g "The format option regulates the type of output file"
						noi di ""
						noi di "The range of value for format option is:"
						noi di "0 to 5"
						noi di "The default value for format option if not specified by user is:"
						noi di "0 (dta dataset)"
						noi di ""
						/*
						0 - dta
						1 - excel
						2 - csv
						3 - sasxport5
						4 - sasxport8
						5 - text
						*/
						noi di "Coding for format option:"
						noi di "0" _col(10) "dta (Stata) file"
						noi di "1" _col(10) "excel file"
						noi di "2" _col(10) "csv file"
						noi di "3" _col(10) "sasxport5 file"
						noi di "4" _col(10) "sasxport8 file"
						noi di "5" _col(10) "text data file with dictionary"
						noi di ""
						noi di as error "Example of use:"
						noi di in g "To output a csv file for the inquiry:"
						noi di "nhanes1, format(2)"
						noi di "To output a sas dataset (sasxport5) for the inquiry:"
						noi di "nhanes1, format(3)"
						noi di ""
						
					}
					
					if (`i' == 6) {
						// help for timeout option
						noi di ""
						noi di as error "Help Page: timeout option"
						noi di in g "The timeout option regulates the time limit for Stata"
						noi di "to initial connection with NHANES database in seconds."
						noi di ""
						noi di "The timeout option takes any positive integers."
						noi di "The default value for timeout option if not specified by user is:"
						noi di "1000"
						noi di ""
						noi di as error "Example of use:"
						noi di in g "To set time limit for initiate connection to be 500 seconds:"
						noi di "nhanes1, timeout(500)"
						noi di "To set time limit for initiate connection to be 2000 seconds:"
						noi di "nhanes1, timeout(2000)"
						noi di ""
					}
					
					if (`i' == 7) {
						// help for mort option
						noi di ""
						noi di as error "Help Page: mort option"
						noi di in g "The mort option allows user to add mortality information to the inquiry."
						noi di "The mort option can be combined with any requested information."
						noi di ""
						noi di as error "Example of use:"
						noi di in g "To request mortality information additional to demographic information"
						noi di "for year between 1999 to 2008"
						noi di "nhanes1, ys(1999) ye(2008) mort ds(2)"
					}
					
				}
				
			}
			
			noi di as error "Do you need further help? (Y/N)", _request(f_help)
			if ((substr(strupper("${f_help}"), 1, 1) == "Y") | (substr("${f_help}", 1, 1) == "1")) {
				noi nhanes1_helper
			} 
			else {
				noi di "The help option is terminated."
				exit
			}
			
		}
	
	end
}
