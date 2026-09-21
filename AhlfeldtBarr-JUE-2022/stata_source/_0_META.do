**********************************************************************
*** Teaching directory for 										   ***
*** Codes for Ahlfeldt & Barr (2022): The economics of skyscrapers ***											
*** (c) Gabriel M Ahlfeldt, 01/2024								   ***
**********************************************************************

*** This teaching directory to Ahlfeldt & Barr (2022): The economics of skyscrapers, Journal of Urban Economics, 129
*** contains the main programs used to solve the model and to invert fundamental production and residential amenities
*** All codes have been commented to enhance accessibility

*** This do file calles other do files used to generate central outputs

* Set working directory
	cd "D:\Dropbox\_HUB_HerreraA\Course\Repository\Replication_Directories\AB2022-toolkit\STATA_WALKTHROUGH"
	capture mkdir "FIGS"
	search grc1leg // install if necessary	
	capture set scheme s1color
	capture set scheme stcolor	// This will only be executed when using Stata 18 or newer
	
* Counterfactuals with stylized city structure
	global large = 0			// If set to 1, generates a larger city. Only useful if you want to generate massive cities with huge population for numerical purposes. 0 is the recommended value. Graphs look better, and model solves faster.
	do "_1_PROGS.do" 			// Load programmes for simulations
	do "_2_ANALYSIS.do"			// Run counterfactuals
	
* Inverting the model to rationalize Chicago height gradient & do over-ID
	do "_3_INVERSION.do"		// Run counterfactuals

* Counterfactuals with inverted fundamentals
	do "_4_INVERTEDCOUNTER.do"	// Generating a secondary CBD
	
************************************************
*** All do files file successfully completed ***	
************************************************