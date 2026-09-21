**********************************************************************
*** Teaching directory for 										   ***
*** Codes for Ahlfeldt & Barr (2022): The economics of skyscrapers ***											
*** (c) Gabriel M Ahlfeldt, 01/2024								   ***
**********************************************************************

*** This do file solves for production and residential amenities that match the 
*** fuzzy height gradient of Chicago and a user-defined target population

*** Keep in mind that the height data are for tall buildings only, but, in reality,
*** many households live in shorter buildings. Hence, the population of Chicago 
*** is not a sensible population target (lower values will be more plausible)

// Load synthetic data
	u "DATA\SIMULATION\INVERTED.dta", clear										// Load synthetic data set with inverted fundamentals
	GHEIGHTB name0																// Plot city structure
// Use equilibrium solver to update endogenous variables	
// Syntax 	theta_C	theta_R	omega_C	omega_R	beta_C	a_bar_C	a_bar_R	tau_C	tau_R	c_C		c_R 	r_a 	S_bar_C	S_bar_R
	FINDEQ 		0.5 	0.55 	0.03 	0.07 	0.030 	2		1 		0.01 	0.005 	1.3 	1.3 	50		999		999	
	scalar list sy sL															// List equilibrium wage and total employment
	
**********************************
*** At a subcenter at 12-14 km ***		
**********************************	
	
// Update fundamentals to reflect a subcenter	
	replace a_rand_C = 1.15*runiform() if x > 12 & x < 15						// Generate positive production amenities with some randomness within 12 and 15km
	egen temp = max(a_rand_C), by(CONVBIN)										// Convert to location bins
	replace a_rand_C = temp														// Update production amenities 

// Use equilibrium solver to update endogenous variables	
// Syntax 	theta_C	theta_R	omega_C	omega_R	beta_C	a_bar_C	a_bar_R	tau_C	tau_R	c_C		c_R 	r_a 	S_bar_C	S_bar_R
	FINDEQ 		0.5 	0.55 	0.03 	0.07 	0.030 	2		1 		0.01 	0.005 	1.3 	1.3 	50		999		999
	GHEIGHTB name0																// Convert to location bins							
	scalar list sy sL															// List equilibrium wage and total employment
	
*** Observe that we have generated a new subcenter with commercial height as expected
*** Notice that the old center has lost height since there is relocation to the new center
*** Yet, overall the city is larger since it is more productive overall. 
*** Wage increases by about 5%
*** Population increases by 30%

**************************************
*** Do file successfully completed ***	
**************************************	
