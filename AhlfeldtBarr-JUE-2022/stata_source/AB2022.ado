**********************************************************************
*** Gradient solver	for 										   ***
*** Codes for Ahlfeldt & Barr (2022): The economics of skyscrapers ***	
*** Version 0.93												   ***											
*** (c) Gabriel M Ahlfeldt, 02/2024								   ***
**********************************************************************



/*************
*** Syntax ***
**************

Once you have copied AB2022.ado file to your ado folder, Stata will generate floor space, height, and land rent gradients under the baseline parameterization. 

Below are a brief description and canonical parameter values from Ahlfeldt & Barr (2022)

theta_C 	Commercial construction cost elasticity of height	0.5
theta_R		Residential construction cost elasticity of height	0.55
omega_C 	Commercial rent elasticity of height				0.03
omega_R 	Residential rent elasticity of height				0.07
beta_C 		Agglomeration elasticity							0.03
a_bar_C 	Fundamental production amenity						2
a_bar_R 	Fundamental residential amenity						1
tau_C 		Commercial amenity decay							0.01
tau_R 		Residential amenity decay							0.005
c_C 		Commerical baseline construction factor				1.4
c_R 		Residential baseline construction factor			1.4
r_a 		Agricultural land rents								150
S_bar_C 	Commercial height limit								999
S_bar_R		Residential height limit							999
x_1_bar		Urban growth boundary (max. radius of urban area)	999

You can change any of the parameters by adding "parameter=value" as an argument. 
For example, to add a 50 floor height limit for commercial buildings, you can enter: 

AB2022 "S_bar_C=50"

You can enter multiple alterations to baseline parameter values. For example, to add   
to add a 50 floor height limit for commercial buildings and a 10-floor height limit for residential buildings, you can enter:

AB2022 "S_bar_R=100"

And so forth...

Enjoy using the toolkit */ 

***********************************************
*** Define programmes for solving the model ***
***********************************************	
	
********************************************************************************
// Solver that computes all endogenous variables, including aggregate labour demand and supply	 
capture program drop SOLVER
program SOLVER // Start defineing SOLVER program

	// Clear any pre-existing value
	qui foreach var of varlist A_tilde_x_C A_tilde_x_R a_x_C a_x_R r_x_C r_x_R U S_x_C S_x_R p_bar_x_C 	p_bar_x_R L_x_C f_bar_x_R n_x {
		replace `var' = .
	}
	
	// Compute Amenity and floor space shifters derived on p. 7
		qui replace A_tilde_x_C = a_bar_C * a_rand_C * L^beta_C * exp(-tau_C*D)	// Commercial production amenity shifter
		qui replace A_tilde_x_R = a_bar_R * a_rand_R * L^beta_R * exp(-tau_R*D)	// Residential amenity shifter
		qui replace a_x_C = A_tilde_x_C^(1/(1-alpha_C))*y^(alpha_C/(alpha_C-1))	// Commercial floor space shifter
		qui replace a_x_R = A_tilde_x_R^(1/(1-alpha_R))*y^(1/(1-alpha_R))		// Residential floor space shifter

	// Land rent according to Eq. (7)
		qui replace r_x_C = a_x_C / (1+omega_C)*(a_x_C / (c_C*(1+omega_C)))^((1+omega_C)/(theta_C-omega_C))-c_C*(a_x_C/(c_C*(1+theta_C)))^((1+theta_C)/(theta_C-omega_C))
		qui replace r_x_R = a_x_R / (1+omega_R)*(a_x_R / (c_R*(1+omega_R)))^((1+omega_R)/(theta_R-omega_R))-c_R*(a_x_R/(c_R*(1+theta_R)))^((1+theta_R)/(theta_R-omega_R))	
		
	// Define land use
		qui replace U = . 														// Clear any pre-existing value
		qui replace U = 3 if r_a > r_x_C & r_a > r_x_R							// Agricultural land rent is highest
		qui replace U = 2 if r_x_R  > r_x_C & r_x_R > r_a						// Residential land rent is highest
		qui replace U = 1 if r_x_C > r_x_R & r_x_C > r_a						// Commercial land rent is highest
	
	// Outer edge of residential zone	
		qui sum x if U < 3 &  x >= 0 											// Take the inner margin of the agricultural zone
		qui scalar x1 = r(max)
		qui local x1 = x1	
		
	// Outer edge of residential zone	
		qui sum x if U != 1 & x >= 0 											// Take the inner margin of any zone that is not commerical
		qui scalar x0 = r(min)
		qui local x0 = x0

	// Profit maximizing height defined just above Eq. (6)
		qui replace S_star_x_C = ( (a_x_C / (c_C*(1+theta_C)))^(1/(theta_C-omega_C)) ) if U == 1
		qui replace S_star_x_R = ( (a_x_R / (c_R*(1+theta_R)))^(1/(theta_R-omega_R)) ) if U == 2
		
	// Compute actual height according to Eq. (6)	
		qui replace S_x_C = min(S_bar_C, (a_x_C / (c_C*(1+theta_C)))^(1/(theta_C-omega_C)) ) if U == 1
		qui replace S_x_R = min(S_bar_R, (a_x_R / (c_R*(1+theta_R)))^(1/(theta_R-omega_R)) ) if U == 2
	
	// Comute bid rents
		qui replace p_bar_x_C = a_x_C * 1/(1-omega_C)*S_x_C^omega_C				// According to Eq. (2)
		qui replace p_bar_x_R = a_x_R * 1/(1-omega_R)*S_x_R^omega_R				// According to Eq. (4)
	
	// Compute labour demand
		qui replace L_x_C = alpha_C/(1-alpha_C)*p_bar_x_C / y * S_x_C if U == 1 // using MRS in Eq. (12)
		qui sum L_x_C															// Read values across locations														
		qui scalar L_hat_demand = r(sum)										// Sum across locations to get aggregate labour demand accordin gto Eq. (14)
		
	// Compute labour supply 
		qui replace f_bar_x_R = (1-alpha_R)/p_bar_x_R*y if U == 2				// Floor space per capita from Marshallian demand function
		qui replace n_x = S_x_R / f_bar_x_R										// Get number of workers from total floor space and floor space per capita
		qui sum n_x																// Read number of workers across locations
		qui return list
		qui scalar L_hat_supply = r(sum)										// Sum across locations to get aggregate labour supply

	// Compute some final statistics	
		qui sum L																// Summarize employment
		qui scalar sL = r(mean)													// Save employment as scalar
		qui sum y																// Summarize wage												
		qui scalar sy = r(mean)													// Save wage as scalar
		qui replace S_x = max(S_x_C,S_x_R)										// Save highest possible building of any potential use at a given location
		scalar list L_hat_demand L_hat_supply sL sy								// Output wage adjustement factor, labour demand & supply, and wage for inspection of the user
	// Land use
		qui replace URBAN = U < 3													// Update indicator for urban use	
		qui replace COM = U == 1													// Update indicator for commercial use
	// Solver ends	
	end
// SOLVER programme completed
********************************************************************************


********************************************************************************
// Algorithm that implmenents iterative procedure to find market-clearing wage	 

capture program drop WAGE
program define WAGE
		qui local obj_int = abs((L_hat_demand+0)/(L_hat_supply+0)-1)			// Define value in the objective function of the internal loop, the percentage difference between demand and supply in the model. 
																				// When this relative difference approaches zero, we have found the WAGE that clears the market 
		while `obj_int' > 0.01 {												// Keep iterating while objective is larger than tolerance level	
			
			// Compute the wage adjustment factor to be used in the iterative procedure to solve for the wage
			qui if L_hat_supply == 0 { // If no labour supply increase wage by 20%
				scalar y_factor = 1.2
			}
				else {
					if L_hat_demand == 0 { // If no labour demand decrease wage by 20%
						scalar y_factor = 0.8
					}
				else {
						scalar  y_factor =  (L_hat_demand / L_hat_supply)^0.01 // Adjust depending on the ration of labour demand to labour supply => if labour demand exceeds supply, increase the wage and vice versa
				}
				}			
			// Now adjust the wage accordingly
			qui replace y = 0.5 * y + 0.5*y*y_factor						// Use the average of the old wage and a new wage, updated by the wage adjustment factor
			// Update the model solutions with the new wage
			qui SOLVER		
			
			// Compute the new internal objective value and report it alongside the current external objective value to the user
			qui local obj_int = abs((L_hat_demand+0.0001)/(L_hat_supply+0.0001)-1)
			display "internal objective " round(`obj_int', 0.001) " external objective " round( `obj_ext', 0.001)
		
		// Close the internal loop, if we pass beyond the point, we have found the market-clearing wage conditional on our guess of total employment
		} 
end
// Programme completed
********************************************************************************


********************************************************************************
// Graph height
capture program drop  GHEIGHT
program GHEIGHT 
		qui sum S_x_C
		qui local top = (int(r(max)/25)+1)*25
		qui replace  SHADE = URBAN*COM*(int(`top'/1)+1)
		qui replace SHADEU = URBAN*(int(`top')+1)	
		qui sum S_x_C
		qui local max = round(r(max),25)
		qui local step = round(`max'/4, 5)
		twoway  (area SHADEU x, color(gs14)) (area SHADE x, color(gs12)) (line S_x_C x , color(red)  	   lpattern(longdash)) (line S_x_R x if x < 0, color(blue) lpattern(shortdash)) (line S_x_R x if 		x > 0, color(blue) lpattern(shortdash)), ///
		graphregion(color(white)) ylabel(0[`step']`max') xlabel(-50[25]50) xtitle("") ytitle("Building height") legend(order(3 "Commercial" 4 "Residential") cols(2) pos(6)) title("Building height") name("`1'", replace) xsize(5) ysize(5)   nodraw  
		end // Solver ends
********************************************************************************

********************************************************************************
// Graph bid rent
capture program drop GBIDRENT
program GBIDRENT 
		*local urban_max = urban_max
		*local urban_min = urban_min
		*local com_max = com_max
		qui sum p_bar_x_C
		qui local top = (int(r(max)/5)+1)*5
		qui replace  SHADE = URBAN*COM*(int(`top'/1))
		qui replace SHADEU = URBAN*(int(`top'))	
		twoway  (area SHADEU x, color(gs14)) (area SHADE x, color(gs12)) (line p_bar_x_C x , color(red) lpattern(longdash)) (line p_bar_x_R x if x > 0, color(blue) lpattern(shortdash)) (line p_bar_x_R x if x < 0, color(blue) lpattern(shortdash))  ///
		,   graphregion(color(white)) xlabel(-50[25]50) xtitle("") ytitle("Floor space rent") legend(order(1 "Urban area" 2 "CBD area" 3 "Commercial" 4 "Residential") cols(4) pos(6)) title("Floor space rent") name(`1', replace)    nodraw  
		end // Solver ends
********************************************************************************


********************************************************************************
// Graph land rent
capture program drop GLANDRENT
program GLANDRENT
	qui sum r_x_C
	qui local top = (int(r(max)/500)+1)*500	
	qui replace SHADE = URBAN*COM*(int(`top')+1)
	qui replace SHADEU = URBAN*(int(`top')+1)
	twoway (area SHADEU x, color(gs14)) (area SHADE x, color(gs12)) (line r_x_C x , color(red)  lpattern(longdash)) (line r_x_R x , color(blue) lpattern(shortdash)) (line r_a x, color(black) lpattern(solid)), ///
	  graphregion(color(white))  xlabel(-50[25]50) xtitle("") ytitle("Land bid rent") legend(order(3 "Commercial" 4 "Residential" 5 "Agricultural" 1 "Urban area" 2 "CBD")  size(vsmall) cols(5) pos(6)) title("Land bid rent") name(`1', replace)   nodraw  
end

********************************************************************************
// Bar height
capture program drop GHEIGHTB
program GHEIGHTB 
		qui sum S_x
			qui gen bin = int(D)
			qui egen mS = mean(S_x), by(bin)
			qui egen sdS = sd(S_x), by(bin)
			qui gen cvS = sdS/mS
			qui sum cvS 		
			qui local CV = round(r(mean), 0.001)
			qui sum cvS if bin <= 5	
			qui local CVcbd = round(r(mean), 0.001)
			qui drop bin mS sdS cvS
		qui sum S_x_C
		qui local top = (int(r(max)/25)+1)*25
		qui replace  SHADE = URBAN*COM*(int(`top'/1)+1)
		qui replace SHADEU = URBAN*(int(`top')+1)	
		qui sum S_x_C
		qui local max = round(r(max),5)
		qui local step = round(`max'/4, 5)
		twoway  (bar S_x_C x if U == 1, color(red) lcolor(none) barwidth(0.1) bargap(1)	   lpattern(solid)) (bar S_x_R x if U == 2, barwidth(0.01) lcolor(none) color(blue) lpattern(solid)) (bar S_x_R x if 	U == 2, barwidth(0.01) color(blue) lcolor(none) lpattern(solid)), ///
		graphregion(color(white)) ylabel(0[`step']`max') xlabel(-50[25]50) xtitle("") ytitle("Building height") legend(order(1 "Commercial" 2 "Residential") cols(2) pos(6))  name("`1'", replace) xsize(10) ysize(5) note("Mean coefficent of variation within one-distance-unit bins = `CV'")   nodraw   // (area SHADEU x, color(gs14)) (area SHADE x, color(gs12)) 
		end // Solver ends
		
********************************************************************************	
// 	Iterative programme that nests SOLVER to find equilibrium conditional on parameter choices

program FINDEQ // Syntax theta_C theta_R omega_C omega_R beta_C a_bar_C a_bar_R tau_C tau_R c_C c_R r_a S_bar_C S_bar_R
	
	// Set scalars based on parameter choices by the user
	scalar theta_C = `1'
	scalar theta_R = `2'
	scalar omega_C = `3'
	scalar omega_R = `4'	
	scalar beta_C = `5'
	scalar a_bar_C = `6'
	scalar a_bar_R = `7'
	scalar tau_C = `8'
	scalar tau_R = `9'
	scalar c_C = `10'
	scalar c_R = `11'
	qui replace r_a = `12'	
	scalar S_bar_C = `13'
	scalar S_bar_R = `14'
	
	// Initialization
	SOLVER 																		// Run initial round of the solver
	local obj_ext = abs(L/(0.5*(L_hat_demand+L_hat_supply))-1)					// Define value in the objective function of the external loop, the percentage difference between total employment L (here the set value) and hte average of current demand and supply in the model
																				// When this relative difference approaches zero, we have have found TOTAL EMPLOYMENT
		
	// Start external loop until we converged
	while 	`obj_ext' > 0.01{													// Keep iterating while objective is larger than tolerance level		
	
		// Initialization for the inner loop  	
		qui SOLVER																// Use SOLVER to solve for endogenous objects
		
		// Solve for market clearing wage										// Use WAGE programme
		WAGE
	// Test if the city is sustainable	
	local test = L_hat_demand+L_hat_supply 										// Check if labour demand and supply converged to zero
		if `test' == 0 { // If so, stop
			display "Error: City does not reach critical size."
			display "Increase productivity to make city more attractive"
			stop
		}
			else {	// Else just continue
				
			}
	// Update the value of the external objective function		
	local obj_ext = abs(L/(0.5*(L_hat_demand+L_hat_supply))-1)					// Percentage difference between the guess of total employment and the average of demand and supply within the model
	
	// Update our guess of total employment
	qui replace L = 0.5 * L + 0.5*0.5*(L_hat_demand + L_hat_supply)				// We use the average of the old guess of total employment and the average of updated demand and supply within the model
	
	// Close outer loop, if we pass this point, we have found a wage that clears the labour market (internal loop) and total employment (outer loop)
	}
	
	// Report key labour market outcomes to reader
	scalar list L_hat_demand L_hat_supply sL sy x0 x1

	// Solver ends	
	end
// Programme completed
********************************************************************************
		
********************************************************************************	
// 	Iterative programme that nests SOLVER to find equilibrium conditional on parameter choices - ADO file version
capture program drop AB2022
program AB2022 // Syntax theta_C theta_R omega_C omega_R beta_C a_bar_C a_bar_R tau_C tau_R c_C c_R r_a S_bar_C S_bar_R
	
	version 17.0 	// The program has not been tested with earlier version, but is likely to run on them
	capture set scheme s1color
	capture set scheme stcolor


************************************
*** Generate artificial data set ***
************************************

	local large = 0						// Reading user choice of city size
	clear								// Drop all data
	program drop _all					// Drop all programs
	
* Residential parameters
	scalar alpha_R = 0.66				// Set expenditure share on non-housing goodes
	scalar beta_R = 0.00				// Set residential agglomeration elasticity to zero
	scalar tau_R = 0.01					// Set residential amenity decay
	scalar omega_R = 0.1				// Set residential height elasticity of rents
	scalar theta_R = 0.5				// Set residential height elasticity of construction cost
	scalar c_R = 1						// Set residential baseline construction cost
	scalar a_bar_R = 1					// Set residential fundamental amenity

* Commercial parameters
	scalar alpha_C = 0.85				// Set input share of floor space in production
	scalar beta_C = 0.04				// Set commercial aggloemration elasticity
	scalar tau_C = 0.01					// Set commercial productivity decay
	scalar omega_C = 0.1				// Set commerical height elasticity of rents
	scalar theta_C = 0.5				// Set commerical height elasticity of construction cost
	scalar c_C = 1						// Set commercial baseline construction cost
	scalar a_bar_C = 1.5				// Set commerical fundamental productivity
	
if `large' == 1 {						// If city is large
	qui set obs 20001						// Generate 20001 observations, one for each location along a line. 1 for the center. 10000 in either direction
	qui gen x = _n -10001					// Generate running variable, 0 at center
	qui replace x = x / 100					// Rescale running variable so that is can be thought of being in km
	qui gen D = abs(x)						// Generate distance from the centre as absolute value of the running variable
}
else {									// city is small 
	qui set obs 10001						// Generate 10001 observations, one for each location along a line. 1 for the center. 5000 in either direction
	qui gen x = _n -5001					// Generate running variable, 0 at center
	qui replace x = x / 100					// Rescale running variable so that is can be thought of being in km
	qui gen D = abs(x)						// Generate distance from the centre as absolute value of the running variable
}

// Starting values
	qui gen y =  2.5 						// Starting value for wage
	qui gen L =  1000000 					// Starting value for employment
	
// Generate placeholders for objects to be set by user 
	qui gen r_a = 100						// Exogenous agricutural land rent
	qui gen a_x_C = . 						// Fudnamental production amenity
	qui gen a_x_R = . 						// Fundamental residential amenity
	
// Generate placeholders for objects to be solved
	qui gen a_rand_C = 1					// Location specific component in production amenity that can be used to generate fuzzy height gradients (not relevant when working with stylized cities)
	qui gen a_rand_R = 1					// Location specific component in residential amenity that can be used to generate fuzzy height gradients (not relevant when working with stylized cities)
	qui gen A_tilde_x_C = . 				// Compositive productivity shifter
	qui gen A_tilde_x_R = . 				// Compositive amenity shifter
	qui gen r_x_C= .						// Commercial land rent
	qui gen r_x_R= . 						// Residential land rent
	qui gen U = . 							// Auxilliary variable for land use category
	qui gen S_star_x_C = .					// Profit-maximizing commerical height
	qui gen S_star_x_R = . 					// Profit-maximizing residential height
	qui gen S_x_C = . 						// Realized commercial height
	qui gen S_x_R = .						// Realized residential height
	qui gen S_x = . 						// Height	
	qui gen p_bar_x_C = .					// Commericial horizontal bid rent
	qui gen p_bar_x_R = .					// Residential horizontal bid rent
	qui gen L_x_C = .						// Labour input
	qui gen f_bar_x_R = .					// Residential floor space demand per capita
	qui gen n_x = .							// Labour supply
	qui gen URBAN = .						// Auxiliary idicator for urban user	
	qui gen COM = .							// Auxiliary indicator for commerical use
	qui gen  SHADE = .						// Auxiliary indicator used in graphs
	qui gen SHADEU = .						// Auxiliary indicator used in graphs	
	
	// Set scalars based on parameter choices by the user
	scalar theta_C 	=	0.5
	scalar theta_R	=	0.55
	scalar omega_C 	=	0.03
	scalar omega_R 	=	0.07
	scalar beta_C 	=	0.03
	scalar a_bar_C 	=	2
	scalar a_bar_R 	=	1
	scalar tau_C 	=	0.01
	scalar tau_R 	=	0.005
	scalar c_C 		=	1.4
	scalar c_R 		=	1.4
	scalar r_a 		=	150
	scalar S_bar_C 	=	999
	scalar S_bar_R	=	999
	scalar x_1_bar 	=	999
	
// Override with user entry
	capture scalar `1'
	capture scalar `2'
	capture scalar `3'
	capture scalar `4'
	capture scalar `5'
	capture scalar `6'
	capture scalar `7'
	capture scalar `8'
	capture scalar `9'
	capture scalar `10'
	capture scalar `11'
	capture scalar `12'
	capture scalar `13'
	capture scalar `14'
	capture scalar `15'
	qui drop r_a
	qui gen r_a = r_a
	
	
// Initialization
	display "<<< solving for the equilibrium >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	qui drop if D > x_1_bar														// Dropping anything outside the user-specified urban growth boundary
	qui SOLVER 																		// Run initial round of the solver
	local obj_ext = abs(L/(0.5*(L_hat_demand+L_hat_supply))-1)					// Define value in the objective function of the external loop, the percentage difference between total employment L (here the set value) and hte average of current demand and supply in the model
																				// When this relative difference approaches zero, we have have found TOTAL EMPLOYMENT
		
	// Start external loop until we converged
	while 	`obj_ext' > 0.01{													// Keep iterating while objective is larger than tolerance level		
	
		// Initialization for the inner loop  	
		qui SOLVER																// Use SOLVER to solve for endogenous objects
		
		// Solve for market clearing wage										// Use WAGE programme
		qui WAGE
	// Test if the city is sustainable	
	local test = L_hat_demand+L_hat_supply 										// Check if labour demand and supply converged to zero
		if `test' == 0 { // If so, stop
			display "Error: City does not reach critical size."
			display "Increase productivity to make city more attractive"
			stop
		}
			else {	// Else just continue
				
			}
	// Update the value of the external objective function		
	local obj_ext = abs(L/(0.5*(L_hat_demand+L_hat_supply))-1)					// Percentage difference between the guess of total employment and the average of demand and supply within the model
	
	// Update our guess of total employment
	qui replace L = 0.5 * L + 0.5*0.5*(L_hat_demand + L_hat_supply)				// We use the average of the old guess of total employment and the average of updated demand and supply within the model
	
	// Close outer loop, if we pass this point, we have found a wage that clears the labour market (internal loop) and total employment (outer loop)
	}
	
	// Report key labour market outcomes to reader
	foreach name in  L_hat_demand L_hat_supply sL sy x0 x1 {
		local `name' = round(`name',0.01)
		}
	display "Labour demand: `L_hat_demand'"
	display "Labour supply: `L_hat_supply'"
	display "Total employment: `sL'"
	display "Wage: `sy'"
	display "CBD radius: `x0'"
	display "Urban radius: `x1'"

* Create graph	
	display "<<< drawing graphs >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	qui GBIDRENT name1
	qui GHEIGHT name0
	qui GLANDRENT name3
	local sL = round(`sL',1)
	grc1leg name1 name0 name3, cols(3) leg(name3) xsize(15) ysize(5) note("Total employment: `sL'" "Wage: `sy'" "CBD radius: `x0'" "Urban area radius: `x1'")	
	display "<<< all done >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"	
	// Solver ends	
	end
// Programme completed
********************************************************************************
