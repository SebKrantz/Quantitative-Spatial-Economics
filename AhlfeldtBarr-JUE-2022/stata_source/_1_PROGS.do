**********************************************************************
*** Teaching directory for 										   ***
*** Codes for Ahlfeldt & Barr (2022): The economics of skyscrapers ***											
*** (c) Gabriel M Ahlfeldt, 01/2024								   ***
**********************************************************************

*** This do file creates a baseline data set and loads relevant programmes

************************************
*** Generate artificial data set ***
************************************

	local large = $large				// Reading user choice of city size
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
	set obs 20001						// Generate 20001 observations, one for each location along a line. 1 for the center. 10000 in either direction
	gen x = _n -10001					// Generate running variable, 0 at center
	replace x = x / 100					// Rescale running variable so that is can be thought of being in km
	gen D = abs(x)						// Generate distance from the centre as absolute value of the running variable
}
else {									// city is small 
	set obs 10001						// Generate 10001 observations, one for each location along a line. 1 for the center. 5000 in either direction
	gen x = _n -5001					// Generate running variable, 0 at center
	replace x = x / 100					// Rescale running variable so that is can be thought of being in km
	gen D = abs(x)						// Generate distance from the centre as absolute value of the running variable
}

// Starting values
	gen y =  2.5 						// Starting value for wage
	gen L =  1000000 					// Starting value for employment
	
// Generate placeholders for objects to be set by user 
	gen r_a = 100						// Exogenous agricutural land rent
	gen a_x_C = . 						// Fudnamental production amenity
	gen a_x_R = . 						// Fundamental residential amenity
	
// Generate placeholders for objects to be solved
	gen a_rand_C = 1					// Location specific component in production amenity that can be used to generate fuzzy height gradients (not relevant when working with stylized cities)
	gen a_rand_R = 1					// Location specific component in residential amenity that can be used to generate fuzzy height gradients (not relevant when working with stylized cities)
	gen A_tilde_x_C = . 				// Compositive productivity shifter
	gen A_tilde_x_R = . 				// Compositive amenity shifter
	gen r_x_C= .						// Commercial land rent
	gen r_x_R= . 						// Residential land rent
	gen U = . 							// Auxilliary variable for land use category
	gen S_star_x_C = .					// Profit-maximizing commerical height
	gen S_star_x_R = . 					// Profit-maximizing residential height
	gen S_x_C = . 						// Realized commercial height
	gen S_x_R = .						// Realized residential height
	gen S_x = . 						// Height	
	gen p_bar_x_C = .					// Commericial horizontal bid rent
	gen p_bar_x_R = .					// Residential horizontal bid rent
	gen L_x_C = .						// Labour input
	gen f_bar_x_R = .					// Residential floor space demand per capita
	gen n_x = .							// Labour supply
	gen URBAN = .						// Auxiliary idicator for urban user	
	gen COM = .							// Auxiliary indicator for commerical use
	gen  SHADE = .						// Auxiliary indicator used in graphs
	gen SHADEU = .						// Auxiliary indicator used in graphs

// Save data set
	save "DATA\SIMULATION\BASE.dta", replace	
	
***********************************************
*** Define programmes for solving the model ***
***********************************************	
	
********************************************************************************
// Solver that computes all endogenous variables, including aggregate labour demand and supply	 

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
		qui replace r_x_C = a_x_C / (1+omega_C)*(a_x_C / (c_C*(1+theta_C)))^((1+omega_C)/(theta_C-omega_C))-c_C*(a_x_C/(c_C*(1+theta_C)))^((1+theta_C)/(theta_C-omega_C))
		qui replace r_x_R = a_x_R / (1+omega_R)*(a_x_R / (c_R*(1+theta_R)))^((1+omega_R)/(theta_R-omega_R))-c_R*(a_x_R/(c_R*(1+theta_R)))^((1+theta_R)/(theta_R-omega_R))	
		
	// Define land use
		qui replace U = . 														// Clear any pre-existing value
		qui replace U = 3 if r_a > r_x_C & r_a > r_x_R							// Agrictultural land rent is highest
		qui replace U = 2 if r_x_R  > r_x_C & r_x_R > r_a						// Commercial land rent is highest
		qui replace U = 1 if r_x_C > r_x_R & r_x_C > r_a						// Residential land rent is highest
	
	// Outer edge of residential zone	
		qui sum x if U == 3 &  x >= 0 											// Take the inner margin of the agricultural zone
		qui scalar x1 = r(min)
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
		qui replace p_bar_x_C = a_x_C * 1/(1+omega_C)*S_x_C^omega_C				// According to Eq. (2)
		qui replace p_bar_x_R = a_x_R * 1/(1+omega_R)*S_x_R^omega_R				// According to Eq. (4)
	
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
		replace URBAN = U < 3													// Update indicator for urban use	
		replace COM = U == 1													// Update indicator for commercial use
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
					if L_hat_demand == 0 { // If no labour demand decrease wage by 10%
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

	
************************************************
*** Define programme for inverting amenities ***	
************************************************

********************************************************************************
// Define program to adjust fundamentals so to bring model heights closer to observed heights

capture program drop CONV 														// Drop any pre-existing program of the same name
	// Start defining program
	program CONV // Syntax: 1 Target population 2 convergence parameter 		
		foreach name in R C {													// Adjust for both land uses
			// If observed height is larger than realized height, we inflate amenities. This is done by using HEIGHT/S as an adjustment factor that inflates amenities a. 
			replace a_rand_`name' = (1-`1') *a_rand_`name' + (`1') * (HEIGHT_`name' / S_star_x_`name') *a_rand_`name'	 // We use the weighted combination of old amenities and new inflated amenities
			replace a_rand_`name' = 0 if a_rand_`name' == .						// If there is no height in data, we set amenity to a theory-consistent zero
		}
	FINDEQ 0.5 0.55 0.03 0.07 0.030 2 1 0.01 0.005 1.4 1.4 30 999 999			// Solve the model under updated fundamentals
	GHEIGHTB name0																// Graph updated spatial structure												
end
// Program ends
********************************************************************************

********************************************************************************
// Define program to adjust fundamentals so to bring model employment closer to target employment
capture program drop EMP
program EMP
foreach name in R C {
	replace a_rand_R = a_rand_R * (`1'/ sL)^0.01								// Adjust residential amenities to reach target population. If target larger than population in model, we inflate residential amenities to make the city more attractive
	FINDEQ 0.5 0.55 0.03 0.07 0.030 2 1 0.01 0.005 1.4 1.4 30 999 999			// Solve the model under updated fundamentals
}
end
// Program ends
********************************************************************************

********************************************************************************
// Iterative programme calling CONV until height has converged to values in data and population has converged to target
capture program drop INVERT 													// Drop any pre-existing program of the same name
program define INVERT
// Initialization	
	local PopGap = abs(sL-`1')													// Initial value of population gap 
	local O = 0 																// Initial value of correlation between obsered and model heights

// Keep iterating
	while `O' <0.999  {															// While correlation is below target
	qui CONV `2'																// Update amenities using program and convergence factor
	qui reg HEIGHT_R S_x_R														// Aux. regression to recover R2 of observed height and model height	
	qui local O =  e(r2) 														// Update correlation
	display("Correlation `O'")
	}
// Past this point we have converged
	while `PopGap' > 1000  {													// While gap is above tolerance
		qui EMP `1'																// Update amenities using program and employment target
		local PopGap = abs(sL-`1')												// Update gap
	display("Employment gap `PopGap'")	
	}
// Past this point we have converged
end
********************************************************************************

***********************************************************************
*** Define programmes for graphical illustration of model solutions ***
***********************************************************************


********************************************************************************
// Graph height
program GHEIGHT 
		sum S_x_C
		local top = (int(r(max)/25)+1)*25
		replace  SHADE = URBAN*COM*(int(`top'/1)+1)
		replace SHADEU = URBAN*(int(`top')+1)	
		qui sum S_x_C
		local max = round(r(max),25)
		local step = round(`max'/4, 5)
		twoway  (area SHADEU x, color(gs14)) (area SHADE x, color(gs12)) (line S_x_C x , color(red)  	   lpattern(longdash)) (line S_x_R x if x < 0, color(blue) lpattern(shortdash)) (line S_x_R x if 		x > 0, color(blue) lpattern(shortdash)), ///
		graphregion(color(white)) ylabel(0[`step']`max') xlabel(-50[10]50) xtitle("") ytitle("Building height") legend(order(3 "Commercial" 4 "Residential") cols(2) pos(6)) title("Building height") name("`1'", replace) xsize(5) ysize(5)
		end // Solver ends
********************************************************************************

********************************************************************************
// Graph bid rent
program GBIDRENT 
		*local urban_max = urban_max
		*local urban_min = urban_min
		*local com_max = com_max
		qui sum p_bar_x_C
		local top = (int(r(max)/5)+1)*5
		replace  SHADE = URBAN*COM*(int(`top'/1))
		replace SHADEU = URBAN*(int(`top'))	
		twoway  (area SHADEU x, color(gs14)) (area SHADE x, color(gs12)) (line p_bar_x_C x , color(red) lpattern(longdash)) (line p_bar_x_R x if x > 0, color(blue) lpattern(shortdash)) (line p_bar_x_R x if x < 0, color(blue) lpattern(shortdash))  ///
		,   graphregion(color(white)) xlabel(-50[10]50) xtitle("") ytitle("Floor space rent") legend(order(1 "Urban area" 2 "CBD area" 3 "Commercial" 4 "Residential") cols(4) pos(6)) title("Floor space rent") name(`1', replace) 
		end // Solver ends
********************************************************************************


********************************************************************************
// Graph land rent
program GLANDRENT
	qui sum r_x_C
	local top = (int(r(max)/500)+1)*500	
	qui replace SHADE = URBAN*COM*(int(`top')+1)
	qui replace SHADEU = URBAN*(int(`top')+1)
	twoway (area SHADEU x, color(gs14)) (area SHADE x, color(gs12)) (line r_x_C x , color(red)  lpattern(longdash)) (line r_x_R x , color(blue) lpattern(shortdash)) (line r_a x, color(black) lpattern(solid)), ///
	  graphregion(color(white))  xlabel(-50[10]50) xtitle("") ytitle("Land bid rent") legend(order(3 "Commercial" 4 "Residential" 5 "Agricultural" 1 "Urban area" 2 "CBD")  size(vsmall) cols(5) pos(6)) title("Land bid rent") name(`1', replace)
end

********************************************************************************
// Bar height
program GHEIGHTB 
		sum S_x
			gen bin = int(D)
			egen mS = mean(S_x), by(bin)
			egen sdS = sd(S_x), by(bin)
			gen cvS = sdS/mS
			sum cvS 		
			local CV = round(r(mean), 0.001)
			sum cvS if bin <= 5	
			local CVcbd = round(r(mean), 0.001)
			drop bin mS sdS cvS
		sum S_x_C
		local top = (int(r(max)/25)+1)*25
		replace  SHADE = URBAN*COM*(int(`top'/1)+1)
		replace SHADEU = URBAN*(int(`top')+1)	
		qui sum S_x_C
		local max = round(r(max),5)
		local step = round(`max'/4, 5)
		twoway  (bar S_x_C x if U == 1, color(red) lcolor(none) barwidth(0.1) bargap(1)	   lpattern(solid)) (bar S_x_R x if U == 2, barwidth(0.01) lcolor(none) color(blue) lpattern(solid)) (bar S_x_R x if 	U == 2, barwidth(0.01) color(blue) lcolor(none) lpattern(solid)), ///
		graphregion(color(white)) ylabel(0[`step']`max') xlabel(-50[10]50) xtitle("") ytitle("Building height") legend(order(1 "Commercial" 2 "Residential") cols(2) pos(6))  name("`1'", replace) xsize(10) ysize(5) note("Mean coefficent of variation within one-distance-unit bins = `CV'") // (area SHADEU x, color(gs14)) (area SHADE x, color(gs12)) 
		end // Solver ends
		
		
******************************
*** All programmes defined ***
******************************
