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
	u "DATA\SIMULATION\BASE.dta", clear

******************************	
*** Merge observed heights ***	
******************************	
	
// Merge actual height and land value data
	gen X = x
	merge 1:1 X using "DATA\EMPIRICAL\CH_skyline.dta"
	
**********************	
*** Initialization ***	
**********************	
	
// Solve the model under the baseline parametrization	
// Syntax 	theta_C	theta_R	omega_C	omega_R	beta_C	a_bar_C	a_bar_R	tau_C	tau_R	c_C		c_R 	r_a 	S_bar_C	S_bar_R
	FINDEQ 		0.5 	0.55 	0.03 	0.07 	0.030 	2		1 		0.01 	0.005 	1.3 	1.3 	50		999		999
	GHEIGHTB name0																// Use porgramme to illustrate spatial structure
	
******************************	
*** Generate location bins ***	
******************************		
	
// Assign heights to use that with greater height in data
	replace HEIGHT_R = . if HEIGHT_C >= HEIGHT_R &  HEIGHT_C != . & HEIGHT_R != .
	replace HEIGHT_C = . if HEIGHT_R > HEIGHT_C &  HEIGHT_C != . & HEIGHT_R != .

// Generate location bins of 0.1 width and assign max heights within these bins
// Impoves visualization and speeds up convergence	
	gen CONVBIN = round(X,0.1)													// Generate bins
	egen BHEIGHT_C = max(HEIGHT_C), by(CONVBIN)									// Generate max height within commercial bins
	replace HEIGHT_C = BHEIGHT_C if HEIGHT_R == . 								// Assign max height
	egen BHEIGHT_R = max(HEIGHT_R), by(CONVBIN)									// Generate max height within residential bins
	replace HEIGHT_R = BHEIGHT_R if HEIGHT_C == . 								// Assign max height
	
***********************************************
*** Define programme for updating amenities ***	
***********************************************
/*	
// Define program to adjust fundamentals so to bring model heights closer to observed heights
capture program drop CONV 														// Drop any pre-existing program of the same name
	// Start defining program
	program CONV // Syntax: 1 Target population 2 convergence parameter 		
		foreach name in R C {													// Adjust for both land uses
			// If observed height is larger than realized height, we inflate amenities. This is done by using HEIGHT/S as an adjustment factor that inflates amenities a. 
			replace a_rand_`name' = (1-`2') *a_rand_`name' + (`2') * (HEIGHT_`name' / S_star_x_`name') *a_rand_`name'	 // We use the weighted combination of old amenities and new inflated amenities
			replace a_rand_`name' = 0 if a_rand_`name' == .						// If there is no height in data, we set amenity to a theory-consistent zero
		}
	FINDEQ 0.5 0.55 0.03 0.07 0.030 2 1 0.01 0.005 1.4 1.4 30 999 999			// Solve the model under updated fundamentals
	GHEIGHTB name0																// Graph updated spatial structure												
	replace a_rand_R = a_rand_R * (`1'/ sL)^0.05								// Adjust residential amenities to reach target population. If target larger than population in model, we inflate residential amenities to make the city more attractive
end
// Program ends

capture program drop INVERT 														// Drop any pre-existing program of the same name
program define INVERT
// Initialization	
	local PopGap = abs(sL-`1')											// Initial value of population gap 
	local O = 0 																// Initial value of correlation between obsered and model heights

// Keep iterating
	while `O' <0.999 | `PopGap'  > 1000 {										// While gap is above tolerance and correlation is below target
	qui CONV `1' `2'													// Update amenities using program and convergence factor
	qui reg HEIGHT_R S_x_R														// Aux. regression to recover R2 of observed height and model height	
	qui local O =  e(r2) 														// Update correlation
	qui local PopGap = abs(sL-`1')												// Update gap
	display "Height correlation: `O', population gap = `PopGap'"				// Displace correlation and gap to user
	}
// Past this point we have converged
end
*/

***************************
*** Invert fundamentals ***	
***************************

	INVERT 1000000 0.05															// Use programme to find fundamentals that rationalize fuzzy height gradient and ensure a population in tall buildings of 1000000 (convergence at convergence parameter 0.05)

************************************
*** Iteratively update amenities ***	
************************************
/*
// set the target population: Select your value
	global emptarget = 1000000

// Initialization	
	local PopGap = abs(sL-$emptarget)											// Initial value of population gap 
	local O = 0 																// Initial value of correlation between obsered and model heights

// Keep iterating
	while `O' <0.999 | `PopGap'  > 1000 {										// While gap is above tolerance and correlation is below target
	qui CONV $emptarget 0.05													// Update amenities using program and convergence factor
	qui reg HEIGHT_R S_x_R														// Aux. regression to recover R2 of observed height and model height	
	qui local O =  e(r2) 														// Update correlation
	qui local PopGap = abs(sL-$emptarget)										// Update gap
	display "Height correlation: `O', population gap = `PopGap'"				// Displace correlation and gap to user
	}
// Past this point we have converged
*/
***************************
*** Inspect the outcome ***	
***************************
	
// Discretize space, aggregate height variables to location bins within 20 km
	drop if D > 20
	egen 	 BS_x = max(S_x), by(CONVBIN)
	egen 	 BS_x_R = max(S_x_R), by(CONVBIN)
	egen 	 BS_x_C = max(S_x_C), by(CONVBIN)
// Plot the skyline
	qui sum BS_x_C																// Recover max height for y-label
	local max = round(r(max),25)
	twoway  (bar BS_x_C x , color(red%50) lcolor(none) barwidth(0.1) bargap(1)	   lpattern(solid)) (bar BS_x_R x , barwidth(0.01) lcolor(none) color(blue%50) lpattern(solid)) (bar BS_x_R x , barwidth(0.01) color(blue%50) lcolor(none) lpattern(solid)), ///
	graphregion(color(white)) ylabel(0[25]`max') xlabel(-20[10]20) xtitle("Distance from centre (km)") ytitle("Building height") legend(order(1 "Commercial" 2 "Residential") cols(2) pos(6))   xsize(10) ysize(5) name(skyline, replace)
		
// now inspect residential floor space price gradient
	GBIDRENT		
		

// Compute coefficient of variation of floor space rent
	sum p_*	
* compute CV of land rent
	gen r_x = max(r_x_C, r_x_R)
	sum r_x
	display r(sd) / r(mean)
	
// Plot the distribution of production and residential amenities	
	twoway (hist a_rand_C if a_rand_C > 0 , color(red%50))(hist a_rand_R if a_rand_R > 0 , color(blue%50))  , legend(order(1 "Commercial" 2 "Residential") cols(2) pos(6)) xtitle("Fundamental amenity") ysize(5) xsize(5) name(amenity, replace)
	sum a_rand_C if a_rand_C > 0 , d
	sum a_rand_R if a_rand_R > 0 , d

// Land value variables
	gen lr_x = ln(r_x)															// Gen log land values
	egen Bllv2000 = max(llv1990), by(CONVBIN)									// Generate log land values in data for locaiton bins
	egen Blr_x = max(lr_x), by(CONVBIN)											// Generate log land values in model for location bins
	foreach var of varlist Bllv2000 Blr_x {
		sum `var'
		replace `var' = `var'-r(mean)
	}

// Plot land rent
	GLANDRENT

// Identify land use
	egen TS_x = sum(S_x), by(CONVBIN)											// Total height by location bin
	egen TS_x_C = sum(S_x_C), by(CONVBIN)										// Commercial height by location bin
	gen comshare_CONVBIN = TS_x_C / TS_x 										// Comercial height share
	gen res = comshare_CONVBIN == 0												// Identify residential bin where commercial height is absent

// Scatter plot of land values in model and data
	twoway 	(scatter Bllv2000 Blr_x if res == 0, mcolor(red%50))  /// 
			(scatter Bllv2000 Blr_x if res == 1, mcolor(blue%50)) ///
			(function y = x , range(-1.25 3) color(black)) ///
			, xlabel(-1[1]3) legend(order(1 "Commercial" 2 "Residential") cols(2) pos(6)) ysize(5) xsize(5) ytitle("Ln land value in data") xtitle("Ln land value in model") name(lv, replace)

// Compile final figure and export
	grc1leg amenity lv, name(tech, replace) scale(1.5)
	grc1leg skyline tech, cols(1) scale(0.8)
	graph export "FIGS\FIG_SkylineCHmodel.png", width(2400) height(1800) replace 	


// Save output
	save "DATA\SIMULATION\INVERTED.dta", replace
	
**************************************
*** Do file successfully completed ***	
**************************************