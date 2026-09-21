{smcl}
{hline}
help for {hi:AB2022}
{hline}

{title:Gradient Solver for "The economics of skyscrapers"}

{title:Description}

AB2022 is a gradient solver developed by Gabriel M Ahlfeldt, designed for the analysis presented in Ahlfeldt & Barr (2022): "The economics of skyscrapers". It is part of the Ahlfeldt-Barr-(2022)-toolkit that also contains a walkthrough replicating some of the results in the paper. This command generates floor space, height, and land rent gradients under the baseline parameterization. It also solves for the land use pattern (commercial vs. residential), total employment, and the wage that clears the labour market.

{title:Syntax}

AB2022 ["parameter=value"]

The following parameters can be adjusted by adding "parameter=value" 
as an argument. Below is a brief description and canonical parameter 
values from Ahlfeldt & Barr (2022):

theta_C 	Commercial construction cost elasticity of height   0.5
theta_R		Residential construction cost elasticity of height  0.55
omega_C 	Commercial rent elasticity of height                0.03
omega_R 	Residential rent elasticity of height               0.07
beta_C 		Agglomeration elasticity                            0.03
a_bar_C 	Fundamental production amenity                      2
a_bar_R 	Fundamental residential amenity                     1
tau_C 		Commercial amenity decay                            0.01
tau_R 		Residential amenity decay                           0.005
c_C 		Commerical baseline construction factor             1.4
c_R 		Residential baseline construction factor            1.4
r_a 		Agricultural land rents                             150
S_bar_C 	Commercial height limit                             999
S_bar_R		Residential height limit                            999
x_1_bar		Urban growth boundary (max. radius urban area)      999

The program will find the equilibrium using the above parameter values if you enter the commend without any argument. 
You can change any of the parameter values by adding "parameter=value" as an argument. You can add as many arguments as there are parameters.


{title:Example 1:} 
Baseline parameterization

AB2022

{title:Example 2}: 
Introduce a height limit of 50 floors for commercial buildings

AB2022 "S_bar_C=50"

{title:Example 3:} 
Introduce a height limit of 50 floors for commercial buildings and height limit of 10 floors for residential buildings
    
AB2022, "S_bar_C=50" "S_bar_R=10"

You can add as many arguments as there are parameters.

{title:Author}

Gabriel M Ahlfeldt, 02/2024
Humboldt University

{title:References}

The complete Ahlfeldt-Barr-(2022)-toolkit is available as a GitHub repository. {browse "https://github.com/Ahlfeldt/AB2022-toolkit":[link]}

The toolkit and the AB2022 ado file build on: 

Ahlfeldt, G., Barr, J. (2022): The Economics Skyscrapers: A synthesis. Journal of Urban Economics, 129. {browse "https://doi.org/10.1016/j.jue.2021.103419":[link]}

When using this programme or the toolkit in your work, please cite Ahlfeldt & Barr (2022).

{title:End of help file}
