# Quantitative Spatial Economics Repository

This repository contains Julia and MATLAB implementations of quantitative models from the spatial economics literature. The code is organized by research paper, with each folder containing implementations of specific theoretical frameworks.

## Contents

### 1. [AhlfeldtBarr-JUE-2022](./AhlfeldtBarr-JUE-2022/)
**The Economics of Skyscrapers: A Synthesis**

- **Authors:** Gabriel M. Ahlfeldt and Jason Barr
- **Journal:** Journal of Urban Economics
- **Year:** 2022
- **Volume:** Vol. 129, 103419
- **Description:** A one-dimensional monocentric city in which developers choose building height against a convex construction cost, commercial and residential uses compete for land through bid rents, and a city-wide wage clears the labour market. The only model in this repository with an endogenous vertical margin. Ported to Julia from the authors' Stata toolkit.
- **Key Features:** Endogenous building height and FAR, bid-rent land-use allocation, perfectly open city, inversion against observed skylines (Chicago), height-limit and subcenter counterfactuals
- **Implementation:** `ahlfeldt_barr_skyscraper_model.jl` (~20 s). Validated against the toolkit's own solved output: all 19 solved columns of `INVERTED.dta` reproduce within a data-derived tolerance, land-use indicators match exactly, Chicago inversion R² = 0.99957

---

### 2. [AhlfeldtReddingSturmWolf-ECMA-2015](./AhlfeldtReddingSturmWolf-ECMA-2015/)
**The Economics of Density: Evidence from the Berlin Wall**

- **Authors:** Gabriel M. Ahlfeldt, Stephen J. Redding, Daniel M. Sturm and Nikolaus Wolf
- **Journal:** Econometrica
- **Year:** 2015
- **Volume:** Vol. 83(6), pp. 2127–2189
- **Description:** Workers with Fréchet-distributed idiosyncratic tastes choose a residence and a workplace separately, trading off wages against commuting costs and floor-space prices. Productivity and amenities decompose into exogenous fundamentals and density-dependent spillovers with distance decay. Ported to Julia from the author's MATLAB toolkit, covering all 16 codebook algorithms.
- **Key Features:** Bilateral commuting with Fréchet shocks, floor-space markets, sequential *and* simultaneous inversion of fundamentals, endogenous agglomeration and amenity spillovers, GMM estimation of the commuting elasticity
- **Implementation:** `ahlfeldt_redding_sturm_wolf_density_model.jl` (~19 s). **Runs on synthetic data** — the real Berlin inputs are large externally hosted files, and `load_berlin_data` is a documented stub. Validation exploits the paper's one-to-one mapping theorem instead: plant known fundamentals, solve forward, invert, recover them to ~1e-15 (66/66 checks)

---

### 3. [AllenArkolakis-RES-2022](./AllenArkolakis-RES-2022/)
**The Welfare Effects of Transportation Infrastructure Improvements**

- **Authors:** Treb Allen and Costas Arkolakis
- **Journal:** The Review of Economic Studies
- **Year:** 2022
- **Description:** Develops a quantitative general equilibrium spatial framework featuring endogenous transportation costs and traffic congestion. The model yields analytical expressions for transportation costs, traffic flows, and the spatial distribution of economic activity. Applications to US highway and Seattle road networks demonstrate highly variable returns to infrastructure investment.
- **Key Features:** Route choice problem, traffic congestion via Leontief inverse, welfare analysis of infrastructure improvements

---

### 4. [BrockhausHinzSerfaty-2026](./BrockhausHinzSerfaty-2026/)
**Navigating Shocks: The Ripple Effects of Shipping Route Closures**

- **Authors:** Carsten Philipp Brockhaus, Julian Hinz and Charles Serfaty
- **Series:** Banque de France Working Paper No. 1057
- **Year:** August 2026
- **Description:** Documents the trade effects of the 2023–24 Red Sea crisis using AIS trajectories, customs data and Turkish firm-level records, then builds a multi-country multi-sector Caliendo–Parro model in which the bilateral trade wedge is endogenous to mode and route choice. Freight costs respond to passage congestion and a global shipping capacity constraint, and canal authorities set revenue-maximizing tolls. A permanent Red Sea closure costs Egypt 3.0% of real income — almost all of it canal rent — while large trading economies lose under 0.02%; a Strait of Hormuz transit fee, a chokepoint without substitutes, instead concentrates losses on the Gulf.
- **Key Features:** Monopoly toll-setting at chokepoints, nested mode-over-route CES transport block, route congestion and global capacity constraints, sector-specific freight cost shares, exact-hat counterfactuals
- **Implementation:** `brockhaus_hinz_serfaty_chokepoint_model.jl` — self-contained 24-country × 8-sector stylized calibration with the paper's real chokepoint geography, monopoly fee calibration to observed Suez revenue, and both counterfactuals with and without chokepoint rents

---

### 5. [FuchsFoongWong-MMN-2026](./FuchsFoongWong-MMN-2026/)
**Multimodal Transport Networks**

- **Authors:** Simon Fuchs and Woan Foong Wong
- **Status:** Working Paper
- **Year:** September 2025
- **Description:** Extends Allen & Arkolakis (2022) to multimodal transport systems incorporating roads, rail, and waterways. Models mode choice via nested CES aggregation and incorporates congestion at intermodal terminals. Calibrated to US freight network to evaluate terminal improvements and policy scenarios.
- **Key Features:** Modal substitution elasticity (η = 1.099), terminal congestion, multimodal Leontief inverse, environmental impact assessment

---

### 6. [MonteReddingRossiHansberg-AER-2018](./MonteReddingRossiHansberg-AER-2018/)
**Commuting, Migration, and Local Employment Elasticities**

- **Authors:** Ferdinando Monte, Stephen J. Redding and Esteban Rossi-Hansberg
- **Journal:** American Economic Review
- **Year:** 2018
- **Volume:** Vol. 108(12), pp. 3855–3890
- **Description:** The elasticity of local employment to a labour demand shock depends on how open the local labour market is to commuting. Combines bilateral commuting with multi-region CES goods trade, so both factor and goods markets link locations. Ported to Julia from the Ahlfeldt–Seidel MATLAB toolkit, with a worked 401-county German calibration.
- **Key Features:** Bilateral commuting *and* goods trade, amenity/productivity inversion, exact-hat counterfactuals, plus a forward equilibrium solver the source toolkit lacks
- **Implementation:** `monte_redding_rossihansberg_commuting_model.jl` (~70 s). Two quantification tracks — the primary one needs only employment by residence and workplace, a wage or rent index, area and a bilateral cost matrix, with **no bilateral commuting flow matrix**, which is what makes the model calibratable to most countries. 21/21 checks; both structural gravity identities hold to machine precision

---

### 7. [QRE-HoRaUE-2025](./QRE-HoRaUE-2025/)
**Quantitative Regional Economics**

- **Authors:** Treb Allen and Costas Arkolakis
- **Journal:** Handbook Chapter
- **Year:** 2025
- **Description:** A comprehensive handbook chapter presenting the major advances in quantitative economic geography over the past decade. Provides a unified quantitative framework incorporating insights from seminal models and details how to combine spatial models with real-world data for estimation, counterfactuals, and welfare analysis.
- **Key Features:** Unified workhorse framework, empirical methodology, extensions and recent developments

---

### 8. [QSE-ARE-2017](./QSE-ARE-2017/)
**Quantitative Spatial Economics**

- **Authors:** Stephen J. Redding and Esteban Rossi-Hansberg
- **Journal:** Annual Review of Economics
- **Year:** 2017
- **Volume:** Vol. 9, pp. 21–58
- **Description:** A survey article synthesizing quantitative advances in spatial economics. Implementations include the Helpman (1998) model with Monte Carlo simulations for countries and regions.
- **Key Features:** Theoretical foundations, empirical applications, connections between spatial economics models

---

### 9. [Redding-JIE-2016](./Redding-JIE-2016/)
**Goods Trade, Factor Mobility and Welfare**

- **Author:** Stephen J. Redding
- **Journal:** Journal of International Economics
- **Year:** 2016
- **Volume:** Vol. 101, pp. 148–167
- **Description:** Develops a quantitative spatial model of trade with imperfect labor mobility. The gravity and topographic friction model (GTFM) combines insights from trade and regional economics to evaluate the welfare effects of infrastructure changes and labor mobility.
- **Key Features:** Gravity-based trade model, quantitative spatial equilibrium, infrastructure analysis tools

---

### 10. [Santamaria-2022](./Santamaria-2022/)
**Reshaping Infrastructure: Evidence from the Division of Germany**

- **Author:** Marta Santamaria
- **Journal:** Journal of the European Economic Association
- **Year:** 2026
- **Description:** Extends Redding (2016) and Allen-Arkolakis frameworks with a benevolent government choosing location-level infrastructure investments on a transport network. Transport costs depend on least-cost paths through the network, and infrastructure has spillover effects. Applies the model to evaluate the impact of German reunification on infrastructure and welfare.
- **Key Features:** Endogenous infrastructure investment, network effects, policy optimization

---

## Also in this repository

Not documented individually above, but present and (where noted) runnable:

- **[National CGE](./National%20CGE/)** — static Cameroon CGE (Dervis–de Melo–Robinson), Julia via JuMP + Ipopt
- **[DSGE_Toolkit](./DSGE_Toolkit/)** — Kleinman, Liu & Redding (2023) dynamic spatial GE with capital accumulation
- **[Optimal Transport Networks](./Optimal%20Transport%20Networks/)** — Fajgelbaum & Schaal (2020) optimal network design, MATLAB
- **[CPRS-USRegTrade-2018](./CPRS-USRegTrade-2018/)** and **[Trade Dynamics](./Trade%20Dynamics/)** — reference paper collections, no code

Each model folder carries its own `README.md` documenting equations, parameters and validation; the ported folders additionally carry a `TRANSLATION_PLAN.md`.

## Attribution for the ported models

Sections 1, 2 and 6 are Julia translations of MIT-licensed MATLAB and Stata toolkits by **Gabriel M. Ahlfeldt** (with Tobias Seidel for the MRRH2018 toolkit). Each folder vendors the original source, the upstream `LICENSE`, and a `TRANSLATION_PLAN.md` recording the equation-by-equation mapping and every decision taken where the source was ambiguous. **Read that plan before modifying the `.jl`** — several deliberate deviations from upstream are documented there and would otherwise look like bugs.

Cite the original papers, and the toolkits where their documentation asks for it.

## Repository Organization

Each folder contains:
- **Model implementations** in Julia (`.jl` files) and/or MATLAB (`.m` files)
- **Markdown summaries** explaining model structure and applications
- **Data files** and supporting materials
- **Graphs** and visualization outputs

## Usage

Each implementation follows self-contained conventions with parameters specified at the top of the file, detailed comments, and appropriate function organization. For Julia implementations, dependencies are listed at the beginning. For MATLAB code, required toolboxes are noted.

## Citation

When using code or implementations from this repository, please cite the relevant papers listed above according to your citation style.

---

*Last updated: September 2026*
