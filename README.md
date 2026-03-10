# Quantitative Spatial Economics Repository

This repository contains Julia and MATLAB implementations of quantitative models from the spatial economics literature. The code is organized by research paper, with each folder containing implementations of specific theoretical frameworks.

## Contents

### 1. [AllenArkolakis-RES-2022](./AllenArkolakis-RES-2022/)
**The Welfare Effects of Transportation Infrastructure Improvements**

- **Authors:** Treb Allen and Costas Arkolakis
- **Journal:** The Review of Economic Studies
- **Year:** 2022
- **Description:** Develops a quantitative general equilibrium spatial framework featuring endogenous transportation costs and traffic congestion. The model yields analytical expressions for transportation costs, traffic flows, and the spatial distribution of economic activity. Applications to US highway and Seattle road networks demonstrate highly variable returns to infrastructure investment.
- **Key Features:** Route choice problem, traffic congestion via Leontief inverse, welfare analysis of infrastructure improvements

---

### 2. [FuchsFoongWong-MMN-2026](./FuchsFoongWong-MMN-2026/)
**Multimodal Transport Networks**

- **Authors:** Simon Fuchs and Woan Foong Wong
- **Status:** Working Paper
- **Year:** September 2025
- **Description:** Extends Allen & Arkolakis (2022) to multimodal transport systems incorporating roads, rail, and waterways. Models mode choice via nested CES aggregation and incorporates congestion at intermodal terminals. Calibrated to US freight network to evaluate terminal improvements and policy scenarios.
- **Key Features:** Modal substitution elasticity (η = 1.099), terminal congestion, multimodal Leontief inverse, environmental impact assessment

---

### 3. [QRE-HoRaUE-2025](./QRE-HoRaUE-2025/)
**Quantitative Regional Economics**

- **Authors:** Treb Allen and Costas Arkolakis
- **Journal:** Handbook Chapter
- **Year:** 2025
- **Description:** A comprehensive handbook chapter presenting the major advances in quantitative economic geography over the past decade. Provides a unified quantitative framework incorporating insights from seminal models and details how to combine spatial models with real-world data for estimation, counterfactuals, and welfare analysis.
- **Key Features:** Unified workhorse framework, empirical methodology, extensions and recent developments

---

### 4. [QSE-ARE-2017](./QSE-ARE-2017/)
**Quantitative Spatial Economics**

- **Authors:** Stephen J. Redding and Esteban Rossi-Hansberg
- **Journal:** Annual Review of Economics
- **Year:** 2017
- **Volume:** Vol. 9, pp. 21–58
- **Description:** A survey article synthesizing quantitative advances in spatial economics. Implementations include the Helpman (1998) model with Monte Carlo simulations for countries and regions.
- **Key Features:** Theoretical foundations, empirical applications, connections between spatial economics models

---

### 5. [Redding-JIE-2016](./Redding-JIE-2016/)
**Goods Trade, Factor Mobility and Welfare**

- **Author:** Stephen J. Redding
- **Journal:** Journal of International Economics
- **Year:** 2016
- **Volume:** Vol. 101, pp. 148–167
- **Description:** Develops a quantitative spatial model of trade with imperfect labor mobility. The gravity and topographic friction model (GTFM) combines insights from trade and regional economics to evaluate the welfare effects of infrastructure changes and labor mobility.
- **Key Features:** Gravity-based trade model, quantitative spatial equilibrium, infrastructure analysis tools

---

### 6. [Santamaria-2022](./Santamaria-2022/)
**Reshaping Infrastructure: Evidence from the Division of Germany**

- **Author:** Marta Santamaria
- **Journal:** Journal of the European Economic Association
- **Year:** 2026
- **Description:** Extends Redding (2016) and Allen-Arkolakis frameworks with a benevolent government choosing location-level infrastructure investments on a transport network. Transport costs depend on least-cost paths through the network, and infrastructure has spillover effects. Applies the model to evaluate the impact of German reunification on infrastructure and welfare.
- **Key Features:** Endogenous infrastructure investment, network effects, policy optimization

---

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

*Last updated: March 2026*
