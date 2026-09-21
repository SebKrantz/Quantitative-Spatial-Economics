# Navigating Shocks: The Ripple Effects of Shipping Route Closures

**Carsten Philipp Brockhaus, Julian Hinz & Charles Serfaty** — Banque de France Working Paper No. 1057, August 2026.

*JEL:* F14, F17, R41, F62 · *Keywords:* shipping routes; chokepoint disruptions; Red Sea crisis; endogenous trade costs; monopoly tolls.

---

## 1. What the Paper Does

The paper asks **who bears the cost when a maritime chokepoint is disrupted**, and shows that the answer depends on two features that are invisible on a map:

1. **Does the chokepoint have a maritime substitute?** Suez and Panama have long but feasible detours; Hormuz places the Gulf economies behind a maritime dead end.
2. **Who prices the passage?** Suez and Panama are *priced infrastructures* whose authorities set revenue-maximizing tolls; Hormuz (currently) is not.

The central result is an **incidence reversal**. Because the Suez Canal Authority already prices to extract most of the surplus created by the shorter route, a permanent Red Sea closure is, in welfare terms, a shock to the *toll collector*, not to shippers: Egypt loses 3.0% of real income (2.76pp of it pure canal rent), while Germany, France, China, Japan, and the US each lose less than 0.02%. A transit fee on the Strait of Hormuz, a chokepoint without substitutes, does the opposite: it falls almost entirely on the Gulf economies themselves (Qatar −3.9%, Kuwait −2.7%, Iraq −2.5%), while importers substitute towards other suppliers at near-zero cost.

The paper combines (i) reduced-form evidence from the 2023–24 Red Sea crisis at three levels of aggregation with (ii) a multi-country, multi-sector Caliendo–Parro model in which the bilateral trade wedge is **endogenous to mode and route choice**, subject to route congestion, a global shipping capacity constraint, and monopoly toll-setting.

---

## 2. Empirical Evidence (Section 3)

### 2.1 Shipping (AIS)

6,032 IMO-registered container vessels, August 2021 – August 2024: 1.2m port calls at 1,593 ports across 30,647 directed port pairs. Port-call boundaries are built by density-based clustering of AIS positions; a voyage is a vessel's sequence of port calls (Appendix A).

| Network statistic | Before crisis | After crisis |
|---|---|---|
| Ports | 1,478 | 1,297 |
| Directed links | 26,413 | 17,495 |
| Density | 0.0121 | 0.0104 |
| Mean degree | 35.7 | 27.0 |
| Mean shortest path (legs) | 5.35 | 6.11 |

Red Sea transits collapsed by two thirds; Cape of Good Hope traffic rose; Panama and Taiwan Strait traffic were flat. Traffic concentrated on **fewer, longer corridors**, and the China–Netherlands travel-duration distribution shifted markedly right. At Hummels & Schaur's (2013) 0.6–2.1% of cargo value per day at sea, these delays are economically meaningful.

### 2.2 Trade (Census + Comext, monthly HS6)

Bilateral maritime distance is computed twice with the SeaRoute algorithm (Gaffuri & Eurostat 2022) — once with all passages open, once with Bab-el-Mandeb closed — and aggregated from ports to countries with distance-decayed IMF PortWatch weights (eq 1, $\alpha = 1$).

Treatment is the interaction $X_{ijs} = D_{ij} \times \pi^{sea}_{ijs}$ of the normalized distance shock $D_{ij} = (d^c_{ij} - d^o_{ij})/9000$ and the 2022 predicted sea share. The event study (eq 2) carries $\alpha_{ijs} + \alpha_{jst} + \alpha_{it}$ fixed effects — so the estimates are *differential*, with general-equilibrium spillovers largely absorbed.

**Result:** sea trade between Suez-dependent pairs fell >10% in December 2023 relative to September 2023 with no comparable air decline, no pre-trend, a single significant December trough, and a return to the pre-crisis path within about half a year.

### 2.3 Firms (Panjiva Turkish customs, weekly)

Turkey is informative because it sat *beside* the northern Suez approach: the Cape reroute moved the Asia–Europe trunk lane off its doorstep. Firm × partner × week event study (eq 3), via-Suez partners vs. westward partners:

| Outcome | Trough | End of window | Pooled DiD |
|---|---|---|---|
| Log maritime export value | −5.1% (wk 3), −5.4% (wk 4) | −0.2% (n.s.) | +0.005 (0.013), n.s. |
| Air share of export value | −1.6pp (wk 2) | **+1.4pp (wk 10)** | **+0.65pp (0.25), p = 0.011** |

**The two-speed story:** *how much* moves recovers within a quarter; *how it moves* shifts persistently toward air. The quantity margin is rerouting; the modal margin is a lasting change in relative transport margins. These are exactly the two margins the model endogenizes.

---

## 3. Model (Section 4)

### 3.1 The toy monopolist (Section 4.1)

A monopolist controls a short route $S$ (cost $d^S_{ij}$) competing with an untolled long alternative $L$ (cost $d^L_{ij}$), with route-substitution elasticity $\mu > 1$:

$$\pi_{ij}(F) = \frac{(d^S_{ij}+F)^{-\mu}}{(d^S_{ij}+F)^{-\mu} + (d^L_{ij})^{-\mu}}, \qquad T(F) = F \cdot \pi_{ij}(F) \cdot T_{ij} \tag{4}$$

With $d^S_{ij} = 0$ the optimal fee is $F^\star = d^L_{ij} (\mu-1)^{-1/\mu}$. As $\mu \to \infty$ the fee converges to the *full* cost difference between the routes — total travel costs are equalized across routes and the monopolist captures the entire surplus of the shortcut. **This is the paper's central mechanism:** when the shortcut's surplus is already levied as a rent, closing the shortcut destroys the rent, not the shippers' surplus.

### 3.2 Freight cost share and endogenous freight wedges (Section 4.2)

The freight cost share $\chi_s \in [0,1]$ is the share of delivered price accounted for by route-level shipping cost. It maps the freight wedge $f$ into the cargo-value iceberg $d$:

$$d_{ij,sr} = f_{ij,sr}^{\chi_s} \tag{5}$$

Route congestion and global capacity tightness are flow aggregates with ton-km capacity weights $\delta$:

$$\Xi_r \equiv \sum_{i,j,s} \delta_{ij,sr} X_{ij,sr}, \qquad \Psi \equiv \sum_{i,j,s,r} \delta_{ij,sr} X_{ij,sr} \tag{6, 7}$$

$$f_{ij,sr} = \bar f_{ij,sr}\cdot(\Xi_r/\bar\Xi_r)^{\lambda_s}\cdot(\Psi/\bar\Psi)^{\gamma_s}\cdot\phi_r
\;\Longrightarrow\;
d_{ij,sr} = \bar d_{ij,sr}\cdot(\Xi_r/\bar\Xi_r)^{\chi_s\lambda_s}\cdot(\Psi/\bar\Psi)^{\chi_s\gamma_s}\cdot\phi_r^{\chi_s} \tag{8, 9}$$

Congestion here accrues on **links** (the alternative passage gets more expensive as traffic reroutes), in contrast to the complementary node/port-terminal margin of Massoni (2025), which is held fixed.

### 3.3 Route choice, then mode choice (Sections 4.3, 4.5)

I.i.d. Gumbel shipper taste shocks over routes give the logit/CES aggregator. With the mode nest of Section 4.5, routes are chosen *within* a mode (elasticity $\eta^R_s$) and modes *across* composites (elasticity $\eta^M_s$):

$$d_{ij,sm} = \Big[\textstyle\sum_{r\in\mathcal M_{ijm}} d_{ij,smr}^{-\eta^R_s}\Big]^{-1/\eta^R_s}, \qquad \pi_{ij,smr} = \frac{d_{ij,smr}^{-\eta^R_s}}{\sum_\rho d_{ij,sm\rho}^{-\eta^R_s}} \tag{14, 15}$$

$$d_{ij,s} = \Big[\textstyle\sum_{m\in\mathcal N_{ij}} d_{ij,sm}^{-\eta^M_s}\Big]^{-1/\eta^M_s}, \qquad \mu_{ij,sm} = \frac{d_{ij,sm}^{-\eta^M_s}}{\sum_\nu d_{ij,s\nu}^{-\eta^M_s}} \tag{16, 17}$$

with route-level costs carrying route congestion, mode capacity, global capacity, and tolls:

$$d_{ij,smr} = \bar d_{ij,smr}\cdot(\Xi_{mr}/\bar\Xi_{mr})^{\chi_s\lambda^R_s}\cdot(\Psi_m/\bar\Psi_m)^{\chi_s\gamma^M_s}\cdot(\Psi/\bar\Psi)^{\chi_s\gamma^G_s}\cdot\phi_{mr}^{\chi_s} \tag{18}$$

Setting $|\mathcal N_{ij}| = 1$ recovers the route-only model of Section 4.3.

> **Note on the effective fee elasticity.** Route shares respond to the *iceberg* $d$ with elasticity $\eta^R_s$, but the toll enters $d$ only as $\phi^{\chi_s}$. The elasticity of route shares with respect to the **fee** is therefore $\eta^R_s \chi_s$ — with the calibrated $\eta^R = 350$ and a trade-weighted $\chi \approx 0.09$, roughly 31. This is why the calibrated route elasticity looks so large in levels.

### 3.4 Embedding in Caliendo–Parro (Section 4.4)

The route-mode composite $d_{ij,s}$ simply **replaces the exogenous bilateral trade wedge** in a multi-sector Eaton–Kortum model with input–output linkages:

$$P_{j,s} = \Big[\textstyle\sum_i T_{i,s}(w_i d_{ij,s})^{-\theta_s}\Big]^{-1/\theta_s}, \quad
\pi_{ij,s} = T_{i,s}(w_i d_{ij,s})^{-\theta_s} P_{j,s}^{\theta_s}, \quad
X_{ij,s} = \pi_{ij,s} E_{j,s}$$

Flows then cascade back down the nest, $X_{ij,sm} = \mu_{ij,sm} X_{ij,s}$ and $X_{ij,smr} = \pi_{ij,smr} X_{ij,sm}$, and feed the aggregates (19)–(20). The framework extends naturally to the general CES production network of Baqaee & Farhi (2024).

### 3.5 Monopoly toll setting (Section 4.4)

Each authority solves $\phi^\star_q \in \arg\max_{\phi_q \ge 1} \Pi_q(\phi_q;\phi_{-q})$ (eq 12), under one of two conventions:

- **Nash pricing** — a simultaneous fee-setting game solved to a fixed point.
- **Pricing in isolation** — each authority optimizes at the others' *reference* wedges. This is the first round of best responses and is what the paper uses to calibrate Suez and Panama.

Where several authorities toll the same route (e.g. Hormuz *then* Suez), wedges **stack multiplicatively**, $\Phi_q = \prod_{a\in A(q)}\phi_a$, and the joint take is divided in proportion to log-contribution:

$$\Pi_a = \sum_{q:\,a\in A(q)}\sum_{i,j,s} \frac{\ln\phi_a}{\ln\Phi_q}\,\frac{\Phi_q^{\chi_s}-1}{\Phi_q^{\chi_s}}\, X_{ij,sq} \tag{13}$$

which collapses to $\sum (\phi_a^{\chi_s}-1)/\phi_a^{\chi_s}\cdot X_{ij,sq}$ for a single authority and to the textbook $(\phi-1)/\phi \cdot \sum X$ when $\chi_s = 1$. Charging each authority's markup separately would over-extract. **Revenue is rebated to the collector country**, which is what makes the incidence result work: the toll is a wedge to shippers but income to Egypt/Panama/Iran/Oman.

> Everything in this subsection — (12), (13), the authorities, the rebate — is precisely what the **no-tolls variant** of §9 deletes. Sections 3.1–3.4 survive there unchanged.

### 3.6 Computation in changes

Exact hat algebra à la Dekle, Eaton & Kortum (2008). The extra block relative to a standard CP implementation is the transport fixed point:

$$\hat d_{ij,smr} = \hat{\bar d}_{ij,smr}\cdot \hat\Xi_{mr}^{\chi_s\lambda^R_s}\cdot\hat\Psi_m^{\chi_s\gamma^M_s}\cdot\hat\Psi^{\chi_s\gamma^G_s}\cdot\hat\phi_{mr}^{\chi_s}$$

together with (14)–(17). Under Nash pricing this nests inside an outer fixed point over $\{\hat\phi_q\}$; under pricing in isolation a single round of best responses replaces it.

Two implementation details matter:

1. **Levels vs. changes.** Toll wedges are levels, hat algebra works in changes. With a baseline embedding fee levels $\phi^{base}_q$, a counterfactual fee enters route costs as $\hat\phi_q = \phi_q/\phi^{base}_q$, while toll revenue (13) is evaluated at the **level** $\phi_q$ — so that re-solving the baseline at its own fees reproduces it exactly, with unchanged route costs and unchanged rents.
2. **Closure.** Nominal trade imbalances are held fixed at baseline values; **world value added is the numeraire**.

### 3.7 Additive alternative (Appendix B)

An economically tighter primitive treats cargo and freight as Leontief complements, $d_{ij,sr} = 1 + \chi_s(f_{ij,sr}-1)$ (eq 3'), with revenue $\Pi_q = \sum \frac{\chi_s(\phi_q-1)}{1+\chi_s(\phi_q-1)} X_{ij,sq}$ (eq 14'). The two coincide at the baseline and to first order; divergence is second-order in $f-1$ and proportional to $\chi(1-\chi)$. At the calibrated $\phi_{Suez}=1.59$ the two differ by <1pp for containerized sectors and <2pp for bulk, with identical sector rankings (Table 7). **All reported results use the multiplicative decomposition.**

---

## 4. Calibration (Section 5)

| Ingredient | Source | Value / detail |
|---|---|---|
| Baseline IO | **GTAP 11** | 160 countries × 65 sectors, reference year 2017 |
| Route shares $\pi_{ij,sr}$ | AIS trajectories + port calls | 4 passages: Suez, Panama, Cape, Direct; common across sectors within a pair; unobserved pairs → Direct |
| Mode shares $\mu_{ij,sm}$ | Predicted HS6 transport modes | sea / air / other; services → other (unexposed) |
| Hormuz exposure | Port geography | $g_i = 1$ (KWT, IRQ, QAT, BHR), 0.95 (IRN), 0.80 (SAU), 0.75 (ARE), 0 otherwise; $e_{ij} = g_i(1-g_j) + g_j(1-g_i)$ |
| $\theta_s$ (trade) | Fontagné et al. (2022) | $\theta_s = \sigma_s - 1$, long-run tariff-based, aggregated to GTAP |
| $\eta^M_s$ (mode) | Ko et al. (2025), IV | 0.39 (food n.e.c.) to 11.36 (pharma); pooled fallback **2.44** (Tolva 2026) |
| $\eta^R_s = \mu^d$ (route) | **Calibrated to Suez revenue** | **350** |
| $\lambda_s$ (route congestion) | $\alpha_\tau \times \beta_\tau$ | $0.38 \times 0.35 \approx$ **0.135** |
| $\gamma^G_s$ (global capacity) | Recursive VAR, impact elasticity | **0.44** (cumulative 0.7–0.9 at 1–3 months; CIs include zero) |
| $\gamma^M_s$ (mode capacity) | — | **0** in the quantification |
| $\chi_s$ (freight share) | Hummels (2007), UNCTAD/OECD | bulk 0.20–0.40; steel/chemicals 0.05–0.10; manufactures 0.02–0.04; electronics/pharma 0.01; trade-weighted mean ≈ **0.09** |

**Congestion elasticity.** $\alpha_\tau \approx 0.38$ from $\log\tau_{Suez,t} = \alpha_\tau\log\Xi_{Suez,t}+\varepsilon_t$ (excluding the Houthi and Ever Given episodes) measures how much longer a passage takes when traffic is higher. Only time-proportional costs (labor, ship utilization) scale with delay, not fuel or (un)loading, so $\alpha_\tau$ is scaled by $\beta_\tau \approx 0.35$, the long-run elasticity of freight cost to distance (Freightos, seven containerized routes, 2022–2025).

**Global capacity elasticity.** Freight prices (Drewry WCI) and a ton-km capacity-use index built from Census/Comext modality data are jointly determined by (21). A recursive VAR (22) with Cholesky ordering $(\log\Xi, \log P^G)$ imposes that within-month capacity use is predetermined, so the innovation to $\log\Xi$ is a pure demand shock and $\gamma_0 = b_{21}/b_{11}$ (23) identifies the inverse supply elasticity. Point estimates are hump-shaped in the horizon (0.44 on impact, 0.7–0.9 cumulative over 1–3 months) but imprecise — every interval contains zero.

**Toll calibration and its validation.** Toll revenues are not in the IO accounts, so Suez and Panama wedges are calibrated in a preliminary equilibrium where each authority sets its revenue-maximizing fee on the baseline data, **pricing in isolation**. $\mu^d$ is pinned by a single revenue moment:

| | Model | Observed |
|---|---|---|
| Suez toll revenue (**targeted**) | $10.5bn | $10.25bn (CY2023) |
| Suez wedge $\phi_{Suez}$ | 1.59 | — |
| Panama toll revenue (**untargeted check**) | $3.62bn | $3.18bn (FY2024, drought-depressed) |
| Panama wedge $\phi_{Panama}$ | 1.16 | — |

---

## 5. Counterfactuals (Section 6)

Solved in changes for 160 countries × 65 sectors. Suez and Panama wedges are held fixed; the Hormuz fee is exogenous.

**Scenarios:** (i) permanent **Red Sea closure** — the Suez route becomes unavailable for all pairs; (ii) **Strait of Hormuz transit fee** $\phi_{Hormuz} \in \{1.05, 1.35, 1.50\}$ on the freight bill of exposed cargo, revenue split equally between Iran and Oman; (iii) the **no-fee counterfactual world** in which Suez and Panama wedges are pinned at $\phi = 1$.

### 5.1 Welfare (% change in real income, Table 8)

| | Red Sea closure | Hormuz $\phi=1.05$ | $\phi=1.35$ | $\phi=1.50$ |
|---|---|---|---|---|
| **Gulf** | | | | |
| Qatar | −0.20 | −0.44 | **−3.86** | −5.42 |
| Kuwait | −0.04 | −0.26 | −2.67 | −3.84 |
| Iraq | −0.02 | −0.25 | −2.45 | −3.52 |
| Bahrain | −0.02 | −0.07 | −0.57 | −0.80 |
| UAE | −0.04 | −0.06 | −0.19 | −0.20 |
| Saudi Arabia | −0.02 | −0.05 | −0.06 | −0.04 |
| **Collectors** | | | | |
| Iran | +0.01 | +0.48 | +1.47 | +1.71 |
| Oman | −0.02 | +4.01 | **+12.25** | +14.24 |
| Egypt | **−3.02** | −0.06 | −0.19 | −0.23 |
| **Rival exporters** | | | | |
| Equatorial Guinea | +0.00 | +0.10 | +0.88 | +1.38 |
| Brunei | −0.02 | +0.11 | +0.82 | +1.14 |
| **Large traders** | | | | |
| Germany / France / Netherlands | −0.01 / −0.00 / −0.01 | ≈0 | −0.01 | −0.01 |
| US / China / Japan / Korea / India | ≈0 | ≈0 | −0.00 / −0.01 / −0.02 / −0.06 / −0.03 | — |

### 5.2 Route reallocation (% change in flow, Table 9)

| Route | Red Sea closure | $\phi=1.05$ | $\phi=1.35$ | $\phi=1.50$ |
|---|---|---|---|---|
| Suez | −100.0 | −0.6 | −2.4 | −2.9 |
| Cape of Good Hope | +18.4 | +0.1 | +0.3 | +0.2 |
| Panama | +19.8 | +0.5 | +1.7 | +1.9 |
| Hormuz transit | −6.6 | −19.0 | **−60.2** | −66.7 |

Hormuz toll revenue: $4.8bn → $14.6bn → $17.0bn. The throughput response is **strongly convex** and revenue flattens between the intermediate and high fee — at the fee under discussion the toll is already near the peak of its Laffer curve.

### 5.3 Rent vs. resource cost (Tables 10 and 6)

The welfare change decomposes into a direct rent component $\Delta\tau$ (change in the country's toll revenue as a share of income) and a residual GE/trade-cost channel:

| Red Sea closure | $\Delta W$ | $\Delta\tau$ | residual |
|---|---|---|---|
| Egypt | −3.02 | −2.76 | −0.27 |
| Panama | +0.15 | +0.19 | −0.04 |
| World | −0.02 | −0.01 | −0.01 |

Re-solving the *same* closure in a world with no chokepoint rents overturns the incidence:

| | With fees | Without fees |
|---|---|---|
| Egypt | −3.02 | **−0.18** |
| Qatar | −0.20 | **−1.09** |
| Haiti | −0.26 | −1.09 |
| Gabon | −0.22 | −0.64 |
| Togo | −0.11 | −0.33 |
| Malta | −0.05 | −0.22 |
| UAE | −0.04 | −0.21 |
| Singapore | −0.01 | −0.07 |

Global real income falls by a **similar amount either way**. The toll does not make the closure cheaper — it moves the incidence off the parties whose trade is disrupted and onto the single collector. *The rent does not insure the collector; it makes the collector the residual risk-bearer.*

### 5.4 Cross-chokepoint spillovers

Chokepoints interact because they are linked through the route structure. Gulf–Mediterranean cargo transits Hormuz *and then* Suez, so the Hormuz fee stacks on the Suez toll: at $\phi=1.35$ the Hormuz–Suez composite falls 54% ($33.9bn, ≈8% of Suez traffic), and Egypt loses ≈0.2% of real income — **an Iranian transit fee taxes the Egyptian canal at one remove**. Conversely, rising Cape traffic under the Hormuz fee is *not* Gulf cargo evading the strait (no such route exists) but the geographic footprint of supplier substitution towards West African hydrocarbon exporters.

### 5.5 Caveats stated by the authors

1. Hormuz exposure is **assigned from port geography**, not estimated; the shares for Saudi Arabia, the UAE, and Iran are judgement calls.
2. The bypass margin for partially exposed countries substitutes at $\mu$, likely **overstating short-run bypass capacity** and understating their losses.
3. The Hormuz fee is **exogenous**; canal fees do not re-optimize across scenarios. Numbers are the incidence of *given* fee levels, not the outcome of a fee-setting game.
4. Route shares and $\eta^R$ come from **container vessels**, whereas Hormuz trade is crude, refined products, and LNG — less elastic and more contract-bound, so the concentration of losses on the Gulf is a lower bound. Relatedly, high long-run substitution elasticities understate short-run importer losses.

---

## 6. Relation to the Other Models in This Repository

| Feature | This paper | AllenArkolakis-RES-2022 | FuchsFoongWong-MMN-2026 | Santamaria-2022 |
|---|---|---|---|---|
| Spatial unit | Countries (GTAP) | Locations on a network | Locations on a network | Locations on a network |
| Route choice | Logit over a **small set of named passages** (data-driven shares) | Fréchet over network links → **Leontief inverse** | Fréchet over links, recursive | **Dijkstra** least-cost paths |
| Mode choice | Nested CES sea/air/other, $\eta^M_s$ | — | Nested CES road/rail/barge, $\eta = 1.099$ | — |
| Congestion | Link (passage) + **global capacity** | Link, $\lambda = 0.07$–0.09 | Link + **terminal** | — |
| Price of passage | **Endogenous monopoly toll** | — | — | Government infrastructure choice |
| Production side | Multi-sector CP with IO | One sector, agglomeration | One sector, agglomeration | Krugman + mobile labor |
| Solution | Exact hat algebra | Exact hat algebra (Prop. 2) | Exact hat algebra (Prop. 1) | Levels + optimization |

The closest sibling is **FuchsFoongWong-MMN-2026**: both use a nested mode-over-route CES transport block with congestion, and both solve in changes. The distinctive contributions here are (a) **the price of passage is set by a revenue-maximizing agent** rather than given, (b) **multi-sector heterogeneity in $\chi_s$** drives differential exposure, and (c) chokepoints are **complements** (Hormuz then Suez), so tolls stack.

The **no-tolls variant** of §9 deletes the "price of passage" row and leaves every other column untouched. That places it squarely between the two traditions in this table: it has FuchsFoongWong's endogenous, congestible, multimodal transport block, but sitting on a multi-sector Caliendo–Parro production side rather than a one-sector economic geography. If you want the transport mechanism without the industrial-organisation layer, that is the file to read.

---

## 7. Contents of This Folder

| File | Description |
|---|---|
| `brockhaus_hinz_serfaty_chokepoint_model.jl` | **The full model** — self-contained Julia implementation (§8) |
| `graphs/` | 10 output figures (PDF) |
| `brockhaus_hinz_serfaty_no_tolls_model.jl` | **The no-tolls variant** — same CP system and transport block, priced-passage layer removed (§9) |
| `graphs_no_tolls/` | 6 output figures (PDF) |
| `Navigating Shocks _ The Ripple Effects of Shipping Route Closures_pdf.pdf` | Working paper (BdF WP 1057) |
| `MinerU_markdown_...md` | OCR'd markdown of the paper |
| `MinerU_latex_.../` | OCR'd LaTeX source + extracted figures |
| `MinerU_....json` | MinerU layout/parse output |
| `README.md` | This summary |
| `Claude_Plan.md` | Implementation plan for the Julia model |

**Known OCR artefacts** in the markdown/LaTeX (the PDF is authoritative):
- Section 4.4, "Computation in changes": `$\phi_q - s0$` should read "toll revenue (13) is evaluated at the level $\phi_q$ — so that re-solving the baseline at its own fees reproduces it exactly".
- Throughout, `Ð`/`±` are mangled en-dashes and `º`/`ª` mangled quotation marks.
- Equations (19)–(20) are merged into a single display block in the markdown.

---

## 8. The Julia Implementation — Full Model

> This folder carries **two** implementations. The full model, with the monopoly toll layer, is described here; the **no-tolls variant** on branch `bhs-no-tolls` is §9. They share the baseline builder, the transport block and the solver architecture, and differ only in whether passages are priced.

```bash
julia BrockhausHinzSerfaty-2026/brockhaus_hinz_serfaty_chokepoint_model.jl
```

One self-contained file (~1,750 lines, no external data, ~5 min, dependencies `LinearAlgebra, Statistics, Random, Printf, Plots`). `BHS_QUICK=1` stops after calibration and verification. Ten figures are written to `graphs/`.

**What it is.** A **stylized 24-country × 8-sector calibration** carrying the paper's real chokepoint geography (Suez, Panama, Cape, Direct, and the four Hormuz composites), its real Gulf-side exposure shares $g_i$, and its calibrated elasticities. GTAP 11, the AIS trajectories and the Panjiva micro-data are not redistributable, so it reproduces the paper's **mechanisms and incidence pattern**, not the magnitudes of Tables 8–10. Loader hooks for a real baseline are in `Claude_Plan.md` §11; nothing below the data-assembly layer would change.

**Baseline it builds:** world VA \$80.8tn, international trade \$14.8tn, seaborne trade \$8.1tn, of which \$1.13tn via Suez, \$648bn via Panama and \$757bn transiting Hormuz; trade-weighted $\bar\chi = 0.112$ (paper ≈ 0.09).

**Calibration it runs.** The full procedure of Section 5.1: each authority's revenue-maximizing fee is found by grid scan plus golden section, priced in isolation, and $\eta^R$ is bisected on the single Suez revenue moment. The equilibrium is then re-baselined at those fees, and the model verifies that re-solving it at its own fees reproduces it exactly.

| | Model | Paper | Observed |
|---|---|---|---|
| Suez toll revenue (**targeted**) | \$10.26bn | \$10.5bn | \$10.25bn |
| $\phi_{Suez}$ | 1.256 | 1.59 | — |
| Panama toll revenue (**untargeted**) | \$5.15bn | \$3.62bn | \$3.18bn |
| $\phi_{Panama}$ | 1.207 | 1.16 | — |
| $\eta^R = \mu^d$ | 140 | 350 | — |

$\eta^R$ differs because it is identified off *this* baseline's Suez traffic and freight-share mix, not GTAP's; what carries over is the effective fee elasticity $\eta^R\bar\chi$ (15.6 here, 31.5 in the paper) and the shape of the identification — the revenue curve peaks at $\phi^\star$ exactly on the \$10.25bn target line (figure 1).

**Results, against the paper.** Signs, orderings and the rent/total split come out right; levels are muted where the stylized baseline understates exposure (Gulf hydrocarbon exports are ~26% of income here against ~45% in reality).

| | Model | Paper |
|---|---|---|
| Red Sea closure, Egypt | **−4.60%** (rank 1 of 24) | −3.02% |
| — of which the rent channel | −4.24pp (**92%**) | −2.76pp (**91%**) |
| Red Sea closure, Egypt with no rents | −0.13% (rank 10) | −0.18% |
| Qatar, with fees → without fees | −0.09% → **−0.21%** | −0.20% → **−1.09%** |
| World cost, with fees / without | −0.025% / −0.026% | ≈ equal |
| Large traders (DEU, CHN, JPN, KOR) | −0.03% to −0.01% | < 0.02% |
| Suez / Cape / Panama flows | −100% / +41% / +9% | −100% / +18% / +20% |
| Hormuz fee 1.35: Gulf losses | −0.66% to −0.74% | −2.45% to −3.86% |
| Hormuz fee 1.35: Oman / Iran | **+7.18% / +0.75%** | +12.25% / +1.47% |
| Hormuz fee 1.35: Egypt | **−0.20%** | −0.19% |
| Hormuz transit flow at 1.05 / 1.35 / 1.50 | −9% / −48% / −57% | −19% / −60% / −67% |
| Hormuz–Suez composite at 1.35 | −41% | −54% |
| Nash vs. pricing in isolation | gap 0.002 in $\phi$ | "closely approximates" |

**Verification.** Every run prints a suite that must pass before the counterfactuals are believable:

- *Identity:* null shock returns all hats = 1 on both baselines (residuals ~1e-10); re-solving the calibrated baseline at its own fees leaves route costs and rents unchanged — the test that the levels-vs-changes handling of $\phi$ is right; market clearing, the income identity, the numeraire, and share normalizations.
- *Analytic:* Table 7 reproduced to 4 decimals under both freight decompositions; the toy monopolist's $F^\star = d^L(\mu-1)^{-1/\mu}$ matched numerically; eq 13's log-share division shown to sum to one and to under-extract relative to charging separately.
- *Limits:* $\chi = 0$ ⇒ a transit fee has exactly zero effect; congestion amplifies the delivered-cost response.
- *Substantive:* 14 checks on the paper's own claims, all passing — including the one that matters, that the with-fees and no-fees loser rankings **cross**.

**Three things worth knowing if you extend it.**

1. **Route shares are data, not model output.** With $\eta^R$ this large, the implied baseline cost gaps are a fraction of a percent, so generating shares from a distance logit would be meaningless. Only *changes* are model-driven — which is all hat algebra needs.
2. **A route closure's cost is $\chi$-independent.** Dropping a route from the CES (14) raises the composite by $\rho_{open}^{-1/\eta^R}$ whatever $\chi_s$ is, because the baseline gaps are inferred from the observed shares. $\chi$ scales the *toll* and *congestion* channels only. The natural-looking check "$\chi=0$ ⇒ a closure has no effect" is false; the correct one is that a *fee* has no effect.
3. **The congestion feedback needs adaptive damping.** A toll drives cargo off a passage, which decongests it, which makes it attractive again; the loop gain is $\eta^R\chi_s\lambda_s(1-\rho_r)$ and exceeds one for the high-freight-share sectors at the calibrated elasticities. Under a fixed step the transport fixed point oscillates and manufactures a spurious branch on which large fees look profitable — which showed up as a badly non-monotone calibration curve before it was fixed.

---

## 9. Variant: No Priced Passages (branch `bhs-no-tolls`)

```bash
julia BrockhausHinzSerfaty-2026/brockhaus_hinz_serfaty_no_tolls_model.jl
```

`brockhaus_hinz_serfaty_no_tolls_model.jl` keeps the Caliendo–Parro system and the whole transport block — freight-share pass-through, the two nests, route congestion, the global capacity constraint — and removes the toll layer entirely. No authorities, no $\phi_q$, no eq (12)–(13), no rent rebate, no re-baselining. Every passage is free; traversing one costs only what the link physically costs, which is still endogenous. ~1,250 lines, ~45 s, 6 figures in `graphs_no_tolls/`.

This is the paper's own "no chokepoint rents" world (Section 6.3) promoted from a robustness check to a standalone framework. Incidence lands on the economies whose cargo travels farther, not on a collector: **correlation(Suez exposure, welfare change) = −0.80**, against a near-zero relationship in the full model. Same plot, opposite finding — that contrast *is* Section 6.3.

**Experiments change character.** With nothing priced, the available shocks are closures and **exogenous cost shocks** — war-risk premia, drought restrictions — via `shock_freight(routes, fhat)`, which raises the freight wedge by `fhat` and the cargo iceberg by `fhat^χ_s`. A *fee* is not available: a fee is a toll.

### What removing the tolls costs you

$\eta^R$ is no longer identified. The full model pins it on a revenue moment — the revenue-maximizing Suez fee must reproduce observed canal income. With no fee there is no such moment, and $\eta^R$ is exactly the parameter the headline number scales against, since closing a route raises the within-mode composite by $(1-\rho_{closed})^{-1/\eta^R}$.

The script sweeps it rather than hiding it, and the result is more reassuring than the closed form alone suggests:

| | |$\eta^R$ = 50| 350 | 800 | $\lvert$loss$\rvert \times \eta^R$ spread |
|---|---|---|---|---|---|
| World loss | congestion **off** | −0.040 % | −0.0062 % | −0.0027 % | **8 %** — exactly $1/\eta^R$ |
| World loss | congestion **on** | −0.045 % | −0.0189 % | −0.0160 % | 473 % |

With congestion off the loss is proportional to $1/\eta^R$ to within 8 %, confirming the closed form. Switching congestion on puts a **floor** under the cost: the volume rerouted onto the Cape barely depends on $\eta^R$, so that component does not shrink as route choice becomes more elastic — at $\eta^R = 800$ it is over 80 % of the total. A 16× range in $\eta^R$ still moves the headline number by a factor of ~2.8.

### Verification specific to this variant

Removing the tolls exposes a clean analytic test that the full model cannot isolate. With congestion off and no exogenous shock, closing route $r_0$ must leave surviving shares exactly renormalized and raise the composite by exactly the closed form:

$$\rho'_r = \frac{\rho_r}{1-\rho_{r_0}}, \qquad \hat d_{ij,sm} = (1-\rho_{r_0})^{-1/\eta^R}$$

Both hold to 4e-16 — an exact check on the route nest (14)–(15). It also makes visible a property that is easy to state wrongly: **a route closure's cost is independent of $\chi_s$**, because the baseline cost gaps are inferred from observed shares. $\chi_s$ scales the congestion and premium channels only.

### What it cannot do

It answers how much a chokepoint disruption costs and who pays when nobody owns the passage. It cannot answer what the passage is worth to its owner, or reproduce the paper's central result that a priced chokepoint's disruption falls on the collector — those need the toll layer.

---

## 10. Key References

- **Caliendo & Parro (2015)**, *REStud* — the multi-sector IO trade framework being extended.
- **Dekle, Eaton & Kortum (2008)** — exact hat algebra.
- **Eaton & Kortum (2002)**, *Econometrica* — Ricardian gravity backbone.
- **Fuchs & Wong (2026)**; **Fuchs et al. (2026)** — multimodal / recursive routing with congestion (see `FuchsFoongWong-MMN-2026/`).
- **Allen & Arkolakis (2022)**, *REStud* — endogenous routing and congestion (see `AllenArkolakis-RES-2022/`).
- **Massoni (2025)** — the complementary *node* (port/terminal) congestion margin.
- **Dunn & Leibovici (2025)** — aggregate shipping-capacity dynamics applied to the Red Sea.
- **Ko et al. (2025)**; **Tolva (2026)** — mode-substitution elasticities.
- **Fontagné et al. (2022)** — sectoral trade elasticities.
- **Hummels & Schaur (2013)** — time in transit as a trade cost.
- **Gaffuri & Eurostat (2022)**; **Halili (2026)** — SeaRoute / `searoute-py` maritime routing.
