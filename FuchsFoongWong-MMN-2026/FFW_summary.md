Here is the summary:

---

## Fuchs, Foong & Wong (2026) — "Multimodal Transport Networks": Technical Summary

### Model Overview

The model extends Allen & Arkolakis (RES 2022) from a single-mode road network to a multimodal transport system with roads (primary), rail, and waterways/barges (secondary modes). It nests AA2022 as a special case when the modal substitution elasticity equals the route elasticity (η = θ) or when only one mode exists.

### Transport Layer

**Route choice.** Shippers choose routes across a network graph where each link *k→l* has a stochastic cost drawn from a Fréchet distribution with shape θ = 8. This yields recursive bilateral transport costs via a Leontief inverse: τ_{ij} = B_{ij}^{−1/θ} where B = (I − A)^{−1} and A_{kl} = t_{kl}^{−θ} is the adjacency matrix (eq 1).

**Mode choice.** On each link, shippers choose among available modes via a nested CES aggregation (eq 5):

$$t_{ik}^{−η} = Σ_m t_{ik,m}^{−η}$$,  yielding mode choice probabilities $$π^m_{ik} = t_{ik,m}^{−η} / t_{ik}^{−η} (eq 6)$$

The key parameter η = 1.099 governs inter-modal substitution. Since η ≪ θ, mode choice is much less elastic than route choice — shippers can easily reroute within a mode but switching modes is costly.

**Secondary mode costs.** For non-road modes, link costs decompose as: t_{ik,m} = s_{ii,m} × ι_{ik,m} × s_{kk,m}, where s are terminal switching costs at intermodal nodes and ι are line-haul iceberg costs. This captures the fact that using rail or barge requires loading/unloading at terminals.

### Congestion

Two distinct congestion mechanisms operate:

1. **Primary network (road) congestion** (eq 11): t_{kl,1} = t̄_{kl,1} · Ξ_{kl,1}^{λ₁} with λ₁ = 0.092, where congestion acts at the link level, scaling free-flow costs by traffic volume.

2. **Terminal congestion** (eq 13): s_{kk,m} = s̄_{kk,m} · [Ξ_{kk,m}]^{λ_m} with λ_m = 0.096, where congestion acts at intermodal terminals — more traffic through a terminal increases dwell time and switching costs for all users. This is estimated using ship-level AIS dwell time data instrumented with a shift-share port trade exposure IV.

### Aggregate Equilibrium

The economic geography layer follows AA2022:
- **Productivity**: A_i = Ā_i · L_i^α (α = 0.1, agglomeration)
- **Amenities**: u_i = ū_i · L_i^β (β = −0.3, congestion)
- **Uniqueness**: α + β ≤ 0 guarantees a unique equilibrium
- **Welfare equalization**: W = w_i · u_i · P_i^{1/θ} (common across locations)

Market access is defined through the Leontief inverse:
- Inward (consumer prices): P_j^{−θ} = (B′ · c^{−θ})_j
- Outward (producer access): Π_i^{−θ} = (B · d)_i where d_j = E_j / P_j^{−θ}

### Traffic Flows

Aggregate and mode-specific traffic emerge as gravity equations:
- **Aggregate** (eq 8): Ξ_{kl} = t_{kl}^{−θ} · P_k^{−θ} · Π_l^{−θ}
- **Mode-specific** (eq 9): Ξ_{kl,m} = t_{kl,m}^{−η} · t_{kl}^{η−θ} · P_k^{−θ} · Π_l^{−θ}

The factor t_{kl}^{η−θ} captures how mode-specific traffic depends on overall link accessibility — a key new term relative to AA2022.

### Counterfactuals (Proposition 1)

Given observed traffic flows {Ξ_{ik}, Ξ_{ik,m}}, economic activity {Y_i, E_j}, and parameters {α, β, θ, λ₁, λ_m, η}, the model evaluates infrastructure changes through exact-hat algebra:
- **Outer loop**: aggregate equilibrium for {P̂_j, Π̂_i} (eq 21)
- **Inner loop**: transport equilibrium for {t̂_{ik}^{−θ}} (eq 22)

Welfare change: Ŵ = χ̂^{−1/θ}. Critically, exogenous fundamentals (Ā_i, ū_i) need not be recovered — the hat algebra cancels them out.

### Key Quantitative Results

The paper applies the model to the US network (228 nodes, 704 edges, 28 ports) and finds:
- **Terminal improvements**: Central terminals (Chicago, Houston, Dallas) yield the largest welfare gains; congestion dampens benefits by 2–3×
- **Rail strike** (20× cost increase): −0.4% welfare with congestion; congestion amplifies losses from losing a mode since displaced traffic worsens road congestion
- **Jones Act repeal** (÷2.7 barge costs): +5% welfare with congestion; congestion dampens gains since cheaper barge attracts traffic that congests terminals
- **Panama Canal disruption**: −0.6% welfare for affected East Coast ports
- **Environmental co-benefits**: Modal substitution from truck to rail/barge reduces GHG emissions (truck: 0.161 vs rail: 0.020 vs barge: 0.015 kg CO₂/ton-mile)
- **Model validation**: 0.75 correlation with CFS bilateral mode-specific trade flows; 0.97 correlation with observed multimodal freight shares by distance

### Parameter Summary

| Parameter | Value | Source |
|-----------|-------|--------|
| θ (route elasticity) | 8.0 | AA2022 literature |
| η (modal substitution) | 1.099 | IV estimation: highways → rail/road ratio |
| α (productivity spillover) | 0.1 | AA2022 |
| β (amenity spillover) | −0.3 | AA2022 |
| λ₁ (road congestion) | 0.092 | AA2022: speed on instrumented traffic |
| λ_m (terminal congestion) | 0.096 | AIS dwell time on shift-share IV port traffic |

---

## Data Requirements for Calibrating to an Arbitrary Country

### Minimum Requirements (Hat Algebra / Counterfactual Analysis)

If you only want to run counterfactuals using Proposition 1 (no need to solve for the full baseline equilibrium), you need:

1. **Multimodal transport network graph**
   - Nodes: cities, junctions, intermodal terminals, ports
   - Edges: mode-specific links (road, rail, waterway) with connectivity
   - Sources: OpenStreetMap (roads, rail), national transport agencies, port authorities

2. **Observed traffic flows on each mode**
   - **Road**: link-level traffic counts (vehicles/day or ton-km). Typically from highway agencies (e.g., AADT counts, toll data, or freight surveys)
   - **Rail**: bilateral rail freight flows (ton-km or carloads by origin-destination). From rail regulators, national freight surveys, or railway operators
   - **Waterway/barge**: port throughput and bilateral waterborne freight. From port authorities, customs data, or army corps equivalents
   - These give you Ξ_{ik} and Ξ_{ik,m} — the observed aggregate and mode-specific traffic matrices

3. **Economic activity by location**
   - Y_i (income/GDP) and E_j (expenditure) at each node. Under balanced trade: Y_i = E_i
   - Sources: national accounts, regional GDP, census income data, nightlights as proxy

4. **Six structural parameters** {θ, η, α, β, λ₁, λ_m}
   - Can be borrowed from the paper if local estimation is infeasible (θ = 8, η ≈ 1.1, α = 0.1, β = −0.3, λ₁ ≈ 0.09, λ_m ≈ 0.10)
   - For country-specific estimates, see below

### Additional Data for Full Baseline Equilibrium Solve

If you want to solve the equilibrium from scratch (as our implementation does with the synthetic network), you also need:

5. **Free-flow transport costs**
   - Road: link-level free-flow travel times or distances (from speed limits × distance, or Google Maps API)
   - Rail: per-km freight rates or iceberg costs by rail corridor
   - Waterway: per-km barge transport costs by waterway segment
   - Terminal switching costs: average dwell times or handling costs at intermodal facilities
   - These give you t̄_{kl,1}, ι_{kl,m}, and s̄_{kk,m}

6. **Location fundamentals**
   - Exogenous productivity Ā_i (can be recovered via the model's inversion: given observed {y_i, l_i} and the equilibrium conditions, back out Ā_i and ū_i)
   - Exogenous amenities ū_i (same inversion)

### Data for Country-Specific Parameter Estimation

7. **For η (modal substitution elasticity)**
   - Cross-sectional variation in road infrastructure across cities/regions (e.g., highway lane-km)
   - Corresponding rail traffic volumes (carloads, ton-km) and road traffic
   - Historical instruments: planned but unbuilt highways, historical rail routes, exploration/colonial routes
   - Minimum: city-level panel of modal traffic shares + road infrastructure variation

8. **For λ_m (terminal congestion elasticity)**
   - Ship-level or train-level dwell time data at ports/terminals
   - Port/terminal throughput volumes (ideally at daily or monthly frequency)
   - Instrument: trade exposure shift-share (import volumes at other ports weighted by historical port shares)
   - Minimum: port-level panel of throughput + average processing times

9. **For λ₁ (road congestion elasticity)**
   - Link-level speed data (from GPS traces, Google traffic, or loop detectors)
   - Link-level traffic counts
   - Instrument: road capacity changes (lane additions, new highways)
   - Can follow AA2022 methodology directly

10. **For θ (route/trade elasticity)**
    - Bilateral trade flows by origin-destination (gravity estimation)
    - Or: observed route choice frequencies + transport cost variation
    - Well-established in the literature; borrowing θ = 8 is standard

### For GHG / Environmental Analysis

11. **Mode-specific emission factors** (kg CO₂ per ton-km)
    - Country-specific factors from national environmental agencies
    - Or IPCC/IEA defaults: truck ≈ 0.06–0.16, rail ≈ 0.01–0.03, barge ≈ 0.01–0.02

### Practical Notes for Developing Countries

- **Data-scarce settings**: The hat-algebra approach (Proposition 1) is attractive because it requires only observed traffic flows and economic activity — not bilateral trade matrices or detailed routing data. Traffic counts are often available even in low-income countries (from road surveys, toll stations, or remote sensing).
- **Network construction**: OpenStreetMap provides global road and rail network geometry. Port data is available from UNCTAD and national maritime authorities.
- **Missing rail/waterway data**: The model nests the unimodal AA2022 case. You can start with road-only and add modes as data becomes available.
- **Spatial resolution**: The model is flexible in granularity — you can use districts, provinces, or cities as nodes. The key requirement is that traffic flow data matches the chosen resolution.