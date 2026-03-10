## Allen & Arkolakis (RES 2022) — Model Summary

### 1. Setup & Externalities

The economy has $N$ locations connected by a transportation network $\bar{\mathbf{T}} = [\bar{t}_{kl}]$ of free-flow link traversal costs. There are two model variants sharing identical mathematical structure:

- **Economic geography (trade):** Agents choose where to live/produce; goods are traded across locations.
- **Urban (commuting):** Agents choose residence and workplace; they commute.

In both variants, productivities and amenities depend on local population (eq 1 / eq 14):

$$A_i = \bar{A}_i L_i^\alpha, \quad u_i = \bar{u}_i L_i^\beta$$

where $\bar{A}_i, \bar{u}_i > 0$ are **fundamental** (exogenous) components, $\alpha$ governs agglomeration in production, and $\beta$ governs amenity externalities. Typically $\alpha > 0$ (agglomeration) and $\beta < 0$ (congestion disamenity).

### 2. Stochastic Route Choice & Transport Costs (Section 3)

Agents choose routes through the network subject to **Fréchet-distributed** idiosyncratic route preferences (shape parameter $\theta > 0$). This yields an analytical formula for bilateral transport costs via the **Leontief inverse** (eq 21):

$$\tau_{ij} = b_{ij}^{-1/\theta}, \quad \mathbf{B} = (\mathbf{I} - \mathbf{A})^{-1}$$

where $\mathbf{A} = [t_{kl}^{-\theta}]$ is the weighted adjacency matrix (with $a_{kk} = 0$). This converges when the spectral radius of $\mathbf{A}$ is less than 1. The key insight: **no Dijkstra needed** — the Fréchet assumption gives a closed-form matrix expression for transport costs across all origin-destination pairs simultaneously.

### 3. Gravity Equations (Sections 2.1, 2.2)

**Trade gravity** (economic geography, eq 6):

$$X_{ij} = \tau_{ij}^{-\theta} \frac{Y_i}{\Pi_i^{-\theta}} \frac{E_j}{P_j^{-\theta}}$$

**Commuting gravity** (urban, eq 16):

$$L_{ij} = \tau_{ij}^{-\theta} \frac{L_i^R}{\Pi_i^{-\theta}} \frac{L_j^F}{P_j^{-\theta}}$$

where $\Pi_i$ and $P_j$ are **outward** and **inward market access** indices:

- Economic geography (eqs 7–8): $\Pi_i^{-\theta} = \sum_j \tau_{ij}^{-\theta} E_j P_j^\theta$ and $P_j^{-\theta} = \sum_i \tau_{ij}^{-\theta} Y_i \Pi_i^\theta$
- Urban (eqs 17–18): $\Pi_i^{-\theta} = \sum_j \tau_{ij}^{-\theta} L_j^F P_j^\theta$ and $P_j^{-\theta} = \sum_i \tau_{ij}^{-\theta} L_i^R \Pi_i^\theta$

These also have closed-form expressions linking them to equilibrium variables (e.g., eq 17: $\Pi_i = u_i (L_i^R)^{-1/\theta} (\bar{L}/\bar{W}^\theta)^{1/(2\theta)}$).

### 4. Traffic Gravity & Congestion (Section 3.2–3.3)

**Link intensity** (eq 23) — the expected number of times link $(k,l)$ is used for origin-destination pair $(i,j)$:

$$\pi_{ij}^{kl} = \left(\frac{\tau_{ij}}{\tau_{ik} \cdot t_{kl} \cdot \tau_{lj}}\right)^\theta$$

**Traffic gravity** (eq 24) — total traffic on link $(k,l)$:

$$\Xi_{kl} = t_{kl}^{-\theta} \times P_k^{-\theta} \times \Pi_l^{-\theta}$$

**Congestion** (eq 25) — link costs depend on traffic:

$$t_{kl} = \bar{t}_{kl} \cdot \Xi_{kl}^\lambda$$

where $\lambda \geq 0$ governs congestion strength. Combining eqs 24 and 25 yields the **congested link cost** (eq 26) and **congested traffic** (eq 27):

$$t_{kl} = \bar{t}_{kl}^{\rho} \cdot P_k^{-\theta\lambda\rho} \cdot \Pi_l^{-\theta\lambda\rho}$$

$$\Xi_{kl} = \bar{t}_{kl}^{-\theta\rho} \cdot P_k^{-\theta\rho} \cdot \Pi_l^{-\theta\rho}$$

where $\rho \equiv 1/(1 + \theta\lambda)$.

### 5. Equilibrium Without Congestion (Section 2)

**Economic geography** (eqs 10–11) — solve for income shares $y_i \equiv Y_i/Y^W$ and labor shares $l_i \equiv L_i/\bar{L}$:

$$\bar{A}_i^{-\theta} y_i^{1+\theta} l_i^{-\theta(1+\alpha)} = \chi \sum_j \tau_{ij}^{-\theta} \bar{u}_j^\theta y_j^{1+\theta} l_j^{\theta(\beta-1)}$$

$$\bar{u}_i^{-\theta} y_i^{-\theta} l_i^{\theta(1-\beta)} = \chi \sum_j \tau_{ji}^{-\theta} \bar{A}_j^\theta y_j^{-\theta} l_j^{\theta(\alpha+1)}$$

where $\chi \equiv (\bar{L}^{(\alpha+\beta)}/\bar{W})^\theta$ is an endogenous scalar (inverse welfare raised to $\theta$), determined by $\sum y_i = 1$, $\sum l_i = 1$.

**Urban** (eqs 19–20) — solve for residential shares $l_i^R$ and workplace shares $l_i^F$:

$$(l_i^R)^{1-\theta\beta} = \chi \sum_j \tau_{ij}^{-\theta} \bar{u}_i^\theta \bar{A}_j^\theta (l_j^F)^{\alpha\theta}$$

$$(l_i^F)^{1-\theta\alpha} = \chi \sum_j \tau_{ji}^{-\theta} \bar{u}_j^\theta \bar{A}_i^\theta (l_j^R)^{\beta\theta}$$

Note the structural similarity: both are $2N$ equations in $2N$ unknowns with the same functional form.

### 6. Equilibrium With Congestion (Section 4.1)

This is the paper's key contribution. By substituting the Leontief inverse (eq 21) into the equilibrium conditions, performing a matrix inversion, and incorporating congestion (eqs 24–27), the equilibrium conditions become self-contained $2N$ systems in $2N$ unknowns that only involve **direct links** $\bar{t}_{ij}$:

**Economic geography** (eqs 28–29):

$$y_i^{(1+\theta)\rho + \theta\lambda\rho} l_i^{-\theta(1+\alpha+(\alpha+\beta)\theta\lambda)\rho} = \chi \bar{u}_i^\theta \bar{A}_i^\theta \cdot y_i^{(1+\theta+\theta\lambda)\rho} l_i^{\theta(\beta-1)\rho} + \chi^{\theta\lambda\rho} \sum_{j} (\bar{L}^\lambda \bar{t}_{ij})^{-\theta\rho} \cdot [\text{fundamentals}]_i \cdot y_j^{(1+\theta)\rho} l_j^{-\theta(1+\alpha)\rho}$$

**Urban** (eqs 30–31):

$$(l_i^R)^{1-\theta\beta} (l_i^F)^{\theta\lambda(1-\alpha\theta)\rho} = \chi \bar{u}_i^\theta \bar{A}_i^\theta (l_i^F)^{\theta(\alpha+\lambda)\rho} + \chi^{\theta\lambda\rho} \sum_j (\bar{L}^\lambda \bar{t}_{ij})^{-\theta\rho} [\text{fund.}]_{ij} (l_j^R)^{(1-\theta\beta)\rho}$$

$$(l_i^R)^{\theta\lambda(1-\beta\theta)\rho} (l_i^F)^{1-\theta\alpha} = \chi \bar{u}_i^\theta \bar{A}_i^\theta (l_i^R)^{\theta(\beta+\lambda)\rho} + \chi^{\theta\lambda\rho} \sum_j (\bar{L}^\lambda \bar{t}_{ji})^{-\theta\rho} [\text{fund.}]_{ij} (l_j^F)^{(1-\alpha\theta)\rho}$$

The crucial structure: each equation splits into a **self-term** (location $i$ only — from the identity in the Leontief inverse) and a **sum over direct links** (from the adjacency matrix). The multi-hop paths are captured implicitly by solving the system self-consistently.

**Key property:** With $\lambda > 0$, equilibrium is **scale dependent** — it depends on $\bar{L}$ (aggregate labor). Increases in $\bar{L}$ are isomorphic to uniform increases in link costs with elasticity $\lambda$.

### 7. Existence & Uniqueness (Proposition 1, Section 4.2)

Existence is guaranteed for all parameter values. Uniqueness holds under:
- **Economic geography** (symmetric $\bar{\mathbf{T}}$): $\alpha + \beta \leq 0$ (eq 32)
- **Urban**: $\alpha \leq \frac{1}{2}(1/\theta - \lambda)$ and $\beta \leq \frac{1}{2}(1/\theta - \lambda)$ (eq 33)

### 8. Traffic-to-Trade/Commuting Prediction (Section 5.1)

An analytical mapping from observed traffic to bilateral flows:

**Trade** (eq 34): $X_{ij} = c_{ij}^X \times Y_i \times E_j$, where $\mathbf{C}^X = (\mathbf{D}^X - \Xi)^{-1}$

**Commuting** (eq 35): $L_{ij} = c_{ij}^L \times L_i^R \times L_j^F$, where $\mathbf{C}^L = (\mathbf{D}^L - \Xi)^{-1}$

with $\mathbf{D}$ a diagonal matrix combining local economic activity and traffic sums. This requires **no model parameters** — only observed traffic and economic activity.

### 9. Exact-Hat Counterfactuals (Proposition 2, Section 5.2)

Given an infrastructure change $\hat{\bar{t}}_{kl}$, the counterfactual system mirrors the level system but uses **observed data as weights**:

**Economic geography** (eqs 36–37): Solve for $(\hat{y}_i, \hat{l}_i, \hat{\chi})$. The weight shares are:
- Outward: $s_{ij}^{out} = \Xi_{ij} / (E_i + \sum_k \Xi_{ik})$, self: $s_i^{own,out} = E_i / (E_i + \sum_k \Xi_{ik})$
- Inward: $s_{ji}^{in} = \Xi_{ji} / (Y_i + \sum_k \Xi_{ki})$, self: $s_i^{own,in} = Y_i / (Y_i + \sum_k \Xi_{ki})$

**Urban** (eqs 38–39): Solve for $(\hat{l}_i^R, \hat{l}_i^F, \hat{\chi})$. Weight shares use $L_i^F$ and $L_i^R$ instead of $E_i$ and $Y_i$.

The exponents in the hat system are **identical** to those in the level system (eqs 28–31). Welfare change: $\hat{W} = \hat{\chi}^{-1/\theta}$.

### 10. Congestion Estimation (Section 5.3)

The congestion parameter decomposes as $\lambda = \delta_0 \cdot \delta_1$ where:
- $\delta_0 = 1/\theta$ is the time elasticity of transport cost (eq 40)
- $\delta_1$ is the congestion elasticity of inverse speed (eq 41), estimated via IV

This yields $\hat{\bar{t}}_{kl} = \widehat{lanes}_{kl}^{-\lambda}$ (eq 42): adding lanes reduces the exogenous cost component.


---

The Fréchet distribution is the linchpin that makes the entire framework tractable. Here's why:

## The Core Problem

Agents choose routes through a network. A route from $i$ to $j$ of length $K$ has cost $\prod_{l=1}^{K} t_{r_{l-1}, r_l}$. In principle, you'd need to solve a shortest-path problem (Dijkstra) for every origin-destination pair — a computational object, not an analytical one. And when transport costs are endogenous (through congestion), you'd need to re-solve all shortest paths at every iteration of the equilibrium solver.

## What Fréchet Buys You

Each agent $\nu$ draws an idiosyncratic route-specific shock $\varepsilon_{ij,r}(\nu)$ from a Fréchet distribution with shape parameter $\theta$. This means agents don't all take the cheapest route — they spread across routes probabilistically, with $\theta$ governing how concentrated choices are on the least-cost route.

**The key mathematical payoff:** The probability of using route $r$ involves products of link costs raised to $-\theta$. When you sum over all possible routes of all lengths, this becomes a geometric matrix series:

$$\sum_{K=0}^{\infty} \mathbf{A}^K = (\mathbf{I} - \mathbf{A})^{-1} = \mathbf{B}$$

where $\mathbf{A} = [t_{kl}^{-\theta}]$. So $\tau_{ij} = b_{ij}^{-1/\theta}$ — a single matrix inversion replaces all shortest-path computations.

This works because the Fréchet is a **max-stable** distribution: the maximum of Fréchet random variables is itself Fréchet. Since agents maximize utility across routes, and each route's utility involves multiplicative link-level shocks, the Fréchet's multiplicative/max-stable structure means the aggregation across routes of all lengths telescopes into the Leontief inverse.

## Five Concrete Consequences

1. **Analytical transport costs (eq 21):** $\tau_{ij} = b_{ij}^{-1/\theta}$. No numerical shortest-path computation needed.

2. **Analytical link intensity (eq 23):** $\pi_{ij}^{kl} = (\tau_{ij}/(\tau_{ik} \cdot t_{kl} \cdot \tau_{lj}))^\theta$. The expected usage of any link for any OD pair is a closed-form ratio of transport costs.

3. **Traffic gravity (eq 24):** $\Xi_{kl} = t_{kl}^{-\theta} P_k^{-\theta} \Pi_l^{-\theta}$. Traffic on a link has a gravity structure with the *same* market access terms as trade/commuting — this is what enables the traffic-to-trade mapping (eqs 34–35).

4. **Dimensionality preservation:** With congestion, the feedback loop (economic activity → traffic → congestion → transport costs → economic activity) could in principle blow up the dimensionality. But because Fréchet gives closed forms at every step, the congested equilibrium (eqs 28–31) remains a $2N$-equation system in $2N$ unknowns — same dimensionality as the standard model without congestion.

5. **Exact-hat algebra with traffic data (eqs 36–39):** The counterfactual system can be written purely in terms of observed traffic flows and economic activity. This is possible because the Fréchet structure creates a tight analytical link between traffic on links and bilateral flows between OD pairs — without it, you couldn't substitute traffic data for the unobserved network.

## The $\theta$ Parameter

$\theta$ does double duty:

- **Dispersion of route choice:** Low $\theta$ = agents spread across many routes (more "noise"); high $\theta$ = agents concentrate on the cheapest route. As $\theta \to \infty$, the model converges to deterministic shortest-path routing (Dijkstra).
- **Trade/commuting elasticity:** The same $\theta$ governs the elasticity of bilateral trade or commuting flows to transport costs. This is a restriction — route choice heterogeneity and origin-destination choice heterogeneity share the same parameter — but it's what makes the traffic gravity equation (eq 24) use the *same* market access terms as the trade/commuting gravity equations (eqs 6, 16), which is essential for the entire framework to close.

In short: the Fréchet turns a computationally intractable network routing problem into linear algebra (matrix inversion), while simultaneously generating the gravity structure needed to connect traffic data to trade theory.

---

## Data Requirements for Calibration (Trade/Economic Geography Version)

### Parameters (4 structural):
| Parameter | Meaning | Typical source |
|-----------|---------|----------------|
| $\theta$ | Trade elasticity / Fréchet shape | Gravity estimation or literature (paper uses $\theta = 8$) |
| $\alpha$ | Productivity agglomeration | Literature or structural estimation ($\alpha = 0.1$) |
| $\beta$ | Amenity externality | Literature or structural estimation ($\beta = -0.3$) |
| $\lambda$ | Congestion strength | IV regression of inverse speed on traffic/lanes ($\lambda = 0.092$) |

### Data for the levels equilibrium (eqs 28–29):
1. **Network topology $\bar{\mathbf{T}}$**: Free-flow traversal costs $\bar{t}_{kl}$ for each link. Derived from road network data: distances, free-flow speeds, and lane counts. Non-connected pairs have $\bar{t}_{kl} = \infty$.
2. **Fundamentals $\bar{A}_i, \bar{u}_i$**: Can be recovered by inverting the equilibrium conditions given observed $y_i, l_i$ and transport costs — standard inversion as in Allen & Arkolakis (2014).
3. **Aggregate labor $\bar{L}$**: Total population/workforce. Only matters when $\lambda > 0$.

### Data for exact-hat counterfactuals (eqs 36–37) — the more practical approach:
1. **Traffic flows $\Xi_{kl}$**: Observed traffic on each link in the network (e.g., from highway traffic counts — AADT data in the US).
2. **Income/expenditure $Y_i, E_i$**: GDP or total income at each node. With balanced trade: $Y_i = E_i$.
3. **Infrastructure change $\hat{\bar{t}}_{kl}$**: The policy counterfactual (e.g., adding lanes: $\hat{\bar{t}}_{kl} = \widehat{lanes}_{kl}^{-\lambda}$).

The major advantage of the exact-hat approach is that it does **not** require:
- Bilateral trade flows $X_{ij}$ (hard to observe)
- The full network topology $\bar{\mathbf{T}}$
- Fundamentals $\bar{A}_i, \bar{u}_i$

It only needs **traffic counts on links** plus **aggregate economic activity at nodes** — both of which are widely available (e.g., US Highway Performance Monitoring System for traffic, BEA for county GDP).

### For estimating $\lambda$ specifically (eq 41):
- Link-level **traffic counts** (vehicle-miles traveled)
- Link-level **speed data** (e.g., from probe vehicle data or speed sensors)
- Link-level **lane counts** (as instrument for capacity)
- Instruments for traffic: the paper uses distance to historical routes (e.g., 1947 highway plan) interacted with current economic activity as IVs, following Duranton & Turner (2011).