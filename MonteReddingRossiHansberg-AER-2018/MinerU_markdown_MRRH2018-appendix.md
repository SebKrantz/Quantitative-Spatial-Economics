# Online Appendix for ìCommuting, Migration and Local Employment Elasticitiesî (Not for Publication)

Ferdinando Montey 

Georgetown University 

Stephen J. Reddingz 

Princeton University 

Esteban Rossi-Hansbergx 

Princeton University 

## A Introduction

Section B of this web appendix contains the proofs of the propositions in the paper, additional technical derivations of results reported in the paper, and further supplementary material for the quantitative analysis of the model. Section C includes additional empirical results and robustness tests. Section D presents further information about the data deÖnitions and sources. 

## B Quantitative Model Appendix

The Örst seven sections of this quantitative part of the web appendix present additional derivations for the main paper. Section B.1 reports the derivations of expected utility and the commuting probabilities. Section B.2 shows how the equilibrium conditions of the model can be used to undertake counterfactuals using the observed values of variables in the initial equilibrium. Section B.3 provides conditions for the existence and uniqueness of the general equilibrium. Section B.4 derives isomorphisms to other trade models with commuting and external economies of scale. Section B.5 shows that unobserved productivity can be uniquely determined from the observed variables and reports additional evidence on gravity in goods trade. Section B.6 shows that unobserved amenities can be uniquely recovered from the observed data and reports additional evidence on gravity in commuting. Section B.7 uses the commuter market clearing condition to show the relationship between di§erent measures of the openness of the local labor market to commuting. 

The remaining sections comprise supplementary material and extensions. Section B.8 reports the derivation of the partial equilibrium local employment elasticities discussed in the main paper. Section B.9 shows that the class of models consistent with a gravity equation for commuting áows implies heterogeneous local employment elasticities. Section B.10 introduces multiple worker types. Section B.11 introduces congestion in commuting. Section B.12 develops an extension of the baseline model to incorporate non traded consumption goods. Section B.13 considers the case where landlords use residential land. Section B.14 generalizes the production technology to incorporate intermediate inputs, commercial land use and capital. Section B.15 introduces heterogeneity in e§ective units of labor. Section B.16 considers the case where commuting costs are incurred in e§ective units of labor rather than in utility. Finally, Section B.17 considers a robustness test in which land is partially-owned locally and partially-owned by a national portfolio, where these ownership shares are chosen to rationalize measured trade deÖcits. 

## B.1 Commuting Decisions

We begin by reporting additional results for the characterization of worker commuting decisions. 

## B.1.1 Distribution of Utility

From all possible pairs of residence and employment locations, each worker chooses the bilateral commute that o§ers the maximum utility. Since the maximum of a sequence of FrÈchet distributed random variables is itself FrÈchet distributed, the distribution of utility across all possible pairs of residence and employment locations is: C 

$$
1 - G (u) = 1 - \prod_ {r = 1} ^ {S} \prod_ {s = 1} ^ {S} e ^ {- \Psi_ {r s} u ^ {- \epsilon}},
$$

where the left-hand side is the probability that a worker has a utility greater than u, and the right-hand side is one minus the probability that the worker has a utility less than u for all possible pairs of residence and employment locations. Therefore we have: 

$$
G (u) = e ^ {- \Phi u ^ {- \epsilon}}, \qquad \Psi = \sum_ {r = 1} ^ {S} \sum_ {s = 1} ^ {S} \Psi_ {r s}.\tag{B.1}
$$

Given this FrÈchet distribution for utility, expected utility is: 

$$
\mathbb {E} [ u ] = \int_ {0} ^ {\infty} \epsilon \Psi u ^ {- \epsilon} e ^ {- \Psi u ^ {- \epsilon}} d u.\tag{B.2}
$$

Now deÖne the following change of variables: 

$$
y = \Phi u ^ {- \epsilon}, \qquad d y = - \epsilon \Psi u ^ {- (\epsilon + 1)} d u.\tag{B.3}
$$

Using this change of variables, expected utility can be written as: 

$$
\mathbb {E} [ u ] = \int_ {0} ^ {\infty} \Psi^ {1 / \epsilon} y ^ {- 1 / \epsilon} e ^ {- y} d y,\tag{B.4}
$$

which can be in turn written as: 

$$
\mathbb {E} [ u ] = \delta \Psi^ {1 / \epsilon}, \quad \delta = \Gamma \left(\frac {\epsilon - 1}{\epsilon}\right),\tag{B.5}
$$

where $\Gamma ( \cdot )$ is the Gamma function. Therefore we have the expression in the paper: 

$$
\mathbb {E} \left[ u \right] = \delta \Psi^ {1 / \epsilon} = \delta \left[ \sum_ {r = 1} ^ {S} \sum_ {s = 1} ^ {S} B _ {r s} \left(\kappa_ {r s} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}\right) ^ {- \epsilon} w _ {s} ^ {\epsilon} \right] ^ {1 / \epsilon}.\tag{B.6}
$$

## B.1.2 Residence and Workplace Choices

Using the distribution of utility for pairs of residence and employment locations, the probability that a worker chooses the bilateral commute from n to i out of all possible bilateral commutes is: 

$$
\begin{array}{l} \pi_ {n i} = \operatorname * {P r} \left[ u _ {n i} \geq \max \{u _ {r s} \}; \forall r, s \right], \\ = \int_ {0} ^ {\infty} \prod_ {s \neq i} G _ {n s} (u) \left[ \prod_ {r \neq n} \prod_ {s} G _ {r s} (u) \right] g _ {n i} (u) d u, \\ = \int_ {0} ^ {\infty} \prod_ {r = 1} ^ {S} \prod_ {s = 1} ^ {S} \epsilon \Psi_ {n i} u ^ {- (\epsilon + 1)} e ^ {- \Psi_ {r s} u ^ {- \epsilon}} d u. \\ = \int_ {0} ^ {\infty} \epsilon \Psi_ {n i} u ^ {- (\epsilon + 1)} e ^ {- \Psi u ^ {- \epsilon}} d u. \end{array}
$$

Note that: 

$$
\frac {d}{d u} \left[ - \frac {1}{\Psi} e ^ {- \Psi u ^ {- \epsilon}} \right] = \epsilon u ^ {- (\epsilon + 1)} e ^ {- \Psi u ^ {- \epsilon}}.\tag{B.7}
$$

Using this result to evaluate the integral above, the probability that the worker chooses to live in location n and commute to work in location i is: 

$$
\lambda_ {n i} = \frac {\Psi_ {n i}}{\Psi} = \frac {B _ {n i} (\kappa_ {n i} P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}) ^ {- \epsilon} (w _ {i}) ^ {\epsilon}}{\sum_ {r = 1} ^ {S} \sum_ {s = 1} ^ {S} B _ {r s} (\kappa_ {r s} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}) ^ {- \epsilon} (w _ {s}) ^ {\epsilon}}.\tag{B.8}
$$

Summing across all possible workplaces s, we obtain the probability that a worker chooses to live in location n out of all possible locations is: 

$$
\lambda_ {n} ^ {R} = \frac {R _ {n}}{\bar {L}} = \frac {\Psi_ {n}}{\Psi} = \frac {\sum_ {s = 1} ^ {S} B _ {n s} (\kappa_ {n s} P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}) ^ {- \epsilon} (w _ {s}) ^ {\epsilon}}{\sum_ {r = 1} ^ {S} \sum_ {s = 1} ^ {S} B _ {r s} (\kappa_ {r s} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}) ^ {- \epsilon} (w _ {s}) ^ {\epsilon}}.\tag{B.9}
$$

Similarly, summing across all possible residence locations $^ { r , }$ we obtain the probability that a worker chooses to work in location i out of all possible locations is: 

$$
\lambda_ {i} ^ {L} = \frac {L _ {i}}{\bar {L}} = \frac {\Psi_ {i}}{\Psi} = \frac {\sum_ {r = 1} ^ {S} B _ {r i} (\kappa_ {r i} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}) ^ {- \epsilon} (w _ {i}) ^ {\epsilon}}{\sum_ {r = 1} ^ {S} \sum_ {s = 1} ^ {S} B _ {r s} (\kappa_ {r s} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}) ^ {- \epsilon} (w _ {s}) ^ {\epsilon}}.\tag{B.10}
$$

For the measure of workers in location i $( L _ { i } )$ , we can evaluate the conditional probability that they commute from location n (conditional on having chosen to work in location i): 

$$
\begin{array}{r l} & {\lambda_ {n i | i} ^ {L} \equiv \frac {\lambda_ {n i}}{\lambda_ {i} ^ {L}} = \operatorname * {P r} \left[ u _ {n i} \geq \max \{u _ {r i} \}; \forall r \right],} \\ & {\qquad = \int_ {0} ^ {\infty} \prod_ {r \neq n} G _ {r i} (u) g _ {n i} (u) d u,} \\ & {\qquad = \int_ {0} ^ {\infty} e ^ {- \Psi_ {i} u ^ {- \epsilon}} \epsilon \Psi_ {n i} u ^ {- (\epsilon + 1)} d u.} \end{array}
$$

Using the result (B.7) to evaluate the integral above, the probability that a worker commutes from location n conditional on having chosen to work in location i is: 

$$
\lambda_ {n i | i} ^ {L} = \frac {B _ {n i} (\kappa_ {n i} P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}) ^ {- \epsilon} (w _ {i}) ^ {\epsilon}}{\sum_ {r = 1} ^ {S} B _ {r i} (\kappa_ {r i} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}) ^ {- \epsilon} (w _ {i}) ^ {\epsilon}},
$$

which simpliÖes to: 

$$
\lambda_ {n i | i} ^ {L} = \frac {B _ {n i} (\kappa_ {n i} P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}) ^ {- \epsilon}}{\sum_ {r = 1} ^ {S} B _ {r i} (\kappa_ {r i} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}) ^ {- \epsilon}}.\tag{B.11}
$$

For the measure of residents of location n $\left( R _ { n } \right)$ , we can evaluate the conditional probability that they commute to location i (conditional on having chosen to live in location n): 

$$
\begin{array}{r l} & {\lambda_ {n i | n} ^ {R} \equiv \frac {\lambda_ {n i}}{\lambda_ {n} ^ {R}} = \operatorname * {P r} [ u _ {n i} \geq \max \{u _ {n s} \}; \forall s ],} \\ & {\qquad = \int_ {0} ^ {\infty} \prod_ {s \neq i} G _ {n s} (u) g _ {n i} (u) d u,} \\ & {\qquad = \int_ {0} ^ {\infty} e ^ {- \Psi_ {n} u ^ {- \epsilon}} \epsilon \Psi_ {n i} u ^ {- (\epsilon + 1)} d u.} \end{array}
$$

Using the result (B.7) to evaluate the integral above, the probability that a worker commutes to location i conditional on having chosen to live in location n is: 

$$
\lambda_ {n i | n} ^ {R} = \frac {\Psi_ {n i}}{\Psi_ {n}} = \frac {B _ {n i} (\kappa_ {n i} P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}) ^ {- \epsilon} (w _ {i}) ^ {\epsilon}}{\sum_ {s = 1} ^ {S} B _ {n s} (\kappa_ {n s} P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}) ^ {- \epsilon} (w _ {s}) ^ {\epsilon}},
$$

which simpliÖes to: 

$$
\lambda_ {n i | n} ^ {R} = \frac {B _ {n i} (w _ {i} / \kappa_ {n i}) ^ {\epsilon}}{\sum_ {s = 1} ^ {S} B _ {n s} (w _ {s} / \kappa_ {n s}) ^ {\epsilon}}.\tag{B.12}
$$

These conditional commuting probabilities provide microeconomic foundations for the reduced-form gravity equations estimated in the empirical literature on commuting patterns.<sup>1</sup> The probability that a resident of location n commutes to location i depends on the wage at i and the amenities and commuting costs from living in n and working in i in the numerator (ìbilateral resistanceî). But it also depends on the wage at all other workplaces s and the amenities and commuting costs from living in n and commuting to all other workplaces s in the denominator (ìmultilateral resistanceî). 

Commuter market clearing requires that the measure of workers employed in each location i $( L _ { i } )$ equals the sum across all locations n of their measures of residents $\left( R _ { n } \right)$ times their conditional probabilities of commuting to $i ~ ( \lambda _ { n i | n } ^ { R } )$ : 

$$
\begin{array}{r l} & L _ {i} = \sum_ {n = 1} ^ {S} \lambda_ {n i | n} ^ {R} R _ {n} \\ & \qquad = \sum_ {n = 1} ^ {S} \frac {B _ {n i} (w _ {i} / \kappa_ {n i}) ^ {\epsilon}}{\sum_ {s = 1} ^ {S} B _ {n s} (w _ {s} / \kappa_ {n s}) ^ {\epsilon}} R _ {n}, \end{array}\tag{B.13}
$$

where, since there is a continuous measure of workers residing in each location, there is no uncertainty in the supply of workers to each employment location. 

Expected worker income conditional on living in location n equals the wages in all possible workplace locations weighted by the probabilities of commuting to those locations conditional on living in n: 

$$
\begin{array}{l} \bar {v} _ {n} = \mathbb {E} [ w | n ] \\ \qquad = \sum_ {i = 1} ^ {S} \lambda_ {n i | n} ^ {R} w _ {i}, \\ \qquad = \sum_ {i = 1} ^ {S} \frac {B _ {n i} (w _ {i} / \kappa_ {n i}) ^ {\epsilon}}{\sum_ {s = 1} ^ {S} B _ {n s} (w _ {s} / \kappa_ {n s}) ^ {\epsilon}} w _ {i}, \end{array}\tag{B.14}
$$

where E denotes the expectations operator and the expectation is taken over the distribution for idiosyncratic amenities. Intuitively, expected worker income is high in locations that have low commuting costs (low $\kappa _ { n s } )$ to high-wage employment locations. 

Finally, another implication of the FrÈchet distribution of utility is that the distribution of utility conditional on residing in location n and commuting to location i is the same across all bilateral pairs of locations with positive residents and employment, and is equal to the distribution of utility for the economy as a whole. To establish this result, note that the distribution of utility conditional on residing in location n and commuting to location i is given by: 

$$
\begin{array}{l} = \frac {1}{\lambda_ {n i}} \int_ {0} ^ {u} \prod_ {s \neq i} G _ {n s} (u) \left[ \prod_ {r \neq n} \prod_ {s} G _ {r s} (u) \right] g _ {n i} (u) d u, \\ = \frac {1}{\lambda_ {n i}} \int_ {0} ^ {u} \left[ \prod_ {r = 1} ^ {S} \prod_ {s = 1} ^ {S} e ^ {- \Psi_ {r s} u ^ {- \epsilon}} \right] \epsilon \Psi_ {n i} u ^ {- (\epsilon + 1)} d u, \\ = \frac {\Psi}{\Psi_ {n i}} \int_ {0} ^ {u} e ^ {- \Psi u ^ {- \epsilon}} \epsilon \Psi_ {n i} u ^ {- (\epsilon + 1)} d u, \\ = e ^ {- \Psi u ^ {\epsilon}}. \end{array}\tag{B.15}
$$

On the one hand, lower land prices in location n or a higher wage in location i raise the utility of a worker with a given realization of idiosyncratic amenities $b ,$ and hence increase the expected utility of residing in n and working in i. On the other hand, lower land prices or a higher wage induce workers with lowe realizations of idiosyncratic amenities b to reside in n and work in i, which reduces the expected utility of residing in n and working in i. With a FrÈchet distribution of utility, these two e§ects exactly o§set one another. Pairs of residence and employment locations with more attractive characteristics attract more commuters on the extensive margin until expected utility is the same across all pairs of residence and employment locations within the economy. 

## B.2 Computing Counterfactuals Using Changes

We now use the structure of the model to solve for a counterfactual equilibrium using the observed values of variables in an initial equilibrium. We denote the value of variables in the counterfactual equilibrium by a prime $( x ^ { \prime } )$ and the relative change of a variable between the initial and the counterfactual equilibrium by a hat $( { \widehat x } = x ^ { \prime } / x )$ . Given the modelís parameters $\{ \alpha , \ \sigma , \ \epsilon , \ \delta , \ \kappa \}$ and counterfactual changes in the bmodelís exogenous variables $\{ \hat { A } _ { n } , \hat { B } _ { n } , \hat { \kappa } _ { n i } , \hat { d } _ { n i } \}$ , we can solve for the counterfactual changes in the modelís endogenous variables $\{ \hat { w } _ { n } , \ \widehat { \bar { v } } _ { n } , \ \hat { Q } _ { n } , \ \hat { \pi } _ { n i } , \ \hat { \lambda } _ { n i } , \ \hat { P } _ { n } , \ \hat { R } _ { n } , \ \hat { L } _ { n } \}$ from the following system of eight equations b(using the iterative algorithm outlined below): 

$$
\widehat {w} _ {i} \widehat {L} _ {i} w _ {i} L _ {i} = \sum_ {n \in N} \pi_ {n i} \widehat {\pi} _ {n i} \widehat {\bar {v}} _ {n} \widehat {R} _ {n} \bar {v} _ {n} R _ {n},\tag{B.16}
$$

$$
\widehat {\overline {{v}}} _ {n} \overline {{v}} _ {n} = \sum_ {i \in N} \frac {\lambda_ {n i} \widehat {B} _ {n i} \left(\widehat {w} _ {i} / \widehat {\kappa} _ {n i}\right) ^ {\epsilon}}{\sum_ {s \in N} \lambda_ {n s} \widehat {B} _ {n s} \left(\widehat {w} _ {s} / \widehat {\kappa} _ {n s}\right) ^ {\epsilon}} \widehat {w} _ {i} w _ {i},\tag{B.17}
$$

$$
\widehat {Q} _ {n} = \widehat {\overline {{v}}} _ {n} \widehat {R} _ {n},\tag{B.18}
$$

$$
\widehat {\pi} _ {n i} \pi_ {n i} = \frac {\pi_ {n i} \widehat {L} _ {i} \left(\widehat {d} _ {n i} \widehat {w} _ {i} / \widehat {A} _ {i}\right) ^ {1 - \sigma}}{\sum_ {k \in N} \pi_ {n k} \widehat {L} _ {k} \left(\widehat {d} _ {n k} \widehat {w} _ {k} / \widehat {A} _ {k}\right) ^ {1 - \sigma}},\tag{B.19}
$$

$$
\widehat {\lambda} _ {n i} \lambda_ {n i} = \frac {\lambda_ {n i} \widehat {B} _ {n i} \left(\widehat {P} _ {n} ^ {\alpha} \widehat {Q} _ {n} ^ {1 - \alpha}\right) ^ {- \epsilon} (\widehat {w} _ {i} / \widehat {\kappa} _ {n i}) ^ {\epsilon}}{\sum_ {r \in N} \sum_ {s \in N} \lambda_ {r s} \widehat {B} _ {r s} \left(\widehat {P} _ {r} ^ {\alpha} \widehat {Q} _ {r} ^ {1 - \alpha}\right) ^ {- \epsilon} (\widehat {w} _ {s} / \widehat {\kappa} _ {r s}) ^ {\epsilon}},\tag{B.20}
$$

$$
\widehat {P} _ {n} = \left(\frac {\widehat {L} _ {n}}{\widehat {\pi} _ {n n}}\right) ^ {\frac {1}{1 - \sigma}} \frac {\widehat {d} _ {n n} \widehat {w} _ {n}}{\widehat {A} _ {n}},\tag{B.21}
$$

$$
\hat {R} _ {n} = \frac {\bar {L}}{R _ {n}} \sum_ {i} \lambda_ {n i} \hat {\lambda} _ {n i},\tag{B.22}
$$

$$
\hat {L} _ {i} = \frac {\bar {L}}{L _ {i}} \sum_ {n} \lambda_ {n i} \hat {\lambda} _ {n i},\tag{B.23}
$$

where these equations correspond to the equality between income and expenditure (B.16), expected worker income (B.17), land market clearing (B.18), trade shares (B.19), commuting probabilities (B.20), price indices (B.21), residential choice probabilities (B.22) and workplace choice probabilities (B.23). 

We solve this system of equations using the following iterative algorithm for the counterfactual equilib rium. Given the modelís parameters $\{ \alpha , \sigma , \epsilon , \delta , \kappa \}$ and changes in the exogenous variables of the model $\{ \hat { A } _ { n } , \hat { B } _ { n } , \hat { \kappa } _ { n i } , \hat { d } _ { n i } \}$ , we can solve for the resulting counterfactual changes in the endogenous variables of the model $\{ \hat { w } _ { n } , \widehat { \bar { v } } _ { n } , { \hat { Q } } _ { n } , \hat { \pi } _ { n i } , \hat { \lambda } _ { n i } , \hat { P } _ { n } , \hat { R } _ { n } , \hat { L } _ { n } \}$ from the system of eight equations (B.16)-(B.23). We solve bthis system of equations using the following iterative algorithm. We Örst conjecture changes in workplace wages and commuting probabilities at iteration t, $\hat { w } _ { i } ^ { ( t ) }$ and $\hat { \lambda } _ { n i } ^ { ( t ) }$ : We next update these conjectures to $\widehat w _ { i } ^ { ( t + 1 ) }$ and $\widehat { \lambda } _ { n i } ^ { ( t + 1 ) }$ using the current guesses and data. We start by computing: 

$$
\widehat {\overline {{v}}} _ {n} ^ {(t)} = \frac {1}{\overline {{v}} _ {n}} \sum_ {i \in N} \frac {\hat {B} _ {n i} \lambda_ {n i} \left(\hat {w} _ {i} ^ {(t)} / \widehat {\kappa} _ {n i}\right) ^ {\epsilon}}{\sum_ {s \in N} \hat {B} _ {n s} \lambda_ {n s} \left(\hat {w} _ {s} ^ {(t)} / \widehat {\kappa} _ {n s}\right) ^ {\epsilon}} \hat {w} _ {i} ^ {(t)} w _ {i},\tag{B.24}
$$

$$
\hat {L} _ {i} ^ {(t)} = \frac {\bar {L}}{L _ {i}} \sum_ {n} \lambda_ {n i} \hat {\lambda} _ {n i} ^ {(t)},\tag{B.25}
$$

$$
\hat {R} _ {n} ^ {(t)} = \frac {\bar {L}}{R _ {n}} \sum_ {i} \lambda_ {n i} \hat {\lambda} _ {n i} ^ {(t)},\tag{B.26}
$$

which are only a function of data and current guesses. We use (B.24) and (B.26) in (B.18) to compute: 

$$
\widehat {Q} _ {n} ^ {(t)} = \widehat {\overline {{v}}} _ {n} ^ {(t)} \widehat {R} _ {n} ^ {(t)}.\tag{B.27}
$$

We use (B.25) and (B.19) to compute: 

$$
\widehat {\pi} _ {n i} ^ {(t)} = \frac {\widehat {L} _ {i} ^ {(t)} (\widehat {d} _ {n i} \widehat {w} _ {i} ^ {(t)} / \widehat {A} _ {i}) ^ {1 - \sigma}}{\sum_ {k \in N} \pi_ {n k} \widehat {L} _ {k} ^ {(t)} (\widehat {d} _ {n k} \widehat {w} _ {k} ^ {(t)} / \widehat {A} _ {k}) ^ {1 - \sigma}}.\tag{B.28}
$$

We use (B.25), (B.28) and (B.21) to compute: 

$$
\widehat {P} _ {n} ^ {(t)} = \binom{\widehat {L} _ {n} ^ {(t)}}{\widehat {\pi} _ {n n} ^ {(t)}} ^ {\frac {1}{1 - \sigma}} \frac {\widehat {w} _ {n} ^ {(t)}}{\widehat {A} _ {n}}.\tag{B.29}
$$

We use (B.24)-(B.29) to rewrite (B.16) and (B.20) as: 

$$
\begin{array}{r c l} \tilde {w} _ {i} ^ {(t + 1)} & = & \frac {1}{Y _ {i} \widehat {L} _ {i} ^ {(t)}} \sum_ {n \in N} \pi_ {n i} \widehat {\pi} _ {n i} ^ {(t)} \widehat {\overline {{v}}} _ {n} ^ {(t)} \widehat {R} _ {n} ^ {(t)} Y _ {n}, \\ \tilde {\lambda} _ {n i} ^ {(t + 1)} & = & \frac {\hat {B} _ {n i} (\widehat {P} _ {n} ^ {(t) \alpha} \widehat {Q} _ {n} ^ {(t) 1 - \alpha}) ^ {- \epsilon} (\widehat {w} _ {i} ^ {(t)} / \widehat {\kappa} _ {n i}) ^ {\epsilon}}{\sum_ {r \in N} \sum_ {s \in N} \hat {B} _ {r s} \lambda_ {r s} (\widehat {P} _ {r} ^ {(t) \alpha} \widehat {Q} _ {r} ^ {(t) 1 - \alpha}) ^ {- \epsilon} (\widehat {w} _ {s} ^ {(t)} / \widehat {\kappa} _ {r s}) ^ {\epsilon}}. \end{array}\tag{B.30}
$$

(B.31) 

Finally, we update our conjectures for wages and commuting probabilities using: 

$$
\begin{array}{r c l} \hat {w} _ {i} ^ {(t + 1)} & = & \zeta \hat {w} _ {i} ^ {(t)} + (1 - \zeta) \tilde {w} _ {i} ^ {(t + 1)}, \\ \hat {\lambda} _ {i} ^ {(t + 1)} & = & \zeta \hat {\lambda} _ {i} ^ {(t)} + (1 - \zeta) \tilde {\lambda} _ {i} ^ {(t + 1)}, \end{array}\tag{B.32}
$$

(B.33) 

where $\zeta \in ( 0 , 1 )$ is an adjustment factor. 

In Section B.3 below, we provide conditions under which the counterfactual equilibrium of this economy is unique. 

## B.3 Existence and Uniqueness

We now provide conditions for the existence and uniqueness of a general equilibrium of this economy. 

## B.3.1 Workplace and Residence Income

From the commuting probabilities in equation (10) in the paper, the labor income received by commuters from residence n to workplace i is: 

$$
w _ {i} \lambda_ {n i} \bar {L} = \left(\frac {\bar {U}}{\delta}\right) ^ {- \epsilon} \bar {L} B _ {n i} \kappa_ {n i} ^ {- \epsilon} \left(P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}\right) ^ {- \epsilon} w _ {i} ^ {1 + \epsilon}.\tag{B.34}
$$

Summing across residences n, total workplace income in location i is: 

$$
Y _ {i} = \sum_ {n \in N} w _ {i} \lambda_ {n i} \bar {L} = \left(\frac {\bar {U}}{\delta}\right) ^ {- \epsilon} \bar {L} w _ {i} ^ {1 + \epsilon} \sum_ {n \in N} B _ {n i} \kappa_ {n i} ^ {- \epsilon} \left(P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}\right) ^ {- \epsilon},\tag{B.35}
$$

Summing across workplaces $i ,$ total residence income in location n (equals total expenditure in residence n) is given by: 

$$
X _ {n} = \sum_ {i \in N} w _ {i} \lambda_ {n i} \bar {L} = \left(\frac {\bar {U}}{\bar {\delta}}\right) ^ {- \epsilon} \bar {L} \left(P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}\right) ^ {- \epsilon} \sum_ {i \in N} B _ {n i} \kappa_ {n i} ^ {- \epsilon} w _ {i} ^ {1 + \epsilon}.\tag{B.36}
$$

Now note that land market clearing in equation (5) in the paper can be written as: 

$$
Q _ {n} = (1 - \alpha) \frac {X _ {n}}{H _ {n}}.\tag{B.37}
$$

Using land market clearing (B.37), total workplace income in location i (B.35) can be re-written as: 

$$
Y _ {i} = \left(\frac {\bar {U}}{\delta}\right) ^ {- \epsilon} \bar {L} (1 - \alpha) ^ {- (1 - \alpha) \epsilon} w _ {i} ^ {1 + \epsilon} \sum_ {n \in N} B _ {n i} \kappa_ {n i} ^ {- \epsilon} P _ {n} ^ {- \alpha \epsilon} X _ {n} ^ {- (1 - \alpha) \epsilon}.\tag{B.38}
$$

Using land market clearing (B.37), total residence income in location n (B.35) can be re-written as: 

$$
X _ {n} = \left(\frac {\bar {U}}{\delta}\right) ^ {- \epsilon} \bar {L} (1 - \alpha) ^ {- (1 - \alpha) \epsilon} H _ {n} ^ {(1 - \alpha) \epsilon} P _ {n} ^ {- \alpha \epsilon} X _ {n} ^ {- \epsilon (1 - \alpha)} \sum_ {i \in N} B _ {n i} \kappa_ {n i} ^ {- \epsilon} w _ {i} ^ {1 + \epsilon}.\tag{B.39}
$$

## B.3.2 Price Index and Goods Market Clearing

Using $Y _ { i } = w _ { i } L _ { i }$ , the price index in equation (8) in the paper can be re-written as: 

$$
P _ {n} ^ {1 - \sigma} = \left(\frac {\sigma}{\sigma - 1}\right) ^ {1 - \sigma} \frac {1}{\sigma F} \left[ \sum_ {i \in N} Y _ {i} \left(\frac {d _ {n i}}{A _ {i}}\right) ^ {1 - \sigma} w _ {i} ^ {- \sigma} \right].\tag{B.40}
$$

Similarly, using $Y _ { i } = w _ { i } L _ { i }$ and $X _ { n } = { \bar { v } } _ { n } R _ { n } .$ , the goods market clearing condition in equation (7) in the paper can be re-written as: 

$$
Y _ {i} = \sum_ {n \in N} \frac {Y _ {i}}{\sigma F} \left(\frac {\sigma}{\sigma - 1}\right) ^ {1 - \sigma} \left(\frac {d _ {n i}}{A _ {i}}\right) ^ {1 - \sigma} w _ {i} ^ {- \sigma} P _ {n} ^ {\sigma - 1} X _ {n},
$$

which simpliÖes to: 

$$
w _ {i} ^ {\sigma} = \sum_ {n \in N} \frac {1}{\sigma F} \left(\frac {\sigma}{\sigma - 1}\right) ^ {1 - \sigma} \left(\frac {d _ {n i}}{A _ {i}}\right) ^ {1 - \sigma} P _ {n} ^ {\sigma - 1} X _ {n}.\tag{B.41}
$$

## B.3.3 System of Equations

Combining workplace income (B.38), residence income (B.39), the price index (B.40), and goods market clearing (B.41), we obtain the following system of equations: 

$$
P _ {n} ^ {1 - \sigma} = \xi^ {P} \sum_ {i \in N} \mathcal {K} _ {n i} ^ {P} Y _ {i} w _ {i} ^ {- \sigma},\tag{B.42}
$$

$$
w _ {n} ^ {\sigma} = \xi^ {w} \sum_ {i \in N} \mathcal {K} _ {n i} ^ {w} P _ {i} ^ {\sigma - 1} X _ {i},\tag{B.43}
$$

$$
Y _ {n} w _ {n} ^ {- (1 + \epsilon)} = \xi^ {Y} \sum_ {i \in N} \mathcal {K} _ {n i} ^ {Y} P _ {i} ^ {- \alpha \epsilon} X _ {i} ^ {- (1 - \alpha) \epsilon},\tag{B.44}
$$

$$
X _ {n} ^ {1 + \epsilon (1 - \alpha)} P _ {n} ^ {\alpha \epsilon} = \xi^ {X} \sum_ {i \in N} \mathcal {K} _ {n i} ^ {X} w _ {i} ^ {1 + \epsilon},\tag{B.45}
$$

where we have assumed symmetric trade costs $( d _ { n i } = d _ { i n } )$ and commuting costs $( B _ { n i } \kappa _ { n i } ^ { - \epsilon } = B _ { i n } \kappa _ { i n } ^ { - \epsilon } )$ ; we have deÖned the following scalars: 

$$
\xi^ {P} \equiv \left(\frac {\sigma}{\sigma - 1}\right) ^ {1 - \sigma} \frac {1}{\sigma F},
$$

$$
\xi^ {w} \equiv \frac {1}{\sigma F} \left(\frac {\sigma}{\sigma - 1}\right) ^ {1 - \sigma},
$$

$$
\xi^ {Y} \equiv \left(\frac {\bar {U}}{\delta}\right) ^ {- \epsilon} \bar {L} (1 - \alpha) ^ {- (1 - \alpha) \epsilon},
$$

$$
\xi^ {X} \equiv \left(\frac {\bar {U}}{\delta}\right) ^ {- \epsilon} \bar {L} (1 - \alpha) ^ {- (1 - \alpha) \epsilon};
$$

and we have deÖned the following kernels: 

$$
\mathcal {K} _ {n i} ^ {P} \equiv \left(\frac {d _ {n i}}{A _ {i}}\right) ^ {1 - \sigma},
$$

$$
\mathcal {K} _ {n i} ^ {w} \equiv \left(\frac {d _ {n i}}{A _ {i}}\right) ^ {1 - \sigma},
$$

$$
\mathcal {K} _ {n i} ^ {Y} \equiv B _ {n i} \kappa_ {n i} ^ {- \epsilon} H _ {n} ^ {(1 - \alpha) \epsilon},
$$

$$
\mathcal {K} _ {n i} ^ {X} \equiv B _ {n i} \kappa_ {n i} ^ {- \epsilon} H _ {n} ^ {(1 - \alpha) \epsilon}.
$$

Note that equations (B.42)-(B.45) take the same form as the class of gravity equation models considered in Allen, Arkolakis and Li (2016). In particular, there are H vectors of endogenous variables $x ^ { h } \in \Re ^ { N } , h =$ $1 , \ldots , H$ , and each vector, $x ^ { h }$ , contains the endogenous variables for the I locations, $x _ { i } ^ { h } \in \Re , i = 1 , \dots , I .$ Using this notation, and denoting the corresponding sets of endogenous variables and locations by $\Omega ^ { H }$ and $\Omega ^ { N }$ respectively, the system of equations (B.42)-(B.45) can be written as: 

$$
\prod_ {h = 1} ^ {H} \left(x _ {i} ^ {h}\right) ^ {\beta_ {k h}} = \xi^ {k} \sum_ {n = 1} ^ {I} \mathcal {K} _ {n i} ^ {k} \left[ \prod_ {h = 1} ^ {H} \left(x _ {n} ^ {h}\right) ^ {\gamma_ {k h}} \right], \qquad i \in \Omega^ {N}, \qquad k, h \in \Omega^ {H},
$$

where the characteristic values $\xi ^ { k } \in \Re$ are endogenous scalars that balance the overall level of the two sides of the equations; the parameters are $\beta _ { k h } , \gamma _ { k h } \in \mathfrak { R } ;$ ; and ${ \kappa _ { n i } ^ { k } }$ is the kernel that regulates interactions across locations, variables and equations. 

We denote B and   as the $H \times H$ matrices, whose elements $( B ) _ { k h } = \beta _ { k h }$ and $( \Gamma ) _ { k h } = \gamma _ { k h }$ are the parameters from the left and right-hand sides of these equations, respectively. From equations (B.42)- (B.45), we have: 

$$
B = \left[ \begin{array}{c c c c} 1 - \sigma & 0 & 0 & 0 \\ 0 & \sigma & 0 & 0 \\ 0 & - (1 + \epsilon) & 1 & 0 \\ \alpha \epsilon & 0 & 0 & 1 + (1 - \alpha) \epsilon \end{array} \right],
$$

$$
\Gamma = \left[ \begin{array}{c c c c} 0 & - \sigma & 1 & 0 \\ \sigma - 1 & 0 & 0 & 1 \\ - \alpha \epsilon & 0 & 0 & - (1 - \alpha) \epsilon \\ 0 & 1 + \epsilon & 0 & 0 \end{array} \right].
$$

Note that all elements of the kernel ${ \kappa _ { n i } ^ { k } }$ are strictly positive. Additionally, both B and   are invertible, and we denote A as the following composite matrix: 

$$
A = \Gamma B ^ {- 1} = \left[ \begin{array}{c c c c} 0 & \frac {1}{\sigma} (\epsilon + 1) - 1 & 1 & 0 \\ \alpha \frac {\epsilon}{(\sigma - 1) (\epsilon - \alpha \epsilon + 1)} - 1 & 0 & 0 & \frac {1}{\epsilon - \alpha \epsilon + 1} \\ \alpha \frac {\epsilon}{\sigma - 1} + \alpha \epsilon^ {2} \frac {\alpha - 1}{(\sigma - 1) (\epsilon - \alpha \epsilon + 1)} & 0 & 0 & \epsilon \frac {\alpha - 1}{\epsilon - \alpha \epsilon + 1} \\ 0 & \frac {1}{\sigma} (\epsilon + 1) & 0 & 0 \end{array} \right].
$$

We also denote A<sup>p</sup> as the matrix whose elements equal the absolute value of the elements of A, such that $\left( A ^ { p } \right) _ { k h } = \vert ( A ) _ { k h } \vert$ , and deÖne $\rho \left( A ^ { p } \right)$ as the largest eigenvalue of $A ^ { p }$ . Applying Theorem 3 of Allen, Arkolakis and Li (2016), a su¢ cient condition for the equilibrium of the economy to be unique is $\rho \left( A ^ { p } \right) \leq 1$ . Having pinned down unique equilibrium values of $\{ P _ { n } , w _ { n } , Y _ { n } , X _ { n } \}$ , all other endogenous variables of the model can be uniquely determined. 

## B.4 Isomorphisms

## B.4.1 New Economic Geography Model with Commuting

We begin by considering our new economic geography model with agglomeration forces through love of variety and increasing returns to scale. The general equilibrium vector $\{ w _ { n } , \ { \bar { v } } _ { n } , \ Q _ { n } , \ L _ { n } , \ R _ { n } , \ P _ { n } \}$ and scalar $\bar { U }$ solve the following system of equations. First, income equals expenditure on goods produced in each location: 

$$
w _ {i} L _ {i} = \sum_ {n \in N} \frac {L _ {i} (d _ {n i} w _ {i} / A _ {i}) ^ {1 - \sigma}}{\sum_ {k \in N} L _ {k} (d _ {n k} w _ {k} / A _ {k}) ^ {1 - \sigma}} \bar {v} _ {n} R _ {n}.\tag{B.46}
$$

Second, expected worker income depends on wages: 

$$
\bar {v} _ {n} = \sum_ {i \in N} \frac {B _ {n i} (w _ {i} / \kappa_ {n i}) ^ {\epsilon}}{\sum_ {s \in N} B _ {n s} (w _ {s} / \kappa_ {n s}) ^ {\epsilon}} w _ {i}.\tag{B.47}
$$

Third, land prices depend on expected worker income and the measure of residents: 

$$
Q _ {n} = (1 - \alpha) \frac {\bar {v} _ {n} R _ {n}}{H _ {n}}.\tag{B.48}
$$

Fourth, workplace choice probabilities solve: 

$$
\frac {L _ {n}}{\bar {L}} = \frac {\sum_ {r \in N} B _ {r n} (\kappa_ {r n} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}) ^ {- \epsilon} w _ {n} ^ {\epsilon}}{\sum_ {r \in N} \sum_ {s \in N} B _ {r s} (\kappa_ {r s} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}) ^ {- \epsilon} w _ {s} ^ {\epsilon}}.\tag{B.49}
$$

Fifth, residential choice probabilities solve: 

$$
\frac {R _ {n}}{\bar {L}} = \frac {\sum_ {s \in N} B _ {n s} (\kappa_ {n s} P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}) ^ {- \epsilon} w _ {s} ^ {\epsilon}}{\sum_ {r \in N} \sum_ {s \in N} B _ {r s} (\kappa_ {r s} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}) ^ {- \epsilon} w _ {s} ^ {\epsilon}}.\tag{B.50}
$$

Sixth, price indices solve: 

$$
P _ {n} = \frac {\sigma}{\sigma - 1} \left(\frac {1}{\sigma F}\right) ^ {\frac {1}{1 - \sigma}} \left[ \sum_ {i \in N} L _ {i} (d _ {n i} w _ {i} / A _ {i}) ^ {1 - \sigma} \right] ^ {\frac {1}{1 - \sigma}}.\tag{B.51}
$$

Seventh, expected utility satisÖes: 

$$
\bar {U} = \delta \left[ \sum_ {r \in N} \sum_ {s \in N} B _ {r s} \left(\kappa_ {r s} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}\right) ^ {- \epsilon} w _ {s} ^ {\epsilon} \right] ^ {\frac {1}{\epsilon}},\tag{B.52}
$$

where $\begin{array} { r } { \delta = \Gamma \left( \frac { \epsilon - 1 } { \epsilon } \right) } \end{array}$ and   ( ) is the Gamma function. 

## B.4.2 Eaton and Kortum (2002) with External Economies of Scale and Commuting

We consider an Eaton and Kortum (2002) with external economies of scale augmented to incorporate heterogeneity in worker preferences over workplace and residence locations. Utility remains as speciÖed in equation (1) in the paper, except that the consumption index $\left( C _ { n } \right)$ is deÖned over a Öxed interval of goods 

$j \in [ 0 , 1 ]$ 

$$
C _ {n} = \left[ \int_ {0} ^ {1} c _ {n} (j) ^ {\rho} d j \right] ^ {\frac {1}{\rho}}.
$$

Productivity for each good j in each location i is drawn from an independent FrÈchet distribution: 

$$
F _ {i} (z) = e ^ {- A _ {i} z ^ {- \theta}}, \qquad A _ {i} = \tilde {A} _ {i} L _ {i} ^ {\eta}, \qquad \theta > 1,
$$

where the scale parameter of this distribution $\left( A _ { i } \right)$ depends on the measure of workers $( L _ { i } )$ and  parameterizes the strength of external economies of scale. The general equilibrium vector $\{ w _ { n } , \bar { v } _ { n } , Q _ { n } , L _ { n }$ 2 $\textstyle R _ { n } , P _ { n } \}$ and scalar $\bar { U }$ solve the following system of equations. First, income equals expenditure on goods produced in each location: 

$$
w _ {i} L _ {i} = \sum_ {n \in N} \frac {\tilde {A} _ {i} L _ {i} ^ {\eta} (d _ {n i} w _ {i}) ^ {- \theta}}{\sum_ {k \in N} \tilde {A} _ {k} L _ {k} ^ {\eta} (d _ {n k} w _ {k}) ^ {- \theta}} \bar {v} _ {n} R _ {n}.\tag{B.53}
$$

Second, expected worker income depends on wages: 

$$
\bar {v} _ {n} = \sum_ {i \in N} \frac {B _ {n i} (w _ {i} / \kappa_ {n i}) ^ {\epsilon}}{\sum_ {s \in N} B _ {n s} (w _ {s} / \kappa_ {n s}) ^ {\epsilon}} w _ {i}.\tag{B.54}
$$

Third, land prices depend on expected worker income and the measure of residents: 

$$
Q _ {n} = (1 - \alpha) \frac {\bar {v} _ {n} R _ {n}}{H _ {n}}.\tag{B.55}
$$

Fourth, workplace choice probabilities solve: 

$$
\frac {L _ {n}}{\bar {L}} = \frac {\sum_ {r \in N} B _ {r n} (\kappa_ {r n} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}) ^ {- \epsilon} w _ {n} ^ {\epsilon}}{\sum_ {r \in N} \sum_ {s \in N} B _ {r s} (\kappa_ {r s} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}) ^ {- \epsilon} w _ {s} ^ {\epsilon}}.\tag{B.56}
$$

Fifth, residential choice probabilities solve: 

$$
\frac {R _ {n}}{\bar {L}} = \frac {\sum_ {s \in N} B _ {n s} (\kappa_ {n s} P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}) ^ {- \epsilon} w _ {s} ^ {\epsilon}}{\sum_ {r \in N} \sum_ {s \in N} B _ {r s} (\kappa_ {r s} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}) ^ {- \epsilon} w _ {s} ^ {\epsilon}}.\tag{B.57}
$$

Sixth, price indices solve: 

$$
P _ {n} = \gamma \left[ \sum_ {i \in N} \tilde {A} _ {i} L _ {i} ^ {\eta} (d _ {n i} w _ {i}) ^ {- \theta} \right] ^ {- \frac {1}{\theta}},\tag{B.58}
$$

where $\begin{array} { r } { \gamma = \left[ \Gamma \left( \frac { \theta - \left( \sigma - 1 \right) } { \theta } \right) \right] ^ { \frac { 1 } { 1 - \sigma } } } \end{array}$ and $\Gamma \left( \cdot \right)$ denotes the Gamma function. Seventh, expected utility satisÖes: 

$$
\bar {U} = \delta \left[ \sum_ {r \in N} \sum_ {s \in N} B _ {r s} \left(\kappa_ {r s} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}\right) ^ {- \epsilon} w _ {s} ^ {\epsilon} \right] ^ {\frac {1}{\epsilon}}.\tag{B.59}
$$

The system of equations (B.53)-(B.59) is isomorphic to the system of equations (B.46)-(B.52) under the following parameter restrictions: 

$$
\begin{array}{r c l} \theta^ {\mathrm{EK}} & = & \sigma^ {\mathrm{NEG}} - 1, \\ \eta^ {\mathrm{EK}} & = & 1, \\ A _ {i} ^ {\mathrm{EK}} & = & \left(A _ {i} ^ {\mathrm{NEG}}\right) ^ {\sigma^ {\mathrm{NEG}} - 1}, \\ \gamma^ {\mathrm{EK}} & = & \frac {\sigma^ {\mathrm{NEG}}}{\sigma^ {\mathrm{NEG}} - 1} \left(\frac {1}{\sigma^ {\mathrm{NEG}} F ^ {\mathrm{NEG}}}\right) ^ {\frac {1}{1 - \sigma^ {\mathrm{NEG}}}}. \end{array}
$$

Under these parameter restrictions, both models generate the same general equilibrium vector $\{ w _ { n } , \bar { v } _ { n }$ 2 $Q _ { n } , L _ { n } , R _ { n } , P _ { n } \}$ and scalar $\bar { U } .$ 

## B.4.3 Armington (1969) with External Economies of Scale and Commuting

We consider an Armington (1969) model with external economies of scale augmented to incorporate heterogeneity in worker preferences over workplace and residence locations. Utility remains as speciÖed in equation (1) in the paper, except that the consumption index $\left( C _ { n } \right)$ is deÖned over goods that are horizontally di§erentiated by location of origin: 

$$
C _ {n} = \left[ \sum_ {i \in N} C _ {i} ^ {\rho} \right] ^ {\frac {1}{\rho}}.
$$

The goods supplied by each location are produced under conditions of perfect competition and externa economies of scale such that the ìcost inclusive of freightî(cif) price of a good produced in location i and consumed in location n is: 

$$
P _ {n i} = \frac {d _ {n i} w _ {i}}{A _ {i}}, \quad A _ {i} = \tilde {A} _ {i} L _ {i} ^ {\eta}.
$$

The general equilibrium vector $\{ w _ { n } , \ { \bar { v } } _ { n } , \ Q _ { n } , \ L _ { n } , \ R _ { n } , \ P _ { n } \}$ and scalar $\bar { U }$ solve the following system of equations. First, income equals expenditure on goods produced in each location: 

$$
w _ {i} L _ {i} = \sum_ {n \in N} \frac {A _ {i} ^ {\sigma - 1} L _ {i} ^ {\eta (\sigma - 1)} (d _ {n i} w _ {i}) ^ {1 - \sigma}}{\sum_ {k \in N} A _ {k} ^ {\sigma - 1} L _ {k} ^ {\eta (\sigma - 1)} (d _ {n k} w _ {k}) ^ {1 - \sigma}} \bar {v} _ {n} R _ {n}.\tag{B.60}
$$

Second, expected worker income depends on wages: 

$$
\bar {v} _ {n} = \sum_ {i \in N} \frac {B _ {n i} (w _ {i} / \kappa_ {n i}) ^ {\epsilon}}{\sum_ {s \in N} B _ {n s} (w _ {s} / \kappa_ {n s}) ^ {\epsilon}} w _ {i}.\tag{B.61}
$$

Third, land prices depend on expected worker income and the measure of residents: 

$$
Q _ {n} = (1 - \alpha) \frac {\bar {v} _ {n} R _ {n}}{H _ {n}}.\tag{B.62}
$$

Fourth, workplace choice probabilities solve: 

$$
\frac {L _ {n}}{\bar {L}} = \frac {\sum_ {r \in N} B _ {r n} (\kappa_ {r n} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}) ^ {- \epsilon} w _ {n} ^ {\epsilon}}{\sum_ {r \in N} \sum_ {s \in N} B _ {r s} (\kappa_ {r s} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}) ^ {- \epsilon} w _ {s} ^ {\epsilon}}.\tag{B.63}
$$

Fifth, residential choice probabilities solve: 

$$
\frac {R _ {n}}{\bar {L}} = \frac {\sum_ {s \in N} B _ {n s} \left(\kappa_ {n s} P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}\right) ^ {- \epsilon} w _ {s} ^ {\epsilon}}{\sum_ {r \in N} \sum_ {s \in N} B _ {r s} \left(\kappa_ {r s} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}\right) ^ {- \epsilon} w _ {s} ^ {\epsilon}}.\tag{B.64}
$$

Sixth, price indices solve: 

$$
P _ {n} = \left[ \sum_ {i \in N} A _ {i} ^ {\sigma - 1} L _ {i} ^ {\eta (\sigma - 1)} (d _ {n i} w _ {i}) ^ {1 - \sigma} \right] ^ {\frac {1}{1 - \sigma}}.\tag{B.65}
$$

Seventh, expected utility satisÖes: 

$$
\bar {U} = \delta \left[ \sum_ {r \in N} \sum_ {s \in N} B _ {r s} \left(\kappa_ {r s} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}\right) ^ {- \epsilon} w _ {s} ^ {\epsilon} \right] ^ {\frac {1}{\epsilon}}.\tag{B.66}
$$

The system of equations (B.60)-(B.66) is isomorphic to the system of equations (B.46)-(B.52) under the following parameter restrictions: 

$$
\begin{array}{r c l} \sigma^ {\mathrm{AR}} & = & \sigma^ {\mathrm{NEG}}, \\ \eta^ {\mathrm{AR}} & = & \frac {1}{\sigma^ {\mathrm{NEG}} - 1}, \\ A _ {i} ^ {\mathrm{AR}} & = & A _ {i} ^ {\mathrm{NEG}}, \\ 1 & = & \frac {\sigma^ {\mathrm{NEG}}}{\sigma^ {\mathrm{NEG}} - 1} \left(\frac {1}{\sigma^ {\mathrm{NEG}} F ^ {\mathrm{NEG}}}\right) ^ {\frac {1}{1 - \sigma^ {\mathrm{NEG}}}}. \end{array}
$$

Under these parameter restrictions, both models generate the same general equilibrium vector $\{ w _ { n } , \bar { v } _ { n }$ 2 $Q _ { n } , L _ { n } , R _ { n } , P _ { n } \}$ and scalar $\bar { U } .$ 

## B.5 Gravity in Goods Trade

As discussed in Section 3.1 of the paper, we use the equality between income and expenditure in equation (7) in the paper to solve for unobserved county productivities $\left( A _ { i } \right)$ 

$$
w _ {i} L _ {i} - \sum_ {n \in N} \frac {L _ {i} (d _ {n i} w _ {i} / A _ {i}) ^ {1 - \sigma}}{\sum_ {k \in N} L _ {k} (d _ {n k} w _ {k} / A _ {k}) ^ {1 - \sigma}} [ \bar {v} _ {n} R _ {n} + D _ {n} ] = 0,\tag{B.67}
$$

where we observe (or have solved for) wages $( w _ { i } )$ , employment $( L _ { i } )$ , average residential income $( { \bar { v } } _ { i } )$ , residents $( R _ { i } )$ and trade deÖcits $( D _ { i } )$ 

Given the elasticity of substitution $( \sigma )$ , our measures for $\left( w _ { i } , L _ { i } , \bar { v } _ { i } , R _ { i } , D _ { i } \right)$ and a parameterization of trade costs $( d _ { n i } ^ { 1 - \sigma } )$ , equation (B.67) provides a system of N equations that can be solved for a unique vector of N unobserved productivities $\left( A _ { i } \right)$ , as summarized in the following proposition. 

Proposition B.1 (Productivity Inversion) Given the elasticity of substitution $( \sigma )$ , our measures of wages, employment, average residential income, residents and trade deÖcits $\{ w _ { i } , ~ L _ { i } , ~ \bar { v } _ { i } , ~ R _ { i } , ~ D _ { i } \}$ , and a parameterization of trade costs $( d _ { n i } ^ { 1 - \sigma } )$ , there exist unique values of the unobserved productivities $\left( A _ { i } \right) f o r$ each location i that are consistent with the data being an equilibrium of the model. 

Proof. Note that the goods market clearing condition (B.67) can be written as the following excess demand system: 

$$
\mathbb {D} _ {i} (\tilde {\mathbf {A}}) = w _ {i} L _ {i} - \sum_ {n \in N} \frac {\tilde {A} _ {i} L _ {i} (d _ {n i} w _ {i}) ^ {1 - \sigma}}{\sum_ {k \in N} \tilde {A} _ {k} L _ {k} (d _ {n k} w _ {k}) ^ {1 - \sigma}} [ \bar {v} _ {n} R _ {n} + D _ {n} ] = 0,\tag{B.68}
$$

where $\tilde { A } _ { i } = A _ { i } ^ { \sigma - 1 } ; \{ w _ { i } , L _ { i } , \bar { v } _ { n } , R _ { n } , d _ { n i } \}$ have already been determined from the observed data or our parameterization of trade costs; and $\textstyle \sum _ { n \in N } D _ { n } = 0$ . This excess demand system exhibits the following properties in ${ \tilde { A } } _ { i } \colon$ 

Property $( \mathbf { i } ) \colon \mathbb { D } ( \tilde { \mathbf { A } } )$ is continuous, as follows immediately from inspection of (B.68). 

Property (ii): $\mathbb { D } ( \tilde { \mathbf { A } } )$ is homogenous of degree zero, as follows immediately from inspection of (B.68). 

Property (iii): $\begin{array} { r } { \sum _ { i \in N } \mathbb { D } _ { i } \left( \tilde { \mathbf { A } } \right) = 0 } \end{array}$ for all $\tilde { \mathbf { A } } \in \Re _ { + } ^ { N }$ . This property can be established by noting: 

$$
\begin{array}{r c l} \sum_ {i \in N} \mathbb {D} _ {i} (\tilde {\mathbf {A}}) & = & \sum_ {i \in N} w _ {i} L _ {i} - \sum_ {n \in N} \frac {\sum_ {i \in N} \tilde {A} _ {i} L _ {i} (d _ {n i} w _ {i}) ^ {1 - \sigma}}{\sum_ {k \in N} \tilde {A} _ {k} L _ {k} (d _ {n k} w _ {k}) ^ {1 - \sigma}} [ \bar {v} _ {n} R _ {n} + D _ {n} ], \\ & = & \sum_ {i \in N} w _ {i} L _ {i} - \sum_ {n \in N} [ \bar {v} _ {n} R _ {n} + D _ {n} ], \\ & = & 0. \end{array}
$$

Property $( \mathbf { i v } ) \colon \mathbb { D } ( \tilde { \mathbf { A } } )$ exhibits gross substitution: 

$$
\begin{array}{r l r l} & {\frac {\partial \mathbb {D} _ {i} (\tilde {\mathbf {A}})}{\partial \tilde {A} _ {r}} > 0} & & {\mathrm{forall} i, r, \neq i, \qquad \qquad \mathrm{forall} \tilde {\mathbf {A}} \in \Re_ {+} ^ {N},} \\ & {\frac {\partial \mathbb {D} _ {i} (\tilde {\mathbf {A}})}{\partial \tilde {A} _ {i}} <   0} & & {\mathrm{forall} i, \qquad \qquad \mathrm{forall} \tilde {\mathbf {A}} \in \Re_ {+} ^ {N}.} \end{array}
$$

This property can be established by noting: 

$$
\frac {\partial \mathbb {D} _ {i} (\tilde {\mathbf {A}})}{\partial \tilde {A} _ {r}} = \sum_ {n \in N} \frac {L _ {r} (d _ {n r} w _ {r}) ^ {1 - \sigma} \tilde {A} _ {i} L _ {i} (d _ {n i} w _ {i}) ^ {1 - \sigma}}{[ \sum_ {k \in N} \tilde {A} _ {k} L _ {k} (d _ {n k} w _ {k}) ^ {1 - \sigma} ] ^ {2}} [ \bar {v} _ {n} R _ {n} + D _ {n} ] > 0.
$$

and using homogeneity of degree zero, which implies: 

$$
\nabla \mathbb {D} (\tilde {\mathbf {A}}) \tilde {\mathbf {A}} = 0,
$$

and hence: 

$$
\frac {\partial \mathbb {D} _ {i} (\tilde {\mathbf {A}})}{\partial \tilde {A} _ {i}} <   0 \quad \text {for all} \tilde {\mathbf {A}} \in \Re_ {+} ^ {N}.
$$

Therefore we have established gross substitution. We now use these Öve properties to establish that the system of equations (B.68) has at most one (normalized) solution. Gross substitution implies that $\mathbb { D } \left( \tilde { \mathbf { A } } \right) = \mathbb { D } \left( \tilde { \mathbf { A } } ^ { \prime } \right)$ cannot occur whenever $\tilde { \mathbf { A } }$ and $\tilde { \mathbf { A } } ^ { \prime }$ are two technology vectors that are not colinear. $\mathrm { B y }$ homogeneity of degree zero, we can assume $\tilde { \mathbf { A } } ^ { \prime } \geq \tilde { \mathbf { A } }$ and ${ \tilde { A } } _ { i } = { \tilde { A } } _ { i } ^ { \prime }$ for some i. Now consider altering the productivity vector $\tilde { \mathbf { A } } ^ { \prime }$ to obtain the productivity vector $\tilde { \mathbf { A } }$ in $N - 1$ steps, lowering (or keeping unaltered) the productivity of all the other $N - 1$ locations $n \neq i$ one at a time. By gross substitution, the excess demand in location i cannot decrease in any step, and because $\tilde { \mathbf { A } } \neq \tilde { \mathbf { A } } ^ { \prime } ,$ , it will actually increase in at least one step. Hence $\mathbb { D } \left( \tilde { \mathbf { A } } \right) > \mathbb { D } \left( \tilde { \mathbf { A } } ^ { \prime } \right)$ and we have a contradiction. 

We next establish that there exists a productivity vector $\tilde { \mathbf { A } } ^ { * } \in \Re _ { + } ^ { N }$ such that $\mathbb { D } ( \tilde { \mathbf { A } } ^ { * } ) = 0$ . By homogeneity of degree zero, we can restrict our search for this productivity vector to the unit simplex $\begin{array} { r } { \Delta = \left\{ \tilde { \mathbf { A } } \in \Re _ { + } ^ { N } : \sum _ { i \in N } \tilde { A } _ { i } = 1 \right\} } \end{array}$ . DeÖne on $\Delta$ the function $\mathbb { D } ^ { + } \left( \cdot \right)$ by $\mathbb { D } _ { i } ^ { + } \left( \tilde { \mathbf { A } } \right) = \operatorname* { m a x } \left\{ \mathbb { D } _ { i } \left( \tilde { \mathbf { A } } \right) , 0 \right\}$ . Note that $\mathbb { D } ^ { + } \left( \cdot \right)$ is continuous. Denote $\begin{array} { r } { \alpha \left( \tilde { \mathbf { A } } \right) = \sum _ { i \in N } \left[ \tilde { A } _ { i } + \mathbb { D } _ { i } ^ { + } \left( \tilde { A } _ { i } \right) \right] } \end{array}$ . We have $\left( \tilde { \mathbf { A } } \right) \geq 1$ for all $\tilde { \mathbf { A } }$ . 

DeÖne a continuous function $f \left( \cdot \right)$ from the closed convex set $\Delta$ into itself by: 

$$
f (\tilde {\mathbf {A}}) = \left[ 1 / \alpha (\tilde {\mathbf {A}}) \right] \left[ \tilde {\mathbf {A}} + \mathbb {D} ^ {+} (\tilde {\mathbf {A}}) \right].
$$

Note that this Öxed-point function tends to increase the productivities of locations with excess demand. By Brouwerís Fixed-point Theorem, there exists $\tilde { \mathbf { A } } ^ { * } \in \Delta$ such that $\tilde { \mathbf { A } } ^ { * } = f \left( \tilde { \mathbf { A } } ^ { * } \right)$ 

Since $\begin{array} { r } { \sum _ { i \in N } \mathbb { D } _ { i } \left( \tilde { \mathbf { A } } \right) = 0 } \end{array}$ , it cannot be the case that $\mathbb { D } _ { i } \left( \tilde { \mathbf { A } } \right) > 0$ for all $i \in N$ or $\mathbb { D } _ { i } \left( \tilde { \mathbf { A } } \right) < 0$ for all $i \in N$ . Additionally, if $\mathbb { D } _ { i } \left( \tilde { \mathbf { A } } \right) > 0$ for some i and $\mathbb { D } _ { r } \left( \tilde { \mathbf { A } } \right) < 0$ for some $r \neq i , { \tilde { \mathbf { A } } } \neq f \left( { \tilde { \mathbf { A } } } \right)$ . It follows that at the Öxed point for productivity, $\tilde { \mathbf { A } } ^ { * } = f \left( \tilde { \mathbf { A } } ^ { * } \right)$ , and $\mathbb { D } _ { i } \left( \tilde { \mathbf { A } } ^ { * } \right) = 0$ for all i. It follows that there exists a unique vector of unobserved productivities $( \tilde { \mathbf { A } } )$ that solves the excess demand system (B.68). 

The resulting solutions for productivities $\left( A _ { i } \right)$ capture characteristics $\left( \mathrm { e . g . } \right.$ . natural resources) that make a location more or less attractive for employment conditional on the observed data and the parameterized values of trade costs. These characteristics include access to international markets. To the extent that such international market access raises employment $( L _ { i } )$ , and international trade áows are not captured in the CFS, this will be reáected in the model in higher productivity $\left( A _ { i } \right)$ to rationalize the higher observed employment. Having recovered these unique unobserved productivities $\left( A _ { i } \right)$ , we can solve for the implied bilateral trade áows between counties $( X _ { n i } )$ using equation (6) and $X _ { n i } = \pi _ { n i } \bar { v } _ { n } R _ { n }$ . We use these solutions for bilateral trade between counties in our counterfactuals for changes in the modelís exogenous variables, as discussed in the paper. 

To parameterize trade costs $( d _ { n i } ^ { 1 - \sigma } )$ , we assume a central value for the elasticity of substitution between varieties from the existing empirical literature of $\sigma = 4$ , which is in line with the estimates of this parameter using price and expenditure data in Broda and Weinstein (2006).<sup>2</sup> We model bilateral trade costs $\left( d _ { n i } \right)$ as a function of distance. For bilateral pairs with positive trade, we assume that bilateral trade costs are a constant elasticity function of distance and a stochastic erro $( d _ { n i } = d i s t _ { n i } ^ { \psi } \tilde { e } _ { n i } )$ . For bilateral pairs with zero trade, the model implies prohibitive trade costs $( d _ { n i } \to \infty ) . ^ { 3 }$ Taking logarithms in the trade share in equation (6) in the paper for pairs with positive trade, the value of bilateral trade between source i and destination n $( X _ { n i } )$ can be expressed as 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/9c8e29c1cb5b8fcc33217959e5d3c68fb327c60ead84869db9d34fcbee3608fa.jpg)



Figure B.1: Gravity in Goods Trade Between CFS Regions


$$
\log X _ {n i} = \zeta_ {n} + \chi_ {i} - (\sigma - 1) \psi \log d i s t _ {n i} + \log e _ {n i},\tag{B.69}
$$

where the source Öxed e§ect $( \chi _ { i } )$ controls for employment, wages and productivity $\left( L _ { i } , \ w _ { i } , \ A _ { i } \right)$ ; the destination Öxed e§ect $\left( \zeta _ { n } \right)$ controls for average income, $\bar { v } _ { n }$ , residents, $R _ { n }$ , and multilateral resistance (as captured in the denominator of equation (6) in the paper); and log $e _ { n i } = ( 1 - \sigma ) \log { \tilde { e } _ { n i } }$ 

Estimating the gravity equation (B.69) for all bilateral pairs with positive trade using OLS, we Önd a regression R-squared of 0.83. In Figure B.1, we display the conditional relationship between the log value of trade and log distance, after removing source and destination Öxed e§ects from both log trade and log distance. Consistent with the existing empirical trade literature, we Önd that the log linear functional form provides a good approximation to the data, with a tight and approximately linear relationship between the two variables. We estimate a coe¢cient on log distance of $- \left( \sigma - 1 \right) \psi = - 1 . 2 9$ . For our assumed value of $\sigma = 4 .$ , this implies an elasticity of trade costs with respect to distance of $\psi = 0 . 4 3$ . The tight linear relationship in Figure B.1, makes us conÖdent in this parametrization of trade costs as $d _ { n i } ^ { 1 - \sigma } = d i s t _ { n i } ^ { - 1 . 2 9 }$ as a way of using equation (B.67) to solve for unobserved productivities $\left( A _ { i } \right)$ 

To provide an alternative check on our speciÖcation, we aggregate the modelís predictions for trade between counties within pairs of CFS regions, and compare these predictions to the data in Figure B.2. The only way in which we used the data on trade between CFS regions was to estimate the distance elasticity $- \left( \sigma - 1 \right) \psi = - 1 . 2 9$ . Given this distance elasticity, we use the goods market clearing condition (B.67) to solve for productivities and generate predictions for bilateral trade between counties and hence CFS regions, as discussed above. Therefore, the modelís predictions and the data can di§er from one another. Nonetheless, we Önd a strong and approximately log linear relationship between the modelís predictions and the data, which is tighter for the larger trade values that account for most of aggregate trade. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/4773bb026827c5c65ecfdc76236b5e2dc1669b911d64e2bf1b69df6045190f4f.jpg)



Figure B.2: Bilateral Trade Shares in the Model and Data


## B.6 Magnitude and Gravity of Commuting Flows

In this subsection of the web appendix, we provide additional evidence on the relevance of commuting as a source of spatial linkages between counties and CZs. In Figure 1 in the paper, we display unweighted kernel densities of the share of residents that work in the same county where they live (the ìresidence own commuting shareî) over time. We focus on these unweighted kernel densities to capture heterogeneity across geographical locations (counties). As a robustness check, Figure B.3 in this web appendix displays analogous kernel densities that are weighted by the number of residents in each county. Therefore these weighted kernel densities capture heterogeneity across residents. As apparent from the two Ögures, we Önd a similar pattern of results whether we use the weighted or unweighted kernel densities. In both cases, we Önd a marked shift in density towards lower values of the residence own commuting share. 

In Table 1 of the paper, we report unweighted descriptive statistics on commuting áows between counties and CZs from 2006-10. As a robustness check, Table B.1 in this web appendix displays reports analogous statistics that are weighted by the number of residents (or workers) in each county. Again the unweighted results capture heterogeneity across counties, while the weighted results capture heterogeneity across people. Whether we use the weighted or unweighted statistics, we Önd that commuting beyond county boundaries is both substantial and heterogeneous. For example, using the unweighted results, we Önd that for the median county around 27 percent of its residents work outside the county and around 20 percent of its workers live outside the county. By comparison, using the weighted results, we Önd that for the median county around 19 percent of its residents work outside the county and around 22 percent of its workers live outside the county. 

In Section 3.2 of the paper, we discuss that these di§erences across counties in openness to commuting generate substantial variation in the ratio of employment to residents $( L _ { i } / R _ { i } )$ . In Table B.2 below, we show that this ratio of employment to residents is not only heterogeneous across counties, but is also hard 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/c4cc69db10e88c6ad97892d5f271bc9580e416e731fea26ce487dbf749383324.jpg)



Figure B.3: Kernel densities of the share of residents that work in the county where they live (weighted by county residents)


<table><tr><td></td><td>Min</td><td>p5</td><td>p10</td><td>p25</td><td>p50</td><td>p75</td><td>p90</td><td>p95</td><td>Max</td><td>Mean</td><td>N</td></tr><tr><td>Commuters from Residence County</td><td>0.00</td><td>0.00</td><td>0.04</td><td>0.08</td><td>0.19</td><td>0.38</td><td>0.51</td><td>0.56</td><td>0.82</td><td>0.24</td><td>3,111</td></tr><tr><td>Commuters to Workplace County</td><td>0.00</td><td>0.01</td><td>0.04</td><td>0.14</td><td>0.22</td><td>0.33</td><td>0.43</td><td>0.52</td><td>0.81</td><td>0.24</td><td>3,111</td></tr><tr><td>County Employment/Residents</td><td>0.26</td><td>0.64</td><td>0.73</td><td>0.88</td><td>1.00</td><td>1.10</td><td>1.21</td><td>1.30</td><td>3.88</td><td>1.00</td><td>3,111</td></tr><tr><td>Commuters from Residence CZ</td><td>0.00</td><td>0.00</td><td>0.00</td><td>0.02</td><td>0.03</td><td>0.07</td><td>0.14</td><td>0.18</td><td>0.49</td><td>0.06</td><td>709</td></tr><tr><td>Commuters to Employment CZ</td><td>0.00</td><td>0.00</td><td>0.00</td><td>0.02</td><td>0.05</td><td>0.08</td><td>0.11</td><td>0.13</td><td>0.25</td><td>0.06</td><td>709</td></tr><tr><td>CZ Employment/Residents</td><td>0.63</td><td>0.91</td><td>0.94</td><td>0.99</td><td>1.00</td><td>1.02</td><td>1.05</td><td>1.07</td><td>1.12</td><td>1.00</td><td>709</td></tr></table>

Tabulations on 3,111 counties and 709 commuting zones. The Örst row shows the fraction of residents that work outside the county. The second row shows the fraction of workers who live outside the county. The third row shows the ratio of county employment to county residents. The fourth row shows the fraction of a CZís residents that work outside the CZ. The Öfth row shows the fraction of a CZís workers that live outside the CZ. The sixth row shows the ratio of CZ employment to CZ residents across all 709 CZs. p5, p10 etc refer to the 5th, 10th etc percentiles of the distribution. Results for commuters from residence are weighted by the number of residents. Results for commuters to workplace are weighted by the number of workers. 

Table B.1: Commuting Across Counties and Commuting Zones (Weighted) 

to explain with the standard empirical controls used in the local labor markets literature (such as various measures of size, area, income and housing supply elasticities). Therefore these results establish that this role of the initial ratio of employment to residents in understanding the e§ects of changes in commuting costs cannot be easily proxied for by these other controls. 

In particular, Table B.2 reports the results of regressing log employment $\left( \log L _ { i } \right)$ , log residents (log R<sub>i</sub>), and the ratio of employment to residents $( L _ { i } / R _ { i } )$ on a number of standard empirical controls from the local labor markets literature. The Örst four columns show that the levels of either employment (log L<sub>i</sub>) or residents (log R ) are strongly related to these standard empirical controls. The Örst column shows tha one can account for most of the variation in county employment using the number of residents and wages. Column (2) shows a similar result for the number of residents and Columns (3) and (4) show that the results are not a§ected when we add land area, developed-land supply elasticities, employment and wages in surrounding counties. In contrast, the remaining four columns demonstrate that it is hard to explain the ratio of employment to residents $( L _ { i } / R _ { i } )$ using these same empirical controls. The level of residents, wages, land area, developed-land supply elasticities, employment, and measures of economic activity in surrounding counties, do a poor job in accounting for the variation in this ratio. None of the R-squaredís in the last four columns of Table B.2 amounts to more than one third. Taken together, these results conÖrm that the ratio of employment to residents $( L _ { i } / R _ { i } )$ cannot be easily proxied for by the standard empirical controls used in the local labor markets literature. 

To examine the extent to which bilateral commuting áows are one-way versus two-way, we use the Grubel and Lloyd (1971) index from the international trade literature. In the context of commuting, this Grubel-Lloyd index captures the extent there is (i) one-way commuting, in which counties either only export or only import commuters, versus (ii) two-way commuting, in which counties simultaneously export and import commuters. SpeciÖcally, the Grubel-Lloyd index for county i is deÖned as 

$$
G L _ {i} = 1 - \frac {\left| \sum_ {n \neq i} L _ {i n} - \sum_ {n \neq i} L _ {n i} \right|}{\sum_ {n \neq i} L _ {i n} + \sum_ {n \neq i} L _ {n i}},\tag{B.70}
$$

where the Örst subscript is the county of residence and the second subscript is the county of workplace. Therefore, $\textstyle \sum _ { n \neq i } L _ { i n }$ is county iís total exports of commuters to workplaces in other counties $n \neq i$ and $\textstyle \sum _ { n \neq i } L _ { n i }$ is county iís total imports of commuters from residences in other counties n $\neq i$ . If there is only one-way commuting, $G L _ { n } = 0$ . In contrast, if there is perfect two-way commuting, with county iís exports of commuters equal to its imports, $G L _ { n } = 1$ 

In Table B.3, we report the mean and percentiles of the distribution of the Grubel-Lloyd index from equation (B.70) across counties. We Önd pervasive two-way commuting, with the mean and median values of the Grubel-Lloyd index closer to perfect two-way commuting than to only one-way commuting. This pattern of results is consistent with the predictions of the model, in which workersíidiosyncratic preferences between pairs of residence and workplace in general induce two-way commuting. As discussed in Subsection 3.2 of the paper, the model rationalizes zero commuting áows from residence n to workplace i in terms of negligible amenities $( B _ { n i } \to 0 )$ and/or prohibitive commuting cost $( \kappa _ { n i }  \infty )$ , which can be used to explain one-way commuting. 

<table><tr><td></td><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td><td>7</td><td>8</td><td>9</td><td>10</td><td>11</td><td>12</td></tr><tr><td>Dep. Variable:</td><td><eq>\log L_i</eq></td><td><eq>\log R_i</eq></td><td><eq>\log L_i</eq></td><td><eq>\log R_i</eq></td><td><eq>L_i/R_i</eq></td><td><eq>L_i/R_i</eq></td><td><eq>L_i/R_i</eq></td><td><eq>L_i/R_i</eq></td><td><eq>L_i/R_i</eq></td><td><eq>L_i/R_i</eq></td><td><eq>L_i/R_i</eq></td><td><eq>L_i/R_i</eq></td></tr><tr><td><eq>\log R_i</eq></td><td>0.974**(0.005)</td><td></td><td>1.001**(0.009)</td><td></td><td></td><td>-0.000(0.010)</td><td></td><td>0.020(0.016)</td><td></td><td>0.064**(0.012)</td><td></td><td>0.050**(0.013)</td></tr><tr><td><eq>\log w_i</eq></td><td>0.460**(0.038)</td><td></td><td>0.480**(0.036)</td><td></td><td>0.341**(0.025)</td><td></td><td>0.331**(0.026)</td><td></td><td>0.468**(0.046)</td><td></td><td>0.479**(0.054)</td><td></td></tr><tr><td><eq>\log L_i</eq></td><td></td><td>0.957**(0.013)</td><td></td><td>0.922**(0.015)</td><td>-0.001(0.007)</td><td></td><td>0.028*(0.012)</td><td></td><td>0.044**(0.006)</td><td></td><td>0.049**(0.006)</td><td></td></tr><tr><td><eq>\log \bar{v}_i</eq></td><td></td><td>0.066(0.051)</td><td></td><td>0.019(0.049)</td><td></td><td>0.171**(0.033)</td><td></td><td>0.239**(0.041)</td><td></td><td>0.287**(0.100)</td><td></td><td>0.273**(0.092)</td></tr><tr><td><eq>\log H_i</eq></td><td></td><td></td><td>0.015(0.010)</td><td>0.037**(0.012)</td><td></td><td></td><td>-0.022(0.011)</td><td>-0.011(0.013)</td><td>-0.055**(0.012)</td><td>-0.067**(0.020)</td><td>-0.058**(0.013)</td><td>-0.059**(0.019)</td></tr><tr><td><eq>\log R_{,-i}</eq></td><td></td><td></td><td>-0.020*(0.008)</td><td></td><td></td><td></td><td>0.389*(0.160)</td><td>0.609**(0.171)</td><td>0.396(0.407)</td><td>0.677(0.524)</td><td>0.391(0.406)</td><td>0.679(0.516)</td></tr><tr><td><eq>\log \bar{w}_{-i}</eq></td><td></td><td></td><td>-0.330**(0.036)</td><td></td><td></td><td></td><td>0.070(0.324)</td><td>0.247(0.364)</td><td>-1.843(1.326)</td><td>-2.619(1.668)</td><td>-1.797(1.308)</td><td>-2.654(1.678)</td></tr><tr><td><eq>\log L_{,-i}</eq></td><td></td><td></td><td></td><td>0.084**(0.011)</td><td></td><td></td><td>-0.435**(0.155)</td><td>-0.654**(0.166)</td><td>-0.405(0.408)</td><td>-0.694(0.525)</td><td>-0.401(0.408)</td><td>-0.691(0.518)</td></tr><tr><td><eq>\log \bar{v}_{-i}</eq></td><td></td><td></td><td></td><td>0.044(0.038)</td><td></td><td></td><td>-0.238(0.315)</td><td>-0.347(0.360)</td><td>1.523(1.310)</td><td>2.410(1.640)</td><td>1.482(1.292)</td><td>2.431(1.648)</td></tr><tr><td>Saiz elasticity</td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td>0.010(0.008)</td><td>-0.022*(0.010)</td></tr><tr><td>Constant</td><td>-4.667**(0.413)</td><td>-0.165(0.431)</td><td>-1.485**(0.323)</td><td>-1.199*(0.473)</td><td>-2.647**(0.282)</td><td>-0.881**(0.285)</td><td>-0.262(0.422)</td><td>-0.057(0.500)</td><td>-0.636(0.560)</td><td>0.000(1.010)</td><td>-0.839(0.588)</td><td>0.373(0.966)</td></tr><tr><td><eq>R^2</eq></td><td>0.98</td><td>0.98</td><td>0.99</td><td>0.98</td><td>0.16</td><td>0.03</td><td>0.30</td><td>0.15</td><td>0.68</td><td>0.53</td><td>0.68</td><td>0.54</td></tr><tr><td>N</td><td>3,111</td><td>3,111</td><td>3,081</td><td>3,081</td><td>3,111</td><td>3,111</td><td>3,081</td><td>3,081</td><td>457</td><td>457</td><td>457</td><td>457</td></tr></table>


Note: $\begin{array} { r } { L , - i \equiv \sum _ { n : d _ { n i } \leq 1 2 0 , n \neq i } L _ { n } } \end{array}$ is the total employment in i neighbors whose centroid is no more than 120km away; $\bar { w } _ { - i } \equiv$ $\begin{array} { r } { \sum _ { n : d _ { n i } \leq 1 2 0 , n \neq i } \frac { L _ { n } } { L , - i } w _ { n } } \end{array}$ is the weighted average of their workplace wage. Analogous deÖnitions apply to $R _ { , - i }$ and $ { \bar { v } } _ { - i } .$ . Columns 1-8 use the whole sample of counties. Columns 9 and 10 repeat the most complete speciÖcations in columns 7 and 8 only for the subsample of counties where we have data on land supply elasticity. Columns 11 and 12 repeat columns 7 and 8 adding the Saiz land supply elasticity as a regressor. Standard errors are clustered by state.  denotes signiÖcance at the 5 percen level;  denotes signiÖcance at the 1 percent level. 



Table B.2: Explaining employment $( L _ { i } )$ , residents $( R _ { i } )$ , and the ratio of employment to residents $( L _ { i } / R _ { i } )$


<table><tr><td>Statistic</td><td>Grubel-Lloyd Index for Commuting</td></tr><tr><td>p5</td><td>0.342</td></tr><tr><td>p10</td><td>0.414</td></tr><tr><td>p25</td><td>0.537</td></tr><tr><td>p50</td><td>0.696</td></tr><tr><td>p75</td><td>0.843</td></tr><tr><td>p90</td><td>0.937</td></tr><tr><td>p95</td><td>0.968</td></tr><tr><td>Mean</td><td>0.681</td></tr></table>


Mean and percentiles of the distribution of the Grubel-Lloyd index from equation (B.70) across counties. 


Table B.3: Grubel-Lloyd Index for Commuting 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/fb2e8b8b4c8efdd2c1b7f52e75dee5a57b1cccb9dc7a2d64da90079b9a5d5d25.jpg)



Figure B.4: Kernel densities of the share of residents that work in the CZ where they live


Additionally, to provide a point of comparison to Figure 1 in the paper for counties, Figure B.4 shows kernel densities of the share of residents that work in the same CZ where they live for 1990 and 2000. We construct these measures for CZs from the matrices of bilateral commuting probabilities between counties, which are only reported in the Population Census from 1990 onwards. We Önd the same pattern of an increase in commuting openness over time, with the increase between 1990 and 2000 for CZs in Figure B.4 larger than the increase over the same period for counties in Figure 1 in the paper. 

As discussed in Section 3.2 of the paper, the gravity equation for the commuting probability in equation (10) in the paper can be written as 

$$
\lambda_ {n i} - \frac {\mathcal {B} _ {n i} \left(\frac {L _ {n}}{\pi_ {n n}}\right) ^ {- \frac {\alpha \epsilon}{\sigma - 1}} A _ {n} ^ {\alpha \epsilon} w _ {n} ^ {- \alpha \epsilon} \bar {v} _ {n} ^ {- \epsilon (1 - \alpha)} \left(\frac {R _ {n}}{H _ {n}}\right) ^ {- \epsilon (1 - \alpha)} w _ {i} ^ {\epsilon}}{\sum_ {r \in N} \sum_ {s \in N} \mathcal {B} _ {r s} \left(\frac {L _ {r}}{\pi_ {r r}}\right) ^ {- \frac {\alpha \epsilon}{\sigma - 1}} A _ {r} ^ {\alpha \epsilon} w _ {r} ^ {- \alpha \epsilon} \bar {v} _ {r} ^ {- \epsilon (1 - \alpha)} \left(\frac {R _ {r}}{H _ {r}}\right) ^ {- \epsilon (1 - \alpha)} w _ {s} ^ {\epsilon}} = 0,\tag{B.71}
$$

where $B _ { n i } \equiv B _ { n i } \kappa _ { n i } ^ { - \epsilon }$ is a composite parameter that captures the ease of commuting. The commuting probabilities (B.71) provide a system of $N \times N$ equations that can be solved for a unique matrix of $N \times N$ values of the ease of commuting $( B _ { n i } )$ , as summarized in the following proposition. 

Proposition B.2 (Amenities Inversion) Given the share of consumption goods in expenditure (), the heterogeneity in location preferences (), the observed data on wages, employment, trade shares, average residential income, residents and land area $\{ w _ { i } , L _ { i } , \pi _ { i i } , \bar { v } _ { i } , R _ { i } , H _ { i } \}$ , there exist unique values of the ease of commuting $( B _ { n i } \equiv B _ { n i } \kappa _ { n i } ^ { - \epsilon } )$ for each pair of locations n and i that are consistent with the data being an equilibrium of the model. 

Proof. Note that the commuting probability (B.71) can be written as the following excess demand system: 

$$
\mathbb {D} _ {i} (\mathcal {B}) = \lambda_ {n i} - \frac {\mathcal {B} _ {n i} \left(\frac {L _ {n}}{\pi_ {n n}}\right) ^ {- \frac {\alpha \epsilon}{\sigma - 1}} A _ {n} ^ {\alpha \epsilon} w _ {n} ^ {- \alpha \epsilon} \bar {v} _ {n} ^ {- \epsilon (1 - \alpha)} \left(\frac {R _ {n}}{H _ {n}}\right) ^ {- \epsilon (1 - \alpha)} w _ {i} ^ {\epsilon}}{\sum_ {r \in N} \sum_ {s \in N} \mathcal {B} _ {r s} \left(\frac {L _ {r}}{\pi_ {r r}}\right) ^ {- \frac {\alpha \epsilon}{\sigma - 1}} A _ {r} ^ {\alpha \epsilon} w _ {r} ^ {- \alpha \epsilon} \bar {v} _ {r} ^ {- \epsilon (1 - \alpha)} \left(\frac {R _ {r}}{H _ {r}}\right) ^ {- \epsilon (1 - \alpha)} w _ {s} ^ {\epsilon}} = 0,\tag{B.72}
$$

where $\{ w _ { i } , L _ { i } , \bar { v } _ { n } , R _ { n } , \pi _ { n n } , A _ { n } , H _ { n } \}$ have already been determined from the observed data or our parameterization of commuting costs. Note that the excess demand system (B.72) exhibits the same properties in $\boldsymbol { B }$ as the excess demand system (B.68) exhibits in $\tilde { \mathbf { A } }$ . It follows that there exists a unique vector of unobserved values of the ease of commuting ( ) that solves the excess demand system (B.72). 

The resulting solutions for the ease of commuting $( B _ { n i } )$ capture all factors that make a pair of residence and workplace locations more $\mathrm { o r }$ less attractive conditional on the observed wages, employment, trade shares, average residential income, residents and land area (e.g. attractive scenery, distance and transport infrastructure). Together productivity $\left( A _ { i } \right)$ and the ease of commuting $( B _ { n i } )$ correspond to structural residuals that ensure that the model exactly replicates the observed data given the parameters. 

To estimate the heterogeneity in location preferences $( \epsilon )$ , we model the determinants of the bilateral ease of commuting. For bilateral pairs with positive commuting áows, we partition the ease of commuting $( B _ { n i } )$ into four components: (i) a residence component $\left( \mathbb { B } _ { n } \right)$ , (ii) a workplace component $\left( { \mathbb { B } } _ { i } \right)$ , (iii) a component that is related to distance $( \mathrm { d i s t } _ { n i } ^ { - \phi } )$ , and (iv) an orthogonal component $\left( { \mathbb { B } } _ { n i } \right)$ 

$$
\log \mathcal {B} _ {n i} \equiv \log (B _ {n i} \kappa_ {n i} ^ {- \epsilon}) = \log \mathbb {B} _ {n} + \log \mathbb {B} _ {i} - \phi \log (\mathrm{dist} _ {n i}) + \log \mathbb {B} _ {n i}.\tag{B.73}
$$

We can always undertake this statistical decomposition of the ease of commuting $\left( \log \boldsymbol { B } _ { n i } \right)$ , where the error term $\left( \log \mathbb { B } _ { n i } \right)$ is orthogonal to distance by construction, because the reduced-form coe¢cient on log distance $\left( - \phi \right)$ captures any correlation of either log bilateral amenities $\left( \log B _ { n i } \right)$ and/or log bilateral commuting costs $\left( \log ( \kappa _ { n i } ^ { - \epsilon } ) \right)$ with log distance. For bilateral pairs with zero commuting, the model implies negligible amenities $( B _ { n i } \to 0 )$ and $1 / \mathrm { o r }$ prohibitive commuting costs $( \kappa _ { n i }  \infty ) .$ .<sup>4</sup> 

In the Örst step of our gravity equation estimation, we use this decomposition (B.73) and our expression for commuting áows (10) to estimate the reduced-form distance coe¢cient $\left( - \phi \right)$ 

$$
\log \lambda_ {n i} = g _ {0} + \eta_ {n} + \mu_ {i} - \phi \log \mathrm{dist} _ {n i} + \log \mathbb {B} _ {n i},\tag{B.74}
$$

where the residence Öxed e§ect $( \eta _ { n } )$ captures the consumption goods price index $\left( P _ { n } \right)$ , the price of residential land $\left( Q _ { n } \right)$ , and the residence component of the ease of commuting $\left( \mathbb { B } _ { n } \right)$ ; the workplace Öxed e§ect $( \mu _ { i } )$ captures the wage $( w _ { i } )$ and the workplace component of the ease of commuting $\left( \mathbb { B } _ { i } \right)$ ; the constant $g _ { 0 }$ captures the denominator of $\lambda _ { n i }$ and is separately identiÖed because we normalize the residence and workplace Öxed e§ects to sum to zero; and the error term $\left( \log \mathbb { B } _ { n i } \right)$ is orthogonal to log distance, because all e§ects of log distance on the composite ease of commuting are captured in the reduced-form distance coe¢cient $( - \phi ) . ^ { 5 }$ 

Estimating the gravity equation (B.74) for all bilateral pairs with positive commuters using OLS, we 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/03456822f1c8029a02e386302636731e23a6a9a289f9032c50142accc7a4bbf8.jpg)



Figure B.5: Gravity in Commuting Between Counties


Önd a regression R-squared of 0.80. In Figure B.5, we display the conditional relationship between log commuters and log distance, after removing residence and workplace Öxed e§ects from both log commuters and log distance. Consistent with the existing empirical literature on commuting, we Önd that the log linear functional form provides a good approximation to the data, with a tight and approximately linear relationship between the two variables, and an estimated coe¢cient on log distance of $- \phi \ : = \ : - 4 . 4 3$ This estimated coe¢cient is substantially larger than the corresponding coe¢cient for trade in goods of $- \left( \sigma - 1 \right) \psi = - 1 . 2 9$ , which is consistent with the view that transporting people is considerably more costly than transporting goods, in line with the substantial opportunity cost of time spent commuting. 

To identify the FrÈchet shape parameter (), the second step of our gravity equation estimation uses additional structure from the model, which implies that the workplace Öxed e§ects $\mu _ { i }$ depend on wages $( w _ { i } )$ and the workplace component of the ease of commuting $\left( \mathbb { B } _ { i } \right)$ : 

$$
\log \lambda_ {n i} = g _ {0} + \eta_ {n} + \epsilon \log w _ {i} - \phi \log \mathrm{dist} _ {n i} + \log u _ {n i},\tag{B.75}
$$

where the error term is given by log $u _ { n i } \equiv \log \mathbb { B } _ { i } + \log \mathbb { B } _ { n i }$ 

We estimate the gravity equation (B.75) imposing $\phi = 4 . 4 3 $ from our estimates above and identify  from the coe¢ cient on wages. Estimating (B.75) using OLS is potentially problematic, because workplace wages (w<sub>i</sub>) depend on the supply of commuters, which in turn depends on amenities that appear in the error term $\left( \log u _ { n i } \right)$ . Therefore we instrument log $w _ { i }$ with the log productivities log $A _ { i }$ that we recovered from the condition (16) equating income and expenditure above, using the fact that the model implies that productivity satisÖes the exclusion restriction of only a§ecting commuting áows through wages. Our Two-Stage-Least-Squares estimate of the FrÈchet shape parameter for the heterogeneity of worker preferences is $\epsilon = 3 . 3 0 . ^ { 6 }$ The tight Öt shown in Figure B.5 makes us conÖdent that our parametrization of the composite ease of commuting in terms of distance Öts the data quite well. 

## B.7 Openness of the Local Labor Market to Commuting

In this section of the web appendix, we use the commuter market clearing condition to derive reduced-form measures of the openness of the local labor market to commuting. We show that the share of residents who work where they live (the ìresidence own commuting share,î $\lambda _ { i i | i } ^ { R } )$ , the share of workers who live where they work (the ìworkplace own commuting share,î $\lambda _ { i i | i } ^ { L } )$ , and the ratio of workers to residents $( L _ { i } / R _ { i } )$ are all closely related to one another through the commuter market clearing condition. 

Re-writing the commuter market clearing condition in equation (13) in the paper, we obtain 

$$
L _ {i} = \lambda_ {i i | i} ^ {R} R _ {i} + \sum_ {n \neq i} \lambda_ {n i | n} ^ {R} R _ {n}.\tag{B.76}
$$

Rearranging this commuter market clearing condition, the importance of commuting from other locations as a source of employment for location i can be written as 

$$
\frac {\sum_ {n \neq i} \lambda_ {n i | n} ^ {R} R _ {n}}{L _ {i}} = 1 - \frac {\lambda_ {i i | i} ^ {R} R _ {i}}{L _ {i}}.\tag{B.77}
$$

We now use the deÖnition of the conditional commuting probabilities $( \lambda _ { n i | i } ^ { L }$ and $\lambda _ { n i \mid n } ^ { R } )$ in equations (B.11) and (B.12), which imply 

$$
\lambda_ {i i | i} ^ {R} = \frac {\lambda_ {i i}}{R _ {i} / \bar {L}} = \frac {L _ {i}}{R _ {i}} \frac {\lambda_ {i i}}{L _ {i} / \bar {L}} = \frac {L _ {i}}{R _ {i}} \lambda_ {i i | i} ^ {L}.\tag{B.78}
$$

Combining equations (B.77) and (B.78), we obtain: 

$$
\frac {\sum_ {n \neq i} \lambda_ {n i | n} ^ {R} R _ {n}}{L _ {i}} = 1 - \lambda_ {i i | i} ^ {L},\tag{B.79}
$$

where higher values of the workplace own commuting share $( \lambda _ { i i | i } ^ { L } )$ imply a local labor market that is more closed to commuting. Alternatively, the commuter market clearing condition can be written equivalently as 

$$
R _ {i} = \lambda_ {i i | i} ^ {L} L _ {i} + \sum_ {n \neq i} \lambda_ {n i | n} ^ {L} L _ {n}.\tag{B.80}
$$

Rearranging this expression, the importance of commuting from other locations as a source of residents for location i can be written as 

$$
\frac {\sum_ {n \neq i} \lambda_ {n i | n} ^ {L} L _ {n}}{R _ {i}} = 1 - \frac {\lambda_ {i i | i} ^ {L} L _ {i}}{R _ {i}}.\tag{B.81}
$$

Combining equations (B.77) and (B.81), we obtain 

$$
\frac {\sum_ {n \neq i} \lambda_ {n i | n} ^ {L} L _ {n}}{R _ {i}} = 1 - \lambda_ {i i | i} ^ {R},\tag{B.82}
$$

where higher values of the residence own commuting share $( \lambda _ { i i | i } ^ { R } )$ again imply a local labor market that is more closed to commuting. 

Together, the residence and workplace own commuting shares $( \lambda _ { i i | i } ^ { R }$ and $\lambda _ { i i | i } ^ { L }$ respectively) are su¢cient to recover the ratio of employment to residents $( L _ { i } / R _ { i } )$ . From equation (B.78), we have: 

$$
\frac {\lambda_ {i i | i} ^ {R}}{\lambda_ {i i | i} ^ {L}} = \frac {L _ {i}}{R _ {i}}.\tag{B.83}
$$

Therefore, knowing whether the minimum value of these two measures is equal to the residence commuting share $( \lambda _ { i i | i } ^ { R } )$ or the workplace commuting share $( \lambda _ { i i | i } ^ { L } )$ reveals whether a location is a net importer or net exporter of commuters: 

$$
\begin{array}{l l} \lambda_ {i i | i} ^ {R} > \lambda_ {i i | i} ^ {L}, & \Leftrightarrow \qquad L _ {i} > R _ {i}, \\ \lambda_ {i i | i} ^ {R} <   \lambda_ {i i | i} ^ {L}, & \Leftrightarrow \qquad L _ {i} <   R _ {i}. \end{array}\tag{B.84}
$$

We Önd that the residence and workplace own commuting shares $( \lambda _ { i i | i } ^ { R }$ and $\lambda _ { i i | i } ^ { L }$ respectively) are strongly positively correlated with one another, with a correlation of 0.60 from 2006-10 that is statistically signiÖcant at the 1 percent level. This positive correlation reáects in part the fact that gross commuting áows are large relative to net commuting áows, as explained by idiosyncratic preference draws in our model. We choose the residence own commuting share $( \lambda _ { i i | i } ^ { R } )$ as our baseline measure, because it is both model consistent and reported in the population census back to 1960. But we show that our results are robust to using either the residence or workplace commuting share or the average or minimum of these two measures. 

## B.8 Partial Equilibrium Elasticities

In this section of the web appendix, we use the model to derive partial equilibrium elasticities that capture the direct e§ect of a productivity shock on wages, employment and residents in the treated location, holding constant all other endogenous variables at their values in the initial equilibrium. Although these partial equilibrium elasticities do not incorporate the full set of interactions between locations that are captured in the general equilibrium elasticities in Figure 2 in the paper, we show in Section 4.1 of the paper that they explain some of the observed variation in these general equilibrium elasticities across locations. 

We now derive these partial equilibrium elasticities of the endogenous variables of the model with respect to a productivity shock. 

Wage Elasticity: Totally di§erentiating the goods market clearing condition in equation (7) in the paper, we have: 

$$
\begin{array}{r} \frac {d w _ {i}}{w _ {i}} w _ {i} L _ {i} + \frac {d L _ {i}}{L _ {i}} w _ {i} L _ {i} = \sum_ {r \in N} (1 - \pi_ {r i}) \pi_ {r i} \bar {v} _ {r} R _ {r} \frac {d L _ {i}}{L _ {i}} - \sum_ {r \in N} \sum_ {s \in N} \pi_ {r s} \pi_ {r i} \bar {v} _ {r} R _ {r} \frac {d L _ {s}}{L _ {s}} \\ - (\sigma - 1) \sum_ {r \in N} (1 - \pi_ {r i}) \pi_ {r i} \bar {v} _ {r} R _ {r} \frac {d w _ {n}}{w _ {n}} + (\sigma - 1) \sum_ {r \in N} \sum_ {s \in N} \pi_ {r s} \pi_ {r n} \bar {v} _ {r} R _ {r} \frac {d w _ {s}}{w _ {s}} \\ + (\sigma - 1) \sum_ {r \in N} (1 - \pi_ {r i}) \pi_ {r i} \bar {v} _ {r} R _ {r} \frac {d A _ {i}}{A _ {i}} - (\sigma - 1) \sum_ {r \in N} \sum_ {s \in N} \pi_ {r s} \pi_ {r i} \bar {v} _ {r} R _ {r} \frac {d A _ {s}}{A _ {s}} \\ + \sum_ {r \in N} \pi_ {r i} \bar {v} _ {r} R _ {r} \frac {d \bar {v} _ {r}}{\bar {v} _ {r}} + \sum_ {r \in N} \pi_ {r i} \bar {v} _ {r} R _ {r} \frac {d R _ {r}}{R _ {r}}. \end{array}
$$

To consider the direct e§ect of a productivity shock in location i on wages, employment and residents in that location, holding constant all other endogenous variables at their values in the initial equilibrium, we set $d A _ { s } = d w _ { s } = d L _ { s } = d R _ { s } = 0$ for $s \neq i$ and $d \bar { v } _ { r } = 0$ for all r, which yields: 

$$
\begin{array}{r} \frac {d w _ {i}}{w _ {i}} w _ {i} L _ {i} + \frac {d L _ {i}}{L _ {i}} w _ {i} L _ {i} = \sum_ {r \in N} (1 - \pi_ {r i}) \pi_ {r i} \bar {v} _ {r} R _ {r} \frac {d L _ {i}}{L _ {i}} - (\sigma - 1) \sum_ {r \in N} (1 - \pi_ {r i}) \pi_ {r i} \bar {v} _ {r} R _ {r} \frac {d w _ {i}}{w _ {i}} \\ + (\sigma - 1) \sum_ {r \in N} (1 - \pi_ {r i}) \pi_ {r i} \bar {v} _ {r} R _ {r} \frac {d A _ {n}}{A _ {n}} + \pi_ {i i} \bar {v} _ {i} R _ {i} \frac {d R _ {i}}{R _ {i}}. \end{array}
$$

This implies: 

$$
\begin{array}{r} \frac {d w _ {i}}{d A _ {i}} \frac {A _ {i}}{w _ {i}} + \frac {d L _ {i}}{d A _ {i}} \frac {A _ {i}}{L _ {i}} = \sum_ {r \in N} (1 - \pi_ {r i}) \frac {\pi_ {r i} \bar {v} _ {r} R _ {r}}{w _ {i} L _ {i}} (\frac {d L _ {i}}{d A _ {i}} \frac {A _ {i}}{L _ {i}}) - (\sigma - 1) \sum_ {r \in N} (1 - \pi_ {r i}) \frac {\pi_ {r i} \bar {v} _ {r} R _ {r}}{w _ {i} L _ {i}} (\frac {d w _ {i}}{d A _ {i}} \frac {A _ {i}}{w _ {i}}) \\ + (\sigma - 1) \sum_ {r \in N} (1 - \pi_ {r i}) \frac {\pi_ {r i} \bar {v} _ {r} R _ {r}}{w _ {i} L _ {i}} + \frac {\pi_ {i i} \bar {v} _ {i} R _ {i}}{w _ {i} L _ {i}} (\frac {d R _ {i}}{d A _ {i}} \frac {A _ {i}}{R _ {i}}), \end{array}
$$

which can be re-written as: 

$$
\begin{array}{r} \frac {d w _ {i}}{d A _ {i}} \frac {A _ {i}}{w _ {i}} + \left(\frac {d L _ {i}}{d w _ {i}} \frac {w _ {i}}{L _ {i}}\right) \left(\frac {d w _ {i}}{d A _ {i}} \frac {A _ {i}}{w _ {i}}\right) = \sum_ {r \in N} (1 - \pi_ {r i}) \frac {\pi_ {r i} \bar {v} _ {r} R _ {r}}{w _ {i} L _ {i}} \left(\frac {d L _ {i}}{d w _ {i}} \frac {w _ {i}}{L _ {i}}\right) \left(\frac {d w _ {i}}{d A _ {i}} \frac {A _ {i}}{w _ {i}}\right) \\ - (\sigma - 1) \sum_ {r \in N} (1 - \pi_ {r i}) \frac {\pi_ {r i} \bar {v} _ {r} R _ {r}}{w _ {i} L _ {i}} \left(\frac {d w _ {i}}{d A _ {i}} \frac {A _ {i}}{w _ {i}}\right) + (\sigma - 1) \sum_ {r \in N} (1 - \pi_ {r i}) \frac {\pi_ {r i} \bar {v} _ {r} R _ {r}}{w _ {i} L _ {i}} \\ + \frac {\pi_ {i i} \bar {v} _ {i} R _ {i}}{w _ {i} L _ {i}} \left(\frac {d R _ {i}}{d w _ {i}} \frac {w _ {i}}{R _ {i}}\right) \left(\frac {d w _ {i}}{d A _ {i}} \frac {A _ {i}}{w _ {i}}\right), \end{array}
$$

where we have used the fact that productivity does not directly enter the commuter market clearing condition in equation (13) in the paper and the residential choice probabilities in equation (11) in the paper, and hence employment and residents only change to the extent that wages change as a result of the productivity shock. Rearranging this expression, we obtain the following partial equilibrium elasticity: 

$$
\frac {\partial w _ {i}}{\partial A _ {i}} \frac {A _ {i}}{w _ {i}} = \frac {(\sigma - 1) \sum_ {r \in N} (1 - \pi_ {r i}) \xi_ {r i}}{\left[ 1 + (\sigma - 1) \sum_ {r \in N} (1 - \pi_ {r i}) \xi_ {r i} \right] + \left[ 1 - \sum_ {r \in N} (1 - \pi_ {r i}) \xi_ {r i} \right] \frac {d L _ {i}}{d w _ {i}} \frac {w _ {i}}{L _ {i}} - \xi_ {i i} \frac {d R _ {i}}{d w _ {i}} \frac {w _ {i}}{R _ {i}}},\tag{B.85}
$$

where $\xi _ { r i } = \pi _ { r i } \bar { v } _ { r } R _ { r } / w _ { i } L _ { i }$ is the share of location iís revenue from market r and we use the partial derivative symbol to clarify that this derivative is not the full general equilibrium one. 

Employment Elasticity: Totally di§erentiating the commuter market clearing condition in equation (13) in the paper, we have: 

$$
\begin{array}{r c l} \frac {d L _ {i}}{L _ {i}} & = & \epsilon \sum_ {r \in N} \left(1 - \lambda_ {r i | r} ^ {R}\right) \frac {d w _ {i}}{w _ {i}} \frac {\lambda_ {r i | r} ^ {R} R _ {r}}{L _ {n}} - \epsilon \sum_ {r \in N} \sum_ {s \neq n} \lambda_ {r s | r} ^ {R} \frac {d w _ {s}}{w _ {s}} \frac {L _ {r i}}{L _ {i}} \\ & & + \sum_ {r} \frac {d R _ {r}}{R _ {r}} \frac {L _ {r i}}{L _ {i}}. \end{array}
$$

To consider the direct e§ect of a productivity shock in location i on its employment and residents through a higher wage in that location, holding constant all other endogenous variables at their values in the initial equilibrium, we set $d w _ { s } = d L _ { s } = d R _ { s } = 0$ for $s \neq i$ , which yields: 

$$
\frac {d L _ {i}}{L _ {i}} = \epsilon \sum_ {r \in N} \left(1 - \lambda_ {r i | r} ^ {R}\right) \frac {\lambda_ {r i | r} ^ {R} R _ {r}}{L _ {i}} \frac {d w _ {i}}{w _ {i}} + \frac {\lambda_ {i i | i} ^ {R} R _ {i}}{L _ {i}} \frac {d R _ {i}}{R _ {i}}.
$$

Rearranging this expression, we obtain the following partial equilibrium elasticity: 

$$
\frac {\partial L _ {i}}{\partial w _ {i}} \frac {w _ {i}}{L _ {i}} = \epsilon \sum_ {r \in N} \left(1 - \lambda_ {r i | r} ^ {R}\right) \vartheta_ {r i} + \vartheta_ {i i} \left(\frac {d R _ {i}}{d w _ {i}} \frac {w _ {i}}{R _ {i}}\right),\tag{B.86}
$$

where $\vartheta _ { r i } = \lambda _ { r i | r } ^ { R } R _ { r } / L _ { i }$ is the share of commuters from residence r in workplace iís employment and we use the partial derivative symbol to clarify that this derivative is not the full general equilibrium one. Residents Elasticity: Totally di§erentiating the residential choice probability $( \lambda _ { n } ^ { R }$ in equation (11) in the paper), we have: 

$$
\begin{array}{r c l} \frac {d R _ {i}}{R _ {i}} \frac {R _ {i}}{\bar {L}} & = & - \epsilon \alpha (1 - \lambda_ {i} ^ {R}) \lambda_ {i} ^ {R} \frac {d P _ {i}}{P _ {i}} + \epsilon \alpha \sum_ {r \neq i} \lambda_ {r} ^ {R} \lambda_ {i} ^ {R} \frac {d P _ {r}}{P _ {r}} \\ & & - \epsilon (1 - \alpha) (1 - \lambda_ {i} ^ {R}) \lambda_ {i} ^ {R} \frac {d Q _ {i}}{Q _ {i}} + \epsilon (1 - \alpha) \sum_ {r \neq i} \lambda_ {r} ^ {R} \lambda_ {i} ^ {R} \frac {d Q _ {r}}{Q _ {r}} \\ & & + \epsilon \lambda_ {i i} \frac {d w _ {i}}{w _ {i}} - \epsilon \lambda_ {i} ^ {L} \lambda_ {i} ^ {R} \frac {d w _ {i}}{w _ {i}} - \epsilon \sum_ {s \neq i} \lambda_ {s} ^ {L} \lambda_ {i} ^ {R} \frac {d w _ {s}}{w _ {s}}. \end{array}
$$

To consider the direct e§ect of a productivity shock in location i on its residents through a higher wage in that location, holding constant all other endogenous variables at their values in the initial equilibrium, we set $\partial P _ { r } = \partial Q _ { r } = 0$ for all r and $\partial w _ { s } = 0$ for $s \neq i$ , which yields: 

$$
\frac {\partial R _ {i}}{R _ {i}} \frac {R _ {i}}{\bar {L}} = \epsilon (\lambda_ {i i} - \lambda_ {i} ^ {L} \lambda_ {i} ^ {R}) \frac {\partial w _ {i}}{w _ {i}}.
$$

This implies the following partial equilibrium elasticity: 

$$
\frac {\partial R _ {i}}{\partial w _ {i}} \frac {w _ {i}}{R _ {i}} = \epsilon \left(\frac {\lambda_ {i i}}{\lambda_ {i} ^ {R}} - \lambda_ {i} ^ {L}\right),\tag{B.87}
$$

where we use the partial derivative symbol to clarify that this derivative is not the full general equilibrium elasticity. Using the residents elasticity (B.87) in the employment elasticity (B.86), and using the residents and employment elasticities ((B.87) and (B.86) respectively) in the wage elasticity (B.85), we obtain the following partial equilibrium elasticities for the productivity shock, 

$$
\begin{array}{l} \frac {\partial R _ {i}}{\partial w _ {i}} \frac {w _ {i}}{R _ {i}} = \epsilon \left(\lambda_ {r s | r} ^ {R} \frac {\lambda_ {i i}}{\lambda_ {i} ^ {R}} - \lambda_ {i} ^ {L}\right), \\ \frac {\partial L _ {i}}{\partial w _ {i}} \frac {w _ {i}}{L _ {i}} = \epsilon \sum_ {r \in N} \left(1 - \lambda_ {r i | r} ^ {R}\right) \vartheta_ {r i} + \vartheta_ {i i} \epsilon \left(\frac {\lambda_ {i i}}{\lambda_ {i} ^ {R}} - \lambda_ {i} ^ {L}\right), \\ \frac {\partial w _ {i}}{\partial A _ {i}} \frac {A _ {i}}{w _ {i}} = \frac {(\sigma - 1) \sum_ {r \in N} (1 - \pi_ {r i}) \xi_ {r i}}{\left[ 1 + (\sigma - 1) \sum_ {r \in N} (1 - \pi_ {r i}) \xi_ {r i} \right] + \left[ 1 - \sum_ {r \in N} (1 - \pi_ {r i}) \xi_ {r i} \right] \left[ \epsilon \sum_ {r \in N} \left(1 - \lambda_ {r i | r} ^ {R}\right) \vartheta_ {r i} + \epsilon \vartheta_ {i i} \left(\frac {\lambda_ {i i}}{\lambda_ {i} ^ {R}} - \lambda_ {i} ^ {L}\right) \right] - \xi_ {i i} \epsilon \left(\frac {\lambda_ {i i}}{\lambda_ {i} ^ {R}} - \lambda_ {i} ^ {L}\right)}. \end{array}
$$

## B.9 Gravity and Local Employment Elasticities

We now show that the class of models consistent with a gravity equation for commuting implies heterogeneous local employment elasticities. Assume that commuting áows satisfy the following gravity equation: 

$$
L _ {n i} = \mathcal {R} _ {n} \mathcal {B} _ {n i} \mathcal {W} _ {i},\tag{B.88}
$$

where $L _ { n i }$ are commuting áows from residence n to workplace $i ; \mathcal { R } _ { n }$ is a residence Öxed e§ect; $\mathcal { W } _ { i }$ is a workplace Öxed e§ect; and $B _ { n i }$ is a measure of the ease of commuting (an inverse measure of bilateral commuting costs). This gravity equation (B.88) implies that the unconditional probability that a worker commutes from residence n to workplace i is: 

$$
\lambda_ {n i} = \frac {L _ {n i}}{\sum_ {r \in N} \sum_ {s \in N} L _ {r s}} = \frac {\mathcal {R} _ {n} \mathcal {B} _ {n i} \mathcal {W} _ {i}}{\sum_ {r \in N} \sum_ {s \in N} \mathcal {R} _ {r} \mathcal {B} _ {r s} \mathcal {W} _ {s}}.\tag{B.89}
$$

The corresponding probability of working in location i is: 

$$
\lambda_ {i} ^ {L} = \frac {\sum_ {r \in N} L _ {r i}}{\sum_ {r \in N} \sum_ {s \in N} L _ {r s}} = \frac {\sum_ {r \in N} \mathcal {R} _ {r} \mathcal {B} _ {r i} \mathcal {W} _ {i}}{\sum_ {r \in N} \sum_ {s \in N} \mathcal {R} _ {r} \mathcal {B} _ {r s} \mathcal {W} _ {s}},\tag{B.90}
$$

and the probability of residing in location n is: 

$$
\lambda_ {n} ^ {R} = \frac {\sum_ {s \in N} L _ {n s}}{\sum_ {r \in N} \sum_ {s \in N} L _ {r s}} = \frac {\sum_ {s \in N} \mathcal {R} _ {n} \mathcal {B} _ {n s} \mathcal {W} _ {s}}{\sum_ {r \in N} \sum_ {s \in N} \mathcal {R} _ {r} \mathcal {B} _ {r s} \mathcal {W} _ {s}},\tag{B.91}
$$

From equations (B.89) and (B.91), the probability of commuting from residence n to workplace i conditional on residing in n is: 

$$
\lambda_ {n i | n} ^ {R} = \frac {\lambda_ {n i}}{\lambda_ {n} ^ {R}} = \frac {\mathcal {R} _ {n} \mathcal {B} _ {n i} \mathcal {W} _ {i}}{\sum_ {s \in N} \mathcal {R} _ {n} \mathcal {B} _ {n s} \mathcal {W} _ {s}} = \frac {\mathcal {B} _ {n i} \mathcal {W} _ {i}}{\sum_ {s \in N} \mathcal {B} _ {n s} \mathcal {W} _ {s}}.\tag{B.92}
$$

Using this conditional probability (B.92), the commuter market clearing condition can be written as: 

$$
L _ {i} = \sum_ {n \in N} \lambda_ {n i | n} ^ {R} R _ {n} = \sum_ {n \in N} \frac {\mathcal {B} _ {n i} \mathcal {W} _ {i}}{\sum_ {s \in N} \mathcal {B} _ {n s} \mathcal {W} _ {s}} R _ {n}.\tag{B.93}
$$

Totally di§erentiating this commuter market clearing condition (B.93) for a given commuting technology $B _ { n i }$ , we have: 

$$
\begin{array}{r c l} \frac {d L _ {i}}{L _ {i}} & = & \sum_ {r \in N} \left(1 - \lambda_ {r i | r} ^ {R}\right) \frac {d \mathcal {W} _ {i}}{\mathcal {W} _ {i}} \frac {\lambda_ {r i | r} ^ {R} R _ {r}}{L _ {i}} - \sum_ {r \in N} \sum_ {s \neq i} \lambda_ {r s | r} ^ {R} \frac {d \mathcal {W} _ {s}}{\mathcal {W} _ {s}} \frac {\lambda_ {r i | r} ^ {R} R _ {r}}{L _ {i}} \\ & & + \sum_ {r \in N} \frac {d R _ {r}}{R _ {r}} \frac {\lambda_ {r i | r} ^ {R} R _ {r}}{L _ {i}}. \end{array}\tag{B.94}
$$

(B.95) 

Now consider the direct e§ect of a shock to the workplace Öxed e§ect for location $i \ ( \partial \mathcal { W } _ { i } \neq 0 )$ evaluated at the values of the variables for all other locations from the initial equilibrium $( \partial \mathcal { W } _ { r } = \partial L _ { r } = \partial R _ { r } = 0$ for $r \neq i )$ 

$$
\frac {\partial L _ {i}}{L _ {i}} = \sum_ {r \in N} \left(1 - \lambda_ {r i | r} ^ {R}\right) \frac {\lambda_ {r i | r} ^ {R} R _ {r}}{L _ {i}} \frac {\partial \mathcal {W} _ {i}}{\mathcal {W} _ {i}} + \frac {\lambda_ {i i | i} ^ {R} R _ {i}}{L _ {i}} \frac {\partial R _ {i}}{R _ {i}}.\tag{B.96}
$$

Rearranging this expression, we obtain the following partial equilibrium local employment elasticity: 

$$
\frac {\partial L _ {i}}{\partial \mathcal {W} _ {i}} \frac {\mathcal {W} _ {i}}{L _ {i}} = \underbrace {\sum_ {r \in N} \left(1 - \lambda_ {r i | r} ^ {R}\right) \vartheta_ {r i}} _ {\text {commuting}} + \underbrace {\vartheta_ {i i} \left(\frac {\partial R _ {i}}{\partial \mathcal {W} _ {i}} \frac {\mathcal {W} _ {i}}{R _ {i}}\right)} _ {\text {migration}},\tag{B.97}
$$

where $\vartheta _ { r i } = \lambda _ { r i | r } ^ { R } R _ { r } / L _ { i }$ is the share of commuters from residence r in workplace iís employment and we use the partial derivative symbol to clarify that this derivative is not the full general equilibrium one. The Örst term on the right-hand side of equation (B.97) captures the impact of the shock to the workplace Öxed e§ect ( <sub>i</sub>) on employment in location i through commuting. The second term on the right-hand side captures its impact on employment in location i through migration. 

This partial equilibrium local employment elasticity (B.97) takes the same form as in the paper (and in the previous section of this web appendix above), where in our model the shock to the workplace Öxed e§ect for location i ( ) corresponds to a shock to the wage at that workplace, which in turn depends on the shock to productivity at that workplace. Therefore our result of a variable local employment elasticit that depends on access to commuters in surrounding locations is a generic feature of the class of models that are consistent with a gravity equation for commuting. We show in the main paper that observed commuting áows are characterized by a strong gravity equation relationship. 

To show empirically that the heterogeneity in local employment elasticities is a generic implication of the gravity equation, we compute the Örst term on the right-hand side of equation (B.97) that captures commuting. This Örst term depends solely on observed variables in the initial equilibrium: (i) the probability of commuting to workplace i conditional on living in residence r and (ii) the share of commuters from residence r in workplace iís employment. In Figure B.6, we show the estimated kernel density of the commuting component of the partial employment elasticity (black line) for counties, and the 95 percent conÖdence intervals (gray shading). As apparent from comparing Figure B.6 to Figure 2 in the paper, the heterogeneity in county local employment elasticities largely reáects the heterogeneity in this Örst commuting term, as conÖrmed in the regressions in Table 2 in the paper. Figure B.7 shows the same commuting component of the partial employment elasticity, but for CZs rather than counties. Comparing Figure B.7 to Figure C.13 in this web appendix, we Önd that the heterogeneity in CZ local employment elasticities also largely reáects the heterogeneity in this Örst commuting term. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/214ac957ce3bef71ee9a5df41d49d42586d45ef39611d16b02627f570ab5b22c.jpg)



Figure B.6: Commuting component of partial employment elasticity for counties


![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/bad367f2b8b34bfd906bd75e0ed6987e0e5bae99228027feeb0304ee031d530f.jpg)



Figure B.7: Commuting component of partial employment elasticity for commuting zones (CZs)


## B.10 Commuting with Multiple Worker Types

In this section of the web appendix, we consider a generalization of our model to allow for multiple worker types, which di§er in their valuation of amenities and the variance of their idiosyncratic preferences. These di§erences in variance in turn imply that the multiple types di§er in the responsiveness of their migration and commuting decisions to economic characteristics of locations (such as wages). This extension of ou FrÈchet model to multiple worker types is analogous to the extension of the logit model to multiple types in the mixed logit model (see for example McFadden and Train 2000), which is in turn closely related to the random coe¢cients model of Berry, Levinsohn and Pakes (1995). We show that our prediction of heterogeneous local employment elasticities across locations is robust to this extension and that there is now an additional source of heterogeneity relative to our baseline speciÖcation. 

In particular, suppose that there are multiple types of workers (e.g. skilled versus unskilled) indexed by $z = 1 , \dots , Z$ . There is a separate labor market and a separate wage for each type of worker z in each workplace i $( w _ { i } ^ { z } )$ . Workers of a given type have idiosyncratic preferences over workplace and residence locations. However, the distributions of these idiosyncratic preferences di§er across types, in terms of both their average preferences for the amenities for each bilateral commute (as determined by $B _ { n i } ^ { z } )$ and the variance of their idiosyncratic preferences across these bilateral commutes (as determined by $\epsilon ^ { z } )$ : 

$$
G _ {n i} ^ {z} (b) = e ^ {- B _ {n i} ^ {z} b ^ {- \epsilon^ {z}}}.\tag{B.98}
$$

## B.10.1 Commuting Decisions for Each Worker Type

Under these assumptions, commuting decisions for each worker type are characterized by a gravity equation, which is analogous to that in our baseline speciÖcation with a single worker type. The probability that workers of type z choose to work in location i conditional on living in location n is: 

$$
\pi_ {n i | n} ^ {z} = \frac {B _ {n i} ^ {z} (w _ {i} ^ {z} / \kappa_ {n i} ^ {z}) ^ {\epsilon^ {z}}}{\sum_ {s \in N} B _ {n s} ^ {z} (w _ {s} ^ {z} / \kappa_ {n s} ^ {z}) ^ {\epsilon^ {z}}}.\tag{B.99}
$$

The corresponding commuter market clearing condition for workers of type z is: 

$$
L _ {i} ^ {z} = \sum_ {r \in N} \frac {B _ {r i} ^ {z} (w _ {i} ^ {z} / \kappa_ {r i} ^ {z}) ^ {\epsilon^ {z}}}{\sum_ {s \in N} B _ {r s} ^ {z} (w _ {s} ^ {z} / \kappa_ {r s} ^ {z}) ^ {\epsilon^ {z}}} R _ {r} ^ {z},\tag{B.100}
$$

which yields a partial elasticity of employment for workers of type z with respect to their wage that takes a similar form as for our baseline speciÖcation with a single worker type: 

$$
\frac {\partial L _ {i} ^ {z}}{\partial w _ {i} ^ {z}} \frac {w _ {i} ^ {z}}{L _ {i} ^ {z}} = \epsilon^ {z} \sum_ {r \in N} \left(1 - \lambda_ {r i | r} ^ {R z}\right) \vartheta_ {r i} ^ {z} + \vartheta_ {i i} ^ {z} \left(\frac {\partial R _ {i} ^ {z}}{\partial w _ {i} ^ {z}} \frac {w _ {i} ^ {z}}{R _ {i} ^ {z}}\right),\tag{B.101}
$$

where $\vartheta _ { r i } ^ { z } = \lambda _ { r i | r } ^ { R z } R _ { r } ^ { z } / L _ { i } ^ { z }$ is the share of commuters from residence r in workplace iís employment for workers of type z. 

## B.10.2 Aggregate Commuting Decisions

Aggregating commuting decisions across worker types, the total number of workers that choose to work in location i is: 

$$
L _ {i} = \sum_ {z = 1} ^ {Z} L _ {i} ^ {z}.\tag{B.102}
$$

Now consider the elasticity of total employment in location $\textit { i } \left( L _ { i } \right)$ with respect to a common increase in the wages of all worker types in that location: 

$$
d w _ {i} ^ {z} = d w _ {i} ^ {k} = d w _ {i} > 0, \quad \forall z, k.\tag{B.103}
$$

Di§erentiating with respect to wages in equation (B.102), we have: 

$$
d L _ {i} = \sum_ {z = 1} ^ {Z} \frac {\partial L _ {i} ^ {z}}{\partial w _ {i} ^ {z}} d w _ {i} ^ {z},\tag{B.104}
$$

which for a common change in wages in equation (B.103) can be re-written as: 

$$
\frac {d L _ {i}}{d w _ {i}} = \sum_ {z = 1} ^ {Z} \frac {\partial L _ {i} ^ {z}}{\partial w _ {i} ^ {z}},\tag{B.105}
$$

which can be further re-written as: 

$$
\frac {d L _ {i}}{d w _ {i}} \frac {w _ {i}}{L _ {i}} = \sum_ {z = 1} ^ {Z} \left(\frac {\partial L _ {i} ^ {z}}{\partial w _ {i} ^ {z}} \frac {w _ {i} ^ {z}}{L _ {i} ^ {z}}\right) \left(\frac {L _ {i} ^ {z} / L _ {i}}{w _ {i} ^ {z} / w _ {i}}\right).\tag{B.106}
$$

Combining equations (B.101) and (B.106), the local employment elasticity for each location is a weighted average of the local employment elasticities for each worker type for that location, where the weights depend on employment shares and relative wages. Therefore, local employment elasticities continue to be heterogeneous across locations in this extension of the model to incorporate multiple worker types, but there is now an additional source of heterogeneity relative to our baseline speciÖcation. First, the local employment elasticity for a given worker type is heterogeneous across locations depending on commuting networks for that worker type (equation (B.101)). This Örst source of heterogeneity is analogous to that in our baseline speciÖcation with a single worker type. Second, the composition of worker types and their relative wages can di§er across locations, which provides an additional source of heterogeneity in local employment elasticities that is not present in our baseline speciÖcation (as in equation (B.106)). Taken together, this extension further reinforces our point that the local employment elasticity is not a structural parameter. 

## B.11 Congestion in Commuting

In this section of the web appendix, we generalize our baseline speciÖcation to allow for congestion in commuting. Assuming that congestion costs are a power function of the volume of commuters, we show that congestion a§ects the interpretation of the estimated parameters in our commuting gravity equation, but leaves the modelís prediction of heterogeneous employment elasticities across locations unchanged. In particular, we assume that each worker draws idiosyncratic preferences for each pair of residence n and workplace i from the following distribution: 

$$
G _ {n i} (b) = e ^ {- B _ {n i} L _ {n i} ^ {\chi} b ^ {- \epsilon}},\tag{B.107}
$$

where the scale parameter of this distribution $( B _ { n i } L _ { n i } ^ { \chi } )$ is a power function of the volume of commuters. Our baseline speciÖcation corresponds to the special case in which $\chi = 0 ; \chi < 0$ corresponds to congestion in commuting decisions, such that the attractiveness of commuting from residence n to workplace i depends negatively on the volume of commuters. Under these assumptions, the probability that a worker commutes from residence n to workplace i is: 

$$
\lambda_ {n i} = \frac {L _ {n i}}{\bar {L}} = \frac {B _ {n i} L _ {n i} ^ {\chi} \left(\kappa_ {n i} P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}\right) ^ {- \epsilon} w _ {i} ^ {\epsilon}}{\sum_ {r \in N} \sum_ {s \in N} B _ {r s} L _ {r s} ^ {\chi} \left(\kappa_ {r s} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}\right) ^ {- \epsilon} w _ {s} ^ {\epsilon}},\tag{B.108}
$$

and expected utility conditional on choosing a given bilateral commute (which is the across all bilateral commutes) is equal to: 

$$
\bar {U} = \mathbb {E} \left[ U _ {n i \omega} \right] = \Gamma \left(\frac {\epsilon - 1}{\epsilon}\right) \left[ \sum_ {r \in N} \sum_ {s \in N} B _ {r s} L _ {r s} ^ {\chi} \left(\kappa_ {r s} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}\right) ^ {- \epsilon} w _ {s} ^ {\epsilon} \right] ^ {\frac {1}{\epsilon}} \mathrm{all} n, i \in N.\tag{B.109}
$$

Combining equations (B.108) and (B.109), the áow of workers that choose to commute from residence n to workplace i can be written as: 

$$
L _ {n i} = \left(\frac {\bar {U}}{\Gamma}\right) ^ {- \epsilon} B _ {n i} L _ {n i} ^ {\chi} \left(\kappa_ {n i} P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}\right) ^ {- \epsilon} w _ {i} ^ {\epsilon} \bar {L},\tag{B.110}
$$

which can be in turn re-written as: 

$$
L _ {n i} = \left(\frac {\bar {U}}{\Gamma}\right) ^ {- \frac {\epsilon}{1 - \chi}} B _ {n i} ^ {\frac {1}{1 - \chi}} \left(\kappa_ {n i} P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}\right) ^ {- \frac {\epsilon}{1 - \chi}} w _ {i} ^ {\frac {\epsilon}{1 - \chi}} \bar {L} ^ {\frac {1}{1 - \chi}}.\tag{B.111}
$$

Dividing equation (B.111) by its sum across all bilateral pairs, the probability that a worker commutes from residence n to workplace i can be equivalently expressed as: 

$$
\lambda_ {n i} = \frac {L _ {n i}}{\sum_ {r \in N} \sum_ {s \in N} L _ {r s}} = \frac {L _ {n i}}{\bar {L}} = \frac {B _ {n i} ^ {\frac {1}{1 - \chi}} \left(\kappa_ {n i} P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}\right) ^ {- \frac {\epsilon}{1 - \chi}} w _ {i} ^ {\frac {\epsilon}{1 - \chi}}}{\sum_ {r \in N} \sum_ {s \in N} B _ {r s} ^ {\frac {1}{1 - \chi}} \left(\kappa_ {r s} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}\right) ^ {- \frac {\epsilon}{1 - \chi}} w _ {s} ^ {\frac {\epsilon}{1 - \chi}}},\tag{B.112}
$$

which takes exactly the same form as in our baseline speciÖcation, except that the exponent on wages, which we interpret as  in our baseline speciÖcation, should be interpreted as $\epsilon / ( 1 - \chi )$ in this extended speciÖcation. Similarly, the exponents on commuting costs $\left( \kappa _ { n i } \right)$ , consumption goods price indices $\left( P _ { n } \right)$ and land prices $\left( Q _ { n } \right)$ are all now multiplied by $1 / ( 1 - \chi )$ . Finally, the values of bilateral amenities implied by this commuting probability, which we interpret as $B _ { n i }$ in our baseline speciÖcation, should be interpreted as $B _ { n i } ^ { 1 / ( 1 - \chi ) }$ in this extended speciÖcation. 

Using the unconditional commuting probabilities (B.112), we can also solve for the probability of commuting to workplace i conditional on living in residence n: 

$$
\lambda_ {n i | n} ^ {R} = \frac {B _ {n i} ^ {\frac {1}{1 - \chi}} (w _ {i} / \kappa_ {n i}) ^ {\frac {\epsilon}{1 - \chi}}}{\sum_ {s \in N} B _ {n s} ^ {\frac {1}{1 - \chi}} (w _ {s} / \kappa_ {n s}) ^ {\frac {\epsilon}{1 - \chi}}}.\tag{B.113}
$$

The corresponding commuter market clearing condition is: 

$$
L _ {i} = \sum_ {r \in N} \frac {B _ {r i} ^ {\frac {1}{1 - \chi}} (w _ {i} / \kappa_ {r i}) ^ {\frac {\epsilon}{1 - \chi}}}{\sum_ {s \in N} B _ {r s} ^ {\frac {1}{1 - \chi}} (w _ {s} / \kappa_ {r s}) ^ {\frac {\epsilon}{1 - \chi}}} R _ {r},\tag{B.114}
$$

which yields a partial elasticity of employment with respect to the wage that takes a similar form as in our baseline speciÖcation: 

$$
\frac {\partial L _ {i}}{\partial w _ {i}} \frac {w _ {i}}{L _ {i}} = \frac {\epsilon}{1 - \chi} \sum_ {r \in N} \left(1 - \lambda_ {r i | r} ^ {R}\right) \vartheta_ {r i} + \vartheta_ {i i} \left(\frac {\partial R _ {i}}{\partial w _ {i}} \frac {w _ {i}}{R _ {i}}\right),\tag{B.115}
$$

where $\vartheta _ { r i } = \lambda _ { r i | r } ^ { R } R _ { r } / L _ { i }$ is the share of commuters from residence r in workplace iís employment. In this extended speciÖcation (B.115), the estimated coe¢cient on the Örst term on the right-hand side is again the exponent on wages from the gravity equation for commuting (B.112), but this estimated coe¢ cient is now interpreted as $\epsilon / ( 1 - \chi )$ rather than as . 

Therefore, taking the results of this section together, the introduction of congestion costs that are a power function of the volume of commuters a§ects the interpretation of the estimated parameters in our gravity equation for commuting, but leaves the modelís prediction of heterogeneous elasticities of employment with respect to wages across locations unchanged. 

## B.12 Non-traded Goods

In the baseline version of the model in the paper, we introduce commuting into a canonical new economic geography model with a single tradable consumption goods sector and land as the only non-traded good. We focus on the implications of introducing commuting into this canonical model for the elasticity of local employment with respect to local labor demand shocks. In this section of the web appendix, we generalize our analysis to incorporate non-traded consumption goods. We show that the commuter market clearing condition and the local elasticity of employment with respect to wages take the same form as in our baseline speciÖcation without non-traded goods. 

The consumption index for worker ! residing at location n and working at location i is now assumed to take the following form: 

$$
U _ {n i \omega} = \frac {b _ {n i \omega}}{\kappa_ {n i}} \left(\frac {C _ {N n \omega}}{\alpha_ {N}}\right) ^ {\alpha_ {N}} \left(\frac {C _ {T n \omega}}{\alpha_ {T}}\right) ^ {\alpha_ {T}} \left(\frac {H _ {n \omega}}{1 - \alpha_ {N} - \alpha_ {T}}\right) ^ {1 - \alpha_ {N} - \alpha_ {T}},\tag{B.116}
$$

$$
\alpha_ {N}, \alpha_ {T} > 0, \qquad 0 <   \alpha_ {N} + \alpha_ {T} <   1,
$$

where $C _ { N n \omega }$ is consumption of the non-traded good; $C _ { T n \omega }$ is consumption of the traded good; and all other terms are deÖned in the same way as in our baseline speciÖcation. As in our baseline speciÖcation, land is owned by immobile landlords, who receive worker expenditure on residential land as income, and consume only goods where they live. Therefore, total expenditure on consumption goods (traded plus non-traded) equals the fraction $\alpha ^ { N } + \alpha ^ { T }$ of the total income of residents plus the entire income of landlords (which equals the fraction $1 - \alpha ^ { N } - \alpha ^ { T } )$ of the total income of residents): 

$$
P _ {n} C _ {n} = (\alpha_ {N} + \alpha_ {T}) \bar {v} _ {n} R _ {n} + (1 - \alpha_ {N} - \alpha_ {T}) \bar {v} _ {n} R _ {n} = \bar {v} _ {n} R _ {n}.\tag{B.117}
$$

Utility maximization implies that a constant fraction $\alpha _ { N } / ( \alpha _ { N } + \alpha _ { T } )$ of total expenditure on consumption goods is allocated to the non-traded sector: 

$$
P _ {N n} C _ {N n} = \frac {\alpha_ {N}}{\alpha_ {N} + \alpha_ {T}} P _ {n} C _ {n} = \frac {\alpha_ {N}}{\alpha_ {N} + \alpha_ {T}} \bar {v} _ {n} R _ {n},\tag{B.118}
$$

and the remaining fraction is allocated to the traded sector: 

$$
P _ {T n} C _ {T n} = \frac {\alpha_ {T}}{\alpha_ {N} + \alpha_ {T}} P _ {n} C _ {n} = \frac {\alpha_ {T}}{\alpha_ {N} + \alpha_ {T}} \bar {v} _ {n} R _ {n},\tag{B.119}
$$

The non-traded good is assumed to be produced under conditions of perfect competition and according to a constant returns to scale technology with a unit labor requirement: 

$$
Y _ {N n} = L _ {N n},\tag{B.120}
$$

where $Y _ { N n }$ is output of the non-traded good in location n and $L _ { N n }$ is employment in the non-traded sector in that location. Perfect competition and constant returns to scale imply that the price of the non-traded good is equal to the wage: 

$$
P _ {N n} = w _ {n}.\tag{B.121}
$$

Combining this result with utility maximization (B.118), and using goods market clearing for the nontraded good $( C _ { N n } ~ = ~ Y _ { N n } )$ and the production technology (B.120), we Önd that the wage bill in the non-traded sector is a constant share of residential income: 

$$
w _ {n} L _ {N n} = \frac {\alpha_ {N}}{\alpha_ {N} + \alpha_ {T}} \bar {v} _ {n} R _ {n}.\tag{B.122}
$$

Using utility maximization and goods market clearing for tradeables, the wage bill in the traded sector is fraction of residential income across all locations: 

$$
w _ {n} L _ {T n} = \frac {\alpha_ {T}}{\alpha_ {N} + \alpha_ {T}} \sum_ {r \in N} \pi_ {r n} \bar {v} _ {r} R _ {r}.\tag{B.123}
$$

Total employment equals the sum of employment in the non-traded and traded sectors: 

$$
L _ {n} = L _ {T n} + L _ {N n} = \frac {\alpha_ {N}}{\alpha_ {N} + \alpha_ {T}} \frac {\bar {v} _ {n} R _ {n}}{w _ {n}} + \frac {\alpha_ {T}}{\alpha_ {N} + \alpha_ {T}} \sum_ {r \in N} \frac {\pi_ {r n} \bar {v} _ {r} R _ {r}}{w _ {n}}.\tag{B.124}
$$

The commuter market clearing condition requires that total employment in each location equals the measure of workers that choose to commute to that location and takes the same form as in our baseline speciÖcation without the non-traded sector: 

$$
L _ {n} = \sum_ {r \in N} \frac {B _ {r n} (w _ {n} / \kappa_ {r n}) ^ {\epsilon}}{\sum_ {s \in N} B _ {r s} (w _ {s} / \kappa_ {r s}) ^ {\epsilon}} R _ {r}.\tag{B.125}
$$

Given the same commuter market clearing condition, the partial elasticity of employment with respect to the wage takes the same form as in our baseline speciÖcation: 

$$
\frac {\partial L _ {n}}{\partial w _ {n}} \frac {w _ {n}}{L _ {n}} = \epsilon \sum_ {r \in N} \left(1 - \lambda_ {r n | r} ^ {R}\right) \vartheta_ {r n} + \vartheta_ {n n} \left(\frac {\partial R _ {n}}{\partial w _ {n}} \frac {w _ {n}}{R _ {n}}\right),\tag{B.126}
$$

where $\vartheta _ { r n } = \lambda _ { r n | r } ^ { R } R _ { r } / L _ { n }$ is the share of commuters from residence r in workplace nís employment. 

Therefore, although the presence of non-traded goods can a§ect the elasticity of wages with respect to productivity, it leaves unchanged the modelís prediction of heterogeneous local employment elasticities with respect to wages. Intuitively, when deciding where to work, workers care about the wage, and not whether this wage is paid in the traded or non-traded sector. Therefore, the gravity equation for commuting takes the same form as in our baseline speciÖcation without the non-traded sector, which in turn implies that the elasticity of local employment with respect to wages takes the same form as in our baseline speciÖcation without the non-traded sector. 

## B.13 Landlords Consume Residential Land

In this subsection of the web appendix, we show that allowing landlords to consume residential land in addition to consumption goods is straightforward, and merely results in less elegant expressions. Under this alternative assumption, consumption goods expenditure that was previously given by equation (4) in the paper is now instead given by: 

$$
P _ {n} C _ {n} = \alpha [ 1 + (1 - \alpha) ] \bar {v} _ {n} R _ {n}.\tag{B.127}
$$

Using this relationship, the equality between income and expenditure that was previously given by equation (7) in the paper is now instead given by: 

$$
w _ {i} L _ {i} = \alpha [ 1 + (1 - \alpha) ] \sum_ {n \in N} \pi_ {n i} \bar {v} _ {n} R _ {n},\tag{B.128}
$$

and the land market clearing condition that was previously given by equation (5) in the paper is now instead given by: 

$$
Q _ {n} = (1 - \alpha) [ 1 + (1 - \alpha) ] \frac {\bar {v} _ {n} R _ {n}}{H _ {n}}.\tag{B.129}
$$

As in our baseline speciÖcation in which landlords consume only consumption goods, the general equilibrium of the model can be referenced by the following vector of six variables $\{ w _ { n } , \ \bar { v } _ { n } , \ Q _ { n } , \ L _ { n } , \ R _ { n } , \ P _ { n } \} _ { n = 1 } ^ { N }$ and a scalar U<sup></sup> . Given this equilibrium vector and scalar, all other endogenous variables of the model can be determined. This equilibrium vector solves the following six sets of equations: income equals expenditure (B.128), land market clearing (B.129), expected labor income (which remains as in equation (14) in the paper), workplace choice probabilities (which continue to equal equation (11) in the paper for $L _ { n } )$ , residence choice probabilities (which are still equal to equation (11) in the paper for ${ \cal R } _ { n } )$ , price indices (again equal to equation (8) in the paper), and the labor market clearing condition (which remains the same as $\begin{array} { r } { \bar { L } = \sum _ { n \in N } R _ { n } = \sum _ { n \in N } L _ { n } ) } \end{array}$ . This system of equations for general equilibrium is exactly the same as in our baseline speciÖcation in which landlords consume only consumption goods, except for the terms in  that appear in equations (B.128) and (B.129). Therefore the properties of this version of the model in which landlords consume residential land as well as consumption goods are similar to those in our baseline speciÖcation. In particular, the model continues to predict heterogeneous local employment elasticities across locations. 

## B.14 Alternative Production Technology

In this subsection of the web appendix, we show how the production technology can be generalized to introduce intermediate inputs, commercial land use and physical capital. We show that the model continues to imply a gravity equation for commuting áows and hence continues to predict heterogeneous local employment elasticities. In our baseline speciÖcation in the paper, we assume the following total cost function for tradeable varieties: 

$$
\Lambda_ {i} (j) = l _ {i} (j) w _ {i} = \left(F + \frac {x _ {i} (j)}{A _ {i}}\right) w _ {i}.\tag{B.130}
$$

We now consider a generalization of this production technology, in which total costs are a Cobb-Douglas function of labor (with wage $w _ { i } )$ , intermediate inputs (with price $P _ { i } )$ , commercial land (with rental rate $Q _ { i } )$ and physical capital (with common rental rate R). We follow Krugman and Venables (1995) and Eaton and Kortum (2002) in assuming that intermediate inputs enter the total cost function through the same CES aggregator as for Önal consumption. Perfect capital mobility ensures that the capital rental rate is the same for all locations $( \mathbb { R } _ { i } = \mathbb { R }$ for all i). Therefore the total cost function now becomes: 

$$
\Lambda_ {i} (j) = \left(F + \frac {x _ {i} (j)}{A _ {i}}\right) w _ {i} ^ {\beta_ {L}} Q _ {i} ^ {\beta_ {Q}} \mathbb {R} ^ {\beta_ {R}} P _ {i} ^ {1 - \beta_ {L} - \beta_ {Q} - \beta_ {R}}.\tag{B.131}
$$

The probability that a worker chooses to live in location n and work in location i remains the same as in equation (10) in the paper, which in turn implies that the commuter market clearing condition takes exactly the same form as in our baseline speciÖcation: 

$$
L _ {n} = \sum_ {r \in N} \frac {B _ {r n} (w _ {n} / \kappa_ {r n}) ^ {\epsilon}}{\sum_ {s \in N} B _ {r s} (w _ {s} / \kappa_ {r s}) ^ {\epsilon}} R _ {r}.\tag{B.132}
$$

Given the same commuter market clearing condition, the partial elasticity of employment with respect to the wage takes the same form as for our baseline speciÖcation: 

$$
\frac {\partial L _ {n}}{\partial w _ {n}} \frac {w _ {n}}{L _ {n}} = \epsilon \sum_ {r \in N} \left(1 - \lambda_ {r n | r} ^ {R}\right) \vartheta_ {r n} + \vartheta_ {n n} \left(\frac {\partial R _ {n}}{\partial w _ {n}} \frac {w _ {n}}{R _ {n}}\right),\tag{B.133}
$$

where $\vartheta _ { r n } = \lambda _ { r n | r } ^ { R } R _ { r } / L _ { n }$ is the share of commuters from residence r in workplace nís employment. 

In general, incorporating additional factors of production a§ects the partial elasticity of wages with respect to productivity, but it leaves the partial elasticity of employment with respect to wages in equation (B.133) unchanged. The reason is that the modelís prediction of heterogeneous local employment elasticities with respect wages is a generic implication of a gravity equation for commuting. 

## B.15 Heterogeneity in E§ective Units of labor

In this section of the web appendix, we consider an alternative speciÖcation of the model, with an idiosyncratic draw to e§ective units of labor instead of amenities. Under this alternative speciÖcation, the idiosyncratic draw $\left( b _ { n i \omega } \right)$ no longer enters the direct utility function, which is now: 

$$
U _ {n i \omega} = \frac {1}{\kappa_ {n i}} \left(\frac {C _ {n \omega}}{\alpha}\right) ^ {\alpha} \left(\frac {H _ {n \omega}}{1 - \alpha}\right) ^ {1 - \alpha}.\tag{B.134}
$$

However, the idiosyncratic draw continues to enter the indirect utility function in exactly the same form as in our baseline speciÖcation, because worker income now depends on the wage per e§ective unit of labor $( w _ { i } )$ times the realization for e§ective units of labor $\left( b _ { n i \omega } \right)$ : 

$$
U _ {n i \omega} = \frac {b _ {n i \omega} w _ {i}}{\kappa_ {n i} P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}}.\tag{B.135}
$$

Therefore the probability that a worker chooses to live in location n and work in location i takes exactly the same form as in our baseline speciÖcation: 

$$
\lambda_ {n i} = \frac {B _ {n i} (\kappa_ {n i} P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}) ^ {- \epsilon} w _ {i} ^ {\epsilon}}{\sum_ {r \in N} \sum_ {s \in N} B _ {r s} (\kappa_ {r s} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}) ^ {- \epsilon} w _ {s} ^ {\epsilon}}.\tag{B.136}
$$

The main di§erence between our baseline speciÖcation and this alternative speciÖcation is the interpretation of wages in the data. In our baseline speciÖcation in terms of amenities, the observed wage for each workplace in the data corresponds directly to the wage in the model, and worker mobility ensures that expected utility is equalized across all workplace-residence pairs (but real wages without taking into account amenities di§er). In contrast, in this alternative speciÖcation in terms of e§ective units of labor, the observed wage for each workplace in the data corresponds to the wage per e§ective unit of labor times average e§ective units of labor conditional on choosing that workplace, and worker mobility ensures that expected real earnings after taking into account average e§ective units of labor are equalized across all workplace-residence pairs. 

## B.16 Commuting Costs in Terms of Labor

In this section of the web appendix, we consider an alternative speciÖcation of the model, in which commuting costs are modeled as a reduction in e§ective units of labor instead of as a reduction in utility. Under this alternative speciÖcation, the iceberg commuting cost $\left( \kappa _ { n i } \right)$ no longer enters the direct utility function, which is now: 

$$
U _ {n i \omega} = b _ {n i \omega} \left(\frac {C _ {n \omega}}{\alpha}\right) ^ {\alpha} \left(\frac {H _ {n \omega}}{1 - \alpha}\right) ^ {1 - \alpha}.\tag{B.137}
$$

However, the iceberg commuting cost continues to enter the indirect utility function in exactly the same form as in our baseline speciÖcation, because worker income now depends on the wage per e§ective unit of labor (w<sub>i</sub>) times e§ective units of labor net of commuting $( 1 / \kappa _ { n i } )$ : 

$$
U _ {n i \omega} = \frac {b _ {n i \omega} w _ {i}}{\kappa_ {n i} P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}}.\tag{B.138}
$$

Therefore the probability that a worker chooses to live in location n and work in location i takes exactly the same form as in our baseline speciÖcation: 

$$
\lambda_ {n i} = \frac {B _ {n i} (\kappa_ {n i} P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}) ^ {- \epsilon} w _ {i} ^ {\epsilon}}{\sum_ {r \in N} \sum_ {s \in N} B _ {r s} (\kappa_ {r s} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}) ^ {- \epsilon} w _ {s} ^ {\epsilon}}.\tag{B.139}
$$

The main di§erence between our baseline speciÖcation and this alternative speciÖcation is whether commuting reduces utility or the labor available for production. One way of interpreting this di§erence is whether workers absorb the commuting cost through reduced leisure or work time. 

## B.17 Partial Local and National Ownership of Land

In our baseline speciÖcation in the paper, we assume that land is owned by immobile landlords, who receive worker expenditure on residential land as income, and consume only goods where they live. This assumption allows us to incorporate general equilibrium e§ects from changes in the value of land, without introducing an externality into workersílocation decisions from the local redistribution of land rents. In this section of the web appendix, we report a robustness test, in which we instead allow for partial local distribution of land rents (as in Caliendo et al. 2014). In particular, we assume that the share $\left( 1 - \iota _ { n } \right)$ of expenditure on residential land is redistributed lump sum to local residents, while the remaining share $\left( \iota _ { n } \right)$ is paid into a national portfolio owned in equal shares by residents throughout the economy. We choose the land ownership share $\left( \iota _ { n } \right)$ to rationalize the trade deÖcit for each county in the data. We show that our Öndings for heterogeneous local employment elasticities are robust to these alternative assumptions about the ownership of land. 

## B.17.1 Expenditure and Income

Let $X _ { n }$ denote the total expenditure of residents in location n. A fraction (1  ) of this expenditure is allocated to land. Of this expenditure on land, we assume that a fraction $\left( 1 - \iota _ { n } \right)$ is redistributed lump sum to local residents, while the remaining fraction $\iota _ { n }$ is paid into a national portfolio owned in equal shares by residents throughout the economy. The per capita return from the national land portfolio is given by 

$$
\xi \equiv \frac {\sum_ {i \in N} \iota_ {i} (1 - \alpha) X _ {i}}{\sum_ {i} R _ {i}}.\tag{B.140}
$$

Using this deÖnition, expenditure in location n can be written as the sum of residential income, nationallyredistributed land rent, and locally-redistributed land rent 

$$
X _ {n} = \bar {v} _ {n} R _ {n} + \xi R _ {n} + (1 - \iota_ {n}) (1 - \alpha) X _ {n},\tag{B.141}
$$

and the trade deÖcit (equal to expenditure minus income) for each location can be expressed as 

$$
D _ {n} \equiv X _ {n} - (\bar {v} _ {n} R _ {n} + Q _ {n} H _ {n}) = \xi R _ {n} - \iota_ {i} (1 - \alpha) X _ {i}.\tag{B.142}
$$

Using equation (B.142) to substitute for $\xi R _ { n }$ in equation (B.141), expenditure in location n can be equivalently written as 

$$
X _ {n} = \frac {\bar {v} _ {n} R _ {n} + D _ {n}}{\alpha}.\tag{B.143}
$$

## B.17.2 Calibrating  to Rationalize the Observed Trade DeÖcits

We calibrate the land ownership shares $\left( \iota _ { n } \right)$ for each location n to rationalize the observed trade deÖcits for each location in the initial equilibrium in the data. Using expenditure (B.141), and denoting the population share of each location in the initial equilibrium by $\textstyle \rho _ { n } \equiv R _ { n } / \sum _ { i } R _ { i }$ , we have 

$$
X _ {n} = \bar {v} _ {n} R _ {n} + \rho_ {n} \sum_ {i \in N} \iota_ {i} (1 - \alpha) X _ {i} + (1 - \iota_ {n}) (1 - \alpha) X _ {n}.\tag{B.144}
$$

Using equations (B.143) and (B.144), we have 

$$
D _ {n} = \alpha X _ {n} - \bar {v} _ {n} R _ {n} = \rho_ {n} \sum_ {i \in N} \iota_ {i} (1 - \alpha) X _ {i} - \iota_ {n} (1 - \alpha) X _ {n},\tag{B.145}
$$

which provides a linear system of equations for each location that can be solved for the unique values of $\iota _ { n }$ that rationalize the observed trade deÖcits as an initial equilibrium of the model. 

## B.17.3 General Equilibrium

We now examine the implications of these alternative assumptions about land ownership for the system of equations that determines general equilibrium. First, workplace-residence choice probabilities $\left( \lambda _ { n i } \right)$ take a similar form as in our baseline speciÖcation in the paper 

$$
\lambda_ {n i} = \frac {B _ {n i} \left(\kappa_ {n i} P _ {n} ^ {\alpha} Q _ {n} ^ {1 - \alpha}\right) ^ {- \epsilon} w _ {i} ^ {\epsilon}}{\sum_ {r \in N} \sum_ {s \in N} B _ {r s} \left(\kappa_ {r s} P _ {r} ^ {\alpha} Q _ {r} ^ {1 - \alpha}\right) ^ {- \epsilon} w _ {s} ^ {\epsilon}}.\tag{B.146}
$$

Therefore the expressions for the number of residents $\left( R _ { n } \right)$ and workers $( L _ { i } )$ in each location take the same form as in our baseline speciÖcation in the paper 

$$
R _ {n} = \bar {L} \sum_ {i \in N} \lambda_ {n i},\tag{B.147}
$$

$$
L _ {i} = \bar {L} \sum_ {n \in N} \lambda_ {n i}.\tag{B.148}
$$

Residential expenditure and income are related through equation (B.141), as reproduced here 

$$
X _ {n} = \bar {v} _ {n} R _ {n} + \xi R _ {n} + (1 - \iota_ {n}) (1 - \alpha) X _ {n},\tag{B.149}
$$

where expected residential income $\left( { \bar { v } } _ { n } \right)$ is given by 

$$
\bar {v} _ {n} = \sum_ {i \in N} \lambda_ {n i | n} ^ {R} w _ {i}.\tag{B.150}
$$

and nationally-redistributed rent per capita () in equation (B.140) can be written as 

$$
\xi = \frac {\sum_ {i \in N} \iota_ {i} (1 - \alpha) X _ {i}}{\bar {L}}.\tag{B.151}
$$

Workplace income equals expenditure on goods produced in that location 

$$
w _ {i} L _ {i} = \sum_ {n \in N} \pi_ {n i} \alpha X _ {n},\tag{B.152}
$$

where the bilateral trade shares $\left( \pi _ { n i } \right)$ are given by 

$$
\pi_ {n i} = \frac {L _ {i} \left(d _ {n i} w _ {i} / A _ {i}\right) ^ {1 - \sigma}}{\sum_ {k \in N} L _ {k} \left(d _ {n k} w _ {k} / A _ {k}\right) ^ {1 - \sigma}}.\tag{B.153}
$$

Finally, the land rent $\left( Q _ { n } \right)$ and price index for tradeables $\left( P _ { n } \right)$ are given by 

$$
Q _ {n} = (1 - \alpha) \frac {X _ {n}}{H _ {n}},
$$

(B.154) 

$$
P _ {n} = \frac {\sigma}{\sigma - 1} \left(\frac {L _ {n}}{\sigma F \pi_ {n n}}\right) ^ {\frac {1}{1 - \sigma}} \frac {d _ {n n} w _ {n}}{A _ {n}}.\tag{B.155}
$$

## B.17.4 Computational Algorithm for Counterfactual Changes

We now discuss the computational algorithm that we use to solve this system of equations for a counterfactual equilibrium given the modelís parameters $\{ \alpha , \sigma , \epsilon , \delta , \kappa \}$ , our calibrated land ownership shares $\iota _ { n } .$ and assumed changes in the exogenous variables of the model $\{ \hat { A } _ { n } , \hat { B } _ { n } , \hat { \kappa } _ { n i } , \hat { d } _ { n i } \}$ . We start with initial guesses for the proportional changes in commuting probabilities, wages and expenditure: $\{ \hat { \lambda } _ { n i } , \ \hat { w } _ { i } , \ \hat { X } _ { i } \}$ Using these initial guesses in the system of equations for general equilibrium, we compute the following proportional changes in the endogenous variables of the model 

$$
\widehat {\overline {{v}}} _ {n} ^ {(t)} = \frac {1}{\overline {{v}} _ {n}} \sum_ {i \in N} \frac {\hat {B} _ {n i} \lambda_ {n i} \left(\hat {w} _ {i} ^ {(t)} / \widehat {\kappa} _ {n i}\right) ^ {\epsilon}}{\sum_ {s \in N} \hat {B} _ {n s} \lambda_ {n s} \left(\hat {w} _ {s} ^ {(t)} / \widehat {\kappa} _ {n s}\right) ^ {\epsilon}} \hat {w} _ {i} ^ {(t)} w _ {i},\tag{B.156}
$$

$$
\hat {L} _ {i} ^ {(t)} = \frac {\bar {L}}{L _ {i}} \sum_ {n} \lambda_ {n i} \hat {\lambda} _ {n i} ^ {(t)},\tag{B.157}
$$

$$
\hat {R} _ {n} ^ {(t)} = \frac {\bar {L}}{R _ {n}} \sum_ {i} \lambda_ {n i} \hat {\lambda} _ {n i} ^ {(t)},\tag{B.158}
$$

which are functions of the observed values of variables in the initial equilibrium and our guesses. From the land market clearing condition (B.154), the proportional change in land rents equals our guess for the proportional change in expenditure 

$$
\widehat {Q} _ {n} ^ {(t)} = \widehat {X} _ {n} ^ {(t)}.\tag{B.159}
$$

Using the proportional change in employment from equation (B.157) and our guess for the proportional change in wages, we can solve for the proportional change in trade shares from equation (B.153) 

$$
\widehat {\pi} _ {n i} ^ {(t)} = \frac {\widehat {L} _ {i} ^ {(t)} \left(\widehat {d} _ {n i} \widehat {w} _ {i} ^ {(t)} / \widehat {A} _ {i}\right) ^ {1 - \sigma}}{\sum_ {k \in N} \pi_ {n k} \widehat {L} _ {k} ^ {(t)} \left(\widehat {d} _ {n k} \widehat {w} _ {k} ^ {(t)} / \widehat {A} _ {k}\right) ^ {1 - \sigma}}.\tag{B.160}
$$

Using the proportional change in employment from equation (B.157), the proportional change in trade shares from equation (B.160) and our guess for the proportional change in wages, we can solve for the proportional change in the tradeables price index from equation (B.155) 

$$
\widehat {P} _ {n} ^ {(t)} = \left(\frac {\widehat {L} _ {n} ^ {(t)}}{\widehat {\pi} _ {n n} ^ {(t)}}\right) ^ {\frac {1}{1 - \sigma}} \frac {\widehat {w} _ {n} ^ {(t)}}{\widehat {A} _ {n}}.\tag{B.161}
$$

Using our guess for the proportional change in expenditure $( \hat { X } _ { i } ^ { ( t ) } )$ , we can also compute the counterfactual change in nationally-redistributed rent per capita from equation (B.151) 

$$
\hat {\xi} ^ {(t)} = \frac {1}{\xi} \frac {\sum_ {i \in N} \iota_ {i} (1 - \alpha) X _ {i} \hat {X} _ {i}}{\bar {L}}.\tag{B.162}
$$

Finally, we use (B.156)-(B.162) in the equality between income and expenditure (B.152), the workplaceresidence choice probabilities (B.146) and expenditure (B.149) to solve for the implied proportional changes in wages, commuting probabilities and expenditure as 

$$
\tilde {w} _ {i} ^ {(t + 1)} = \frac {1}{w _ {i} L _ {i} \widehat {L} _ {i} ^ {(t)}} \sum_ {n \in N} \alpha \pi_ {n i} \widehat {\pi} _ {n i} ^ {(t)} X _ {n} \hat {X} _ {n} ^ {(t)},\tag{B.163}
$$

$$
\tilde {\lambda} _ {n i} ^ {(t + 1)} = \frac {\hat {B} _ {n i} (\widehat {P} _ {n} ^ {(t) \alpha} \widehat {Q} _ {n} ^ {(t) 1 - \alpha}) ^ {- \epsilon} (\widehat {w} _ {i} ^ {(t)} / \widehat {\kappa} _ {n i}) ^ {\epsilon}}{\sum_ {r \in N} \sum_ {s \in N} \hat {B} _ {r s} \lambda_ {r s} (\widehat {P} _ {r} ^ {(t) \alpha} \widehat {Q} _ {r} ^ {(t) 1 - \alpha}) ^ {- \epsilon} (\widehat {w} _ {s} ^ {(t)} / \widehat {\kappa} _ {r s}) ^ {\epsilon}},\tag{B.164}
$$

$$
\tilde {X} _ {n} ^ {(t + 1)} = \frac {\bar {v} _ {n} R _ {n} \hat {v} _ {n} \hat {R} _ {n} + \xi R _ {n} \hat {\xi} ^ {(t)} \hat {R} _ {n} ^ {(t)} + (1 - \iota_ {n}) (1 - \alpha) \hat {X} _ {n} ^ {(t)}}{X _ {n}}.\tag{B.165}
$$

Using these solutions, we update our guesses for wages, commuting probabilities, and expenditures as 

$$
\hat {w} _ {i} ^ {(t + 1)} = \zeta \hat {w} _ {i} ^ {(t)} + (1 - \zeta) \tilde {w} _ {i} ^ {(t + 1)},\tag{B.166}
$$

$$
\hat {\lambda} _ {i} ^ {(t + 1)} = \zeta \hat {\lambda} _ {i} ^ {(t)} + (1 - \zeta) \tilde {\lambda} _ {i} ^ {(t + 1)},\tag{B.167}
$$

$$
\hat {X} _ {i} ^ {(t + 1)} = \zeta \hat {X} _ {i} ^ {(t)} + (1 - \zeta) \tilde {X} _ {n} ^ {(t + 1)},\tag{B.168}
$$

where $\zeta \in ( 0 , 1 )$ is an adjustment factor. 

## B.17.5 Local Employment Elasticities

As in Section 4 in the paper, we compute 3,111 counterfactual exercises where we shock each county with a 5 percent productivity shock (holding productivity in all other counties and holding all other exogenous variables constant). Figure B.8 shows the estimated kernel densities for the distributions of the general equilibrium elasticities of employment (solid blue line) and residents (dashed red line) with respect to the productivity shock across these treated counties. We also show the 95 percent conÖdence intervals around these estimated kernel densities (gray shading). This Ögure is analogous to Figure 2 in the paper, but reports results for this robustness speciÖcation, in which the local rents from land are partially redistributed locally and partially contributed to a global portfolio. We continue to Önd substantial heterogeneity in local employment elasticities that is around the same magnitude as in our baseline speciÖcation. This pattern of results is consistent with the heterogeneity in local employment elasticities being a generic prediction of a gravity equation for commuting áows. As a result, our Öndings of heterogeneous local employment elasticities are robust across di§erent assumptions about the ownership of land. 

## C Additional Empirical Results

In this part of the web appendix, we report additional empirical results and robustness tests. Subsection C.1 shows that the modelís predictions for land prices are strongly positively correlated with median house prices in the data. Subsection C.2 reports standardized coe¢cients for the regressions examining the determinants of the local employment elasticity in Table 2 in the paper. 

Subsection C.3 reports additional results from estimating ìdi§erence-in-di§erencesî regressions using the counterfactuals from the model, as discussed in Subsection 4.1 of the paper. We show that the modelsuggested controls are more successful in explaining the heterogeneous treatment e§ects than the standard empirical controls from the local labor markets literature. Subsection C.4 shows that the heterogeneity in local employment elasticities remains if we shock counties with spatially-correlated shocks reproducing the industrial composition of the U.S. economy. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/e4aac14fe1ee10907dc64712b582b76046c70a1edf970da5701c226d1bd895d0.jpg)



Figure B.8: Kernel density for the distribution of employment and residents elasticities in response to a productivity shock across counties (partial local and national ownership of land)


Subsection C.5 reports additional results from the extension of the model to incorporate heterogeneous positive supply elasticities for developed land following Saiz (2010), as considered in Subsection 4.2 of the paper. Subsection C.6 provides further evidence on the role of commuting in generating heterogeneity in local employment elasticities in our quantitative model. We show that there is substantially less heterogeneity in these elasticities in a counterfactual world with no commuting between counties. 

Subsection C.7 reports additional results for the MDP experiment from Section 5 of the paper. Subsection C.8 reports the shift-share decompositions of cross-section and time-series variation in employment discussed in Section 5 of the paper. 

Subsection C.10 reports counterfactuals for a 20 percent reduction in the costs of trading costs, both starting from the initial equilibrium in the data with commuting, and starting from a counterfactual equilibrium with no commuting. Subsection C.11 shows that we continue to Önd substantial heterogeneity in local employment elasticities when we replicate our entire quantitative analysis for commuting zones (CZs) rather than for counties. 


Dashed line: linear fit; slope: 2.04


## C.1 Land Prices

In this subsection of the web appendix, we show that the modelís predictions for land prices are strongly positively correlated with observed median house prices. In our baseline speciÖcation, we assume Cobb-Douglas utility and interpret land area as geographical land area. In Figure C.1, we show the predictions for land prices from this baseline speciÖcation against median house prices in the data. We Önd a strong and approximately log linear relationship, with a regression slope coe¢cient of 2.04 and R-squared of 0.26. Therefore, although our model is necessarily an abstraction, and there are a number of potential sources of di§erences between land prices in the model and median house prices in the data, we Önd that the mode has strong predictive power. In Section 4.2 of the paper, we generalize this baseline speciÖcation to allow for a positive supply elasticity for developed land that is heterogeneous across locations. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/e7772fd0a91a7d2dcd92931e41ba094e132f2d762542a9fe2fa5669ed03e0fbd.jpg)



Figure C.1: Land Prices in the Model and House Prices in the Data


## C.2 Standardized Employment Elasticities Regression


Table C.1 reports the estimated coe¢ cients from the same set of regressions presented in Table 2 in the paper, after standardizing all variables to make their means zero and standard deviations one. Hence, all coe¢ cients can be interpreted as the fraction of standard deviations by which the dependent variable changes with a one standard deviation change in each independent variable.


<table><tr><td></td><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td><td>7</td><td>8</td><td>9</td></tr><tr><td colspan="10">Dependent Variable: Elasticity of Employment</td></tr><tr><td><eq>\log L_i</eq></td><td></td><td>-0.012(0.056)</td><td>0.036(0.046)</td><td>-0.217**(0.025)</td><td></td><td></td><td></td><td>0.147**(0.018)</td><td>0.132**(0.017)</td></tr><tr><td><eq>\log w_i</eq></td><td></td><td></td><td>-0.126**(0.037)</td><td>-0.100**(0.024)</td><td></td><td></td><td></td><td>-0.162**(0.010)</td><td>-0.166**(0.010)</td></tr><tr><td><eq>\log H_i</eq></td><td></td><td></td><td>-0.621**(0.045)</td><td>-0.372**(0.033)</td><td></td><td></td><td></td><td>0.007(0.020)</td><td>0.020(0.020)</td></tr><tr><td><eq>\log L_{,-i}</eq></td><td></td><td></td><td></td><td>0.429**(0.061)</td><td></td><td></td><td></td><td>-0.097**(0.032)</td><td>-0.097**(0.033)</td></tr><tr><td><eq>\log \bar{w}_{-i}</eq></td><td></td><td></td><td></td><td>0.090*(0.037)</td><td></td><td></td><td></td><td>0.072**(0.016)</td><td>0.091**(0.017)</td></tr><tr><td><eq>\lambda_{ii|i}^R</eq></td><td></td><td></td><td></td><td></td><td>-0.945**(0.019)</td><td></td><td></td><td></td><td></td></tr><tr><td><eq>\sum_{n\in N}(1-\lambda_{Rni})\vartheta_{ni}</eq></td><td></td><td></td><td></td><td></td><td></td><td>1.462**(0.101)</td><td></td><td>1.343**(0.093)</td><td></td></tr><tr><td><eq>\vartheta_{ii}\left(\frac{\lambda_{ii}}{\lambda_{Ri}}-\lambda_{Li}\right)</eq></td><td></td><td></td><td></td><td></td><td></td><td>0.487**(0.112)</td><td></td><td>0.322**(0.093)</td><td></td></tr><tr><td><eq>\frac{\partial w_i}{\partial A_i}\frac{A_i}{w_i}</eq></td><td></td><td></td><td></td><td></td><td></td><td>-0.110**(0.013)</td><td></td><td>-0.090**(0.016)</td><td></td></tr><tr><td><eq>\frac{\partial w_i}{\partial A_i}\frac{A_i}{w_i}\cdot\sum_{r\in N}(1-\lambda_{rn|r})\vartheta_{rn}</eq></td><td></td><td></td><td></td><td></td><td></td><td></td><td>0.544**(0.047)</td><td></td><td>0.576**(0.048)</td></tr><tr><td><eq>\frac{\partial w_i}{\partial A_i}\frac{A_i}{w_i}\cdot\vartheta_{ii}\left(\frac{\lambda_{ii}}{\lambda_{Ri}}-\lambda_{Li}\right)</eq></td><td></td><td></td><td></td><td></td><td></td><td></td><td>-0.428**(0.051)</td><td></td><td>-0.444**(0.048)</td></tr><tr><td>Constant</td><td>-0.000(0.090)</td><td>-0.000(0.090)</td><td>0.000(0.046)</td><td>0.000(0.036)</td><td>-0.000(0.031)</td><td>0.000(0.028)</td><td>-0.000(0.029)</td><td>-0.006(0.026)</td><td>-0.006(0.026)</td></tr><tr><td><eq>R^2</eq></td><td>0.00</td><td>0.00</td><td>0.40</td><td>0.51</td><td>0.89</td><td>0.93</td><td>0.93</td><td>0.95</td><td>0.95</td></tr><tr><td>N</td><td>3,111</td><td>3,111</td><td>3,111</td><td>3,081</td><td>3,111</td><td>3,111</td><td>3,111</td><td>3,081</td><td>3,081</td></tr></table>


Note: $\begin{array} { r } { L , \ L _ { - n } \equiv \sum _ { r : d _ { r n } \leq 1 2 0 , r \not = n } L , } \end{array}$ is the total employment in n neighbors whose centroid is no more than 120km away; $\begin{array} { r } { \bar { w } _ { - n } \equiv \sum _ { r : d _ { r n } \leq 1 2 0 , r \neq n } \frac { L _ { r } } { L _ { , - n } } w _ { r } } \end{array}$ is the weighted average of their workplace wage. All variables are standardized. Standard errors are clustered by state.  denotes signiÖcance at the 5 percent level;  denotes signiÖcance at the 1 percent level. 


Table C.1: Explaining the general equilibrium local employment elasticities to a 5 percent productivity shock (standardized regression) 

## C.3 Additional Treatment Heterogeneity Results

In this subsection of the web appendix, we supplement the results reported in Subsection 4.1 of the paper, and provide further evidence that the model-suggested controls are more successful in explaining the heterogeneity in treatment e§ects in our quantitative model than the standard empirical controls from the local labor markets literature. We compute the deviation between the general equilibrium elasticity in the model and the predicted elasticity from the reduced-form regression for each of the control groups (i)-(v) discussed in the paper 

$$
\beta_ {i} = \frac {a _ {1} + a _ {3} X _ {i t}}{0 . 0 5} - \frac {d L _ {i}}{d A _ {i}} \frac {A _ {i}}{L _ {i}},
$$

where we scale the regression estimates by size of the productivity shock. 

In Figure C.2, we show that this deviation between the general equilibrium elasticity and the ìdi§erencesin-di§erencesî prediction is systematically related to the size of the general equilibrium employment elasticity in the model. For the speciÖcations using reduced-form controls (left panel) and model-generated controls (right panel), we display the results of locally-linear weighted least squares regressions of the deviation term $\beta _ { i }$ against the general equilibrium employment elasticity ${ \frac { d L _ { i } } { d A _ { i } } } { \frac { A _ { i } } { L _ { i } } }$ , along with 95% conÖdence intervals. In each panel, we show the results of these regressions for each group of control counties, where the results using random county ((i) above), non-neighbors ((iv) above) and all counties ((v) above) are visually indistinguishable. 

Using reduced-form controls (left panel) and all deÖnitions of the control group except for the closest county (red line), we Önd that low elasticities are substantially over-estimated, while high elasticities are substantially under-estimated. This pattern of results is intuitive: low and high elasticities occur where commuting linkages are weak and strong respectively. A reduced-form speciÖcation that ignores commuting linkages cannot capture this variation and hence tends to overpredict for low elasticities and underpredic for high elasticities. This e§ect is still present for the closest county control group (red line), as reáected in the downward-sloping relationship between the deviation term and the general equilibrium elasticity. However, the closest county tends to be negatively a§ected by the productivity shock, which shifts the distribution of predicted treatment e§ects (and hence the distribution of the deviation term) upwards. 

Using model-suggested controls (right panel) and all deÖnitions of the control group except for the closest county (red line), we Önd that the deviation term for the ìdi§erence-in-di§erencesî predictions is close to zero and has a much weaker downward-sloping relationship with the general equilibrium elasticity in the model. The exception is the deviation term using the closest-county as a control, which has an upward-sloping relationship with the general equilibrium elasticity in the model and becomes large for high values of this elasticity. The reason is that the productivity shock to treated counties has larger negative e§ects on the closest county for higher values of the general equilibrium elasticity in the model, which leads to a larger upward shift in the distribution of the deviation term. This pattern of results again highlights the potentially large discrepancies from the general equilibrium elasticity from using contiguous locations as controls in the presence of spatial linkages in goods and factor markets. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/91ae9aa3a893802f210d62fd9e72dd2c707c862fcd8e645956291cddbdbff33c.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/1a05bf3db13894ff63c90d06391621c62241b8b21a3b14fdce4edc20553a5089.jpg)



Figure C.2: Average deviation term $\beta _ { i }$ vs. general equilibrium employment elasticity


## C.4 Spatially Correlated Productivity Shocks

In this section of the web appendix, we show that the heterogeneity in local employment elasticities remains if we shock counties with spatially-correlated productivity shocks reproducing the industrial composition of the US economy. We construct these spatially-correlated shocks using aggregate productivity growth in manufacturing and non-manufacturing and the observed shares of these sectors within each countyís employment. In particular, we proceed as follows. Data from BLS shows that between 2004 and 2010 TFP grew 6.2% for the manufacturing sector and 3.4% for the overall private business sector. Given a U.S. employment share in manufacturing of about 11% in 2007 (computed from County Business Patterns; see Data Appendix below), we infer a growth in the non-manufacturing sectorís TFP of 3.1%. We use the County Business Patterns 2007 data to also compute the share of each countyís manufacturing employment over total employment. Figure C.3 shows a map of these shares across the United States. 

We Örst show the consequences of a spatially correlated shock to manufacturing. We compute the equilibrium change in employment and residents in a single counterfactual exercise where each countyís productivity is changed by 6.1% times the share of manufacturing employment in that county: hence, the spatial correlation in manufacturing shares induces a spatial correlation in productivity shocks. Figure C.4 shows the resulting distribution of elasticities of employment and residents. 

Figure C.5 shows an analogous exercise for a shock to the non-manufacturing sector. Finally, Figure C.6 shows the same elasticities when both sectors are shocked: in this case, each countyís shock is a weighted average of the national increase in TFP in the manufacturing and non-manufacturing sectors, where the weights are the corresponding employment shares in the county. Across all of these speciÖcations, we continue to Önd substantial heterogeneity in local employment elasticities. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/69ceeadfe87e96883cf5bf96927d42e1c16733e8a1e945fb3129800c3e0ab00c.jpg)



Figure C.3: U.S. countiesíshare of employment in manufacturing, 2007.


![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/25d09a237375adf9f04ca7d3fc702d2b7d56a32a0ed6e7917230e5c712b01559.jpg)



Figure C.4: Kernel density for the distribution of employment and residents elasticities in response to a spatially correlated productivity shock in the manufacturing sector


![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/3fc73d63e94632f8f71c9d858531c148e270c0681ccdaa84bbbaf37848f2e532.jpg)



Figure C.5: Kernel density for the distribution of employment and residents elasticities in response to a spatially correlated productivity shock in the non-manufacturing sector


![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/856d0d7d506c90098894cf5d861b780848bbea10d7468cc778f5ed5b76a72512.jpg)



Figure C.6: Kernel density for the distribution of employment and residents elasticities in response to a spatially correlated shock in both sectors


## C.5 Positive Developed Land Elasticities

In Subsection 4.2 of the paper, we develop an extension of the model in which we interpret the non-traded amenity as developed land and allow for a positive developed land supply elasticity that can di§er across locations. In this subsection of the web appendix, we provide further details on this robustness check. We introduce a positive developed land supply elasticity by following Saiz (2010) in assuming that the supply of land $\left( H _ { n } \right)$ for each residence n depends on the endogenous price of land $\left( Q _ { n } \right)$ as well as on the exogenous characteristics of locations $\left( { \bar { H } } _ { n } \right)$ 

$$
H _ {n} = \bar {H} _ {n} Q _ {n} ^ {\eta_ {n}},\tag{C.1}
$$

where $\eta _ { n } \geq 0$ is the developed land supply elasticity, which we allow to vary across locations; $\eta _ { n } = 0$ is our baseline speciÖcation of a perfectly inelastic land supply; and $\eta _ { n } \to \infty$ is the special case of a perfectly elastic land supply. 

Introducing a positive and heterogeneous developed land supply elasticity only a§ects one of the conditions for general equilibrium in the model, namely, the land market clearing condition. The rest of the model remains identical. Using the supply function for land (C.1) in the land market clearing condition (5), we obtain the following generalization of our earlier expression for the equilibrium price of land $\left( Q _ { n } \right)$ 

$$
Q _ {n} = \left((1 - \alpha) \frac {\bar {v} _ {n} R _ {n}}{\bar {H} _ {n}}\right) ^ {\frac {1}{1 + \eta_ {n}}}.\tag{C.2}
$$

We now show that our Önding that commuting linkages are important for explaining di§erences in local employment elasticities is robust to controlling for variable housing supply elasticities. As in the main body of the paper, we focus on the subset of counties for which an estimate of the housing supply elasticity is available from Saiz (2010) and no imputation is required. We undertake counterfactuals for productivity shocks for these counties and undertake a horse race, in which we regress the general equilibrium employment elasticities in the model on our measures of commuting linkages, the Saiz housing supply elasticities, and other controls. 

In Table C.2, we report the results from these regressions. In Columns 1-6, we begin by replicating the speciÖcations from Table 2 in the paper for the subsample of counties for which Saiz housing supply elasticities are available. We Önd a similar pattern of results as for the full sample of counties in Table 2 in the paper. In particular, Column 2 shows that the residence own commuting share $( \lambda _ { i i | i } ^ { R } )$ alone explains 63 percent of the variation in local employment elasticities (compared with 89 percent for the full sample). Columns 7-12 of Table C.2 augment the speciÖcations in Columns 1-6 with the Saiz housing supply elasticity. We Önd that both the estimated coe¢ cient and statistical signiÖcance of our commuting measure are robust to the inclusion of the Saiz housing supply elasticity. The Saiz housing supply elasticity is statistically signiÖcant across speciÖcations but contributes only to a minority of the explanatory power of the regression. This pattern of results in consistent with Figure 3 in the paper, where we show that introducing di§erences in housing supply elasticities increases the heterogeneity in the elasticity of residents with respect to the productivity shock, but has relatively little impact on the heterogeneity in the elasticity of employment with respect to the productivity shock. 

This pattern of results is also consistent with existing research on housing supply elasticities. This existing research has typically not distinguished between employment and residents (often focusing on population) and has typically been concerned on metropolitan statistical areas (MSAs) rather than counties. Therefore, the housing supply elasticity can be important for the response of the overall population of metropolitan areas to local labor demand shocks, but there can be considerable variation in the response of employment relative to residents across counties within these metropolitan areas. An important implication is that improvements in commuting technologies provide an alternative approach to relaxing housing supply elasticities in enabling individuals to access high productivity locations. While this possibility has been informally discussed in the existing literature on housing supply elasticities (as for example in Hsieh and Moretti 2017), our paper is the Örst study of which we are aware to provide quantitative empirical evidence on the relevance of commuting for local employment elasticities. 

<table><tr><td></td><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td><td>7</td><td>8</td><td>9</td><td>10</td><td>11</td><td>12</td></tr><tr><td rowspan="2">Dependent Variable:</td><td colspan="12">Elasticity of Employment</td></tr><tr><td colspan="6">without housing supply elasticity</td><td colspan="6">with housing supply elasticity</td></tr><tr><td><eq>\log L_i</eq></td><td>-0.138**(0.021)</td><td></td><td></td><td></td><td>-0.042*(0.020)</td><td>-0.052*(0.022)</td><td>-0.058*(0.022)</td><td></td><td></td><td></td><td>0.054*(0.023)</td><td>0.052*(0.024)</td></tr><tr><td><eq>\log w_i</eq></td><td>-0.318*(0.129)</td><td></td><td></td><td></td><td>-0.359**(0.084)</td><td>-0.417**(0.103)</td><td>-0.183*(0.088)</td><td></td><td></td><td></td><td>-0.217**(0.053)</td><td>-0.249**(0.065)</td></tr><tr><td><eq>\log H_i</eq></td><td>-0.078**(0.019)</td><td></td><td></td><td></td><td>0.022(0.025)</td><td>0.037(0.024)</td><td>-0.117**(0.019)</td><td></td><td></td><td></td><td>-0.029(0.026)</td><td>-0.021(0.025)</td></tr><tr><td><eq>\log L_{, - i}</eq></td><td>-0.009(0.036)</td><td></td><td></td><td></td><td>-0.031(0.032)</td><td>-0.046(0.037)</td><td>-0.034(0.035)</td><td></td><td></td><td></td><td>-0.072*(0.033)</td><td>-0.085*(0.037)</td></tr><tr><td><eq>\log \bar{w}_{-i}</eq></td><td>0.516**(0.147)</td><td></td><td></td><td></td><td>0.137(0.112)</td><td>0.384**(0.128)</td><td>0.667**(0.146)</td><td></td><td></td><td></td><td>0.324**(0.117)</td><td>0.507**(0.120)</td></tr><tr><td><eq>\lambda_{ii|i}^R</eq></td><td></td><td>-1.738**(0.098)</td><td></td><td></td><td></td><td></td><td></td><td>-1.380**(0.108)</td><td></td><td></td><td></td><td></td></tr><tr><td><eq>\sum_{n\in N}(1-\lambda_{Rni})\vartheta_{ni}</eq></td><td></td><td></td><td>5.500**(0.693)</td><td></td><td>3.653**(0.791)</td><td></td><td></td><td></td><td>3.266**(0.511)</td><td></td><td>3.227**(0.720)</td><td></td></tr><tr><td><eq>\vartheta_{ii}\left(\frac{\lambda_{ii}}{\lambda_{Ri}}-\lambda_{Li}\right)</eq></td><td></td><td></td><td>4.014**(0.746)</td><td></td><td>2.181*(0.829)</td><td></td><td></td><td></td><td>1.932**(0.543)</td><td></td><td>1.822*(0.747)</td><td></td></tr><tr><td><eq>\frac{\partial w_i}{\partial A_i}\frac{A_i}{w_i}</eq></td><td></td><td></td><td>-0.448(0.304)</td><td></td><td>-1.275**(0.283)</td><td></td><td></td><td></td><td>-0.743**(0.273)</td><td></td><td>-0.527*(0.239)</td><td></td></tr><tr><td><eq>\frac{\partial w_i}{\partial A_i}\frac{A_i}{w_i}\cdot\sum_{r\in N}(1-\lambda_{rn|r})\vartheta_{rn}</eq></td><td></td><td></td><td></td><td>1.862**(0.216)</td><td></td><td>0.424(0.283)</td><td></td><td></td><td></td><td>0.892**(0.206)</td><td></td><td>0.982**(0.185)</td></tr><tr><td><eq>\frac{\partial w_i}{\partial A_i}\frac{A_i}{w_i}\cdot\vartheta_{ii}\left(\frac{\lambda_{ii}}{\lambda_{Ri}}-\lambda_{Li}\right)</eq></td><td></td><td></td><td></td><td>0.459(0.251)</td><td></td><td>-0.992**(0.317)</td><td></td><td></td><td></td><td>-0.406*(0.189)</td><td></td><td>-0.405(0.236)</td></tr><tr><td>Saiz elasticity</td><td></td><td></td><td></td><td></td><td></td><td></td><td>0.158**(0.025)</td><td>0.140**(0.016)</td><td>0.148**(0.019)</td><td>0.154**(0.021)</td><td>0.152**(0.022)</td><td>0.158**(0.023)</td></tr><tr><td>Constant</td><td>2.088(1.491)</td><td>3.042**(0.049)</td><td>-2.283**(0.469)</td><td>0.828**(0.207)</td><td>3.416**(1.227)</td><td>3.482**(1.218)</td><td>-1.418(1.245)</td><td>2.678**(0.071)</td><td>-0.052(0.340)</td><td>1.526**(0.170)</td><td>-0.741(0.945)</td><td>-0.542(0.980)</td></tr><tr><td><eq>R^2</eq></td><td>0.51</td><td>0.63</td><td>0.64</td><td>0.60</td><td>0.70</td><td>0.68</td><td>0.67</td><td>0.82</td><td>0.82</td><td>0.81</td><td>0.84</td><td>0.83</td></tr><tr><td>N</td><td>457</td><td>460</td><td>460</td><td>460</td><td>457</td><td>457</td><td>457</td><td>460</td><td>460</td><td>460</td><td>457</td><td>457</td></tr></table>


<sub>the</sub> <sub>general</sub> <sub>equilibrium</sub> <sub>local</sub> <sub>e</sub>m<sup>ployment</sup> <sup>elasticities</sup> <sup>to</sup> <sup>a</sup> <sup>5</sup> <sup>percent</sup> <sup>productivity</sup> <sup>shock</sup>


## C.6 Additional Results with No Commuting Between Counties

In this subsection of the web appendix, we provide further evidence that the heterogeneity in local employment elasticities is driven by commuting, by reporting local employment elasticities for a counterfactual world with no commuting between counties. As in our counterfactuals in Section 4 in the paper, we start with the initial equilibrium in the observed data. We Örst undertake a counterfactual for prohibitive com muting costs between counties $( \kappa _ { n i }  \infty { \mathrm { f o r } } n \neq i )$ and solve for the new spatial equilibrium distribution of economic activity. Starting from this counterfactual world with no commuting between counties, we next compute 3,111 counterfactual exercises where we shock each county with a 5 percent productivity shock (holding productivity in all other counties and holding all other exogenous variables constant). 

Figure C.7 shows the estimated kernel density for the distribution of the general equilibrium elasticity of employment with respect to the productivity shock across the treated counties (red dashed line). In this counterfactual world with no commuting, the employment and residents elasticity are equal to one another. To provide a point of comparison, the Ögure also displays the estimated kernel density for the general equilibrium employment elasticity from our baseline speciÖcation in the paper with commuting between counties (blue solid line). Even in the absence of commuting between counties, we expect local employment elasticities to be heterogeneous, because counties di§er substantially from one another in terms of their initial shares of U.S. employment. Consistent with this, we Önd that local employment elasticities in the world with no commuting between counties range from around 0.5 to 1. However, this variation is substantially less than in our baseline speciÖcation with commuting between counties, where the local employment elasticities range from around 0.5 to 2.5. Therefore, these results provide further evidence that commuting indeed plays a central role in generating the heterogeneity in local employment elasticities. Comparing the two speciÖcations in Figure C.7, local employment elasticities are also larger on average with commuting than in the counterfactual world without commuting. This pattern of results is consistent with commuting weakening congestion forces in the model. As a county experiences an increase in productivity, commuting enables it to increase employment by drawing residents from surrounding counties, thereby bidding up land prices less than otherwise would be the case in a world without commuting between counties. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/61c6bb6d1cdd6416bd5d6a8e32d312c5f5d6ea2e048671544b444d5921be144e.jpg)



Figure C.7: Kernel density for the distribution of employment and resident elasticities in response to a productivity shock across counties (with and without commuting between counties)


## C.7 Million Dollar Plants Natural Experiment

In this section of the web appendix, we report additional results for the MDP experiment from Section 5 of the paper.<sup>7</sup> First, we report a balance table that compares the observed characteristics of winner and runner-up counties before a MDP announcement for the full set of 82 cases. Second, we consider the non-parametric speciÖcation that estimates separate treatment e§ects for each MDP case from equation (22) in the paper. 

First, we compare the observed characteristics of winner and runner-up counties before a MDP an nouncement. Table C.3 reports the mean and standard error of the mean for employment, wages, land area and population density Öve years before a MDP announcement for these two groups of counties for the full set of 82 MDPs. We also report the same statistics for workplace and residence own commuting shares in 1990 as the closest Census year.<sup>8</sup> We Önd that winner counties have somewhat lower prior values of levels of employment, wages, population and population density than runner-up counties. We also Önd that they have somewhat more open local labor markets in terms of workplace and residence own commuting shares. Despite these di§erences in individual observed characteristics, the fact that the Örms selected these counties as winners and runners-up suggests that they have similar implied proÖtability for plant location. As a check on the identifying assumption that the losers form a valid counterfactual for the winners, we report an event-study speciÖcation following GHM in Section 5 of the paper. 

<table><tr><td>Variable</td><td>Winner</td><td>Runners-up</td></tr><tr><td>Log employment</td><td>11.122(0.176)</td><td>11.660(0.116)</td></tr><tr><td>Log wages</td><td>2.758(0.032)</td><td>2.802(0.023)</td></tr><tr><td>Log land area</td><td>14.213(0.085)</td><td>14.152(0.063)</td></tr><tr><td>Log Population</td><td>11.999(0.153)</td><td>12.446(0.100)</td></tr><tr><td>Log population density</td><td>-2.214(0.153)</td><td>-1.706(0.115)</td></tr><tr><td>Workplace own commuting share</td><td>0.742(0.015)</td><td>0.764(0.011)</td></tr><tr><td>Residence own commuting share</td><td>0.737(0.020)</td><td>0.786(0.015)</td></tr></table>


Means and standard errors of the mean of observed characteristics; standard errors of the means are in parentheses; employment, wages, land area, population and population density for winner and runner-up counties in each case are measured Öve years before the MDP announcement; workplace and residence commuting shares are measured in 1990.


## Table C.3: Characteristics of Winner and Runner-up Counties Before a MDP Annoucement

Second, we turn to the heterogeneous treatment e§ects speciÖcation from equation (22) in the paper, which estimates a separate treatment e§ect for each of the 82 MDP cases. These heterogeneous treatment e§ects are identiÖed as the mean change in employment in winner counties relative to control counties for each case (the excluded category is the runner-up counties for each case). In Figure C.8 below, we display these estimated treatment e§ects for each case. As apparent from the Ögure, we Önd substantial heterogeneity in these estimated treatment e§ects, which range from less than zero to just below one. 

Therefore, although the average estimated treatment e§ect is positive, there is substantial variation around this average. We reject the null hypothesis that these estimated treatment e§ects all take the same value at conventional levels of statistical signiÖcance (p-value 0.000). 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/3130c6d00250de1a8ee29229f201d2933dec0801498d564be9e3699a1c60fd5a.jpg)



Note: Heterogeneous treatment e§ects for 82 cases. SpeciÖcation includes county, case and year Öxed e§ects, a post-MDP announcement dummy, interaction terms between the dummy for winner county and treatment year, interaction terms between commuting openness and treatment year, and three-way interaction terms between the winner dummy, commuting openness and treatment year (equation (22) in the paper).



Figure C.8: Heterogeneous Treatment E§ects Across MDPs


## C.8 Shift-Share Decomposition

In this subsection of the web appendix, we provide further evidence on the importance of commuting for employment changes using shift-share decompositions, as discussed in the paper. We undertake these decompositions for both cross-section and time-series variation in employment. 

## C.8.1 Cross-section Decomposition

We begin with our cross-section decomposition. We use the accounting identity provided by the commuter market clearing condition, which requires that employment in each county i equals the sum of commuting áows from all counties: 

$$
L _ {i t} = \sum_ {n \in N} \lambda_ {n i | n t} ^ {R} R _ {n t}.\tag{C.3}
$$

Separating these commuting áows into those from the own county and those from other counties, this commuter market clearing condition can be re-written as: 

$$
L _ {i t} = \underbrace {\lambda_ {i i | i t} ^ {R} R _ {i t}} _ {\text {(a) own residents}} + \underbrace {\sum_ {n \neq i} \lambda_ {n i | n t} ^ {R} R _ {n t}} _ {\text {(b) commuters}}.\tag{C.4}
$$

The same accounting also holds for the county with the median level of employment m: 

$$
L _ {m t} = \underbrace {\lambda_ {m m | m t} ^ {R} R _ {m t}} _ {\text {(a) own residents}} + \underbrace {\sum_ {n \neq m} \lambda_ {n m | n t} ^ {R} R _ {n t}} _ {\text {(b) commuters}}.\tag{C.5}
$$

Taking di§erences between equations (C.4) and (C.5), we obtain: 

$$
\Delta^ {I} L _ {i t} = \left[ \lambda_ {i i | i t} ^ {R} R _ {i t} - \lambda_ {m m | m t} ^ {R} R _ {m t} \right] + \left[ \sum_ {n \neq i} \lambda_ {n i | n t} ^ {R} R _ {n t} - \sum_ {n \neq m} \lambda_ {n m | n t} ^ {R} R _ {n t} \right],\tag{C.6}
$$

where $\Delta ^ { I }$ is the cross-section di§erence operator between an individual county i and the county with the median level of employment m (such that $\Delta ^ { I } L _ { i t } = L _ { i t } - L _ { m t } )$ . Subtracting and adding $\lambda _ { i i | i t } ^ { R } R _ { m t }$ from the Örst term in square parentheses, and subtracting and adding $\begin{array} { r } { \sum _ { n \neq m } \lambda _ { n i | n t } ^ { R } R _ { n t } } \end{array}$ from the second term in square parentheses, we have: 

$$
\begin{array}{r l} & {\Delta^ {I} L _ {i t} = \lambda_ {i i | i t} ^ {R} R _ {i t} - \lambda_ {i i | i t} ^ {R} R _ {m t} - \lambda_ {i i | m t} ^ {R} R _ {m t} + \lambda_ {i i | i t} ^ {R} R _ {m t}} \\ & {\qquad + \sum_ {n \neq i} \lambda_ {n i | n t} ^ {R} R _ {n t} - \sum_ {n \neq m} \lambda_ {n i | n t} ^ {R} R _ {n t} - \sum_ {n \neq m} \lambda_ {n m | n t} ^ {R} R _ {n t} + \sum_ {n \neq m} \lambda_ {n i | n t} ^ {R} R _ {n t}.} \end{array}\tag{C.7}
$$

which can be re-written as: 

$$
\Delta L _ {i t} = \underbrace {\lambda_ {i i | i t} ^ {R} \Delta^ {I} R _ {i t}} _ {\mathrm{(i)ownresidents}} + \underbrace {\underbrace {R _ {m t} \Delta^ {I} \lambda_ {i i | i t} ^ {R}} _ {\mathrm{(ii)owncommutingshares}}} + \underbrace {\left(\sum_ {n \neq i} \lambda_ {n i | n t} ^ {R} R _ {n t} - \sum_ {n \neq m} \lambda_ {n i | n t} ^ {R} R _ {n t}\right)} _ {\mathrm{(iii)otherresidents}} + \underbrace {\sum_ {n \neq m} \left(\lambda_ {n i | n t} ^ {R} - \lambda_ {n m | n t} ^ {R}\right) R _ {n t}} _ {\mathrm{(iv)othercommutingshares}}.\tag{C.8}
$$

We thus obtain a decomposition of cross-section di§erences in employment between counties into the following four contributions: (i) di§erences in own residents holding own commuting shares constant; (ii) di§erences in own commuting shares holding own residents constant; (iii) di§erences in other residents holding other commuting shares constant; and (iv) di§erences in other commuting shares holding other residents constant. In the other residents term (iii), the only thing that varies between the two components of the term is the lower limit of the summation, which captures di§erences in the sets of other counties n $\neq i$ and $n \neq m$ . In the other commuting term (iv), the only thing that varies between the two components of the term is the commuting shares with other counties: $\lambda _ { n i | n t } ^ { R } \neq \lambda _ { n m | n t } ^ { R }$ for $i \neq m$ . All four terms are equal to zero for the county m with median employment, we report the distribution of results for all other counties $i \neq m$ 

In interpreting this decomposition, several points are worth bearing in mind. First, we interpret any term involving workers commuting across county borders as capturing commuting, which implies that we view terms (ii), (iii) and (iv) as capturing commuting. In the special case of no commuting between counties, these Önal three terms are all necessarily equal to zero, because in this special case $\lambda _ { i i | i t } ^ { R } = 1$ and $\lambda _ { n i | n t } ^ { R } = \lambda _ { n m | n t } ^ { R } = 0$ for $n \neq i$ and $n \neq m .$ . Nevertheless, we acknowledge that other interpretations are possible, such as only viewing terms involving variation in commuting shares $( \lambda _ { n i | n t } ^ { R } )$ as capturing commuting (terms (ii) and (iv)). Second, we note that counties with similar total residents to the median county might have very di§erent commuting links with other counties. We intentionally capture this variation in our decomposition, but we acknowledge that the importance of these commuting links can be a§ected by many idiosyncratic factors, such as the drawing of county boundaries. 

Third, in principle, the relative importance of the di§erent terms in the cross-section decomposition in equation (C.8) can vary depending on which county is chosen as the base. However, in practice, we Önd a similar qualitative and quantitative pattern of results for alternative choices of the base county other than the county with median employment. Fourth, we recognize that total employment and total residents are still strongly positively correlated across counties, because of third factors that a§ect both employment and residents (e.g. productivity and climate). 

In Table C.4, we report the results of this cross-section decomposition (C.8) using our bilateral commuting data for 2006-10. As all four terms are equal to zero for the county with median employment, we report the distribution of results across all other counties. For each individual county, the four terms add up to the total di§erence in employment, which in turn implies that the mean of the four terms adds up to the mean total di§erence in employment (bottom row), because the mean is a linear operator. The same need not be true for the percentiles of the distribution of each contribution (other rows), because the county at a given percentile for one contribution may be di§erent from the county at the same percentile for another contribution. Each individual term in the decomposition can be positive or negative, as reáected in the negative values in a number of cells in the table. 

As apparent from the table, we Önd substantial contributions from all four terms of the decomposition. On average, we Önd a di§erence in employment from the median county of 45,211, where the fact that this di§erence is positive reáects the fact that the distribution of employment across counties is skewed. Of this 45,211, we Önd that own residents contribute 32,712 (bottom row, second column), own commuting is responsible for 2,712 (bottom row, third column), other residents make a negative contribution of -44,119 (bottom row, fourth column), and other commuting accounts for the remaining 53,906 (bottom row, Öfth column). We also Önd substantial heterogeneity across counties in the relative importance of these four terms. Between the 10th and 90th percentiles, these contributions range from -10,048 to 73,406 for own residents (second column), -1,145 to 6,350 for own commuting (third column), -86,445 to -2,063 for other residents (fourth column), and -1,339 to 108,003 for other commuting (Öfth column). 


Although the individual terms in these shift-share decompositions can be interpreted in di§erent ways, we view these results as supporting the idea that commuting patterns are a quantitatively relevant margin for accounting for cross-sectional di§erences in employment across counties. In Figure 1 and Table 1 in the paper, we report results using a simpler and more intuitive measure of the relevance of commuting, given by the share of residents who work in the county where they live.


<table><tr><td>2006-10</td><td>(i) Changes Own Residents, Constant Own Commuting</td><td>(ii) Changes Own Commuting, Constant Own Residents</td><td>(iii) Changes Other Residents, Constant Other Commuting</td><td>(iv) Changes Other Commuting, Constant Other Residents</td><td>Sum (i)-(iv)</td></tr><tr><td>10th percentile</td><td>-10,048</td><td>-1,145</td><td>-86,445</td><td>-1,339</td><td>-</td></tr><tr><td>25th percentile</td><td>-6,242</td><td>613</td><td>-26,593</td><td>1,372</td><td>-</td></tr><tr><td>50th percentile</td><td>-1,077</td><td>2,917</td><td>-9,822.9</td><td>8,693</td><td>-</td></tr><tr><td>75th percentile</td><td>14,744</td><td>4,999</td><td>-4,157.7</td><td>30,616</td><td>-</td></tr><tr><td>90th percentile</td><td>73,406</td><td>6,350</td><td>-2,063</td><td>108,003</td><td>-</td></tr><tr><td>Mean</td><td>32,712</td><td>2,712</td><td>-44,119</td><td>53,906</td><td>45,211</td></tr></table>


Mean and percentiles of the distribution of the contributions to cross-section di§erences in employment between each county and the median county for 2006-10. The four terms are di§erences in (i) the number of residents holding own commuting shares constant; (ii) own commuting shares holding own residents constant; (iii) other residents holding other commuting shares constant; and (iv) other commuting shares holding other residents constant. 


Table C.4: Cross-section Decomposition of Employment Di§erences across Counties for 2006-10 

## C.8.2 Time-series Decomposition

We next consider our time-series decomposition. Taking di§erences between equation (C.4) for time t and the analogous equation for time t  1, we obtain: 

$$
\Delta^ {T} L _ {i t} = \lambda_ {i i | i t} ^ {R} R _ {i t} - \lambda_ {i i | i t - 1} ^ {R} R _ {i t - 1} + \sum_ {n \neq i} \lambda_ {n i | n t} ^ {R} R _ {n t} - \sum_ {n \neq i} \lambda_ {n i | n t - 1} ^ {R} R _ {n t - 1},\tag{C.9}
$$

where $\Delta ^ { T }$ is the time-series di§erence operator such that $\Delta ^ { T } L _ { i t } = L _ { i t } - L _ { i t - 1 }$ . Subtracting and adding $\lambda _ { i i | i t } ^ { R } R _ { i t - 1 }$ from the Örst term in parentheses, and subtracting and adding $\begin{array} { r } { \sum _ { n \neq i } \lambda _ { n i | n t } ^ { R } R _ { n t - 1 } } \end{array}$ from the second 

term in parentheses, we obtain: 

$$
\begin{array}{r l} & {\Delta^ {T} L _ {i t} = \lambda_ {i i | i t} ^ {R} R _ {i t} - \lambda_ {i i | i t} ^ {R} R _ {i t - 1} - \lambda_ {i i | i t - 1} ^ {R} R _ {i t - 1} + \lambda_ {i i | i t} ^ {R} R _ {i t - 1}} \\ & {\qquad + \sum_ {n \neq i} \lambda_ {n i | n t} ^ {R} R _ {n t} - \sum_ {n \neq i} \lambda_ {n i | n t} ^ {R} R _ {n t - 1} - \sum_ {n \neq i} \lambda_ {n i | n t - 1} ^ {R} R _ {n t - 1} + \sum_ {n \neq i} \lambda_ {n i | n t} ^ {R} R _ {n t - 1}.} \end{array}\tag{C.10}
$$

which can be re-written as: 

$$
\Delta^ {T} L _ {i t} = \underbrace {\lambda_ {i i | i t} ^ {R} \Delta^ {T} R _ {i t}} _ {\text {(i) own residents}} + \underbrace {R _ {i t - 1} \Delta^ {T} \lambda_ {i i | i t} ^ {R}} _ {\text {(ii) own commuting shares}} + \underbrace {\sum_ {n \neq i} \lambda_ {n i | n t} ^ {R} \Delta^ {T} R _ {n t}} _ {\text {(iii) other residents}} + \underbrace {\sum_ {n \neq i} R _ {n t - 1} \Delta^ {T} \lambda_ {n i | n t} ^ {R}} _ {\text {(iv) other commuting shares}}.\tag{C.11}
$$

In interpreting this decomposition, we again view any term involving workers commuting across county borders as capturing commuting (terms (ii), (iii) and (iv)). In the special case of no commuting between counties, the Örst term for changes in own residents $( \lambda _ { i i | i t } ^ { R } \Delta ^ { T } R _ { i t } )$ is the only source of employment changes, because in this special case $\begin{array} { r } { \lambda _ { i i | i t } ^ { R } = 1 , \lambda _ { n i | n t } ^ { R } = 0 } \end{array}$ for $n \neq i ,$ , and $\Delta ^ { T } \lambda _ { n i | n t } ^ { R } = 0$ for all $n , i ,$ which implies that the Önal three terms are necessarily all equal to zero. But we acknowledge that other interpretations are possible, such as only viewing terms that involve variation in commuting shares $( \lambda _ { n i | n t } ^ { R } )$ as capturing commuting (terms (ii) and (iv)). 

In Table C.5, we report the results of this time-series decomposition (C.11) using the change in the bilateral commuting probabilities between 1990 and 2006-2010. As for our cross-section decomposition above, the four terms add up to the total change in employment for an individual county. However, the same need not be true for the percentiles of the distribution of each contribution (other rows), because the county at a given percentile for one contribution may be di§erent from the county at the same percentile for another contribution. Each individual term in the decomposition can be positive or negative, as reáected in the negative values in a number of cells in the table. 

<table><tr><td>1990 to 2006-10</td><td>(i) Changes Own Residents, Constant Own Commuting</td><td>(ii) Changes Own Commuting, Constant Own Residents</td><td>(iii) Changes Other Residents, Constant Other Commuting</td><td>(iv) Changes Other Commuting, Constant Other Residents</td></tr><tr><td>10th percentile</td><td>23</td><td>-2,619</td><td>28</td><td>-184</td></tr><tr><td>25th percentile</td><td>474</td><td>-1,155</td><td>130</td><td>32</td></tr><tr><td>50th percentile</td><td>1,728</td><td>-457</td><td>447</td><td>335</td></tr><tr><td>75th percentile</td><td>6,094</td><td>-82</td><td>1,517</td><td>1,109</td></tr><tr><td>90th percentile</td><td>21,170</td><td>181</td><td>5,626</td><td>3,268</td></tr></table>


Mean and percentiles of the distribution of the contributions to time-series changes in employment between 1990 and 2006-10 from (i) the number of residents holding own commuting shares constant; (ii) own commuting shares holding own residents constant; (iii) other residents holding other commuting shares constant; and (iv) other commuting shares holding other residents constant.



Table C.5: Time-series Decomposition of County Employment Changes between 1990 and 2006-10


We again Önd quantitatively relevant contributions from all four terms in the decomposition. For the median county, we Önd a change in employment of 1,981, of which own residents contribute 1,728 (fourth row, second column), own commuting is responsible for -457 (fourth row, third column), other residents make a contribution of 447 (fourth row, fourth column), and other commuting accounts for the remaining 335 (fourth row, Öfth column). We also Önd substantial heterogeneity across counties in the relative importance of these four terms. Between the 10th and 90th percentiles, these contributions range from 23 to 21,170 for own residents (second column), -2,619 to 181 for own commuting (third column), 28 to 5,626 for other residents (fourth column), and -184 to 3,268 for other commuting (Öfth column). 

In both the cross-section and over time, variation in county employment is ultimately driven by variation in productivity and other county characteristics. Therefore, notwithstanding the caveats discussed above, we view the Öndings of these cross-section and time-series decompositions as supporting the idea that response of employment to such county characteristics is shaped by heterogeneous patterns of commuting áows. 

## C.9 Changes in Commuting Costs

Figure C.9 presents the changes in local employment against the initial labor to resident ratio $( L _ { i } / R _ { i } )$ for the counterfactual in which we reduce commuting cost by the median change between 1990 and 2010, as discussed in Section 6 of the paper. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/efe27c0baee235679185d90b245a696710563ac59c519e67d43f3d5594645eca.jpg)



Employment/Residents Ratio (log scale)



Figure C.9: Counterfactual relative change in county employment $( \hat { L } )$ for median decrease in commuting costs throughout U.S. against initial employment to residents ratio (L=R).


## C.10 Interaction Between Trade and Commuting Costs

In this subsection of the web appendix, we examine the extent to which trade and commuting costs interact in the model, as discussed in Section 6 of the paper. To provide evidence on this interaction, we compare the e§ects of reductions in trade costs, both with and without commuting between counties. To do so, we Örst undertake a counterfactual for a 20 percent reduction in trade costs between locations $( \hat { d } _ { n i } = 0 . 8$ for $n \neq i$ and $\hat { d } _ { n n } = 1 )$ starting from the observed initial equilibrium with commuting between counties (using the observed bilateral commuting shares to implicitly reveal the magnitude of bilateral commuting costs). We next undertake a counterfactual for the same 20 percent reduction in trade costs between locations from a counterfactual equilibrium with no commuting between counties. That is, starting from the observed equilibrium, we Örst undertake a counterfactual for prohibitive commuting costs between counties $( \kappa _ { n i }  \infty$ for $n \neq i )$ , before then undertaking the counterfactual for the reduction in trade costs. 

We Önd that commuting between counties has a relatively small impact on the welfare gains from trade cost reductions. Starting from the observed equilibrium, we Önd aggregate welfare gains from the trade cost reduction of 11.66 percent. In contrast, starting from the counterfactual equilibrium without commuting between counties, we Önd aggregate welfare gains from the same trade cost reduction of 11.56 percent. However, we Önd that commuting between counties plays a major role in ináuencing the impact of trade cost reductions on the spatial distribution of economic activity. Figure C.10 shows the relative change in employment from a 20 percent reduction in trade costs in the New York region (without commuting in the left panel and with commuting in the right panel). In general, reductions in trade costs lead to a more dispersed spatial distribution of economic activity in the model. But this dispersal is smaller with commuting between counties than without it. As trade costs fall, commuting increases the ability of the most productive locations to serve the national market by drawing workers from a suburban hinterland, without bidding up land prices as much as would otherwise occur. 

Intuitively, lower trade costs and higher commuting costs are both forces for the dispersion of economic activity in the model. On the one hand, lower trade costs weaken agglomeration forces by reducing the incentive for Örms and workers to locate close to one another. On the other hand, higher commuting costs increase congestion forces by forcing workers to live where they work, thereby bidding up land prices in congested locations. These two sets of forces interact with one another, so that the impact of a reduction in trade costs depends on the level of commuting costs. While lower trade costs necessarily redistribute employment away from the most congested locations, this redistribution is smaller with commuting between counties than without it. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/c5823ca426d9c558f157e832545c7d2146694af9841b80837a56e535d3264c39.jpg)



Figure C.10: Relative change in employment (L<sup>^</sup>) from a 20 percent reduction in trade costs (with and without commuting between counties) in the New York area


This exercise also illustrates more generally the role of commuting linkages in shaping the consequences of a reduction in trade costs. Figure C.11 shows changes in county employment and real income following a reduction in trade costs in an economy without commuting (vertical axis) and with commuting (horizontal axis), alongside a 45-degree line. We Önd a relatively low correlation between changes in employment with and without commuting between counties. In particular, commuting and trade tend to be complements in expanding areas: whenever employment increases with the reduction in trade costs, the commuting technology allows a larger expansion because it alleviates the increase in congestion (employment changes are below the diagonal in the left panel of Figure C.11). But trade and commuting tend to be local substitutes from the perspective of real income: whenever real income increases with trade, the increase is larger without commuting because production is more spatially dispersed without commuting (real income changes are above the diagonal in the right panel of Figure C.11). These results further underscore the prominence of commuting linkages in shaping the equilibrium spatial distribution of economic activity, and the necessity of incorporating them in models of economic geography. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/a5aa6207d9d58d47922cd71058c3ebfe8e87eb5e7b638665748c21b700f00758.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/b89b72f5c92a54756981ebdefb59559231bb4937f1966751f4c3e36c61ed656f.jpg)



Figure C.11: Relative change in employment (L<sup>^</sup>) and real income $( \widehat { \bar { v } } _ { n } / \left( \hat { P } _ { n } ^ { \alpha } \hat { Q } _ { n } ^ { 1 - \alpha } \right) )$ from a 20 percent breduction in trade costs across all counties (with and without commuting between counties)


## C.11 Commuting Zones (CZs)

As discussed in the paper, previous research has often worked at relatively high levels of spatial aggre gation (e.g. commuting zones (CZs)) to reduce commuting áows. In contrast, we explicitly model the spatial interactions between locations in goods and commuting markets, thereby providing a framework for examining the local impact of labor demand shocks at alternative spatial scales, including those Öner than CZs. 

In our baseline speciÖcation in the paper, we report results for counties, because this is the Önest level of geographical detail at which commuting data are reported for the entire United States in the American Community Survey (ACS) and Census of Population, and a number of ináuential papers in the local labor markets literature have used county data (such as Greenstone, Hornbeck and Moretti 2010). In this section of the web appendix, we report the results of a robustness check, in which we replicate our entire analysis for Commuting Zones (CZs) (aggregations of counties). This replication involves undertaking the full quantitative analysis of the model at this higher level of spatial aggregation. First, we aggregate our employment and wage to the CZ level. Second, we aggregate our bilateral commuting data between pairs of counties to construct bilateral commuting áows between pairs of CZs. Third, we use our data on bilateral trade between CFS regions to solve for implied CZ productivity $\left( A _ { i } \right)$ and bilateral trade between CZs $\left( \pi _ { n i } \right)$ , using the same approach as for counties in our baseline speciÖcation in Section 3.1 of the paper and Section B.5 of this web appendix. Fourth, we use our data on bilateral commuting between pairs of CZs to solve for implied bilateral amenities $\left( B _ { n i } \right)$ , using the same approach as in Section 3.2 of the paper and Section B.6 of this web appendix. In Figure C.12, we show the conditional relationship between the log value of commuting áows and log distance between pairs of CZs, after removing workplace and residence Öxed e§ects. This Ögure is analogous to Figure B.5 in this web appendix, but uses CZs rather than counties. Again we Önd that the gravity equation provides a good approximation to the data, with a tight and approximately log linear relationship between the two variables. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/15b3993d69ffb172f0e4e860e4839f204d0aafd5b76ba6b652238cbfb8e21465.jpg)



Figure C.12: Gravity in Commuting Between Commuting Zones (CZs)


Having calibrated the model to match the initial equilibrium in the observed data at the CZ level, we next shock each of the 709 CZs with a 5 percent productivity shock, following the same approach as for counties in Section 4 of the paper. Figure C.13 shows the estimated kernel density for the general equilibrium elasticities of employment and residents with respect to the productivity shock across the treated CZs (blue solid and red dashed lines). We also show the 95 percent conÖdence intervals around these estimated kernel densities (gray shading). As CZs are aggregations of counties, there is necessarily less commuting between pairs of CZs than between pairs of counties. Nonetheless, CZs di§er substantially in the extent to which their boundaries capture commuting linkages. Therefore we Önd that there is su¢ cient variation in the importance of commuting networks across CZs to generate substantial heterogeneity in the local employment elasticity, which ranges from just above 0.5 to just over 2.5, a similar range as for the employment elasticity distribution across counties. Again we Önd substantial di§erences between the employment and residents elasticities, with the residents elasticity having less dispersion. Since employment and residents can only di§er through commuting, these Öndings reinforce the importance of commuting in understanding the local response to local economic shocks, even at the more aggregated level of CZs. 

In Table C.6, we provide further evidence on the role of commuting linkages in explaining the heterogeneity in employment elasticities across CZs. This table is analogous to Table 2 in the paper, but reports results for CZs rather than for counties. In Columns (1)-(4), we regress the local employment elasticity on standard empirical controls from the local labor markets literature. Although some of these controls are statistically signiÖcant, we Önd that they are not particularly successful in explaining the variation in employment elasticities. Adding a constant and all these controls yields an R-squared of only just over one quarter in Column (4). Therefore there is considerable variation in local employment elasticities not explained by these standard empirical controls. In contrast, when we include the share of workers that work in i conditional on living in i $( \lambda _ { i i | i } ^ { R } )$ in Column (5) as as summary statistic for openness to commuting, we Önd that this variable is highly statistically signiÖcant, and results in an R-squared of over one half. Including the partial equilibrium elasticities that capture commuting linkages in the model further increases the R-squared to around 0.60, more than double that using the standard controls in Column (4). In the last two columns, we combine these partial equilibrium elasticities with the standard controls used in the Örst four columns. Although some of these standard controls are statistically signiÖcant, we Önd that they add little once we control for the partial equilibrium elasticities. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/b61e26864fb2d89ba79d3904146f4373d7233b83e94c541c1cbe8bb475a07132.jpg)



Figure C.13: Kernel density for the distribution of employment and residents elasticities in response to a productivity shock across CZs


Taken together, these results conÖrm that the use of CZs is an imperfect control for commuting. There remains substantial heterogeneity in employment elasticities across CZs, because they di§er in the extent to which their boundaries are successful in capturing commuting patterns. This heterogeneity in employment elasticities across CZs is not well explained by standard controls from the local labor markets literature. In contrast, consistent with our results for counties above, we Önd that adding a summary statistic of commuting, or the partial equilibrium elasticities from the model, can go a long way in explaining the heterogeneous responses of CZs to productivity shocks. 

We next examine the impact of reductions in the costs of commuting between CZs on the spatial distribution of economic activity. We undertake a counterfactual in which we reduce commuting costs between CZs by the same proportional amount as for counties in our central exercise in Section 6 of the 

<table><tr><td></td><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td><td>7</td><td>8</td><td>9</td></tr><tr><td>Dependent Variable:</td><td colspan="9">Elasticity of Employment</td></tr><tr><td><eq>\log L_i</eq></td><td></td><td>0.025*(0.011)</td><td>0.044(0.022)</td><td>0.002(0.018)</td><td></td><td></td><td></td><td>0.057**(0.017)</td><td>0.055**(0.017)</td></tr><tr><td><eq>\log w_i</eq></td><td></td><td></td><td>-0.037(0.176)</td><td>-0.168(0.136)</td><td></td><td></td><td></td><td>-0.002(0.088)</td><td>-0.020(0.089)</td></tr><tr><td><eq>\log H_i</eq></td><td></td><td></td><td>-0.166**(0.042)</td><td>-0.087(0.049)</td><td></td><td></td><td></td><td>-0.010(0.023)</td><td>-0.011(0.023)</td></tr><tr><td><eq>\log L_{,-i}</eq></td><td></td><td></td><td></td><td>0.081**(0.023)</td><td></td><td></td><td></td><td>-0.038*(0.015)</td><td>-0.040*(0.016)</td></tr><tr><td><eq>\log \bar{w}_{-i}</eq></td><td></td><td></td><td></td><td>0.107(0.146)</td><td></td><td></td><td></td><td>0.012(0.107)</td><td>0.036(0.108)</td></tr><tr><td><eq>\lambda_{ii|i}^R</eq></td><td></td><td></td><td></td><td></td><td>-3.434**(0.216)</td><td></td><td></td><td></td><td></td></tr><tr><td><eq>\sum_{n\in N}(1-\lambda_{Rni})\vartheta_{ni}</eq></td><td></td><td></td><td></td><td></td><td></td><td>8.815**(2.887)</td><td></td><td>9.936*(3.868)</td><td></td></tr><tr><td><eq>\vartheta_{ii}\left(\frac{\lambda_{ii}}{\lambda_{Ri}}-\lambda_{Li}\right)</eq></td><td></td><td></td><td></td><td></td><td></td><td>6.044*(2.916)</td><td></td><td>6.670(3.836)</td><td></td></tr><tr><td><eq>\frac{\partial w_i}{\partial A_i}\frac{A_i}{w_i}</eq></td><td></td><td></td><td></td><td></td><td></td><td>-1.624**(0.194)</td><td></td><td>-0.997**(0.333)</td><td></td></tr><tr><td><eq>\frac{\partial w_i}{\partial A_i}\frac{A_i}{w_i}\cdot\sum_{r\in N}(1-\lambda_{rn|r})\vartheta_{rn}</eq></td><td></td><td></td><td></td><td></td><td></td><td></td><td>1.345**(0.210)</td><td></td><td>2.546**(0.323)</td></tr><tr><td><eq>\frac{\partial w_i}{\partial A_i}\frac{A_i}{w_i}\cdot\vartheta_{ii}\left(\frac{\lambda_{ii}}{\lambda_{Ri}}-\lambda_{Li}\right)</eq></td><td></td><td></td><td></td><td></td><td></td><td></td><td>-1.391**(0.196)</td><td></td><td>-0.680*(0.307)</td></tr><tr><td>Constant</td><td>1.376**(0.031)</td><td>1.098**(0.142)</td><td>2.779(1.810)</td><td>1.747(1.698)</td><td>4.522**(0.194)</td><td>-3.459(2.824)</td><td>2.347**(0.185)</td><td>-4.959(3.920)</td><td>1.387(1.321)</td></tr><tr><td><eq>R^2</eq></td><td>0.00</td><td>0.01</td><td>0.15</td><td>0.27</td><td>0.54</td><td>0.60</td><td>0.59</td><td>0.69</td><td>0.68</td></tr><tr><td>N</td><td>709</td><td>709</td><td>709</td><td>636</td><td>709</td><td>709</td><td>709</td><td>636</td><td>636</td></tr></table>


Note: $\begin{array} { r } { L , - i \equiv \sum _ { n : d _ { n i } \leq 1 2 0 , n \neq i } L _ { n } } \end{array}$ is the total employment in i neighbors whose centroid is no more than 120km away; $\bar { w } _ { - i } \equiv$ Pn:dni120;n6=i $\frac { L _ { n } } { L , - i } w _ { n }$ is the weighted average of their workplace wage. Standard errors are clustered by state; when a CZ overlaps di§erent states, the state that accounts for most of the CZ population is assigned.  denotes signiÖcance at the 5 percent level;  denotes signiÖcance at the 1 percent level. 


Table C.6: Explaining the general equilibrium local employment elasticities to a 5 percent productivity shock for commuting zones (CZs) 

paper $( \hat { B } _ { n i } = 0 . 8 8 )$ . In Figure C.14, we show the proportional change in employment for each CZ against its initial commuting intensity $( L _ { i } / R _ { i } )$ , where $L _ { i } / R _ { i } > 1$ implies that a CZ is a net importer of commuters and $L _ { i } / R _ { i } < 1$ implies that a CZ is a net exporter of commuters. We Önd substantial changes in employment for individual CZs, which range from increases of 10 percent to reductions of 20 percent. Furthermore, these changes in the distribution of employment across CZs are well explained by initial commuting intensity In contrast, in Figure C.15, we show the same proportionate change in employment for each CZ against its initial employment size. We Önd little relationship between the impact of the reduction in commuting costs on employment and initial CZ size. Therefore, these results conÖrm our Öndings for counties that the importance of commuting is by no means restricted to large cities. 

More generally, in Table C.7, we show that it is not easy to proxy for CZ commuting intensity $( L _ { i } / R _ { i } )$ using standard empirical controls from the local labor markets literature. This table is analogous to Table 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/086617281e92114f36ca3a6d59f5b628525336ed1288db514070e2b6d56f9c3c.jpg)



Figure C.14: Counterfactual relative change in commuting zone (CZ) employment $( \hat { L } )$ from median proportional reduction in commuting costs $( \hat { B } _ { n i } = 0 . 8 8 )$ and initial dependence on commuting


![image](https://cdn-mineru.openxlab.org.cn/result/2026-09-21/d6fa5ef8-8ced-435c-a40c-bb395d90ffb9/f4adec98b3d9c459d71542c1fd58e6f1f2583892b55ca87adc0b9edf50966a0e.jpg)



Figure C.15: Counterfactual relative change in commuting zone (CZ) employment (L<sup>^</sup>) from median proportional reduction in commuting costs $( \hat { B } _ { n i } = 0 . 8 8 )$ and initial employment size


B.2 earlier in this web appendix, but reports results for CZs rather than for counties. The Örst four columns show that the levels of either employment $\left( \log L _ { i } \right)$ or residents $( \log { R _ { i } } )$ are strongly related to these standard empirical controls. The Örst column shows that one can account for most of the variation in CZ employment using the number of residents and wages. Column (2) shows a similar result for the number of residents and Columns (3) and (4) show that the results are not a§ected when we add land area, developed-land supply elasticities, employment and wages in surrounding CZs. In contrast, the remaining four columns demonstrate that it is hard to explain the ratio of employment to residents $( L _ { i } / R _ { i } )$ using these same empirical controls. The level of residents, wages, land area, developed-land supply elasticities, employment, and measures of economic activity in surrounding $\mathrm { C Z s , }$ do a poor job in accounting for the 


variation in this ratio. None of the R-squaredís in the last four columns of Table C.7 amounts to more than one third. Therefore, as with our earlier results for counties, we Önd that there is substantial additional information in patterns of commuting that is not captured by the standard empirical controls from the local labor markets literature.


<table><tr><td></td><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td><td>7</td><td>8</td><td>9</td><td>10</td><td>11</td><td>12</td></tr><tr><td>Dep. Variable:</td><td><eq>\log L_i</eq></td><td><eq>\log R_i</eq></td><td><eq>\log L_i</eq></td><td><eq>\log R_i</eq></td><td><eq>L_i/R_i</eq></td><td><eq>L_i/R_i</eq></td><td><eq>L_i/R_i</eq></td><td><eq>L_i/R_i</eq></td><td><eq>L_i/R_i</eq></td><td><eq>L_i/R_i</eq></td><td><eq>L_i/R_i</eq></td><td><eq>L_i/R_i</eq></td></tr><tr><td><eq>\log R_i</eq></td><td>0.991**(0.002)</td><td></td><td>0.992**(0.003)</td><td></td><td></td><td>-0.004**(0.001)</td><td></td><td>0.001(0.002)</td><td></td><td>0.049**(0.010)</td><td></td><td>0.048**(0.015)</td></tr><tr><td><eq>\log w_i</eq></td><td>0.116**(0.024)</td><td></td><td>0.161**(0.027)</td><td></td><td>0.099**(0.020)</td><td></td><td>0.123**(0.022)</td><td></td><td>0.195*(0.074)</td><td></td><td>0.196*(0.074)</td><td></td></tr><tr><td><eq>\log L_i</eq></td><td></td><td>1.001**(0.001)</td><td></td><td>0.993**(0.003)</td><td>-0.006**(0.002)</td><td></td><td>-0.002(0.002)</td><td></td><td>0.025**(0.008)</td><td></td><td>0.029**(0.009)</td><td></td></tr><tr><td><eq>\log \bar{v}_i</eq></td><td></td><td>-0.042**(0.015)</td><td></td><td>-0.057*(0.028)</td><td></td><td>0.057**(0.014)</td><td></td><td>0.078**(0.022)</td><td></td><td>0.014(0.058)</td><td></td><td>0.014(0.058)</td></tr><tr><td><eq>\log H_i</eq></td><td></td><td></td><td>0.006(0.007)</td><td>0.005(0.007)</td><td></td><td></td><td>-0.001(0.006)</td><td>0.001(0.006)</td><td>-0.005(0.011)</td><td>-0.014(0.015)</td><td>-0.009(0.013)</td><td>-0.014(0.018)</td></tr><tr><td><eq>\log R_{,-i}</eq></td><td></td><td></td><td>-0.001(0.003)</td><td></td><td></td><td></td><td>0.638**(0.107)</td><td>0.693**(0.113)</td><td>0.428*(0.168)</td><td>0.499*(0.196)</td><td>0.443*(0.188)</td><td>0.497*(0.195)</td></tr><tr><td><eq>\log \bar{w}_{-i}</eq></td><td></td><td></td><td>-0.149**(0.028)</td><td></td><td></td><td></td><td>1.007**(0.219)</td><td>1.017**(0.222)</td><td>0.934(0.502)</td><td>0.890(0.489)</td><td>0.963(0.559)</td><td>0.887(0.520)</td></tr><tr><td><eq>\log L_{,-i}</eq></td><td></td><td></td><td></td><td>0.010**(0.003)</td><td></td><td></td><td>-0.644**(0.107)</td><td>-0.700**(0.114)</td><td>-0.414*(0.168)</td><td>-0.489*(0.195)</td><td>-0.431*(0.190)</td><td>-0.488*(0.193)</td></tr><tr><td><eq>\log \bar{v}_{-i}</eq></td><td></td><td></td><td></td><td>0.106**(0.029)</td><td></td><td></td><td>-1.119**(0.226)</td><td>-1.111**(0.230)</td><td>-1.176*(0.517)</td><td>-1.097*(0.500)</td><td>-1.207*(0.574)</td><td>-1.093*(0.535)</td></tr><tr><td>Saiz elasticity</td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td>0.005(0.008)</td><td>-0.000(0.008)</td></tr><tr><td>Constant</td><td>-1.146**(0.237)</td><td>0.450**(0.154)</td><td>-0.095(0.159)</td><td>-0.591**(0.183)</td><td>0.007(0.202)</td><td>0.425**(0.147)</td><td>0.974**(0.191)</td><td>1.227**(0.199)</td><td>0.993(0.637)</td><td>2.389**(0.529)</td><td>0.993(0.646)</td><td>2.388**(0.532)</td></tr><tr><td><eq>R^2</eq></td><td>1.00</td><td>1.00</td><td>1.00</td><td>1.00</td><td>0.07</td><td>0.02</td><td>0.31</td><td>0.26</td><td>0.61</td><td>0.52</td><td>0.61</td><td>0.52</td></tr><tr><td>N</td><td>709</td><td>709</td><td>636</td><td>636</td><td>709</td><td>709</td><td>636</td><td>636</td><td>110</td><td>110</td><td>110</td><td>110</td></tr></table>


Note: $\begin{array} { r } { L , - i \equiv \sum _ { n : d _ { n i } \leq 1 2 0 , n \neq i } L _ { n } } \end{array}$ is the total employment in i neighbors whose centroid is no more than 120km away; $\bar { w } _ { - i } \equiv$ $\begin{array} { r } { \sum _ { n : d _ { n i } \leq 1 2 0 , n \neq i } \frac { L _ { n } } { L , - i } w _ { n } } \end{array}$ is the weighted average of their workplace wage. Analogous deÖnitions apply to $R _ { , - i }$ and $\bar { v } _ { - i }$ . Columns 1-8 are unweighted regressions. Columns 9 and 10 repeat the most complete speciÖcations in columns 7 and 8 giving to each CZ a weight proportional to the fraction of the CZís population living in counties where we have data on land supply elasticity; this process excludes from the regressions CZ for which no county has data on land supply elasticity. Columns 11 and 12 then repeat columns 9 and 10 adding the Saiz land supply elasticity as a regressor. Land supply elasticity for a CZ is the populationweighted average of its countiesíland supply elasticities. Standard errors are clustered by state.  denotes signiÖcance at the 5 percent level;  denotes signiÖcance at the 1 percent level. 


## Table C.7: Explaining employment levels and commuting intensity for commuting zones (CZs)

Taking the results of this section as a whole, we Önd that the heterogeneity in commuting linkages across commuting zones (CZs) is su¢cient to generate substantial heterogeneity in local employment elasticities, in response to either productivity shocks or reductions in commuting costs. This heterogeneity is hard to explain with the standard empirical controls from the local labor markets literature, but is well explained by measures of commuting linkages, highlighting the importance of incorporating this commuting information into the analysis of regional economies. 

## D Data Appendix

This section of the web appendix contains further information on the data sources and deÖnitions, as well as additional details about the construction of Ögures and tables. In Section D.1, we discuss the data used for the quantitative analysis of the model in Sections 3-4 of the paper. In Section D.2, we discuss the data used to provide independent evidence in support of the modelís predictions in Section 5 of the paper and Sections C.7 and C.8 of this web appendix. 

## D.1 Quantitative Analysis of the Model (Sections 3-4 of the paper)

## D.1.1 Data Sources and DeÖnitions

In what follows we list the sources and the variable deÖnitions that we use. We consider them understood in the following section on data processing. 

Earnings by Place of Work. This data is taken from the Bureau of Economic Analysis (BEA) website, under Regional Data, Economic ProÖles for all U.S. counties. The BEA deÖnes this variable as "the sum of Wages and Salaries, supplements to wages and salaries and proprietorsí income. [...] Proprietorís income [...] is the current-production income (including income in kind) of sole proprietorships and partnerships and of tax-exempt cooperatives. Corporate directorsí fees are included in proprietors income, but the imputed net rental income of owner-occupants of all dwellings is included in rental income of persons. Proprietorsíincome excludes dividends and monetary interest received by nonÖnancial business and rental incomes received by persons not primarily engaged in the real estate business." The BEA states that earnings by place of work "can be used in the analyses of regional economies as a proxy for the income that is generated from participation in current production". We use the year 2007. 

Total Full-Time and Part-Time Employment (Number of Jobs). This data is taken from the BEA website, under Regional Data, Economic ProÖles for all U.S. counties. The BEA deÖnes this series as an estimate "of the number of jobs, full-time plus part-time, by place of work. Full-time and part-time jobs are counted at equal weight. Employees, sole proprietors, and active partners are included, but unpaid family workers and volunteers are not included. Proprietors employment consists of the number of sole proprietorships and the number of partners in partnerships. [...] The proprietors employment portion of the series [...] is more nearly by place of residence because, for nonfarm sole proprietorships, the estimates are based on IRS tax data that reáect the address from which the proprietorís individual tax return is Öled, which is usually the proprietorís residence. The nonfarm partnership portion of the proprietors employment series reáects the tax-Öling address of the partnership, which may be either the residence of one of the partners or the business address of the partnership." We use the year 2007. 

County-to-County Worker Flows. This data contains county-level tabulations of the workforce "residence-to-workplace" commuting áows from the American Community Survey (ACS) 2006-2010 5-year Öle. The ACS asks respondents in the workforce about their principal workplace location during the reference week. People who worked at more than one location are asked to report the location at which they worked the greatest number of hours. We use data for all the 50 States and the District of Columbia. 

County Land Area, County Centroids. This data comes from the 2010 Census Gazetteer Files. Land area is geographical land area. When we need to aggregate counties (see below), the geographical land area is the sum of that for the aggregated counties, and the centroid of the new county formed by the aggregation is computed using spatial analysis software. In Subsection 4.2 of the paper, we develop an extension to allow for a heterogeneous positive supply elasticity for developed land following Saiz (2010). 

County Median Housing Values. This data reports the countyís median value of owner-occupied housing units from the American Community Survey 2009-2013 5-year Öle. 

Commodity Flows among CFS Area. We use the 2007 Origin-Destination Files of the Commodity Flow Survey for internal trade áows of all merchandise among the 123 Commodity Flow Survey areas in the United States. 

Share of county employment in manufacturing. We use the County Business Pattern Öle for the year 2007. We use the information on total employment, and employment in manufacturing only. For some counties, employment is suppressed to preserve non-disclosure of individual information, and employment is only reported as a range. In those cases, we proceed as follow. We Örst use the information on the Örm-size distribution, reported for all cases, to narrow the plausible employment range in the cell. We run these regressions separately for employment in manufacturing and total employment. We then use this estimated relationship to predict the employment level where the data only reports information on the Örm size-distribution. Whenever the predicted employment lies outside the range identiÖed above, we use the employment at the relevant corner of the range. 

## D.1.2 Initial Data Processing

We start by assigning to each workplace county in the County-to-County Worker Flows data, information on the Earnings by Place of Work and the Number of Jobs. Note that the commuting data contains 3,143 counties while the BEA data contains 3,111 counties. This happens because, for example, some independen cities in Virginia for which we have separate data on commuting are included in the surrounding county in the BEA data. We make the two sources consistent by aggregating the relevant commuting áows by origin-destination, and so we always work with 3,111 counties. 

The ACS data reports some unrealistically long commutes, which arise for example for itinerant professions. We call these áows "business trips" and we remove them as follow. We measure the distance between counties as the distance between their centroids computed using the Haversine formula. We start by assuming that no commute can be longer than 120km: hence, áows with distances longer than 120km are assumed to only be business trips, while áows with distances less than or equal to 120km are a mix business trips and actual commuting. We choose the 120km threshold based on a change in slope of the relationship between log commuters and log distance at this distance threshold. To split total travellers into commuters and business travellers, we write the identit $\tilde { \lambda } _ { i j } = \psi _ { i j } ^ { B } \tilde { \lambda } _ { i j } ^ { B }$ ;where $\tilde { \lambda } _ { i j }$ is total travellers, $\tilde { \lambda } _ { i j } ^ { B }$ is business travellers, $\tilde { \lambda } _ { i j } ^ { C }$ is commuters, and $\psi _ { i j }$ is deÖned as an identity as the ratio of total travellers to business travellers: 

$$
\psi_ {i j} = \frac {\tilde {\lambda} _ {i j} ^ {C} + \tilde {\lambda} _ {i j} ^ {B}}{\tilde {\lambda} _ {i j} ^ {B}}.
$$

We assume that business travel follows the gravity equation $\tilde { \lambda } _ { i j } ^ { B } = S _ { i } M _ { j } \mathrm { d i s t } _ { i j } ^ { \delta ^ { B } } u _ { i j }$ ;where $S _ { i }$ is a residence Öxed e§ect, $M _ { j }$ is a workplace Öxed e§ect, dis $\mathfrak { t } _ { i j }$ is bilateral distance, and $u _ { i j }$ is a stochastic error. We 

assume that $\psi _ { i j }$ takes the following form: 

$$
\psi_ {i j} = \left\{ \begin{array}{c c} 1 & \mathrm{dist} _ {i j} > \bar {d} \\ \gamma \mathrm{dist} _ {i j} ^ {\delta^ {C}} & \mathrm{dist} _ {i j} \leq \bar {d} \end{array} \right.,
$$

where we expect $\gamma > 1$ and $\delta _ { C } < 0$ . Therefore we have the following gravity equation for total travellers: 

$$
\ln \tilde {\lambda} _ {i j} = \ln S _ {i} + \ln M _ {j} + \gamma \mathbb {I} _ {i j} + (\delta_ {B} + \delta_ {C} \mathbb {I} _ {i j}) \ln \mathrm{dist} _ {i j} + u _ {i j},\tag{D.1}
$$

where $\mathbb { I } _ { i j }$ is an indicator variable that is one if dis $\bar { \iota } _ { i j } \leq \bar { d }$ and zero otherwise. Estimating the above equation for total travellers, we can generate the predicted share of commuters as: 

$$
\hat {s} _ {i j} ^ {C} = 1 - \frac {\widehat {\hat {\lambda}} _ {i j} ^ {B}}{\widehat {\hat {\lambda}} _ {i j}} = 1 - \frac {\hat {S} _ {i} \hat {M} _ {j} \mathrm{dist} _ {i j} ^ {\hat {\delta} _ {B}}}{\widehat {\hat {\lambda}} _ {i j}},
$$

where $\widehat { \tilde { \lambda } } _ { i j } = \exp \left( \ln \widehat { \tilde { \lambda } } _ { i j } \right)$ are the Ötted values from gravity (D.1). Note that this predicted share satisÖes the requirements that (a) commuters are zero beyond the threshold <sup></sup>d, (b) the predicted share of commuters always lies in between zero and one, (c) commuters, business travellers and total travellers all satisfy gravity. Note also that since the regression cannot be run on áows internal to a county $\tilde { \lambda } _ { i i } ,$ we set $\hat { s } _ { i i } ^ { C } = 1 ~ ( \mathrm { i . e . }$ áows of agents who live and work in the same county are assumed to contain no business trips). Therefore we can construct commuting áows as: a 

$$
\hat {\tilde {\lambda}} _ {i j} ^ {C} = \hat {s} _ {i j} ^ {C} \tilde {\lambda} _ {i j}.
$$

The total business trips originating from residence i are then $\begin{array} { r } { \sum _ { j } \left( 1 - \hat { s } _ { i j } ^ { C } \right) \tilde { \lambda } _ { i j } } \end{array}$ . For any residence i, we reimpute these business trips across destinations $j$ in proportion to the estimated workplace composition of the residence $i , \hat { \tilde { \lambda } } _ { i j } ^ { C } / \sum _ { i } \hat { \tilde { \lambda } } _ { i j } ^ { C }$ . The total employment (and average wage) in a county in the initial equilibrium is taken from the BEA, while total residents (and average residential income) in a county are reconstructed using the estimated residence composition of each workplace. Table 1, Figure 3, and all the results in the paper are based on these ìcleanedî commuting áows and initial equilibrium values. 

Whenever necessary, we allow for expenditure imbalances across counties. We compute these imbalances as follows. We start from the CFS trade áows. The total sales of a CFS area anywhere must correspond, in a model with only labor (such as the one in this paper), to total payments to workers employed in the area. We rescale the total sales from a CFS area to the value of the total wage bill from the BEA data.<sup>9</sup> For any origin CFS, we keep the destination composition of sales as implied by the CFS bilateral áows. This procedure gives us, for any CFS, total expenditures and total sales consistent with the total labor payments in the economy. We compute the deÖcit of any CFS area by subtracting total sales from total expenditure. We apportion this deÖcit across all the counties in the CFS in proportion to the total residential income of the county, as computed above. The total expenditure of the county in the initial equilibrium is always total residential income plus deÖcit. In any counterfactual equilibrium, the dollar value of the deÖcit is kept Öxed. 

## D.1.3 Further Information on Figures and Tables

We now report additional technical details related to the data sources and manipulation for some of the tables and Ögures in the paper or this web appendix. 

Table 1. The table reports statistics on the out-degree distribution (Örst and third row) and in-degree distribution of the fraction of commuters across counties. Commuting áows are cleaned with the procedure described above. The correspondence between counties and commuting zones is taken from the Economic Research Service of the United States Department of Agriculture.<sup>10</sup> 

Figure 1. This Ögure reports kernel densities of the distribution of the share of a countyís residents working in their county of residence for 4 decades. Data on the share of residents working in the county are constructed from the ICPSR Study $7 7 3 6 ^ { 1 1 }$ , and the 1983, 1994 and 2000 editions of the ìCounty and City Data Bookîpublished by the U.S. Department of Commerce. 

Figure B.1. This Ögure reports a scatterplot of the log trade áows among CFS areas against log distance between these areas, after removing origin and destination Öxed e§ects. The distance between CFS areas is the average distance travelled by shipments, computed dividing the total ton-miles travelled by the total tons shipped, as reported in the CFS data. Whenever this distance cannot be computed (in about 1/3 of the áows) we supplement it with an estimated distance as follows. We compute the centroids of CFS areas using the Freight Analysis Framework Regions shape-Öles provided by the Bureau of Transportation Statistics<sup>12</sup> and bilateral distances among these centroids using the Haversine formula. We then regress the actual distance shipped on these centroid-based distances, in logs, and Önd strong predictive power (slope of 1.012, $R ^ { 2 } = 0 . 9 5 )$ . We use the predicted distances from this regression for áows where the average distance shipped cannot be computed. If we restrict our sample to only áows for which the distance can be computed directly, we Önd a slope of -1.23, and $R ^ { 2 }$ of 0.82 (similar to the ones used in the paper of -1.29 and 0.83, respectively). 

Figure B.2. This Ögure reports a scatterplot of expenditure shares across CFS areas in the data and the model-implied expenditure shares after recovering the productivity of each county, with the procedure described in Section 3.1 of the paper. Both the estimated productivities and the implied trade shares are calculated using the expenditure of a county allowing for deÖcits computed as above. 

Figure B.3. This Ögure reports kernel densities analogous to Figure 1 that are weighted by the number of residents in each county, and it shares with that Figure the data source. 

Figure B.5. This Ögure reports a scatterplot of log commuting áows against log distance between countyís centroids after removing residence and workplace Öxed e§ects. The commuting áows used in the regression are cleaned of the business trips as described above. 

Figure C.1. This Ögure reports a scatterplot of log of land price, as computed from the model, and the County Median Housing Value from the ACS. To compute the price of land in the model we use residents expenditure allowing for trade deÖcits. For counties that are aggregated at the BEA level (see above), we compute the population weighted average of the median values. 

## D.2 Additional Empirical Evidence

We now discuss the data sources and deÖnitions for the independent evidence in support of the predictions of the model in Section 5 of the paper and Sections C.7 and C.8 of this web appendix. 

Commuting Data (Section 5 of the paper and Sections C.7 and C.8 of this web appendix). We construct three bilateral commuting matrices for 1990, 2000 and 2006-2010. We use these matrices for both Section 5 of the paper and C.8 of this web appendix. Our bilateral commuting data comes from the County-to-County Worker Flows tabulation Öles based on the U.S. Census (for years 1990 and 2000) and American Community Survey (for 2006-2010). We construct commuting áows following the same procedure indicated in Section D.1.2 for the contiguous United States. We compute distances between county centroids using the coordinates provided in the corresponding years of the Census Gazetteer Öles. To construct a balanced panel of counties over time, some aggregation of counties is needed, and we end up with a cross-section of 3,108 spatial units for all three years. 

Million Dollar Plants (Section 5 of the paper and Section C.7 of this web appendix). We use the full list of 82 plants openings gathered by Greenstone and Moretti (2004) from the Journal Site Selection. For each county, yearly workplace employment is taken from the Bureau of Economic Analysis, County Economic ProÖles (Table CA30). In particular, we use the measure of Wage and Salary Employment (data line 250). This measure includes ìAll jobs for which wages and salaries are paid are countedî, which cover all industries covered by Unemployment Insurance, plus adjustments for industries not fully covered by Unemployment Insurance as detailed in the ìLocal Area Personal Income Methodologyî(November 2016) from the BEA. In weighted regressions, the population at the beginning of the sample for each county also comes from the same BEA source (line 100). For each county, the measured own commuting share is for the closest available year to the plant opening date from the commuting data discussed at the beginning of this subsection. For all 82 plant openings, the closest available year is 1990. To control for industry-year Öxed e§ects, we assign industries to cases using the reported industry for each case from Appendix Table 2 in Greenstone and Moretti (2004). Cases are classiÖed into 5 broad industries: Manufacturing (63 cases), Financial (1 case), Services (6 cases), Trade (4 cases), Transportation and Utilities (8 cases). 

Shift-share Decompositions (Section C.8 of this web appendix). We use the bilateral matrix for 2006-2010 for the cross-section decomposition and the bilateral matrices for 1990 and 2006-2010 for the time-series decomposition. 

## References



Allen, Treb, Costas Arkolakis and Xiangliang Li (2016) ìOn the Existence and Uniqueness of Trade Equilibria,îYale University, mimeograph. 





Armington, Paul S. (1969) ìA Theory of Demand for Products Distinguished by Place of Production, IMF Sta§ Papers, 16(1), 159-178. 





Berry, Steven, James Levinsohn and Ariel Pakes (1995) ìAutomobile Prices in Market Equilibrium, Econometrica, 63(4), 841-890. 





Caliendo, Lorenzo, Fernando Parro, Esteban Rossi-Hansberg, and Pierre-David Sarte (2014) ìThe Impact of Regional and Sectoral Productivity Changes on the U.S. Economy,îNBER Working Paper, 20168. 





Eaton, Jonathan and Kortum, Samuel (2002) ìTechnology, Geography, and Trade,îEconometrica, 70(5), 1741-1779. 





Erlander, Stewart and Neil F. Stewart (1990) The Gravity Model in Transportation Analysis: Theory and Extensions, Utrecht: VSP. 





Grubel, Herbert G. and Peter Lloyd (1975) Intra-Industry Trade: The Theory and Measurement of International Trade in Di§erentiated Products, London: MacMillan. 





Hsieh, Chang-Tai and Enrico Moretti (2017) ìHousing Constraints and Spatial Misallocation,îUniversity of California, Berkeley, mimeogarph. 





Krugman, Paul and Anthony J. Venables (1995) ìGlobalization and the Inequality of Nations,îQuarterly Journal of Economics, 110(4), 857-880. 





McFadden, Daniel (1974) ìThe Measurement of Urban Travel Demand,î Journal of Public Economics, 3(4), 303-328. 





McFadden, Daniel and Kenneth Train (2000) ìMixed MNL Models of Discrete Response,î Journal of Applied Econometrics, 15, 447-470. 





Saiz, Albert (2010) ìThe Geographic Determinants of Housing Supply,îQuarterly Journal of Economics, 125(3), 1253-1296. 





Sen, Ashish K. and Tony E. Smith (1995) Gravity Models of Spatial Interaction Behavior, New York: Springer Verlag. 

