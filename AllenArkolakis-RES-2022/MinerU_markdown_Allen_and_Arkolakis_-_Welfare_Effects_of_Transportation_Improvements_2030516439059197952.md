# The Welfare Effects of Transportation Infrastructure Improvements∗

Treb Allen 

Costas Arkolakis 

Dartmouth and NBER 

Yale and NBER 

First draft: August 2016 

This version: October 2021 

# Abstract

Each year in the U.S., hundreds of billions of dollars are spent on transportation infrastructure and billions of hours are lost in traffic. We develop a quantitative general equilibrium spatial framework featuring endogenous transportation costs and traffic congestion and apply it to evaluate the welfare impact of transportation infrastructure improvements. Our approach yields analytical expressions for transportation costs between any two locations, the traffic along each link of the transportation network, and the equilibrium distribution of economic activity across the economy, each as a function of the underlying quality of infrastructure and the strength of traffic congestion. We characterize the properties of such an equilibrium and show how the framework can be combined with traffic data to evaluate the impact of improving any segment of the infrastructure network. Applying our framework to both the U.S. highway network and the Seattle road network, we find highly variable returns to investment across different links in the respective transportation networks, highlighting the importance of well-targeted infrastructure investment. 

# 1 Introduction

More than a trillion dollars is spent on transportation infrastructure across the world each year (Lefevre, Leipziger, and Raifman, 2014). In the U.S. alone – where annual spending on highways exceeds $150 billion – the average driver spends an average of 42 hours a year in traffic, generating economic losses exceeding these direct costs (ASCE, 2017). Evaluating the impact of infrastructure investments in the presence of such traffic congestion is difficult. On the one hand, improvements to one part of the infrastructure network causes drivers to alter their routes, changing traffic patterns and congestion throughout the network. On the other hand, changes in traffic patterns affects the spatial distribution of economic activity, as individuals re-optimize where to live, work, and/or consume. But as the spatial distribution of economic activity determines the underlying traffic patterns, these two hands are intricately intertwined, resulting in a complex feedback loop between routing, traffic, congestion, and the spatial distribution of economic activity. 

We develop a new tractable spatial framework featuring endogenous transportation costs and traffic congestion and apply it to evaluate the welfare impact of transportation infrastructure improvements. We embed a route choice problem into two spatial models where the cost of traversing a particular link depends on the equilibrium amount of traffic on that link. Our approach yields analytical expressions for transportation costs between any two locations, the traffic along each link of the transportation network, and the equilibrium distribution of economic activity across the economy. We characterize the properties of such an equilibrium, highlighting how the presence of traffic congestion shapes those properties. We then show how the framework can be combined with readily available traffic data to evaluate the welfare impact of improving any segment of the infrastructure network. Finally, we evaluate the welfare impact in two settings: (1) the U.S. highway network; and (2) the Seattle road network. In both cases we find on average positive but highly variable returns to investment, showing the importance of well-targeted infrastructure investment. 

Our framework begins with a modest departure from two widely used quantitative general equilibrium models: an economic geography model where agents choose a location to live (as in Allen and Arkolakis (2014)) and engage in trade between locations (as in Eaton and Kortum (2002)), and an urban model where agents choose where to live and where to work within a city (as in Ahlfeldt, Redding, Sturm, and Wolf (2015)). In Eaton and Kortum (2002), it is assumed that while each location has a idiosyncratic productivity for producing each type of good, the transportation technology is identical for all goods. Similarly, in Ahlfeldt, Redding, Sturm, and Wolf (2015), while it is assumed that each individual has idiosyncratic preferences for each home-work pair of locations, all individuals incur the same 

transportation costs when commuting from home to work. In our framework, we allow for transportation costs in both models to also be subject to idiosyncrasies at the route-level. As a result, simultaneous to their choice of where to purchase goods (in the economic geography model) or where to live and work (in the urban model), agents also choose an optimal route through the transportation network. 

This departure allows us to derive an analytical expression for the endogenous transportation costs between all pairs of locations as a function of the transportation network. It also allows us to derive an analytical expression for the equilibrium traffic along a link. This expression takes an appealing “gravity” form, where traffic depends only on the cost of travel along the link and the economic conditions at the beginning and end of the link. Those economic conditions turn out to be the familiar market access terms (see e.g. Anderson and Van Wincoop (2003); Redding and Venables (2004)) – the “inward” market access at the start of the link and the “outward” market access at the end – highlighting the close relationship between equilibrium traffic flows and the equilibrium distribution of economic activity. It is this close relationship that allows us to tractably introduce traffic congestion, which we do so in the spirit Vickrey (1967), by assuming transportation costs of traversing a link depend on both the underlying infrastructure and amount of traffic along the link. 

Ultimately, we can express the equilibrium distribution of economic activity solely as a function of the underlying infrastructure matrix, the geographic fundamentals of each location, and four model elasticities, one of which is new (the traffic congestion elasticity) and three of which are not (a trade/commuting elasticity, a productivity externality, and an amenity externality). While the mathematical structure the equilibrium system takes has to our knowledge not been studied before, we prove an equilibrium will exist and provide conditions under which it will be unique. The new mathematical structure also yields new implications: most notably, the presence of traffic congestion implies that the equilibrium is no longer scale invariant. Increasing the size of an economy results in disproportionate changes in bilateral transportation costs due to changes in traffic congestion, reshaping the equilibrium distribution of economic activity. 

We then turn to the question of how to apply our framework empirically. We begin by developing a few new tools. First, we derive an analytical relationship between traffic flows along a network and bilateral trade/commuting flows between an origin and destination; in contexts such as our own where we observe both, this serves as a model validation check. Second, we show that the “exact-hat” approach of conducting counterfactuals (see Dekle, Eaton, and Kortum (2008); Costinot and Rodr´ıguez-Clare (2014); Redding and Rossi-Hansberg (2017)) can be applied to our framework, albeit using (readily available) traffic data rather than harder to observe bilateral trade/commuting data. Third, we provide conditions 

under which one can recover the necessary traffic congestion elasticity from a regression of speed of travel on traffic, where the traffic gravity equation provides guidance in the search for an appropriate instrument for traffic. 

Finally, we calculate the welfare impact of transportation infrastructure improvements in two settings: (1) the U.S. highway network (using the economic geography variant of the framework); and (2) the Seattle road network (using the urban variant). In both cases, we begin by showing that the observed network of traffic flows, appropriately inverted through the lens of the model, does a good job predicting the observed matrix of trade and commuting flows, respectively. We then estimate the strength of traffic congestion, finding in both cases substantial traffic congestion. We proceed by estimating the welfare elasticity of improving each link on each road network. We find highly variable elasticities across different links, with the greatest gains in the densest areas of economic activity and at choke-points in the network. Here, traffic congestion plays a particularly important role, as there is only a modest positive correlation between these welfare elasticities and those that one would calculate in a standard model ignoring congestion forces. 

Finally, we combine our welfare elasticities with detailed cost estimates of improving each link (which depends on the number of lane-miles needed to be added as well as the geographic topography and the density of economic activity along the link) to construct an estimate of the return on investment for each link. For the U.S. highway network, we estimate an average annual return on investment of 108%; for the Seattle road network that figure is 16%. Both averages, however, belie substantial heterogeneity across links. For the U.S. highway network, the returns on investment for a handful of highways serving as connectors just outside major metropolitan areas exceed $4 0 0 \%$ ; in Seattle, a number of links surrounding downtown have annualized returns exceeding $6 0 \%$ . Conversely, a substantial fraction of U.S. highway links (mainly through the mountain west) and nearly half the links in Seattle are estimated to have a negative return on investment. Taken together, these results highlight the importance of targeting infrastructure improvements to the appropriate locations in the infrastructure network. 

The primary contribution of the paper is to develop a quantitative general equilibrium spatial framework that incorporates traffic congestion and can be applied to empirically evaluate the welfare impact of transportation infrastructure improvements. In doing so, we seek to connect two related – but thus far distinct – literatures. 

The first literature seeks to understand the impacts of infrastructure improvements on the distribution of economic activity. This literature is mostly the domain of spatial economists; early examples include Fogel (1962, 1964); recent quantitative work on the subject that incorporates rich geographies and general equilibrium linkages across locations include Don-

aldson (2018), Allen and Arkolakis (2014), Donaldson and Hornbeck (2016) in an inter-city context Ahlfeldt, Redding, Sturm, and Wolf (2015), Tsivanidis (2018), Heblich, Redding, and Sturm (2020) in an intra-city context, and Monte, Redding, and Rossi-Hansberg (2018) combining intra-city and inter-city analyses; Redding and Turner (2015) and Redding and Rossi-Hansberg (2017) offer excellent reviews. While the details of these models vary, a unifying characteristic is that the transportation costs are treated as exogenous model parameters (usually determined by the least cost route, as computed using Dijkstra’s algorithm or the “Fast Marching Method” pioneered by Osher and Sethian (1988) and Tsitsiklis (1995)). As a result, this literature abstracts from the effect of infrastructure improvements on how changes in the use of the transportation network affects the transportation costs themselves through traffic congestion. 

Relative to this literature, we make two contributions: first, we provide an analytical relationship between the transportation network and the bilateral costs of travel through the network, obviating the need to rely on computational methods. Second (and more importantly), we allow the transportation costs to respond endogenously through traffic congestion to changes in the distribution of economic activity. This force has been identified as empirically relevant (see Duranton and Turner (2011)) but thus far has been absent in such quantitative modeling. Our analysis retains the key analytical benefits of that previous work but also provides a comprehensive framework to analyze the effects of traffic both theoretically and empirically. 

The second literature seeks to understand the impacts of infrastructure improvements on the transportation network. This literature is mostly the domain of transportation economics; early examples include Beckmann, McGuire, and Winsten (1955) and seminal textbook of Sheffi (1985); recent work on the subject includes Bell (1995), Akamatsu (1996), De Palma, Kilani, and Lindsey (2005), Eluru, Pinjari, Guo, Sener, Srinivasan, Copperman, and Bhat (2008), Mattsson, Weibull, and Lindberg (2014); Galichon (2016) provides a comprehensive theoretical treatment and Chapter 10 of De Palma, Lindsey, Quinet, and Vickerman (2011) provides an excellent review. While the details of these models vary, a unifying characteristic is that the economic activity at each node in the network is taken as given, so the literature abstracts from how changes in the transportation costs affects this distribution of economic activity. 

Relative to this literature, we also make two contributions: first, we provide an analytical solution for the equilibrium traffic along each link in the network that highlights the close relationship between traffic and the equilibrium distribution of economic activity. Second (and more importantly), we allow infrastructure improvements to affect traffic not only through changing route choices (and congestion) on the network, but also through the 

resulting equilibrium changes in the distribution of economic activity. 

Most closely related to this paper is parallel work by Fajgelbaum and Schaal (2020), who characterize the optimal transportation network in a similarly rich geography and also in the presence of traffic congestion. In that important work, the focus is on an efficient equilibrium of a flexible spatial model, as it is assumed that the social planner can implement optimal Pigouvian taxes to offset the externalities created by traffic congestion. Our focus, instead, is on the competitive equilibrium of constant elasticity quantitative spatial models where the presence of productivity and amenity externalities and/or traffic congestion given the absence of congestion tolls means the equilibrium is (generically) inefficient. Relative to Fajgelbaum and Schaal (2020), a separate contribution is that the analytical tractability of the framework developed here facilitates the use of many of the tools developed previously by the quantitative spatial literature, such as the ability to evaluate the welfare impact of infrastructure improvements using readily available traffic data and the use of “exact hat algebra” methodology to compute counterfactuals.1 

The remainder of the paper proceeds as follows. In the next section, we incorporate the routing choice of agents in economic geography and urban variants of the framework. In Section 3, we provide analytical expressions for the endogenous transportation costs and traffic flows in the presence of traffic congestion. In Section 4, we combine the results of the previous sections to characterize the equilibrium distribution of economic activity and traffic. In Section 5, we develop a set of tools for applying the framework empirically. In Section 6, we implement these tools to examine the welfare impacts of improvements to the U.S. highway network and the Seattle road network. Section 7 concludes. 

# 2 Optimal Routing in Two Spatial Models

In this section, we embed an optimal routing problem into two quantitative spatial models: an economic geography model (where goods are traded between locations subject to trade costs) and an urban model (where workers commute between locations subject to commuting costs). We show that both models yield identical expressions for the endogenous transportation costs, and mathematically identical equilibrium conditions as a function of these costs. This allows us to derive analytical expressions for costs, traffic, and congestion in both frameworks, 

a task we undertake in Section 3; we refer the interested reader to Online Appendix B for detailed derivations of the results that follow in this section. 

For both models, we posit the following geography. Suppose the economy consists of a finite number of locations $i \in \{ 1 , . . . , N \} \equiv \mathcal { N }$ arrayed on a network and inhabited by $L$ individuals. Mathematically, this network is represented by an $N \times N$ matrix $\mathbf { T } = [ t _ { k l } \ge 1 ]$ , where $t _ { k l }$ indicates the (ad valorem) cost incurred from moving directly from $k$ to $l$ along a link (if no link between $k$ and $l$ exists, then $t _ { k l } = \infty$ ).2 We refer to $\mathbf { T }$ as the transportation network and emphasize that it is endogenous and will depend on the equilibrium traffic congestion. 

Moving goods (in the economic geography model) or people (in the commuting model) from an origin $i$ to a destination $j$ entails taking a route $r$ through the network. Mathematically, $r$ is a sequence of locations beginning with location $i$ and ending with location $j$ , i.e. $r \equiv \{ i = r _ { 0 } , r _ { 1 } , . . . , r _ { K } = j \}$ , where $K$ is the number of links crossed on the route, i.e. the length of route $r$ . Because iceberg costs are multiplicative, the total costs incurred from moving from $i$ to $j$ along route $r$ of length $K$ is then $\textstyle \prod _ { k = 1 } ^ { K } t _ { r _ { k - 1 } , r _ { k } }$ .3 Let $\Re _ { i j }$ denote the set of all the (countably infinite) possible routes from $i$ to $j$ . 

# 2.1 An Economic Geography Model with Optimal Routing

We first embed a routing framework into an economic geography model where goods are traded across locations and labor is mobile, as in Allen and Arkolakis (2014) and Redding (2016). 

# 2.1.1 Setup

An individual residing in location $i$ supplies her endowed unit of labor inelastically for the production and shipment of goods, for which she receives a wage $w _ { i }$ and from which she purchases quantities of a continuum of consumption goods $\nu \in \left[ 0 , 1 \right]$ with constant elasticity of substitution (CES) preferences with elasticity of substitution $\sigma \geq 0$ . Labor is the only factor used in the production and shipment of goods. Let $Y ^ { W }$ and $L$ denote the total income and total labor endowment in the economy, respectively; in what follows, we choose average 

per-capita income as our numeraire, i.e. $Y ^ { W } / L = 1$ , which implies that the value of trade is measured in average units of labor. 

Each location $i \in \mathcal N$ is endowed with a constant returns to scale technology for producing and shipping each good $\nu \in [ 0 , 1 ]$ to each destination $j \in \mathcal N$ along each route $r ~ \in ~ \Re _ { i j }$ , which is subject to idiosyncratic productivity shocks $\varepsilon _ { i j , r } \left( \nu \right)$ , meant to capture the various uncertainties that production and shipping are subject to. Under perfect competition the price of good $\nu$ in destination $j \in \mathcal N$ from origin $i \in \mathcal N$ along route $r \in \Re _ { i j }$ is 

$$
p _ {i j, r} (\nu) = w _ {i} \frac {\prod_ {k = 1} ^ {K} t _ {r _ {k - 1} , r _ {k}}}{\varepsilon_ {i j , r} (\nu)}.
$$

Individuals in destination $j$ then purchase each good $\nu \in \left[ 0 , 1 \right]$ from the cheapest source (i.e. location-route). Following Eaton and Kortum (2002), we assume $\varepsilon _ { i j , r } \left( \nu \right)$ is independently and identically Frechet distributed across routes and goods distributed with scale parameter $1 / A _ { i }$ , where $A _ { i }$ captures an origin-specific efficiency, and shape parameter $\theta$ , which regulates the (inverse of) shock dispersion.4 

The main innovation in our setup is that individuals choose both a location and route to source each good (rather than just a location). But why would a consumer not simply choose to purchase the goods from the cheapest source along the least cost route? Some of the value of this choice of modeling arises from the great tractability it yields below. Yet this added “noise” is also is plausible in the presence of traffic congestion, as there will be many alternative routes that yield approximately the same costs.5 If all consumers were to use the least cost route, then infinitesimal deviations from Mogridge’s hypothesis would result in large changes in agents’ route choice; empirically, an infinitely elastic route choice is unrealistic; theoretically, it would lead to a nightmare of corner solutions. Avoiding corner solutions by adding such noise is cited as the original impetus for the Eaton and Kortum (2002) framework and the Frechet assumption allows us to further retain the tractability and extend the analytical solutions of that framework in the presence of traffic.6 

A related concern is with the assumption that agents simultaneously choose the location that sources the good and the route over which it is supplied. Should agents not first choose where to purchase a good and then decide how to ship it? It turns out the timing assumption is not crucial: one can construct a model with just such a timing assumption that is formally isomorphic to the framework presented here (see Online Appendix D.2). Instead, what is enormously helpful (and which the simultaneous choice over locations and routes ensures) is that agents’ elasticities of substitution among locations and among routes are the same. Deviations from this assumption – while computationally straightforward – come so at the loss of substantial analytical tractability and ensuing economic insight.7 

We further allow for the possibility that productivities and amenities potentially depend on the measure of workers in a given location as follows: 

$$
A _ {i} = \bar {A} _ {i} L _ {i} ^ {\alpha}, u _ {i} = \bar {u} _ {i} L _ {i} ^ {\beta}, \qquad (1)
$$

where $A _ { i } > 0$ and $\bar { u } _ { i } > 0$ are the local geography of productivity and amenities, and $\alpha , \beta \in \mathbb { R }$ govern the strength of the productivity and amenity externalities, respectively. As noted in Allen and Arkolakis (2014), the presence of productivity and amenity spillovers create formal isomorphisms between a large set of economic geography models and also play an important role in determining the qualitative and quantitative implications of the model. For example the parameter $\alpha$ can be considered as capturing entry externalities as in Krugman (1991), which lead to more concentration of economic activity, and the parameter $\beta$ negative amenity spillovers or the presence of a housing market, which lead to dispersion of economic activity. We will contrast the implications of these (now standard) spillovers to the (new) traffic congestion spillovers below. 

# 2.1.2 An analytical expression for transportation costs

We now characterize the fraction and value of goods shipped on each route between each origin and destination. Given the Frechet assumption, the probability that $j \in \mathcal N$ purchases 

good $\nu \in [ 0 , 1 ]$ from $i \in \mathcal N$ along route $r \in \Re _ { i j }$ , $\pi _ { i j , r }$ , can be written as: 

$$
\pi_ {i j, r} = \frac {(w _ {i} / A _ {i}) ^ {- \theta} \left(\prod_ {l = 1} ^ {K} t _ {r _ {l - 1} , r _ {l}} ^ {- \theta}\right)}{\sum_ {k \in \mathcal {N}} (w _ {k} / A _ {k}) ^ {- \theta} \sum_ {r ^ {\prime} \in \Re_ {k j}} \prod_ {l = 1} ^ {K} t _ {r _ {l - 1} ^ {\prime} , r _ {l} ^ {\prime}} ^ {- \theta}}. \tag {2}
$$

To determine the total value of goods shipped from $i \in \mathcal N$ to $j \in \mathcal N$ , $X _ { i j }$ , we sum across all routes, recalling from Eaton and Kortum (2002) that the expenditure shares are equal to the probability of purchasing a good: 

$$
X _ {i j} = \sum_ {r \in \Re_ {i j}} \pi_ {i j, r} E _ {j} = \frac {\tau_ {i j} ^ {- \theta} (w _ {i} / A _ {i}) ^ {- \theta}}{\sum_ {k \in \mathcal {N}} \tau_ {k j} ^ {- \theta} (w _ {k} / A _ {k}) ^ {- \theta}} E _ {j}, \tag {3}
$$

where: 

$$
\tau_ {i j} \equiv \left(\sum_ {r \in \Re_ {i j}} \left(\prod_ {l = 1} ^ {K} t _ {r _ {l - 1}, r _ {l}} ^ {- \theta}\right)\right) ^ {- \frac {1}{\theta}} \tag {4}
$$

is the transportation costs from $i$ to $j$ . Note that expression (3) is identical to that of Eaton and Kortum (2002); however, rather than the transportation cost $\tau _ { i j }$ being taken as given, here it is determined by the least cost routing problem through the (endogenous) transportation network. 

# 2.1.3 Market Access and Gravity

While (3) provides an analytical expression for the value of bilateral trade flows, it turns out it is convenient for what follows to express it in market access terms, as in Anderson and Van Wincoop (2003) and Redding and Venables (2004). To do so, we first impose two equilibrium market clearing conditions: (1) total income $Y _ { i }$ in each location is is equal to its total sales; and (2) total expenditure $E _ { i }$ in each location is equal to its total purchases: 

$$
Y _ {i} = \sum_ {j = 1} ^ {N} X _ {i j}, E _ {i} = \sum_ {j = 1} ^ {N} X _ {j i}. \tag {5}
$$

We can re-write the gravity equation (3) as follows: 

$$
X _ {i j} = \tau_ {i j} ^ {- \theta} \times \frac {Y _ {i}}{\Pi_ {i} ^ {- \theta}} \times \frac {E _ {j}}{P _ {j} ^ {- \theta}}, (6)
$$

where $\Pi _ { i }$ is a producer price index capturing the (inverse) of producer market access: 

$$
\Pi_ {i} \equiv \left(\sum_ {j = 1} ^ {N} \tau_ {i j} ^ {- \theta} E _ {j} P _ {j} ^ {\theta}\right) ^ {- \frac {1}{\theta}} = A _ {i} L _ {i} Y _ {i} ^ {- \frac {\theta + 1}{\theta}}, \tag {7}
$$

and $P _ { j }$ is the consumer price index capturing the (inverse) of consumer market access: 

$$
P _ {j} = \left(\sum_ {i = 1} ^ {N} \tau_ {i j} ^ {- \theta} Y _ {i} \Pi_ {i} ^ {\theta}\right) ^ {- \frac {1}{\theta}}. \tag {8}
$$

A lower value of $P _ { j }$ indicates that consumers in location $i$ have greater access to producers in other markets, and a lower value of $\Pi _ { i }$ indicates that producers have greater access to consumers in other markets. 

# 2.1.4 Equilibrium

Finally, we calculate the equilibrium distribution of population and economic output across space. Following Allen and Arkolakis (2014), we write the welfare of residents in location $j \in \mathcal N$ , $W _ { j }$ , as: 

$$
W _ {j} = \frac {w _ {j}}{P _ {j}} u _ {j}, \tag {9}
$$

where $u _ { j }$ is an amenity value of living in location $j \in \mathcal N$ . We assume that there is free labor mobility across locations and we focus in equilibria where welfare equalizes across locations, $W _ { j } = W$ , and every location is populated.8 

Combining the definitions in (1), equation (6), the market clearing conditions (5), imposing balanced trade (i.e. $E _ { i } = Y _ { i }$ ) and welfare equalization (i.e. condition (9)), we obtain the following equilibrium conditions: 

$$
\bar {A} _ {i} ^ {- \theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {- \theta (1 + \alpha)} = \chi \sum_ {j = 1} ^ {N} \tau_ {i j} ^ {- \theta} \bar {u} _ {j} ^ {\theta} y _ {j} ^ {1 + \theta} l _ {j} ^ {\theta (\beta - 1)} \tag {10}
$$

$$
\bar {u} _ {i} ^ {- \theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (1 - \beta)} = \chi \sum_ {j = 1} ^ {N} \tau_ {j i} ^ {- \theta} \bar {A} _ {j} ^ {\theta} y _ {j} ^ {- \theta} l _ {j} ^ {\theta (\alpha + 1)} \tag {11}
$$

where $y _ { i } \equiv Y _ { i } / Y ^ { W }$ and $l _ { i } \equiv L _ { i } / L$ are the share of total income and total labor in location $i \in \mathcal N$ , respectively, and $\begin{array} { r } { \chi \equiv \left( \frac { \bar { L } ^ { ( \alpha + \beta ) } } { \bar { W } } \right) } \end{array}$  L¯(α+β)  θis an endogenous scalar capturing the (inverse) of the 

equilibrium welfare of the system.9 Conditional on $\tau _ { i j }$ this equilibrium system is identical to the one in Allen and Arkolakis (2014) and Redding (2016). In particular, given productivities $\left\{ A _ { i } \right\}$ , amenities $\{ \bar { u } _ { i } \}$ , and transportation costs $\{ \tau _ { i j } \}$ , the $2 N$ equations (10) and (11) can be solved for the $2 N$ equilibrium shares of income $\{ y _ { i } \}$ and labor $\{ l _ { i } \}$ in all locations. However, it is essential to note that the transportation costs themselves are endogenous and – through traffic congestion – will respond to the equilibrium distribution of economic activity; hence, these conditions only provide part of the story. We address the remainder of the story below in Section 3. First, however, we turn to another spatial model. 

# 2.2 An Urban Model with Optimal Routing

We next embed a routing framework in an urban model where agents commute between their place of residence and their place of work, as in Ahlfeldt, Redding, Sturm, and Wolf (2015). 

# 2.2.1 Setup

An individual $\nu \in [ 0 , 1 ]$ residing in city block $i \in \mathcal N$ who works in city block $j \in \mathcal N$ and commutes via route $r$ of length $K$ to work receives a payoff $V _ { i j , r } \left( \nu \right)$ that depends on the wage in the workplace, $w _ { j }$ ; the amenity value of residence, $u _ { i }$ ; the time spent commuting; and an idiosyncratic (Frechet distributed with shape parameter $\theta$ ) route-, origin-, and destinationspecific term, $\varepsilon _ { i j , r } \left( \boldsymbol { v } \right)$ : 

$$
V _ {i j, r} (\nu) = \left(u _ {i} w _ {j} / \prod_ {l = 1} ^ {K} t _ {r _ {l - 1}, r _ {l}}\right) \times \varepsilon_ {i j, r} (v).
$$

Individual $\nu$ chooses where to live, work, and which route to take in order to maximize $V _ { i j , r } \left( \nu \right)$ . That is, we extend the framework of Ahlfeldt, Redding, Sturm, and Wolf (2015) to introduce heterogeneity across individuals in their preference not only of where to live and work but also of what route to take when commuting between the two. Like in the economic geography framework above, this additional “noise” both substantially increases the tractability and generates an empirically plausible finite elasticity to the costs of different routes between home and work. And as above, the assumption that the three choices of 

where to live, where to work, and what route to take share the same elasticity – while straightforward to relax – greatly facilitate the tractability of the derivations and ensuing economic insight that follows. 

We assume each location $j$ produces a homogeneous and costlessly traded good with a constant returns to scale production function where labor is the only factor of production with productivity $A _ { j }$ . Taking the price of the good as the numeraire, this implies that the equilibrium real wage is the marginal product of labor $w _ { j } = A _ { j }$ . 

# 2.2.2 An analytical expression for transportation costs

The probability a worker chooses to live in $i$ , work in $j$ , and commute via route $r$ can be written as: 

$$
\pi_ {i j, r} = \frac {\prod_ {l = 1} ^ {K} t _ {r _ {l - 1} , r _ {l}} ^ {- \theta} \times u _ {i} ^ {\theta} \times w _ {j} ^ {\theta}}{\sum_ {i , j} \prod_ {l = 1} ^ {K} t _ {r _ {l - 1} , r _ {l}} ^ {- \theta} \times u _ {i} ^ {\theta} \times w _ {j} ^ {\theta}}, \tag {12}
$$

where we re-use the notation from the economic geography model for reasons that will become apparent below. This implies that the total number of workers residing in $i$ and working in $j$ , $L _ { i j }$ , can then be determined by simply summing across all routes and multiplying by the aggregate population $L$ , yielding for all $i \in \mathcal N$ and $j \in \mathcal N$ : 

$$
L _ {i j} = \sum_ {r \in \Re_ {i j}} L _ {i j, r} = \tau_ {i j} ^ {- \theta} \times u _ {i} ^ {\theta} \times w _ {j} ^ {\theta} \times \frac {\bar {L}}{\bar {W} ^ {\theta}}, \tag {13}
$$

where transportation costs $\tau _ { i j }$ are given again by (4) and $W \equiv E \left[ \operatorname* { m a x } _ { i , j , r } V _ { i j , r } \left( \nu \right) \right] =$ $\begin{array} { r } { \left( \sum _ { i j } \tau _ { i j } ^ { - \theta } \times u _ { i } ^ { \theta } \times w _ { j } ^ { \theta } \right) ^ { \frac { 1 } { \theta } } } \end{array}$ is the expected welfare of a resident in the city. 

As in the economic geography model, we assume that productivities and amenities are affected by commercial and residential population, respectively, as follows: 

$$
A _ {i} = \bar {A} _ {i} \left(L _ {i} ^ {F}\right) ^ {\alpha}, u _ {i} = \bar {u} _ {i} \left(L _ {i} ^ {R}\right) ^ {\beta}, \tag {14}
$$

where $A _ { i } > 0$ and $\bar { u } _ { i } > 0$ are again the fundamental components of productivity and amenities and $\alpha , \beta$ the respective elasticities. 

# 2.2.3 Market Access and Gravity

We can now express the gravity commuting equation (13) in market access terms. To do so, we impose the following two market clearing conditions: (1) we require that the total number of residents in $i$ , $L _ { i } ^ { R }$ , is equal to the commuting flow to all workplaces; and (2) we require that the total number of workers in $j$ , $L _ { j } ^ { F ^ { \prime } }$ , is equal to the commuting flow from all 

residences: 

$$
L _ {i} ^ {R} \equiv \sum_ {j} L _ {i j}, L _ {j} ^ {F} \equiv \sum_ {i} L _ {i j}. \tag {15}
$$

We can write the gravity commuting equation (13) as follows: 

$$
L _ {i j} = \tau_ {i j} ^ {- \theta} \times \frac {L _ {i} ^ {R}}{\Pi_ {i} ^ {- \theta}} \times \frac {L _ {j} ^ {F}}{P _ {j} ^ {- \theta}}, \tag {16}
$$

where $\Pi _ { i }$ is a resident price index capturing the (inverse of) the commuting market access residents in $i$ have to firms in all locations: 

$$
\Pi_ {i} = \left(\sum_ {j} \tau_ {i j} ^ {- \theta} L _ {j} ^ {F} P _ {j} ^ {\theta}\right) ^ {- \frac {1}{\theta}} = u _ {i} \left(L _ {i} ^ {R}\right) ^ {- \frac {1}{\theta}} \left(\frac {\bar {L}}{\bar {W} ^ {\theta}}\right) ^ {\frac {1}{2 \theta}}, \tag {17}
$$

and $P _ { i }$ is a firm price index capturing the (inverse of) the commuting market access firms in $j$ have to residents in all locations: 

$$
P _ {j} = \left(\sum_ {j} \tau_ {i j} ^ {- \theta} L _ {i} ^ {R} \Pi_ {i} ^ {\theta}\right) ^ {- \frac {1}{\theta}} = w _ {j} \left(L _ {j} ^ {F}\right) ^ {- \frac {1}{\theta}} \left(\frac {\bar {L}}{\bar {W} ^ {\theta}}\right) ^ {\frac {1}{2 \theta}}. \tag {18}
$$

Note that we re-use the notation from the economic geography framework above: both models $\Pi _ { i } ^ { - \theta }$ captures the “outward” market access and $P _ { j } ^ { - \theta }$ captures the “inward” market access with respect to the flows from $i$ to $j$ . 

# 2.2.4 Equilibrium

Substituting equations (14) into the commuting gravity equation (13) and imposing the equilibrium market clearing conditions (15) yields the following system of equations: 

$$
\left(l _ {i} ^ {R}\right) ^ {- \theta \beta + 1} = \chi \sum_ {j} \tau_ {i j} ^ {- \theta} \bar {u} _ {i} ^ {\theta} \bar {A} _ {j} ^ {\theta} \left(l _ {j} ^ {F}\right) ^ {\alpha \theta} \tag {19}
$$

$$
\left(l _ {i} ^ {F}\right) ^ {- \theta \alpha + 1} = \chi \sum_ {j} \tau_ {j i} ^ {- \theta} \bar {u} _ {j} ^ {\theta} \bar {A} _ {i} ^ {\theta} \left(l _ {j} ^ {R}\right) ^ {\beta \theta}, \tag {20}
$$

where $l _ { i } ^ { R } \equiv L _ { i } ^ { R } / L$ and $l _ { i } ^ { F } \equiv L _ { i } ^ { F } / L$ are the share of workers living and working, respectively, in location $i$ and $\begin{array} { r } { \chi \equiv \left( \frac { \bar { L } ^ { ( \alpha + \beta ) } } { \bar { W } } \right) ^ { \theta } } \end{array}$ θ is again the (inverse) of the equilibrium welfare of the system. As in the trade model above, given transportation costs $\{ \tau _ { i j } \}$ , productivities $\left\{ A _ { i } \right\}$ , and amenities $\{ u _ { i } \}$ , equations (19) and (20) can be solved to determine the equilibrium 

distribution of where people live $\left\{ l _ { i } ^ { R } \right\}$ and where they work $\left\{ l _ { i } ^ { F } \right\}$ . Once again, however, the transportation costs themselves are endogenously determined and will respond to the distribution of economic activity through traffic congestion. 

# 2.3 Taking Stock: Gravity and Optimal Routing on the Network

We now compare the aggregate outcomes of the economic geography and urban models. As is evident, the two setups are very similar, sharing (1) identical expressions for the (endogenous) bilateral trade/commuting costs (summarized in equation (4)); (2) identical gravity expressions for the bilateral flow of goods / commuters as a function of bilateral costs and market access (summarized in equations (6) and (13), respectively); and (3) mathematically equivalent equilibrium conditions (summarized in equations (10) and (11) for the economic geography model and equations (19) and (20) for the urban model). Indeed, the only distinction between the two models is the particular log linear relationship between market access variables $\Pi _ { i } ^ { - \theta }$ and $P _ { j } ^ { - \theta }$ and the equilibrium economic activity in the origin ( $Y _ { i }$ and $L _ { i } ^ { R }$ , respectively) and the destination ( $E _ { j }$ and $L _ { j } ^ { F ^ { \prime } }$ , respectively): the equilibrium conditions in both models as functions of the market access variables and economic activities are identical.10 These similarities allow us to introduce endogenous transportation costs through equilibrium traffic congestion in both frameworks using a unified set of tools we develop, which we turn to next. 

# 3 Transportation Costs, Traffic, and Congestion

In this section, we provide analytical solutions for the equilibrium transportation costs, traffic, and congestion throughout the infrastructure network. We refer the interested reader to Appendix A for detailed derivations of the results that follow in this section. 

# 3.1 Transportation Costs

Both the economic geography and urban models yield transportation costs of the form given in equation (4). By explicitly enumerating all possible routes, equation (4) can be written 

in matrix notation as follows:11 

$$
\tau_ {i j} ^ {- \theta} = \sum_ {K = 0} ^ {\infty} A _ {i j} ^ {K},
$$

where $\mathbf { A } \equiv [ t _ { i j } ^ { - \theta } ]$ , i.e. A is an $N \times N$ matrix with $( i , j )$ element $t _ { i j } ^ { - \theta }$ (not to be confused with the vector of productivities) and $\mathbf { A } ^ { K } = \lfloor A _ { i j } ^ { K } \rfloor$ , i.e. $A _ { i j } ^ { K }$ is the $( i , j )$ element of the matrix A to the matrix power $K$ .12 As in Bell (1995), as long as the spectral radius of $\mathbf { A }$ is less than one, the geometric sum can be expressed as:13 

$$
\sum_ {K = 0} ^ {\infty} \mathbf {A} ^ {K} = (\mathbf {I} - \mathbf {A}) ^ {- 1} \equiv \mathbf {B},
$$

where $\mathbf { B } = \lfloor b _ { i j } \rfloor$ is simply the Leontief inverse of the weighted adjacency matrix. As a result, the transportation cost from $i$ to $j$ can be written as a simple function of the infrastructure matrix: 

$$
\tau_ {i j} = b _ {i j} ^ {- \frac {1}{\theta}}. \tag {21}
$$

Equation (21) provides an analytical relationship between the transportation network $\mathbf { T } \equiv$ $[ t _ { k l } ]$ and the resulting transportation costs $\left\{ \tau _ { i j } \right\} _ { i , j \in \mathcal { N } ^ { 2 } }$ , accounting for the choice of the least cost route. 

Notice that in the limit case of no heterogeneity ( $\theta \to \infty$ ), the transportation costs converge to those of the least cost route, which is typically solved computationally using the Dijkstra algorithm (see e.g. Donaldson (2018)). Our formulation results in an analytical solution by extending the idiosyncratic heterogeneity already assumed in spatial models to also incorporate heterogeneity over the route chosen. In doing so, our setup bears resemblance to stochastic path-assignment methods used in transportation and computer science literature (c.f. Bell (1995); Akamatsu (1996)); here, however, the endogenous transportation costs arise from –and are determined simultaneously with– a larger general equilibrium spatial model.14 

# 3.2 Traffic Flows

We next characterize traffic along a particular link in the infrastructure matrix. This will allow us to introduce traffic congestion into the framework and relate it to observed measures of economic activity.15 

To begin, we characterize the expected number of times in which link $( k , l )$ is used in trade between $( i , j )$ , $\pi _ { i j } ^ { k l }$ , which we refer to as the link intensity. We sum across all routes from $i$ to $j$ the product of the probability a particular route is used (conditional on purchasing a product from $i$ to $j$ ) and the number of times that route passes through link $( k , l )$ , $n _ { r } ^ { k l }$ (a s some routes may use a link more than once): 

$$
\pi_ {i j} ^ {k l} \equiv \sum_ {r \in \Re_ {i j}} \left(\frac {\pi_ {i j , r}}{\sum_ {r ^ {\prime} \in \Re_ {i j}} \pi_ {i j , r ^ {\prime}}}\right) n _ {r} ^ {k l}. \tag {22}
$$

Note that for any route $r$ of length $K$ that travels through link $( k , l )$ at least once, there must exist some length $B \in [ 1 , 2 , . . . , K - 1 ]$ at which the route arrives at link $( k , l )$ . As a result, we can calculate $\pi _ { i j } ^ { k l }$ by explicitly enumerating all possible routes from $i$ to $k$ of length $B$ and all possible routes from $l$ to $j$ of length $K - B - 1$ , which can be expressed as elements of matrix powers of A. With some matrix calculus, we obtain: 

$$
\pi_ {i j} ^ {k l} = \left(\frac {\tau_ {i j}}{\tau_ {i k} t _ {k l} \tau_ {l j}}\right) ^ {\theta}. \tag {23}
$$

This expression – which resembles the one of Akamatsu (1996) derived using an exponential distribution – has a simple intuition: the more “out of the way” the transportation link $( k , l )$ is from the optimal path between $i$ and $j$ (and hence the greater the cost of traveling through link $( k , l )$ along the way from $i$ to $j$ relative to the unconstrained cost of traveling from $i$ to $j$ ) the less frequently that link is used. 

We now use the above derivation to characterize equilibrium traffic flows along each link of the network. Let $\Xi _ { k l }$ be the total traffic over link $( k , l )$ , by which we mean the total value of goods shipped (in the economic geography model) or the total number of commuters (in the urban model) over the link $( k , l )$ . To calculate $\Xi _ { k l }$ , we sum across all origins, destinations, 

and routes which travel over link $k l$ , which can be written as: 

$$
\Xi_ {k l} \equiv \sum_ {i \in \mathcal {N}} \sum_ {j \in \mathcal {N}} \sum_ {r \in \Re_ {i j}} \pi_ {i j, r} n _ {r} ^ {k l} E _ {j} = \sum_ {i \in \mathcal {N}} \sum_ {j \in \mathcal {N}} \pi_ {i j} ^ {k l} X _ {i j},
$$

$$
\Xi_ {k l} \equiv \sum_ {i \in \mathcal {N}} \sum_ {j \in \mathcal {N}} \sum_ {r \in \Re_ {i j}} \pi_ {i j, r} n _ {r} ^ {k l} \bar {L} = \sum_ {i \in \mathcal {N}} \sum_ {j \in \mathcal {N}} \pi_ {i j} ^ {k l} L _ {i j},
$$

in the economic geography and urban models, respectively. In either case, combining the market access gravity equation ((6) in the economic geography model or (16) in the urban model) with the link intensity equation (23), we obtain the following expression for equilibrium traffic flows: 

$$
\Xi_ {k l} = t _ {k l} ^ {- \theta} \times P _ {k} ^ {- \theta} \times \Pi_ {l} ^ {- \theta}. \tag {24}
$$

Equation (24) offers a gravity equation for traffic , where all determinants of the flow of traffic along link $( k , l )$ are fully summarized by the cost to travel along the link (tkl) and the economic conditions at the beginning and end of the link. It shows the tight connection between the gravity equation for traffic and trade/commuting flows, as the variables summarizing the economic conditions for the traffic gravity equation are the same market access terms $P _ { k }$ and $\Pi _ { k }$ that shape the economic conditions in the origin and destination in the economic geography and urban models. The intuition for the role that the market access terms play in the traffic gravity equation is straightforward: the greater the inward market access $\left( P _ { k } ^ { - \theta } \right)$ , the more traffic that flows into a link $k$ , and the greater the outward market access $\left( \Pi _ { l } ^ { - \theta } \right)$ , the more traffic that flows out of link $l$ .16 

Equation (24) takes the cost of traveling along a link $t _ { k l }$ as given – we now introduce traffic congestion by a parametric relationship between this cost and the traffic along the link. 

# 3.3 Traffic Congestion

To complete our modeling of traffic flows, we now suppose that the direct cost of traveling over a particular link depends in part on the total traffic flowing over that link through traffic congestion. In particular, we assume that the direct cost of traveling over a link, $t _ { k l }$ , depends in part on the amount of traffic over that link $\Xi _ { k l }$ through the following simple functional 

form: 

$$
t _ {k l} = \bar {t} _ {k l} \left(\Xi_ {k l}\right) ^ {\lambda}, \tag {25}
$$

where $\lambda > 0$ governs the strength of traffic congestion and $\mathbf { T } \equiv [ \bar { t } _ { k l } ]$ is the infrastructure network. Intuitively, if $\lambda > 0$ , the greater the fraction of total economic activity that passes through a link, the more costly traversing that link is. Like the amenity and productivity externalities in equations (1) and (14), the choice of the functional form of equation (25) succinctly allows for transportation costs to depend on an exogenous component (the infrastructure network) and an endogenous component (traffic), with a single structural parameter (λ) governing the relative strength of the two. And like with the amenity and productivity externalities, it has the unattractive feature that the transportation costs is equal to zero when the endogenous component (traffic) is equal to zero. Just as with the amenity and productivity externalities, however, this never occurs in equilibrium, as all agents’ idiosyncratic preferences over routes ensures there will be strictly positive traffic on all links. An additional attractive feature of equation (25) is that can be derived from a simple microfoundation (presented in Section 5.3) where transportation costs are log-linear functions of travel time and speed is a log-linear function of traffic congestion. 

It is important to note that the measure of traffic – and hence traffic congestion – is in the same units that we measure bilateral flows, i.e. in the economic geography model, traffic is measured in the value of goods flowing over a link, whereas in the urban model, traffic is measured in the quantity of commuters flowing over a link. There are several advantages to this approach. First, by measuring traffic in the same units that we measure bilateral flows, we generate a close connection between the (new) gravity equation (24) for traffic on a link and the (traditional) gravity equation for flows between an origin and destination (i.e. equations (6) and (16) for the economic geography and urban models, respectively). Second, as we will see below, this close connection allows us to derive analytical equilibrium conditions for the distribution of economic activity solely as a function of the model fundamentals by solving the same number of equations for the same number of unknowns despite the additional complicated feed-back loop that the presence of traffic congestion generates. Third, retaining the same units for traffic and bilateral flows – along with the assumed log-linear congestion relationship in equation equation (25) – ensures that the transportation costs between origin and destination remain ad-valorem in the presence of traffic congestion, i.e. our framework follows the large literature focusing on so-called iceberg transportation costs. 

These advantages notwithstanding, however, a reasonable objection that applies to the economic geography framework is that traffic congestion actually is increasing in the quantity rather than the value of trade: e.g., a truck carrying cheap apples generates the same traffic congestion as one carrying expensive apples. In the Online Appendix D.4, we show how 

the economic geography framework can be easily altered to to assume instead that traffic (and traffic congestion) are measured in the quantity of labor used to produce the goods and Online Appendix D.5 for the case where traffic is measured in the quantity of goods. In both cases, we show that equilibrium traffic flows also follow a gravity equation nearly identical to that of equation (24), differing only in that the market access measures are quantity-based rather than value-based. However, the need to simultaneously consider both quantity- and value-based market access measures increases the complexity of the equilibrium system, e.g. increasing the set of endogenous variables (and systems to solve) from $2 N$ to $3 N$ when traffic congestion depends on the quantity of labor used to produce the goods. 

Combining equation (25) with the gravity equation for $\Xi _ { k l }$ from equation (24) we immediately obtain: 

$$
t _ {k l} = \bar {t} _ {k l} ^ {\frac {1}{1 + \theta \lambda}} \times P _ {k} ^ {- \frac {\theta \lambda}{1 + \theta \lambda}} \times \Pi_ {l} ^ {- \frac {\theta \lambda}{1 + \theta \lambda}}, \tag {26}
$$

$$
\Xi_ {k l} = \bar {t} _ {k l} ^ {- \frac {\theta}{1 + \theta \lambda}} \times P _ {k} ^ {- \frac {\theta}{1 + \theta \lambda}} \times \Pi_ {l} ^ {- \frac {\theta}{1 + \theta \lambda}}. \tag {27}
$$

Equation (26) shows how the distribution of economic activity affects transportation costs through traffic congestion. It says that the cost of transiting a link $t _ { k l }$ is higher the better the inward market access (lower $P _ { k }$ ) at the beginning of the link and/or the better the outward market access (lower $\Pi _ { l }$ at the end of the link), as both increase traffic along the link, with $\lambda$ governing the strength of the forces. Equation (27) – which provides the basis for estimating the strength of traffic congestion below – shows traffic flows retain a gravity structure in the presence of traffic congestion. It also highlights that improvements in infrastructure quality endogenously increases the traffic demand for the infrastructure with an elasticity $\begin{array} { r } { \frac { \partial \ln \Xi _ { k l } } { \partial \ln \bar { t } _ { k l } } = - \frac { \theta } { 1 + \theta \lambda } } \end{array}$ , a fact highlighted by Duranton and Turner (2011), and a point we return to in Section 5.3. 

# 4 Traffic Congestion in the Spatial Economy

In Section 2, we characterized the equilibrium distribution of economic activity given transportation costs. In Section 3, we characterized the equilibrium transportation costs given the distribution of economic activity. In this section, we characterize both simultaneously as a function of the fundamental infrastructure network. 

# 4.1 General Equilibrium with Traffic

We begin by formally defining our equilibrium: Given a local geography $\left\{ A _ { i } , \bar { u } _ { i } \right\} _ { i \in \mathcal { N } }$ , an aggregate labor endowment $L$ , an infrastructure network $\bar { \bf T } \equiv [ \bar { t } _ { k l } ]$ , and model parameters $\{ \alpha , \beta , \theta , \lambda \}$ , we define an equilibrium to be a distribution of economic activity $\{ y _ { i } , l _ { i } \} _ { i \in \mathcal { N } }$ in the economic geography model and $\big \{ l _ { i } ^ { F } , l _ { i } ^ { R } \big \} _ { i \in \mathcal { N } }$ in the urban model and an aggregate (inverse) welfare $\chi > 0$ such that: 

1. Given equilibrium transportation costs $\left\{ \tau _ { i j } \right\} _ { i , j \in \mathcal { N } ^ { 2 } }$ , the equilibrium distribution of economic activity ensures markets clear, i.e. equations (10) and (11) hold in the economic geography model and equations (19) and (20) hold in the urban model; 

2. Given the equilibrium transportation network $\mathbf { T } \equiv [ t _ { k l } ]$ , agents optimally choose their routes through the network, i.e. equilibrium transportation costs are determined by equation (21); and 

3. Given the equilibrium distribution of economic activity, the infrastructure network $\mathbf { T } \equiv [ \bar { t } _ { k l } ]$ , and agents’ optimal route choice, the equilibrium transportation network $\mathbf { T } \equiv [ t _ { k l } ]$ is determined by the equilibrium levels of traffic congestion, i.e. equation (26) holds. 

We further define a strictly positive equilibrium to be one where the distribution of economic activity is strictly greater than zero in all locations, i.e. $y _ { i } > 0$ and $l _ { i } > 0$ for all $i \in \mathcal N$ in an economic geography model and $l _ { i } ^ { F } > 0$ and $l _ { i } ^ { R } > 0$ for all $i \in \mathcal N$ in an urban model. While the first equilibrium condition – market clearing given transportation costs – is standard to all general equilibrium spatial models, the second and third conditions are new, introducing optimal routing on the part of agents and endogenous traffic congestion, respectively. Despite the added complexity of the system, however, it turns out that the equilibrium of the system remains surprisingly tractable. 

Before deriving the new equilibrium system, two remarks are in order. First, in the absence of traffic congestion – i.e. $\lambda = 0$ – then conditional on the equilibrium transportation costs $\{ \tau _ { i j } \}$ that arise from agents optimal routing decision, the equilibrium is equivalent to the standard spatial setup upon which it is based, i.e. our framework tractably nests the standard no-congestion case. Second, with traffic congestion – i.e. $\lambda > 0$ – the equilibrium will differ from the no-congestion case, as the level of economic activity across space determines the cost of shipping in each link through traffic congestion, differences which we discuss further below. This also implies that the counterfactual predictions of our new setup with traffic congestion cannot be determined by substituting unobserved transportation costs with observed data following the “exact hat” approach of Dekle, Eaton, and Kortum (2008), as 

now $\tau _ { i j }$ is endogenous and depends on the entire network of connections through traffic and not just on bilateral flows. We nevertheless devise a new procedure, same in spirit to their exercise, but which instead replaces the need of knowledge of the entire network of connections with the use of traffic data. We discuss this in subsections 5.1 and 5.2. 

Consider first the economic geography model. Recall that equations (10) and (11) characterize the equilibrium distribution of population and income as a function of the endogenous transportation costs $\{ \tau _ { i j } \}$ , i.e. they satisfy equilibrium condition 1. To satisfy equilibrium condition 2, we substitute in equation (21) for the endogenous transportation costs and perform a matrix inversion to re-write the equilibrium conditions as a functions of the infrastructure network rather than the transportation costs and then substitute the endogenous transport costs using equations (24),(26), (27), yielding:17 

$$
y _ {i} ^ {\frac {1 + \theta + \theta \lambda}{1 + \theta \lambda}} l _ {i} ^ {- \frac {\theta (1 + \alpha + (\alpha + \beta) \theta \lambda)}{1 + \theta \lambda}} = \chi \bar {u} _ {i} ^ {\theta} \bar {A} _ {i} ^ {\theta} y _ {i} ^ {\frac {1 + \theta + \theta \lambda}{1 + \theta \lambda}} l _ {i} ^ {\frac {\theta (\beta - 1)}{1 + \theta \lambda}} + \chi^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j = 1} ^ {N} \left(\bar {L} ^ {\lambda} \bar {t} _ {i j}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \bar {A} _ {i} ^ {\theta} \bar {u} _ {i} ^ {\theta} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \bar {A} _ {j} ^ {- \frac {\theta}{1 + \theta \lambda}} y _ {j} ^ {\frac {1 + \theta}{1 + \theta \lambda}} l _ {j} ^ {- \frac {\theta (1 + \alpha)}{1 + \theta \lambda}} \tag {28}
$$

$$
y _ {i} ^ {\frac {- \theta (1 - \lambda)}{1 + \theta \lambda}} l _ {i} ^ {\frac {\theta (1 - \beta - \theta \lambda (\alpha + \beta))}{1 + \theta \lambda}} = \chi \bar {A} _ {i} ^ {\theta} \bar {u} _ {i} ^ {\theta} y _ {i} ^ {\frac {- \theta (1 - \lambda)}{1 + \theta \lambda}} l _ {i} ^ {\frac {\theta (\alpha + 1)}{1 + \theta \lambda}} + \chi^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j = 1} ^ {N} \left(\bar {L} ^ {\lambda} \bar {t} _ {j i}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \bar {A} _ {i} ^ {\frac {\theta \lambda}{1 + \theta \lambda} \theta} \bar {u} _ {i} ^ {\theta} \bar {u} _ {j} ^ {- \frac {\theta}{1 + \theta \lambda}} y _ {j} ^ {- \frac {\theta}{1 + \theta \lambda}} l _ {j} ^ {\frac {\theta (1 - \beta)}{1 + \theta \lambda}}. \tag {29}
$$

An identical process for the urban model – starting from equilibrium conditions (19) and (20), substituting in equation (21) for the endogenous transportation costs, performing a matrix inversion, and incorporating endogenous traffic congestion from equation (24),(26), (27), – yields: 

$$
\left(l _ {i} ^ {R}\right) ^ {1 - \theta \beta} \left(l _ {i} ^ {F}\right) ^ {\frac {\theta \lambda (1 - \alpha \theta)}{1 + \theta \lambda}} = \chi \bar {A} _ {i} ^ {\theta} \bar {u} _ {i} ^ {\theta} \left(l _ {i} ^ {F}\right) ^ {\frac {\theta (\alpha + \lambda)}{1 + \theta \lambda}} + \chi^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j = 1} ^ {N} \left(\bar {L} ^ {\lambda} \bar {t} _ {i j}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \bar {u} _ {i} ^ {\theta} \bar {A} _ {i} ^ {\theta} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \bar {u} _ {j} ^ {- \frac {\theta}{1 + \theta \lambda}} \left(l _ {j} ^ {R}\right) ^ {\frac {1 - \theta \beta}{1 + \theta \lambda}} \tag {30}
$$

$$
\left(l _ {i} ^ {R}\right) ^ {\frac {\theta \lambda (1 - \beta \theta)}{1 + \theta \lambda}} \left(l _ {i} ^ {F}\right) ^ {1 - \theta \alpha} = \chi \bar {u} _ {i} ^ {\theta} \bar {A} _ {i} ^ {\theta} \left(l _ {i} ^ {R}\right) ^ {\frac {\theta (\beta + \lambda)}{1 + \theta \lambda}} + \chi^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j = 1} ^ {N} \left(\bar {L} ^ {\lambda} \bar {t} _ {j i}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \bar {A} _ {i} ^ {\theta} \bar {u} _ {i} ^ {\theta} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \bar {A} _ {j} ^ {- \frac {\theta}{1 + \theta \lambda}} \left(l _ {j} ^ {F}\right) ^ {\frac {1 - \theta \alpha}{1 + \theta \lambda}}. \tag {31}
$$

Equations (28) and (29) for the economic geography model and equations (30) and (31) for the urban model determine the equilibrium distribution of economic activity $\{ y _ { i } , l _ { i } \}$ or $\{ l _ { i } ^ { F } , l _ { i } ^ { R } \}$ as a function of the model elasticities $\{ \alpha , \beta , \theta , \lambda \}$ , geography $\left\{ { { \bar { A } } _ { i } } , { { \bar { u } } _ { i } } \right\}$ , and fundamental infrastructure matrix $\bar { \bf T } \equiv [ \bar { t } _ { k l } ]$ , accounting for both the (standard) effect of transportation costs on the distribution of economic activity and the (new) effect of the distribution of economic activity on agents’ optimal routing choice, the resulting traffic congestion, and the 

equilibrium transportation costs. 

Despite the complicated feedback loop between the two effects and the necessity of solving the resulting fixed point, the dimensionality of resulting equilibrium system is not larger than the typical system treating transportation costs as exogenous, as the number of equations and number of unknowns remains the same. That allows us to make some progress in characterizing their positive properties (existence and uniqueness), which we turn to next. 

# 4.2 Existence and Uniqueness of Equilibrium with Traffic

While systems of equations with a structure as in (28) and (29) have, to our knowledge, not been studied previously, it turns out that the tools developed in Allen, Arkolakis, and Takahashi (2020) can be extended to analyze the properties of such an equilibrium.18 We first make an additional assumption on the infrastructure matrix: 

Assumption 1. The infrastructure matrix $\mathbf { T }$ is strongly connected, i.e. there exists a path with finite costs between any two locations i and $j$ , where $i \neq j$ . 

Given Assumption 1, we now provide conditions regarding existence and uniqueness in the following proposition. 

Proposition 1. For any strictly positive local geography $\left\{ \bar { A } _ { i } > 0 , \bar { u } _ { i } > 0 \right\} _ { i \in \mathcal { N } }$ , aggregate labor endowment $\bar { L } > 0$ , strongly connected infrastructure network $\bar { \bf T } \equiv [ \bar { t } _ { k l } ]$ , and model parameters $\{ \alpha \in \mathbb { R } , \beta \in \mathbb { R } , \theta > 0 , \lambda \ge 0 \}$ , then: 

1. (Existence): There exists a strictly positive equilibrium. 

2. (Uniqueness): For any $\alpha \in [ - 1 , 1 ]$ and $\beta \in [ - 1 , 1 ]$ : 

(a) In an economic geography model with a symmetric infrastructure matrix, i.e. $\bar { t } _ { k l } =$ $\bar { t } _ { l k }$ for all $l \in \mathcal N$ and $k \in \mathcal N$ , the equilibrium is unique if: 

$$
\alpha + \beta \leq 0. \tag {32}
$$

(b) In an urban model, the equilibrium is unique if: 

$$
\alpha \leq \frac {1}{2} \left(\frac {1}{\theta} - \lambda\right) a n d \beta \leq \frac {1}{2} \left(\frac {1}{\theta} - \lambda\right) \tag {33}
$$

Part 1 of Proposition 1 relies on showing that the equilibrium system defined by Equations (28) and (29) for the economic geography model and equations (30) and (31) for the urban model can be transformed into a continuous operator on a compact space so that Brouwer’s fixed point theorem applies; whereas Part 2 uses a bounding argument in the spirit of Karlin and Nirenberg (1967) and Allen, Arkolakis, and Takahashi (2020) to show that a (different) transformation of the respective systems would generate a contradiction under the reported parameter constellations. 

Despite the added complexity of endogenous traffic congestion (and the involved nature of the proofs), the sufficient conditions for uniqueness in the economic geography model provided in part (a) of the Proposition are identical to those of an economic geography model with exogenous transportation costs, provided by Allen and Arkolakis (2014): the sum of the productivity and amenity externalities must be (weakly) negative to ensure a unique equilibrium i.e. on net the forces that cause dispersion need to dominate the forces that cause concentration. In the urban model, we achieve a similar result but since we do not impose symmetry the productivity and amenity spillovers must satisfy a related condition individually, rather than combined).19 Unlike in the economic geography model, however, the strength of traffic congestion ( $\lambda$ ) does play a role in ensuring uniqueness: the stronger the traffic congestion, the lower the values of the productivity and amenity externalities must be to satisfy these sufficient conditions for uniqueness. Unlike productivity and amenity externalities where the forces occur within a location, traffic congestion forces arise on flows between locations; loosely speaking, stronger traffic congestion forces can induce greater economic concentration by reducing the flows of goods or people between locations. 

# 4.3 Traffic Congestion and Scale Dependence

In the absence of traffic congestion, equilibrium of the economic geography and urban models do not depend on the size of the aggregate labor endowment $L$ , i.e. both (standard) spatial models are scale invariant.20 In the presence of traffic congestion, however, the equilibrium 

distribution of economic activity does depend on the size of the aggregate labor endowment $L$ , i.e. the equilibrium is scale dependent. As is evident from equations (28) and (29) (in the economic geography model) and equations (30) and (31) (in the urban model), increases in $L$ are isomorphic to increases in costs of travel through the infrastructure network $t _ { i j }$ , with an elasticity equal to the strength of the traffic congestion $\lambda$ . Intuitively, the greater the aggregate labor endowment, the greater the traffic flowing through the network, and the greater the resulting traffic congestion. While the increases in the cost of travel through the infrastructure network are uniform, the impact on equilibrium transportation costs is not. To see this, we ask how a small uniform increase in the cost of travel through the entire infrastructure matrix by a factor of $c > 1$ , i.e. suppose $t _ { k l }$ increases to $c t _ { k l }$ , changes equilibrium transportation costs (holding constant traffic congestion fixed). Differentiating equation (21) around $c = 1$ yields:21 

$$
\frac {\partial \ln \tau_ {i j} (c)}{\partial \ln c} | _ {c = 1} = \sum_ {k = 1} ^ {N} \sum_ {l = 1} ^ {N} \pi_ {i j} ^ {k l},
$$

i.e. a uniform increase in the cost of travel results in a non-uniform increase in bilateral transportation costs, where origins and destinations whose link intensity across the entire network is greater face the largest increases. These disproportionate changes in transportation costs alter the equilibrium distribution of economic activity, as the following example highlights. 

# 4.4 Example

Consider a city comprising 25 locations arranged in a $5 \times 5$ grid, where, apart from their location in the grid, all locations are identical. Panel (a) of Figure 1 depicts the equilibrium distribution of economic activity in the absence of congestion forces (i.e. $\lambda = 0$ ). Locations in the center of the grid with better market access enjoy greater equilibrium economic activity (as indicated by taller “buildings”), and links in the center of the grid experience greater traffic (as indicated by their color), as they are more heavily used to travel through the network. 

In panel (b), we introduce traffic congestion, setting $\lambda = 0 . 0 5$ , but holding everything else constant. Traffic congestion disproportionately increases the cost of traversing the more heavily traveled central network segments. This disproportionately reduces the amount of traffic on those segments, causing relatively greater declines in central locations’ market access and resulting in a fall in economic activity falls in the center of the city and rises in 

the outskirts: i.e. traffic congestion forces agents out of the center of the city and into the suburbs. 

In panels (c) and (d), we increase the size of the economy from $L = 1 0 0$ to $L = 1 0 0 0$ (panel c) and $L = 1 0 0 0 0$ (panel d). As discussed above, this would have no effect on the distribution of economic activity in the absence of traffic congestion, but in the presence of traffic congestion, scale matters. Increasing the aggregate population increases traffic everywhere, but the center of city is the worse affected: the resulting gridlock induces a reallocation of economic activity away from the center and toward the edges, further amplifying the move to the suburbs. 

# 5 From Theory to Data

We now turn to applying our framework to evaluate the welfare impact of transportation infrastructure improvements. To do so, we begin by developing three helpful empirical tools: (1) we derive an equilibrium relationship between traffic flows on the one hand and trade (in the economic geography model) or commuting (in the urban model) on the other; (2) we show how to re-write the equilibrium conditions in terms of “exact hat” changes that depend only on observed traffic flows and economic activity and model parameters (e.g. the strength of traffic congestion); and (3) we present a procedure for estimating the strength of traffic congestion. 

# 5.1 Traffic, Trade, and Commuting Flows

As we discussed in Section 3.2, there is a close link between the gravity equations for trade/commuting flows (equations 6 and 16, respectively) and the gravity equation for traffic (27). It turns out that this close link admits an analytical relationship between trade/commuting flows and traffic. Combining the two gravity equations (along with the definitions of the respective market access terms), one can express equilibrium trade flows in the economic geography model as:22 

$$
X _ {i j} = c _ {i j} ^ {X} \times Y _ {i} \times E _ {j}, \tag {34}
$$

where $c _ { i j } ^ { X }$ is the $( i , j ) ^ { t h }$ element of the matrix ${ \bf C } ^ { X } \equiv \left( { \bf D } ^ { X } - \Xi \right) ^ { - 1 }$ , $\mathbf { D } ^ { X }$ is a diagonal matrix with $i ^ { t h }$ element $\begin{array} { r } { d _ { i } \equiv \frac { 1 } { 2 } \left( Y _ { i } + E _ { i } \right) + \frac { 1 } { 2 } \left( \sum _ { j = 1 } ^ { N } \left( \Xi _ { j i } + \Xi _ { i j } \right) \right) } \end{array}$ and $\Xi \equiv [ \Xi _ { i j } ]$ . 

Similarly, one can express equilibrium commuting flows in the urban model as: 

$$
L _ {i j} = c _ {i j} ^ {L} \times L _ {i} ^ {R} \times L _ {j} ^ {F}, \tag {35}
$$

where $c _ { i j } ^ { L }$ is the $( i , j ) ^ { t h }$ element of the matrix ${ \bf C } ^ { L } \equiv \left( { \bf D } ^ { L } - \Xi \right) ^ { - 1 }$ , $\mathbf { D } ^ { L }$ is a diagonal matrix with $i ^ { t h }$ element $\begin{array} { r } { d _ { i } \equiv \frac { 1 } { 2 } \left( L _ { i } ^ { R } + L _ { i } ^ { F } \right) + \frac { 1 } { 2 } \left( \sum _ { j = 1 } ^ { N } \left( \Xi _ { j i } + \Xi _ { i j } \right) \right) } \end{array}$ and $\Xi \equiv [ \Xi _ { i j } ]$ . 

Equations (34) and (35) show that in both the economic geography and urban models, the equilibrium flows from origin to destination can be written only in terms of the economic activity in the origin ( $Y _ { i }$ and $L _ { i } ^ { R }$ , respectively), economic activity in the destination ( $E _ { i }$ and $L _ { i } ^ { F }$ , respectively), and the matrix of traffic flows through the network, $\Xi$ .23 In particular, equations (34) and (35), show that trade and commuting flows can be expressed as (an appropriately scaled) Leontief inverse of the traffic flows. Note that the expression depends only on available data and hence can be accomplished without knowledge of the underlying model elasticities. This result had two advantages, depending on the empirical availability of trade/commuting flows. In settings where both traffic flows and commuting / trade flows are observed (such as our empirical contexts discussed below), it provides an out-of-sample test of the model predictions about traffic flows. In addition, if trade/commuting data are not available, but traffic data is (e.g. much of the developing world), it still enables one to evaluate the welfare impacts of infrastructure improvements, a point we turn to next. 

# 5.2 Counterfactuals

To evaluate the welfare impact of transportation infrastructure improvements in the presence of traffic congestion we next analyze how to conduct counterfactuals. To do so we follow the “exact hat algebra” approach pioneered by Dekle, Eaton, and Kortum (2008), where we denote with hats the change in variables, $\begin{array} { r } { \hat { \gamma } _ { i } \equiv \frac { \gamma _ { i } ^ { \prime } } { \gamma _ { i } } } \end{array}$ , where we denote with prime the counterfactual outcome. We summarize the result in the following proposition. 

Proposition 2. Suppose an observed economy has infrastructure network $\mathbf { T } \equiv [ \bar { t } _ { k l } ]$ and is in equilibrium. Consider any change in the underlying infrastructure network denoted by $\hat { \bar { t } } _ { k l }$ . Given observed traffic flows, $\left[ \Xi _ { i j } \right]$ , economic activity in the geography $( Y _ { i } , E _ { j } )$ or urban model $\left( L _ { i } ^ { R } , L _ { j } ^ { F } \right)$ and parameters $\{ \alpha , \beta , \theta , \lambda \}$ , the equilibrium change in economic outcomes 

$\left( \hat { y } _ { i } , \hat { l } _ { i } , \hat { \chi } \right)$ is the solution the following system of equations: 

$$
\begin{array}{l} \hat {y} _ {i} ^ {\frac {1 + \theta \lambda + \theta}{1 + \theta \lambda}} \hat {l} _ {i} ^ {- \frac {\theta (1 + \alpha + \theta \lambda (\beta + \alpha))}{1 + \theta \lambda}} = \hat {\chi} \left(\frac {E _ {i}}{E _ {i} + \sum_ {k} \Xi_ {i k}}\right) \hat {y} _ {i} ^ {\frac {1 + \theta \lambda + \theta}{1 + \theta \lambda}} \hat {l} _ {i} ^ {\frac {\theta (\beta - 1)}{1 + \theta \lambda}} + \hat {\chi} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j} \left(\frac {\Xi_ {i j}}{E _ {i} + \sum_ {k} \Xi_ {i k}}\right) \hat {\bar {t}} _ {i j} ^ {- \frac {\theta}{1 + \theta \lambda}} \hat {y} _ {j} ^ {\frac {1 + \theta}{1 + \theta \lambda}} \hat {l} _ {j} ^ {- \frac {\theta (1 + \alpha)}{1 + \theta \lambda}} \\ \hat {y} _ {i} ^ {\frac {- \theta (1 - \lambda)}{1 + \theta \lambda}} \hat {l} _ {i} ^ {\frac {\theta (1 - \beta - \theta \lambda (\alpha + \beta))}{1 + \theta \lambda}} = \hat {\chi} \left(\frac {Y _ {i}}{Y _ {i} + \sum_ {k} \Xi_ {k i}}\right) \hat {y} _ {i} ^ {\frac {- \theta (1 - \lambda)}{1 + \theta \lambda}} \hat {l} _ {i} ^ {\frac {\theta (\alpha + 1)}{1 + \theta \lambda}} + \hat {\chi} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j = 1} ^ {N} \left(\frac {\Xi_ {j i}}{Y _ {i} + \sum_ {k} \Xi_ {k i}}\right) \hat {\tilde {t}} _ {j i} ^ {- \frac {\theta}{1 + \theta \lambda}} \hat {y} _ {j} ^ {- \frac {\theta}{1 + \theta \lambda}} \hat {l} _ {j} ^ {\frac {\theta (1 - \beta)}{1 + \theta \lambda}} \\ \end{array}
$$

for the economic geography model and as: 

$$
\left(\hat {l} _ {i} ^ {R}\right) ^ {1 - \theta \beta} \left(\hat {l} _ {i} ^ {F}\right) ^ {\frac {\theta \lambda (1 - \alpha \theta)}{1 + \theta \lambda}} = \hat {\chi} \left(\frac {L _ {i} ^ {F}}{L _ {i} ^ {F} + \sum_ {k} \Xi_ {i k}}\right) \left(\hat {l} _ {i} ^ {F}\right) ^ {\frac {\theta (\alpha + \lambda)}{1 + \theta \lambda}} + \hat {\chi} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j = 1} ^ {N} \left(\frac {\Xi_ {i j}}{L _ {i} ^ {F} + \sum_ {k} \Xi_ {i k}}\right) \hat {\tilde {t}} _ {i j} ^ {- \frac {\theta}{1 + \theta \lambda}} \left(\hat {l} _ {j} ^ {R}\right) ^ {\frac {1 - \theta \beta}{1 + \theta \lambda}} \tag {38}
$$

$$
\left(\hat {l} _ {i} ^ {R}\right) ^ {\frac {\theta \lambda (1 - \beta \theta)}{1 + \theta \lambda}} \left(\hat {l} _ {i} ^ {F}\right) ^ {1 - \theta \alpha} = \hat {\chi} \left(\frac {L _ {i} ^ {R}}{L _ {i} ^ {R} + \sum_ {k} \Xi_ {k i}}\right) \left(\hat {l} _ {i} ^ {R}\right) ^ {\frac {\theta (\beta + \lambda)}{1 + \theta \lambda}} + \hat {\chi} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j = 1} ^ {N} \left(\frac {\Xi_ {j i}}{L _ {i} ^ {R} + \sum_ {k} \Xi_ {k i}}\right) \hat {\tilde {t}} _ {j i} ^ {- \frac {\theta}{1 + \theta \lambda}} \left(\hat {l} _ {j} ^ {F}\right) ^ {\frac {1 - \theta \alpha}{1 + \theta \lambda}}
$$

for the urban model. Moreover, existence and uniqueness of the counterfactuals are given by the same conditions as in Proposition (1). 

Proof. See Online Appendix C.2. 

Proposition 2 says that given observed traffic flows and the observed distribution of economic activity – and knowledge of the model parameters $\{ \theta , \alpha , \beta , \lambda \} - \mathrm { i t }$ is possible to evaluate the impact of any transportation infrastructure improvements $\left\{ \hat { \bar { t } } _ { i j } \right\}$ on the equilibrium distribution of economic activity and aggregate welfare.24 Note that equations (36)-(39) all say that some log linear combination of endogenous changes in location $i$ depend on a weighted average of a (different) log linear combination of endogenous changes in location $i$ and a (third) log linear combination of endogenous changes in all location $j$ , where the weights are determined by the relative size of observed local economic activity and traffic flows. Loosely speaking, this locations with large amounts of traffic flows to $i$ will play a greater role in determining the counterfactual outcomes in $i$ , all the more so if these traffic flows are large relative to the economic activity in $i$ . It is worth emphasizing that conducting counterfactuals using this result requires easily observed traffic flows along links in the network, instead of potentially harder to observe bilateral trade or commuting flows between origins and destinations upon which traditional implementations of the Dekle, Eaton, and Kortum (2008) “exact hat” algorithm rely (see e.g. Redding (2016); Caliendo, Parro, Rossi-Hansberg, and Sarte (2018); Adao, Arkolakis, and Esposito (2021)). 

The second part of Proposition 2 says that existence and the sufficient conditions for uniqueness for the counterfactuals are the same as for the the system in levels. This result arises from the fact that the systems of equations that determine the counterfactual outcomes in changes are mathematically equivalent to their level variants above, where the local geography and infrastructure matrix are simply replaced with shares that depend only on observed traffic flows and the observed distribution of economic activity. 

While the first three parameters $\{ \theta , \alpha , \beta \}$ are familiar ingredients in spatial models (and we will be calibrating their values to those of the literature below), the strength of traffic congestion $\lambda$ is new to our framework. We turn now to its estimation. 

# 5.3 Estimating the Strength of Traffic Congestion

To derive a straightforward estimating equation for the strength of the endogenous traffic congestion, we make two additional assumptions. First, we follow an extensive literature on trade cost estimation, and assume that transportation costs $t _ { k l }$ are a log-linear function of travel time.25 As a result, we can write $t _ { k l }$ as a function of the distance of the link and the speed of travel on the link: 

$$
t _ {k l} = \left(d i s t a n c e _ {k l} \times s p e e d _ {k l} ^ {- 1}\right) ^ {\delta_ {0}}, \tag {40}
$$

where $\delta _ { 0 }$ is the time elasticity of the transportation cost. In our preferred results below, we set $\delta _ { 0 } = 1 / \theta$ to imply a “distance elasticity” of negative one, which is consistent with a large gravity literature, see e.g. Disdier and Head (2008) and Chaney (2018).26 

Our second assumption is that time per unit distance (inverse speed) is a log-linear function of traffic congestion (measured as total vehicle miles traveled per lane-miles, or equivalently, traffic per average lanes) as follows: 

$$
s p e e d _ {k l} ^ {- 1} = m _ {0} \times \left(\frac {\Xi_ {k l}}{l a n e s _ {k l}}\right) ^ {\delta_ {1}} \times \varepsilon_ {k l} \tag {41}
$$

where $\delta _ { 1 }$ is the congestion elasticity of inverse speed, $m _ { 0 }$ is the average rate of flow without 

congestion, $\boldsymbol { l a n e s _ { k l } }$ are the average number of lanes on a link, and $\varepsilon _ { k l }$ is a segment specific idiosyncratic free rate of flow. The log-linear specification was first posited by Vickrey (1967), and while simple, has a number of advantages in our setting.27 First, combined with equations (40), and (41) immediately implies: 

$$
t _ {k l} = \bar {t} _ {k l} \times \left(\Xi_ {k l}\right) ^ {\lambda},
$$

where $\bar { t } _ { k l } \equiv l a n e s _ { k l } ^ { - \delta _ { 0 } \delta _ { 1 } } \times \left( d i s t a n c e _ { k l } \times m _ { 0 } \times \varepsilon _ { k l } \right) ^ { \delta _ { 0 } }$ −δ0δ1kl × (distancekl × m0 × εkl)δ0 and λ ≡ δ0δ1. That is, this simple setup $\lambda \equiv \delta _ { 0 } \delta _ { 1 }$ offers a micro-foundation for the traffic congestion formulation (25) posited in Section 4. Second, treating distance and the free rate of flow as segment specific time-invariant characteristics, equation (41) provides a simple relationship between infrastructure improvements and the change in the infrastructure matrix: 

$$
\hat {\bar {t}} _ {k l} = \hat {l a n e s} _ {k l} ^ {- \lambda}. \tag {42}
$$

As additional lane-miles are added to a segment, congestion on the segment falls, reducing the exogenous component of transportation costs with an elasticity of $\lambda$ . This is intuitive: the greater the strength of traffic congestion, the larger the impact of adding additional lanes. However, it is important to (re-)emphasize that improvements in the infrastructure matrix will also result in an endogenous increase in traffic demand. Indeed, combining equation (42) with (27), we see that the elasticity of traffic to lanes is ∂ ln laneskl $\begin{array} { r } { \frac { \partial \ln \Xi _ { k l } } { \partial \ln l a n e s _ { k l } } = \frac { \lambda \theta } { 1 + \lambda \theta } } \end{array}$ ∂ ln Ξkl , i.e. the limiting case as traffic congestion becomes infinitely large is that traffic increases proportionately with the adding of additional lanes, as in “the fundamental law of road congestion” identified by Duranton and Turner (2011).28 

The final advantage of this setup is that it delivers a straightforward estimating equation and, combined with the traffic gravity equation (27), an appropriate identification strategy. Taking logs of equation (41) yields: 

$$
\ln s p e e d _ {k l} ^ {- 1} = \ln m _ {0} + \delta_ {1} \ln \left(\frac {\Xi_ {k l}}{l a n e s _ {k l}}\right) + \ln \varepsilon_ {k l}, \tag {43}
$$

i.e. a regression of inverse speed on traffic congestion can in principal identify the congestion 

elasticity of inverse speed $\delta _ { 1 }$ . An ordinary least squares regression is inappropriate in this case, as the residual is the free rate of flow on the segment $k l$ , which enters into $t _ { k l }$ and so, from the traffic gravity equation (27) is negatively correlated with traffic $\Xi _ { k l }$ , biasing the estimate of $\dot { \delta _ { 1 } }$ downwards. Instead, we propose to use an instrumental variables strategy, instrumenting for traffic $\Xi _ { k l }$ with observables that affect traffic demand for a segment but are uncorrelated with the free rate of flow on the segment.29 From the traffic gravity equation (27), conditional on $k$ and $l$ fixed effects, any component of $t _ { k l }$ that does not affect the free rate of flow is a suitable instrument. Intuitively, we can use observables that shift the traffic gravity (demand) equation to identify the slope of the traffic congestion (supply) equation. We describe such instruments in the next section, where we apply our procedure to determine the welfare impact of transportation infrastructure improvements in two different settings. 

# 6 The welfare impact of transportation infrastructure improvements

We first apply the economic geography variant of our framework to evaluate the welfare impact (and, given cost estimates, the return on investment) of small improvements to every single segment of the U.S. Interstate Highway network. We then apply the urban variant of our framework to do the same for each segment of the road network in Seattle, WA.30 

# 6.1 Traffic across the Country: The U.S. Highway Network

The U.S. National Highway System is the largest highway system in the world. The main backbone of the National Highway System – the Interstate Highway System – is one of the world’s largest infrastructure projects in history (Kaszynski, 2000), taking more than thirty five years to construct at an estimated cost $650 billion (in 2014 dollars), and total annual maintenance costs are approximately $70 billion (CBO, 1982; FHA, 2008; NSTIFC, 2009; ASCE, 2017). However, little is known about the relative importance of different segments of the highway system in terms of how each affects the welfare of the U.S. population. Such knowledge is crucial for appropriately targeting future infrastructure investments. 

Our strategy to estimate the welfare impact of improvements to the U.S. Highway System 

is straightforward: for each segments of the network, we will use equations (36) and (37) from Proposition 2 for the economic geography variant of our approach to estimate the aggregate welfare impact ( $\hat { \bar { W } } = \hat { \chi } ^ { - \frac { 1 } { \theta } }$ ) of a small (1%) improvement to the infrastructure network. We then use equation (42) to calculate how many lane-miles must be added in order to achieve a $1 \%$ improvement in order to estimate such an infrastructure cost. Given costs and benefits, we can then identify the highway segments with the greatest return on investment. This procedure requires just two ingredients: (1) data on traffic $\{ \Xi _ { k l } \}$ and income $\{ Y _ { i } = E _ { i } \}$ ; and (2) knowledge of the four model parameters $\{ \theta , \alpha , \beta , \lambda \}$ . We discuss the source of these ingredients in turn. 

# 6.1.1 Data

We briefly summarize the data used here; see Online Appendix F.1 for more details. The primary source of data we use to construct the infrastructure network is the 2012 Highway Performance Monitoring System (HPMS) dataset by the Federal Highway Administration. This dataset comprises the length, location, number of lanes, and average annual daily traffic (AADT) over 330,021 segments of the U.S. highway system.31 

To create the infrastructure network, we begin by placing nodes at each endpoint and intersection between two different Interstate highways and collapsing all nodes within the same core-based statistical area (CBSA) to a single CBSA point. This resulting 228 locations and 704 links between adjacent nodes, where for each link we construct a length-weighted average of AADT and number of lanes. Panel (a) of Figure 3 depicts the actual highway network and the resulting infrastructure network. 

To this network, we append four additional data sources. First, to estimate the strength of congestion, we recover the time of travel $( t i m e _ { k l } )$ across each link from the HERE API using the georoute Stata command by Weber and P´eclat (2017). Second, we calculate the population and income at each node by summing the population and averaging the median income of all cities from Edwards (2017) (which is itself based on the U.S. Census and American Community Survey) within 25 miles of the node. Third, we estimate the cost of improving each link based on the topography of its constituent segments. To do so, we classify each segment of the Interstate Highway System into one of seven categories from the Federal Highway Administration’s Highway Economic Requirements System (HERS) Federal Highway Administration (2015), each of which is associated with an estimated cost of adding 

one lane-mile.32 To determine the average cost of adding one lane-mile to a link, we construct a distance-weighted average of the cost of improving each of its constituent segments. Fourth, we rely on the 2012 Commodity Flow Survey (CFS) to construct measures of the value of bilateral trade flows between each CBSA; for CFS areas comprising more than one CBSA, we allocate observed CFS area flows to CBSAs proportionally to their share of the CFS area’s total income. 

# 6.1.2 Predicted versus observed trade flows

As a first check of the validity of the framework developed above, we compare the observed value of bilateral trade flows between CBSAs from the CFS to the backed out bilateral trade flows using equation (34) and the the observed traffic flows. To do so, we assume that each element of the matrix of traffic flows $\Xi \equiv [ \Xi _ { k l } ]$ is equal to the observed AADT along the highway segment, which is equivalent to assuming that each car is carrying a value of trade equal to the average value of a single individual’s labor. This of course abstracts from many nuances of traffic flows, including shipments via truck (where the trade value exceeds this average) as well as traffic for non-trade purposes such as commuting and shopping (where the trade value falls below this average). Given these abstractions, it is all the more remarkable how well traffic across the interstates is able to predict actual trade between CBSAs. Panel (a) of Figure 2 shows the scatter plot between observed and predicted (log) trade flows, conditional on origin and destination fixed effects (so the only variation arises from the bilateral flows and not e.g. income in the origin or destination). As is evident, there is a strong positive correlation of 0.60, indicating the traffic matrix – through the lens of the theory and despite obvious measurement issues – does a good job of predicting trade flows.33 

Finally, panel (a) of Figure 4 provides an example of intensity of usage of different links for a specific origin and destination pair, Los Angeles, California to New York, New York. As expected, the links that are on very direct routes, such as for example segments of the I-95 interstate near New York, are very intensively used to serve that pair, whereas more indirect links such as highway segments in California north of Los Angeles, have negligible usage. 

# 6.1.3 Estimation

We now discuss our choice of the four model parameters $\{ \theta , \alpha , \beta , \lambda \}$ . As the first three model parameters – the trade elasticity $\theta$ , productivity externality $\alpha$ , and amenity externality $\beta$ – are standard in the economic geography literature, we choose central values from the literature. We set $\theta = 8$ to match previous estimates of the trade elasticity.34 We also choose $\alpha = 0 . 1$ , and $\beta = - 0 . 3$ , which corresponds to the estimated scale economies found in the literature, as e.g. summarized in Rosenthal and Strange (2004) and Combes and Gobillon (2015) and the share of consumption allocated to housing, see e.g. Allen and Arkolakis (2014).35 From Proposition 1, this choice of parameter values guarantees the existence of a unique equilibrium. 

To estimate the strength of traffic congestion, we follow the estimation procedure described in Section 5.3, regressing observed inverse speed on (appropriately instrumented) traffic congestion as in equation (43). As implied by the traffic gravity equation (27), recall that an appropriate instrument would be something that – conditional on start-location and end-location fixed effects – affects the cost of travel $t _ { k l }$ but is uncorrelated with the free-flow speed of travel on the link. In the context of the U.S. highway system, we propose that the distance along the link is such an appropriate instance. Distance clearly affects the cost of travel (and so is relevant), and given the relative homogeneity of U.S. highways in terms of speed limits, lanes, limited access, etc., we have no reason to believe that longer or shorter links have different free flow rates of speed (so it is plausibly excludable). 

Panel (a) of Table 1 presents the results. Columns (1) and (2) show using OLS that there is a positive, but small, correlation between inverse speed of travel and congestion. Column (3) presents the first stage regression of traffic on distance: as expected, conditional on start-location and end-location fixed effects, distance is strongly negatively correlated with traffic. Column (4) presents the IV regression: Consistent with OLS exhibiting downward bias due to traffic demand being lower on slower links, the IV is substantially larger, finding a coefficient $\delta _ { 1 } = . 7 3 9$ (with standard error of .181). Recall from above that we set $\delta _ { 0 } = 1 / \theta$ to match the unit distance elasticity, so this implies $\lambda = \dot { \delta } _ { 1 } \delta _ { 0 } = 0 . 0 9 2$ , i.e. a 10% increase in 

traffic flows is associated with a $7 . 4 \%$ increase in travel time, resulting in a $0 . 9 \%$ increase in the transportation cost.36 

# 6.1.4 Results

Given the observed traffic data and estimated parameters, we calculate the aggregate welfare elasticity to a 1% reduction in iceberg transportation costs on every link (in both directions of travel) of the U.S. Highway System using equations (36) and (37) of Proposition 2, i.e. $\begin{array} { r } { \frac { 1 } { 2 } \left( \frac { \partial \ln \bar { W } } { \partial \ln \bar { t } _ { k l } } + \frac { \partial \ln \bar { W } } { \partial \ln \bar { t } _ { l k } } \right) } \end{array}$ Panel (a) of Figure 5 presents our results. While all highway segments have positive welfare elasticities, the elasticities are largest on short segments connecting CB-SAs in densely populated areas, e.g. along I-95 between Boston and Philadelphia and on I-5 between Los Angeles and San Diego. Welfare elasticities are also large along longer highway segments that do not directly connect large urban areas but that are major thoroughfares for trade, e.g. in the interstates passing through Indiana (“the crossroads of America”). Conversely, highway segments that neither connect major urban areas nor are used intensively for trade – such as I-90 through Montana – have the lowest positive impact on aggregate welfare. 

How much does incorporating endogenous traffic congestion affect our welfare elasticity estimates? Panel (a) of Figure 6 compares the welfare elasticity for each segment with and without congestion. From the scatter plot on the right, it is clear that in the absence of traffic congestion, the welfare gains from reducing transportation costs are greater. What is surprising, however, is that there is substantial variation in welfare gains with and without congestion across segments. From the map on the left, we see that ignoring traffic congestion overstates the welfare gains from infrastructure improvements the most along highly trafficked segments of the highway system such as I-5 between Los Angeles and San Diego, California as well as along highway segments around important hubs for intrastate shipping such as those surrounding Atlanta, Georgia, highlighting the fact that traffic congestion plays an important role in determining which segments would achieve the greatest welfare gains. 

The benefit of improving a link, of course, is only half of the story. To calculate a return on investment, we pursue a cost-benefit approach. On the benefit side, we translate the welfare elasticity into a dollar amount use a compensating variation approach, asking how much the annual U.S. real GDP (of $\$ 19$ trillion) would have to increase (in millions of chained 2012 US dollars) to bring about the same welfare increase we estimate. On the cost side, we first use equation (42) to calculate how many additional lane-miles would need to be 

added to the route to achieve a $1 \%$ reduction in transportation costs. We then multiply this number of lane-miles by the cost per lane-mile to get a total construction cost. We assume a 20 year depreciation schedule (as in Appendix C of Office of the State Auditor (2002)), a 5% annual maintenance cost, and a 3% borrowing cost, which together imply 10% of the construction cost is incurred each year.37 

Panel (a) of Figure 7 reports the annual return on investment (RoI) for each segment of the U.S. highway system. On average, infrastructure improvement return are well-worth the investment, with a mean RoI of just over 108%. However, there is also huge variance in returns, with some segments offering negative RoI (such as I-90 through Montana) and others offering much higher than average. Panel (a) of Table 2 presents the ten links with the highest RoI (each of which exceed 400%). All ten are for links outside the largest cities, where reducing transportation costs is less costly. This does not mean that returns are entirely driven by costs: the links with the highest returns are those on the periphery of densely populated areas with high welfare elasticities, reflecting the importance of trade between these regions. 

# 6.2 Traffic in the City: The Seattle Road Network

We now analyze the urban variant of our framework to examine the welfare impacts of transportation infrastructure improvement in Seattle, WA. Seattle provides an ideal testcase for our framework for several reasons, notably: (1) it has some of the worst traffic in the U.S.; (2) with limited (non-bus) public transit options, its road network plays a critical role in commuting; and (3) its road network is particularly interesting, with multiple natural choke points created by the waterways which intersect the city.38 

Our strategy for estimating the welfare impacts of improvements to the Seattle road network proceeds analogously to the U.S. highway system above: for each link in the road network, we estimate the change in the aggregate welfare $\left( \hat { W } = \hat { \chi } ^ { - \frac { 1 } { \theta } } \right)$ from a small (1%) improvement using equations (38) and (39). Doing so requires just two ingredients: (1) data on traffic $( \Xi _ { k l } )$ , residential population $\left( L _ { i } ^ { R } \right)$ , and workplace population $\left( L _ { i } ^ { F } \right)$ and (2) values 

for the model parameters $\{ \theta , \alpha , \beta , \lambda \}$ . We discuss the source of both ingredients in turn. 

# 6.2.1 Data

We briefly summarize the data used here; see Online Appendix F.2 for more details. Data on the location, functional system (i.e., interstate, arterial road, local road, etc.), ownership, AADT, lane width, and possibility for lane expansion of the 9,188 road segments within the municipal boundaries of Seattle were taken from the 2016 HPMS release for the state of Washington.39 To construct our adjacency matrix of Seattle, we divide Seattle into $\tilde { \mathbf { \Gamma } } ^ { 1 }$ sq. mi. grids, place the center point of each of these grids as a node into ArcGIS Network Analyst, and find the least-cost path between each of these nodes.40 This gives us a total of 217 nodes, with 1,384 links between adjacent nodes, 1,338 for which we observe traffic.41 Panel (b) of Figure 3 depicts the actual Seattle road network and the resulting infrastructure network. 

We append to this network five additional sources of data. First, we calculate the time of travel between each link from the HERE API using the georoute Stata command by Weber and P´eclat (2017). Second, we observe the labor force and residential population density at the census block group level from the 2017 Longitudinal Employer-Household Dynamics Origin-Destination Employment Statistics (LODES), which we aggregate to our constructed grids (allocating population from block groups intersected by our grids proportional to the area of the block group within each grid). Third, the LODES data also provide bilateral commuting flows between census block groups, which we aggregate to bilateral grid cell pairs using a similar procedure. Fourth, we estimate the cost of adding an additional lane-mile to each link in the network. To do so, we classify each Seattle’s road sections into the major urbanized road type based on the population of the Seattle urban area (as defined by the Census Bureau’s 2012 Urban Area data) and additionally indicate if the section is “restricted” if the HPMS indicate that additional lanes cannot be added. Then, based on a road section’s functional system classification, its major urbanized classification, and whether it is a high 

39Traffic data on a road segment is reported without regard to the direction of travel. As such, we evaluate simultaneous improvements to each link in the Seattle road network in both directions of travel. This has the added advantage of reconciling our urban framework – where traffic is modeled as flowing from an agents’ residence to her workplace – to the (presumed) empirical reality that the agent returns home after work. 

40This approach is necessary because, at this level, typical units of observation like census blocks and block grounds are endogenous to the road structure of Seattle; this leaves us with concerns that census blocks which are larger are in a less dense area of Seattle with less traffic. 

$^ { 4 1 }$ Unlike the interstates, where we observe all segments of the highway system, our analysis does not cover every road in Seattle, just those along the least-cost path between adjacent nodes. We do, however, observe the entirety of the Seattle road network in our dataset. We assume the route along the least-cost path between nodes reasonably captures the costs of moving across similar paths, on different roads, between the same nodes. 

cost road to improve or not, we code each road section with the cost of adding a lane-mile to it, as estimated by the FHA’s HERS from Federal Highway Administration (2015).42 Fifth, for the construction of our instrument, we calculate the number of intersections and turns along each link of the network using the ArcGIS network analyst. 

# 6.2.2 Predicted versus observed commuting flows

As a first pass of the validity of the urban variant of our framework to the data, we compare the observed bilateral commuting flows from LODES to backed out from equation (35) using the observed traffic flows using equation (35). To do so, we assume that each element of the matrix of traffic flows $\Xi \equiv [ \Xi _ { k l } ]$ is equal to the observed AADT along that road segment. This assumes every vehicle carries one commuter. As with the interstates, this introduces obvious measurement error: some vehicles contain many commuters (e.g. buses), whereas other vehicles contain none (e.g. when driving to go shopping). And like with interstates, it is remarkable how well observed traffic flows are able to predict commuting flows, as panel (b) of Figure 2 illustrates. Even conditional on origin and destination fixed effects, there is a positive correlation between predicted and observed commuting flows of 0.43, indicating that the urban model with traffic congestion is able to successfully predict observed commuting flows.43 

Finally, panel (b) of Figure 4 shows the intensity of usage of different links for an example commute from Safeco Field to the University of Washington, both on opposite sides of the the city center. As with the interstate highway system, links along the most direct routes are most intensively used. The figure also highlights the fact that different links through the city center are quite substitutable with each other, with no one link being used more than about half the time, whereas the natural choke-points – e.g. the bridges over Lake Union – are traversed essentially on all routes. In contrast, routes not along the direct route are used negligibly. 

# 6.2.3 Estimation

We now discuss our choice of model parameters $\{ \theta , \alpha , \beta , \lambda \}$ . As the first three model parameters are standard in the quantitative urban literature, for our preferred estimates presented 

here we set them equal to the values of estimated in the seminal work of Ahlfeldt, Redding, Sturm, and Wolf (2015), with $\theta = 6 . 8 3$ , $\alpha = - 0 . 1 2$ , and $\beta = - 0 . 1$ .44 From Proposition 1, this choice of parameter values guarantees the existence of a unique equilibrium. 

To estimate the strength of traffic congestion, we again proceed as discussed in Section (39), regressing the observed inverse speed of travel over a link on the traffic congestion, appropriately instrumented by a demand shifter uncorrelated with the free-flow rate of speed over the link. Unfortunately, the instrument used for the U.S. highway system – distance – is inappropriate in a city setting. There exists enormous variation in the types of roads and speed of travel within Seattle (e.g. surface streets with stop signs, larger streets with major intersections, highways, etc.), so it is likely that the distance of a segment is correlated with its free-flow rate of speed (e.g. a link which travels along a highway might be longer but faster). As an alternative, we propose that the complexity of a route is a suitable instrument: conditional on the free-flow rate of speed, drivers would prefer to take routes that are less complex. To measure complexity, we use the number of turns along the route as our instrument, conditioning on the number of intersections.45 Intuitively, intersections reduce the free-flow rate of speed of travel regardless if one turns or not, while turns themselves present an additional inconvenience to drivers. 

Panel (b) of Table 1 presents the results. Column (1) shows that there is actually a small negative correlation between inverse speed and traffic, consistent with substantial downward bias due to the heterogeneity in free-flow speed across links (e.g. faster links on highways also have higher traffic). Column (2) presents the first stage results; as expected, the greater the number of turns along a route (conditional on the number of intersections), the lower the traffic along that link. Column (3) presents the IV results, where we estimate $\delta _ { 1 } =$ 0.118 (with a standard error of 0.048). One potential concern with the instrument is that controlling for the number of intersections alone may not be sufficient to allay the concern that more complex routes are more likely to travel over smaller (and slower) roads. In Columns (4) and (5) present the first and second stage results where we nonparametrically control for the share of the route that travels over arterial and local roads.46 Such a procedure compares links with similar road compositions, mitigating the concern that route complexity 

is correlated with unobserved speed of travel. Adding these controls increases our estimate of $\delta _ { 1 } = 0 . 4 8 8$ (with standard error of 0.278). Combined with the maintained assumption that $\delta _ { 0 } = 1 / \theta$ (to generate a unit distance elasticity), this implies a traffic congestion parameter of $\lambda = \delta _ { 1 } \delta _ { 0 } = 0 . 0 7 1$ , i.e. a 10% increase in traffic flows is associated with a 4.9% increase in travel time, resulting in a 0.7% increase in the transportation cost. It is interesting to note that while the elasticity of travel time to congestion is smaller in Seattle than U.S. highways – perhaps due to the lower free-flow rates of speed within a city – the impact of traffic congestion on transportation costs in both settings is quite similar. 

# 6.2.4 Results

For each link in the road network, we simulate a 1% reduction in transportation costs in each direction and calculate the change in aggregate welfare elasticity $\begin{array} { r } { \frac { 1 } { 2 } \left( \frac { \partial \ln \bar { W } } { \partial \ln \bar { t } _ { k l } } + \frac { \partial \ln \bar { W } } { \partial \ln \bar { t } _ { l k } } \right) } \end{array}$ ∂ ln t¯lk Panel (b) of Figure 5 presents our findings. While a reduction in transportation costs on all links are welfare improving, the largest welfare elasticities are greatest in the center of the city (downtown). Welfare elasticities are also higher for the various choke-points in the road network (oftentimes corresponding to bridges over water). 

Panel (b) of Figure 5 compares these estimated welfare elasticities to those estimated without traffic congestion. As with the U.S. highway system, ignoring congestion would not just result in overestimates of the welfare elasticities, it would also substantially change which links one would identify as having the largest welfare effects. The left figure shows the variation across links in the degree to which one would overestimate welfare gains by ignoring congestion. As is evident, heavily trafficked links near the city center and along interstate I-5 whose gains fall the most in the presence of traffic congestion. For example, ignoring traffic congestion would cause one to identify a stretch along interstate I-5 as the one whose improvement would yield the greatest welfare gains for the city. Accounting for the endogenous change in traffic congestion throughout the whole network, the aggregate welfare elasticity to improving this link is not even in the top fifty of links. 

Finally, we combine these welfare elasticities with estimated costs of construction to estimate a return on investment for each link of the Seattle road network. We proceed analogously to the U.S. highway system case, first calculating the necessary lane-miles to achieve a $1 \%$ reduction in transportation costs, assuming 10% of construction costs are incurred each year, and then using a compensating variation approach to assign a dollar value to the aggregate welfare gains.47 We find that improving the average link in Seattle 

yields an annual return of 16.8% for the residents of the city.48 Like with the U.S. highway system, however, there is substantial heterogeneity, with returns varying from less than 25% to more than 250%. Panel (b) of Figure (7) shows the RoI for each segment; the highest returns are concentrated in the center of the city. Panel (b) of Table 2 lists the top 20 links in terms of their RoI; half of the list are either entirely within downtown Seattle or between downtown Seattle and another part of the city. Other locations with high returns on infrastructure improvement include the area around the University of Washington campus and Lake City Way in the neighborhood of North Seattle. On the other hand, we estimate that nearly half (331 of 692) links in the Seattle road network would generate negative returns of investment, highlighting the importance of well-targeted infrastructure improvements. 

# 7 Conclusion

This paper proposes a new spatial framework that incorporates traffic congestion and uses it to evaluate the welfare impact of transportation infrastructure improvements. In doing so, it combines the rich geography and general equilibrium structure of existing quantitative spatial models with the endogenous routing and traffic congestion of transportation models, but where both the distribution of economic activity and the resulting traffic patterns are determined jointly in equilibrium. 

The approach generates analytical expressions for transportation costs between any two locations, the traffic along each link of the transportation network, and the equilibrium spatial distribution of economic activity. This tractability not only allows us to characterize the equilibrium properties of the framework, but it also facilitates applying the framework to evaluate the welfare impacts of transportation infrastructure improvements empirically. Using readily available traffic data we show that for both the U.S. highway network and the Seattle road network, congestion matters for, where you improve the road network. 

The goal of this paper has been to provide a tractable framework that bridges the gap between the quantitative spatial and transportation economics literatures. Qe hope it can facilitate the answering of a number interesting and unresolved research questions, including: How does traffic congestion impact urban land use? What is the best way to design congestion tolls? How does the presence of multiple uses of transportation infrastructure (e.g. 

trade, commuting, consumption) interact in determining the spatial distribution of economic activity? We look forward to fruitful future research on these topics. 

# References



Adao, R., C. Arkolakis, and F. Esposito (2021): “General Equilibrium Effects in Space: Theory and Measurement,” mimeo. 





Ahlfeldt, G. M., S. J. Redding, D. M. Sturm, and N. Wolf (2015): “The economics of density: Evidence from the Berlin Wall,” Econometrica, 83(6), 2127–2189. 





Akamatsu, T. (1996): “Cyclic flows, Markov process and stochastic traffic assignment,” Transportation Research Part B: Methodological, 30(5), 369–386. 





Allen, T., and C. Arkolakis (2014): “Trade and the Topography of the Spatial Economy,” The Quarterly Journal of Economics. 





Allen, T., C. Arkolakis, and X. Li (2020): “On the Equilibrium Properties of Network Models with Heterogeneous Agents,” NBER Working Paper, (w27837). 





Allen, T., C. Arkolakis, and Y. Takahashi (2020): “Universal gravity,” Journal of Political Economy, 128(2), 393–433. 





American Association of State Highway and Transportation Officials (2016): “A Policy on Design Standards - Interstate System,” Discussion paper. 





Anderson, J. E., and E. Van Wincoop (2003): “Gravity with Gravitas: A Solution to the Border Puzzle,” American Economic Review, 93(1), 170–192. 





(2004): “Trade Costs,” Journal of Economic Literature, 42(3), 691–751. 





ASCE (2017): “2017 infrastructure report card,” . 





Beckmann, M. J., C. B. McGuire, and C. B. Winsten (1955): “Studies in the Economics of Transportation,” . 





Bell, M. G. (1995): “Alternatives to Dial’s logit assignment algorithm,” Transportation Research Part B: Methodological, 29(4), 287–295. 





Caliendo, L., F. Parro, E. Rossi-Hansberg, and P.-D. Sarte (2018): “The impact of regional and sectoral productivity changes on the US economy,” The Review of economic studies, 85(4), 2042–2096. 





CBO (1982): “The interstate highway system: Issues and options,” Congressional Budget Office of the United States. 





CFS (2012): “2012 Commodity Flow Survey Public Use Microdata,” Bureau of Transportation Statistics. 





CGIAR-CSI (2017): “Shuttle Radar Topography Mission,” Discussion paper, NASA. 





Chaney, T. (2018): “The gravity equation in international trade: An explanation,” Journal of Political Economy, 126(1), 150–177. 





Chartrand, G. (1977): Introductory graph theory. Courier Corporation. 





Combes, P.-P., and L. Gobillon (2015): “The empirics of agglomeration economies,” in Handbook of regional and urban economics, vol. 5, pp. 247–348. Elsevier. 





Costinot, A., and A. Rodr´ıguez-Clare (2014): “Trade theory with numbers: Quantifying the consequences of globalization,” in Handbook of international economics, vol. 4, pp. 197–261. Elsevier. 





Couture, V., G. Duranton, and M. A. Turner (2018): “Speed,” Review of Economics and Statistics, 100(4), 725–739. 





De Palma, A., M. Kilani, and R. Lindsey (2005): “Congestion pricing on a road network: A study using the dynamic equilibrium simulator METROPOLIS,” Transportation Research Part A: Policy and Practice, 39(7-9), 588–611. 





De Palma, A., R. Lindsey, E. Quinet, and R. Vickerman (2011): A handbook of transport economics. Edward Elgar Publishing. 





Dekle, R., J. Eaton, and S. Kortum (2008): “Global Rebalancing with Gravity: Measuring the Burden of Adjustment,” IMF Staff Papers, 55(3), 511–540. 





Disdier, A.-C., and K. Head (2008): “The puzzling persistence of the distance effect on bilateral trade,” The Review of Economics and statistics, 90(1), 37–48. 





Donaldson, D. (2015): “The gains from market integration,” economics, 7(1), 619–647. 





(2018): “Railroads of the Raj: Estimating the impact of transportation infrastructure,” American Economic Review, 108(4-5), 899–934. 





Donaldson, D., and R. Hornbeck (2016): “Railroads and American economic growth: A ”market access” approach,” The Quarterly Journal of Economics, 131(2), 799–858. 





Ducruet, C., R. Juhasz, D. K. Nagy, C. Steinwender, et al. ´ (2020): “All aboard: The effects of port development,” Discussion paper. 





Duranton, G., and M. A. Turner (2011): “The fundamental law of road congestion: Evidence from US cities,” American Economic Review, 101(6), 2616–52. 





Eaton, J., and S. Kortum (2002): “Technology, Geography and Trade,” Econometrica, 70(5), 1741–1779. 





Edwards, A. (2017): “U.S. Cities List,” Discussion paper. 





Eluru, N., A. Pinjari, J. Guo, I. Sener, S. Srinivasan, R. Copperman, and C. Bhat (2008): “Population updating system structures and models embedded in the comprehensive econometric microsimulator for urban systems,” Transportation Research Record: Journal of the Transportation Research Board, (2076), 171–182. 





EMC Research (2016): “2016 Center City Commuter Mode Split Survey,” Discussion paper, Commute Seattle. 





Fajgelbaum, P. D., and E. Schaal (2020): “Optimal transport networks in spatial equilibrium,” Econometrica, 88(4), 1411–1452. 





Fan, J., Y. Lu, and W. Luo (2019): “Valuing Domestic Transport Infrastructure: A View from the Route Choice of Exporters,” . 





Fan, J., and W. Luo (2020): “A Tractable Model of Transshipment,” . 





Federal Highway Administration (2015): “2015 Status of the Nation’s Highways, Bridges, and Transit: Conditions and Performance,” Discussion paper, US Department of Transportation. 





(2016): “Highway Performance Monitoring System Field Manual,” Discussion paper, Department of Transportation. 





Feigenbaum, B., M. G. Fields, and S. Purnell (2020): “25th Annual Highway Report,” Reason Foundation. 





Feyrer, J. (2019): “Trade and income: Exploiting time series in geography,” American Economic Journal: Applied Economics, 11(4), 1–35. 





FHA (2008): “2008 status of the nation’s highways, bridges, and transit: Conditions and perfomance,” Federal Highway Administration and Federal Transit Administration. 





Fogel, R. (1962): “A Quantitative Approach to the Study of Railroads in American Economic Growth: a Report of Some Preliminary Findings,” Journal of Economic History, 22(2), 163–197. 





Fogel, R. W. (1964): Railroads and American economic growth: essays in econometric history, vol. 296. Johns Hopkins Press Baltimore. 





Galichon, A. (2016): Optimal transport methods in economics. Princeton University Press. 





Ganapati, S., W. F. Wong, and O. Ziv (2020): “Entrepˆot: Hubs, Scale, and Trade Costs,” Discussion paper. 





Head, K., and J. Ries (2001): “Increasing Returns versus National Product Differentiation as an Explanation for the Pattern of U.S.-Canada Trade,” American Economic Review, 91(4), 858–876. 





Heblich, S., S. J. Redding, and D. M. Sturm (2020): “The making of the modern metropolis: evidence from London,” The Quarterly Journal of Economics, 135(4), 2059– 2133. 





Hillberry, R., and D. Hummels (2008): “Trade responses to geographic frictions: A decomposition using micro-data,” European Economic Review, 52(3), 527–550. 





HMPS (2016): “Highway Performance Monitoring System,” Discussion paper, United States Department of Transportation. 





Holden, D. J. (1989): “Wardrop’s third principle: urban traffic congestion and traffic policy,” Journal of Transport Economics and Policy, pp. 239–262. 





Hummels, D. L., and G. Schaur (2013): “Time as a trade barrier,” American Economic Review, 103(7), 2935–59. 





Karlin, S., and L. Nirenberg (1967): “On a theorem of P. Nowosad,” Journal of Mathematical Analysis and Applications, 17(1), 61–67. 





Kaszynski, W. (2000): The American Highway: The History and Culture of Roads in the United States. McFarland. 





Krugman, P. (1991): “Increasing Returns and Economic Geography,” The Journal of Political Economy, 99(3), 483–499. 





Lefevre, B., D. Leipziger, and M. Raifman (2014): “The Trillion Dollar Question: Tracking public and private investment in transport,” Washington, D. C: World Resources Institute. 





Lind, N., and N. Ramondo (2018): “Trade with correlation,” Discussion paper, National Bureau of Economic Research. 





LODES (2017): “Longitudinal Employer-Household Dynamics,” Discussion paper, United States Census Bureau. 





Mattsson, L.-G., J. W. Weibull, and P. O. Lindberg (2014): “Extreme values, invariance and choice probabilities,” Transportation Research Part B: Methodological, 59, 81–95. 





Monte, F., S. J. Redding, and E. Rossi-Hansberg (2018): “Commuting, migration, and local employment elasticities,” American Economic Review, 108(12), 3855–90. 





MPC (2011): “National Historical Information System: Version 2.0,” Minnesota Population Center, http://www.nhgis.org. 





NSTIFC (2009): “Paying our way: A new framework for transportation finance,” National Surface Transportation Infrastructure Financing Commission. 





Office of the State Auditor, D. o. T. A. (2002): Capitalization and Depreciation of Infrastructure. Mississippi Association of Governmental Purchasing and Property Agents. 





Osher, S., and J. Sethian (1988): “Fronts Propagating with Curvature-dependent Speed: Algorithms based on Hamilton-Jacobi Formulations,” Journal of computational physics, 79(1), 12–49. 





Pascali, L. (2017): “The wind of change: Maritime technology, trade, and economic development,” American Economic Review, 107(9), 2821–54. 





Polyanin, A., and A. Manzhirov (2008): Handbook of Integral Equations. Chapman & Hall/CRC. 





Ramondo, N., and A. Rodr´ıguez-Clare (2013): “Trade, Multinational Production, and the Gains from Openness,” forthcoming, Journal of Political Economy. 





Redding, S., and A. J. Venables (2004): “Economic geography and international inequality,” Journal of international Economics, 62(1), 53–82. 





Redding, S. J. (2016): “Goods trade, factor mobility and welfare,” Journal of International Economics, 101, 148–167. 





Redding, S. J., and E. Rossi-Hansberg (2017): “Quantitative spatial economics,” Annual Review of Economics, 9, 21–58. 





Redding, S. J., and M. A. Turner (2015): “Transportation costs and the spatial organization of economic activity,” in Handbook of regional and urban economics, vol. 5, pp. 1339–1398. Elsevier. 





Rosenthal, S. S., and W. C. Strange (2004): “Evidence on the nature and sources of agglomeration economies,” Handbook of regional and urban economics, 4, 2119–2171. 





Seattle (2017): “Seattle GeoData,” Discussion paper, City of Seattle. 





Sheffi, Y. (1985): Urban transportation networks, vol. 6. Prentice-Hall, Englewood Cliffs, NJ. 





Szabo, F. (2015): The linear algebra survival guide: illustrated with Mathematica. Academic Press. 





Tsitsiklis, J. (1995): “Efficient Algorithms for Globally Optimal Trajectories,” Automatic Control, IEEE Transactions on, 40(9), 1528–1538. 





Tsivanidis, N. (2018): “The aggregate and distributional effects of urban transit infrastructure: Evidence from bogot´a’s transmilenio,” . 





Vickrey, W. (1967): “Optimization of traffic and facilities,” Journal of Transport Economics and Policy, pp. 123–136. 





Washington State Legislature (1965): “RCW 46.61.400,” Electronic. 





Weber, H. J., and G. B. Arfken (2003): Essential Mathematical Methods for Physicists. Academic Press. 





Weber, S., and M. Peclat ´ (2017): “GEOROUTE: Stata module to calculate travel distance and travel time between two addresses or two geographical points,” . 



# Tables and Figures


Figure 1: Traffic Congestion and the Distribution of Economic Activity


(a) No congestion ( $\lambda = 0$ , $L = 1 0 0$ ) 

(b) Congestion, low scale ( $\lambda = 0 . 0 5$ , $L = 1 0 0$ ) 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/d3f56580acfbaf046b0f3ef56839e3c6e5187cf7761cae7e431128c242c80a97.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/24b3f83aeba5f64bc0fbeb4a1ec3aaeb2c38b60c9ee004aaefd99be155379255.jpg)



(c) Congestion, medium scale ( $\lambda ~ = ~ 0 . 0 5$ , $\bar { L } \ =$ 1000) （



(d) Congestion, large scale ( $\lambda = 0 . 0 5$ , $L = 1 0 0 0 0$ )


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/65add94030094e98b6cb1855e35783195a7c15709166673c6bdaf387cc3ca7d2.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/a1d03ad5b90cb244ccfb699486212e8782675c20fb1457bc6e2b1e0c61a95eb0.jpg)


Notes: This figure shows how traffic congestion ( $\lambda$ ) and the scale of the economy   $L$  shapes the distribution of economic activity within an example 5x5 grid network using the urban model. The height of the buildings (and the rooftop colors, associated with the color bar on the right) indicate the equilibrium residential population $( L _ { i } ^ { R } )$ at each location in the city, and the color of each link (associated with the color bar on the left) indicates the equilibrium traffic along the link. Throughout, $\alpha = \beta = 0$ , $\theta = 4$ , and $t _ { k l } = 1 . 5$ for connected links. 


Figure 2: Predicting Flows Using Traffic


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/8c83726f9c87545829f7a7514bda59c8c6fe4662528aba1773ae7dba502a70a8.jpg)



(a) Trade flows in an economic geography model


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/edc0783e3f9c4dc19ff744cb919c54fccb131169970508fd489c1fa25559e1ca.jpg)



(b) Commuting flows in an urban model


Notes: This figure compares the observed bilateral origin to destination flows to those predicted from the observed traffic along the transportation network. In panel (a), we compare the predicted (log) trade flows on the x-axis to the observed (log) trade flows between metropolitan areas from the Commodity Flow Survey (CFS) data on the y-axis using the economic geography model. In panel (b), we compare the predicted (log) commuting flows on the x-axis to the observed (log) commuting flows from the Longitudinal Employer-Household Dynamics Origin-Destination Employment Statistics (LODES) between grid cells within Seattle. In both figures, the predicted and observed flows are residualized using origin and destination fixed effects, so the observed correlation only arises through similarity at the pair level. 50 50 

Figure 3: Transportation Systems and their Network Representations 

Panel A: The Interstate Highway System 

(a) U.S. Highway Network 

Traffic (AADT) 

≤11580 

≤17260 

≤20780 

≤26470 

≤35650 

≤50490 

≤74460 

≤113200 

≤175800 

≤276850 


Node Population


≤511901 

≤1308587 

·≤2252276 

•≤4532390 

≤14745610 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/3047e52b93477d9b189e96a155e19d882721929be39cb8ba267bed98e1a4421d.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/38b31eb8bb6b8b2cf078d9a1ef16a89d4bf71aee2df713c49ccad36cb10d8606.jpg)



(b) Seattle Road Network



Traffic (AADT)



≤4696



≤5453



≤10030



≤37680



≤204800



Node Population


≤2168 

≤4260 

≤7061 

≤11470 

≤18695 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/cb36ab918dea20b044340262cbc4eb95fec755ef5ba562ab6f8ede54b0c53045.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/15f3440e2b19b29d68f8783ced489411b5860a34b8a73e3997ae860476ff0764.jpg)


Notes: This figure presents the observed transportation network (on the top) and the constructed infrastructure matrix (on the bottom) for the U.S. highway network (panel a) and the observed transportation network (on the right) and the constructed infrastructure matrix (on the right) for the Seattle road network (panel b). In both panels, the size of each node reflects its population and the color of each link reflects the amount of traffic with red (blue) indicating high (low) levels of traffic. The gray roads in panel (b) are roads not on the least cost route between grid centers. 


Figure 4: Example Link Intensities


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/95d15ff184e497caec58a052c1c3db3560806242b7694ec4eabf2cbe81869173.jpg)



(a) Economic geography: Los Angeles, CA to New York, NY


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/566844df8d353052e93c7626b8abba52aeedfd1a8332a72fd6f0205e97e533bd.jpg)



(b) Urban model: Safeco Field to the University of Washington


Notes: This figure shows an example of the link intensity $\pi _ { i j } ^ { k l } -$ i.e. the expected number of traverses New York in panel (a) and from Safeco Field to the University of Washington in panel (b). These link intensities are calculated using only observed data on traffic flows and the economic activity in each location (i.e. no assumptions on model parameters are necessary); see Online Appendix B.4 for details. 


Figure 5: Welfare Elasticities of Infrastructure Improvement



(a) U.S. Highway Network


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/f4d9e3c744c5759b78e01c361689f29e94ada0f4f88d82aaa2d63ce31dcb5544.jpg)



(b) Seattle Road Network


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/f419bb63952328eab85cacfe502b237fc034f7d9f0a409ad1d028a453a391af7.jpg)


Notes: This figure presents the elasticity of aggregate welfare to improving each link in the U.S. Highway Network (panel A) and the Seattle road network (Panel B). The color ramp goes from blue (lower welfare elasticity) to red (higher welfare elasticity). Nodes in the network are marked by the black circles, which are increasing the population size of the node. 


Figure 6: Comparing Welfare Elasticities With and Without Congestion



(a) U.S. Highway Network


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/5edafcfec610cc436fbf3df7f42df0a577d8c2d6b97ce2fed383df79d48658f3.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/dbf837a9a558e156e3597cb7e3fbc0ca0a4e12a9763687c610c4c06b32d49b59.jpg)



(b) Seattle Road Network


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/5802764b9d1001cb75e74c8c4948ffca76c31f43609743fd34e4a8ab2b67f4c7.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/28989a057fdf5bfba5a1372bada5c9d7cc3a61f756f68e0b1394768454f16d57.jpg)


Notes: This figure compares the welfare elasticity calculated allowing for traffic congestion (given the estimated strength of congestion $\lambda$ ) to the welfare elasticity that would be calculated if traffic congestion were ignored (i.e. if $\lambda = 0$ ), as in a standard spatial model for each link in the U.S. highway network (panel a) and the Seattle road network (panel b). The left figure in the panel shows the difference in welfare with and without congestion across links, whereas the right figure panel shows a scatter plot of the two estimated elasticities. 


Figure 7: Returns On Investment of Infrastructure Improvement



(a) U.S. Highway Network


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/7662e5848de3129a0f585de391f75353b64e969896daa543eaab3eeba096556e.jpg)



(b) Seattle Road Network


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/e3f119dd319b1fddf891654ecc874d20c960329540d7ad9d1b3d0e257ddfad3b.jpg)


Notes: This figure presents the return on investment of improving links in the Interstate Highway System (Panel A) and the Seattle road network (Panel B). Return on investment is annual and in decimals of the initial investment (i.e. 0.75 means a 75% return on initial investment per annum). The color ramps goes from blue (negative returns) to red (high positive returns). Nodes in the network are marked by the black circles, which are increasing the population size of the node. 


stimating the strength of traffic co


<table><tr><td></td><td colspan="5">Panel A: Interstate Highway System</td></tr><tr><td></td><td>(1) OLS</td><td>(2) OLS</td><td>(3) IV: 1st stage</td><td colspan="2">(4) IV: 2nd stage</td></tr><tr><td>Log congestion</td><td>0.109***</td><td>0.050***</td><td></td><td colspan="2">0.739***</td></tr><tr><td></td><td>(0.010)</td><td>(0.012)</td><td></td><td colspan="2">(0.181)</td></tr><tr><td>Log distance</td><td></td><td></td><td>-0.156***</td><td></td><td></td></tr><tr><td></td><td></td><td></td><td>(0.033)</td><td></td><td></td></tr><tr><td>Start-location FE</td><td>No</td><td>Yes</td><td>Yes</td><td colspan="2">Yes</td></tr><tr><td>End-location FE</td><td>No</td><td>Yes</td><td>Yes</td><td colspan="2">Yes</td></tr><tr><td>F-statistic</td><td>120.994</td><td>17.929</td><td>22.748</td><td colspan="2">16.723</td></tr><tr><td>Observations (excl. singletons)</td><td>630</td><td>630</td><td>630</td><td colspan="2">630</td></tr><tr><td>Observations (incl. singletons)</td><td>704</td><td>704</td><td>704</td><td colspan="2">704</td></tr><tr><td colspan="6">Panel B: Seattle Road Network</td></tr><tr><td></td><td>(1) OLS</td><td>(2) IV: 1st stage</td><td>(3) IV</td><td>(4) IV: 1st stage</td><td>(5) IV</td></tr><tr><td>AADT per Lane</td><td>-0.048***</td><td></td><td>0.118**</td><td></td><td>0.488*</td></tr><tr><td></td><td>(0.007)</td><td></td><td>(0.048)</td><td></td><td>(0.278)</td></tr><tr><td>Turns along Route</td><td></td><td>-0.252***</td><td></td><td>-0.091**</td><td></td></tr><tr><td></td><td></td><td>(0.049)</td><td></td><td>(0.039)</td><td></td></tr><tr><td>Start-location FE</td><td>Yes</td><td>Yes</td><td>Yes</td><td>Yes</td><td>Yes</td></tr><tr><td>End-location FE</td><td>Yes</td><td>Yes</td><td>Yes</td><td>Yes</td><td>Yes</td></tr><tr><td>No. of Intersections</td><td>No</td><td>Yes</td><td>Yes</td><td>Yes</td><td>Yes</td></tr><tr><td>Bilateral Route Quality</td><td>No</td><td>No</td><td>No</td><td>Yes</td><td>Yes</td></tr><tr><td>F-statistic</td><td>41.546</td><td>26.347</td><td>6.191</td><td>5.336</td><td>3.084</td></tr><tr><td>Observations</td><td>1338</td><td>1338</td><td>1338</td><td>1338</td><td>1338</td></tr></table>


end-location are reported in parentheses. Stars indicate statistical significance: * p<.10 cations (functional system are taken from the HPMS. For both panels, standard errors which are generated by binning each segment into deciles based on the shares of arterial ) (.     ,   he (log) number of turns, conditional on number of intersections traversed. In column5 nefromthehighwayperformancemonitoringsystemHPMSIncolumns3and5wein e is the (log time of travel per unit distance, calculated using the HERE API, and the in ) ,     (      (    eter estimates for Seattle, and each observation is a segment of the Seattle’s Network. I ) ,       (       olumn3weinstrumentforthelog)trafficperlaneusingthelog)lengthofthesegme , ,  ,      (   EREAPIandtheindependentvariableisthelog)AADTperlanefromthehighwaype waynetwork.Incolumns12and4thedependentvariableisthelogtimeoftra ts the congestion parameter estimates for the Interstate Highway System, and each obs ) 



Return on investment ra


<table><tr><td colspan="7">Panel A: Interstate Highway System</td></tr><tr><td></td><td>Origin</td><td>Destination</td><td>Rt. Number</td><td>RoI</td><td>Benefit ($m)</td><td>Cost ($m)</td></tr><tr><td>1</td><td>Kingsport-Bristol-Bristol, TN-VA</td><td>Johnson City, TN</td><td>I-26</td><td>10.43</td><td>55.492</td><td>5.271</td></tr><tr><td>2</td><td>Greensboro-High Point, NC</td><td>Winston-Salem, NC</td><td>I-40</td><td>9.54</td><td>158.548</td><td>16.455</td></tr><tr><td>3</td><td>Rochester, NY</td><td>Batavia, NY</td><td>I-90</td><td>7.31</td><td>117.779</td><td>15.896</td></tr><tr><td>4</td><td>Springfield, OH</td><td>Dayton, OH</td><td>I-70</td><td>6.38</td><td>175.536</td><td>27.07</td></tr><tr><td>5</td><td>Durham-Chapel Hill, NC</td><td>Raleigh-Cary, NC</td><td>I-40</td><td>5.95</td><td>274.062</td><td>45.313</td></tr><tr><td>6</td><td>Grand Rapids-Wyoming, MI</td><td>Muskegon-Norton Shores, MI</td><td>I-96</td><td>5.32</td><td>108.719</td><td>20.075</td></tr><tr><td>7</td><td>Toledo, OH</td><td>Monroe, MI</td><td>I-75</td><td>4.97</td><td>159.681</td><td>31.504</td></tr><tr><td>8</td><td>Vallejo-Fairfield, CA</td><td>Sacramento-Arden-Arcade-Roseville, CA</td><td>I-80</td><td>4.51</td><td>236.69</td><td>51.33</td></tr><tr><td>9</td><td>Breezewood, PA</td><td>Bedford, PA</td><td>I-76</td><td>4.27</td><td>41.902</td><td>9.598</td></tr><tr><td>10</td><td>Rockford, IL</td><td>Rochelle, IL</td><td>I-39</td><td>4.02</td><td>75.877</td><td>18.411</td></tr></table>


l B: Seattle Road Net


<table><tr><td></td><td>Home Neighborhood</td><td>Home Address</td><td>Work Neighborhood</td><td>Work Address</td><td>RoI</td><td>Benefit ($m)</td><td>Cost ($m)</td></tr><tr><td>1</td><td>Downtown Seattle</td><td>805 S Jackson St</td><td>Downtown Seattle</td><td>83 Spring St</td><td>2.3</td><td>24.587</td><td>10.26</td></tr><tr><td>2</td><td>Downtown Seattle</td><td>83 Spring St</td><td>Downtown Seattle</td><td>1001 Columbia St</td><td>1.3</td><td>5.316</td><td>3.794</td></tr><tr><td>3</td><td>Downtown Seattle</td><td>805 S Jackson St</td><td>Downtown Seattle</td><td>1001 Columbia St</td><td>1.29</td><td>21.189</td><td>15.201</td></tr><tr><td>4</td><td>Capitol Hill</td><td>1698 Yale Ave</td><td>Capitol Hill</td><td>107A 13th Ave E</td><td>1.03</td><td>2.598</td><td>2.296</td></tr><tr><td>5</td><td>Downtown Seattle</td><td>2020 4th Ave</td><td>Capitol Hill</td><td>1698 Yale Ave</td><td>1.01</td><td>4.684</td><td>4.207</td></tr><tr><td>6</td><td>Downtown Seattle</td><td>2020 4th Ave</td><td>South Lake Union</td><td>401 Aurora Ave N</td><td>.91</td><td>3.474</td><td>3.44</td></tr><tr><td>7</td><td>Downtown Seattle</td><td>2020 4th Ave</td><td>Lower Queen Anne</td><td>201 Queen Anne Ave N</td><td>.89</td><td>7.973</td><td>8.024</td></tr><tr><td>8</td><td>Capitol Hill</td><td>1698 Yale Ave</td><td>South Lake Union</td><td>1214 E Roy St</td><td>.87</td><td>6.596</td><td>6.823</td></tr><tr><td>9</td><td>SoDo</td><td>562 1st Ave S</td><td>Downtown Seattle</td><td>83 Spring St</td><td>.79</td><td>3.021</td><td>3.404</td></tr><tr><td>10</td><td>Downtown Seattle</td><td>83 Spring St</td><td>Downtown Seattle</td><td>2020 4th Ave</td><td>.76</td><td>5.619</td><td>6.518</td></tr><tr><td>11</td><td>South Lake Union</td><td>401 Aurora Ave N</td><td>South Lake Union</td><td>1214 E Roy St</td><td>.75</td><td>4.738</td><td>5.545</td></tr><tr><td>12</td><td>North Seattle</td><td>2517 NE 120th St</td><td>Cedar Park</td><td>13504 WA-522</td><td>.75</td><td>2.026</td><td>2.381</td></tr><tr><td>13</td><td>Northeast Seattle</td><td>4200 Mary Gates Memorial Dr NE N207</td><td>Northeast Seattle</td><td>4700 45th Ave NE</td><td>.74</td><td>1.261</td><td>1.5</td></tr><tr><td>14</td><td>Downtown Seattle</td><td>1001 Columbia St</td><td>Downtown Seattle</td><td>2020 4th Ave</td><td>.73</td><td>5.001</td><td>6.052</td></tr><tr><td>15</td><td>Downtown Seattle</td><td>83 Spring St</td><td>Capitol Hill</td><td>1698 Yale Ave</td><td>.71</td><td>4.261</td><td>5.259</td></tr><tr><td>16</td><td>Greater Duwamish</td><td>1150 S Atlantic St</td><td>Downtown Seattle</td><td>805 S Jackson St</td><td>.67</td><td>9.416</td><td>12.306</td></tr><tr><td>17</td><td>SoDo</td><td>1705 4th Ave S</td><td>SoDo</td><td>562 1st Ave S</td><td>.66</td><td>2.696</td><td>3.567</td></tr><tr><td>18</td><td>Downtown Seattle</td><td>1001 Columbia St</td><td>Capitol Hill</td><td>1698 Yale Ave</td><td>.65</td><td>3.802</td><td>5.066</td></tr><tr><td>19</td><td>Cedar Park</td><td>13504 WA-522</td><td>Olympic Hills</td><td>14322 19th Ave NE</td><td>.65</td><td>2.373</td><td>3.177</td></tr><tr><td>20</td><td>Downtown Seattle</td><td>2020 4th Ave</td><td>South Lake Union</td><td>1214 E Roy St</td><td>.64</td><td>4.735</td><td>6.4</td></tr></table>


in Panel A are reported in millions of chains 2012 USD, and benefits in Panel B are reported in millions of 2016 USD. means an annual 120% return. Maps Geocoding API. RoIs are reported as annualized decimal returns (i.e 1.2 means an annual 120% return). Benefits re sourced from the Go ) Pane and from Google Maps. Names and addresses of the home and work locations in Panel B are sourced from the Google . . ionandareaofCB he taon sourcedfromNH A ane Names of the origin and destinati 


# A Appendix: Derivations

The appendix presents the derivations of the results in Sections 3.1, 3.2, and 4.1; derivations for other sections are presented in Online Appendix B. 

# A.1 Section 3.1: Transportation Costs

Define the $N \times N$ matrix ${ \bf A } = \left[ a _ { i j } \equiv t _ { i j } ^ { - \theta } \right]$ . We can write $\tau _ { i j }$ from equation 4 by explicitly summing across all possible routes of all possible lengths. To do so, we sum across all locations that are traveled through all the possible paths as follows: 

$$
\tau_ {i j} ^ {- \theta} = \sum_ {K = 0} ^ {\infty} \left(\sum_ {k _ {1} = 1} ^ {N} \sum_ {k _ {2} = 1} ^ {N} \dots \sum_ {k _ {K - 1} = 1} ^ {N} a _ {i, k _ {1}} \times a _ {k _ {1}, k _ {2}} \times \dots \times a _ {k _ {K - 2}, k _ {K - 1}} \times a _ {k _ {K - 1}, j}\right),
$$

where $k _ { n }$ is the sub-index for the $n ^ { t h }$ location arrived at on a particular route. Note that pairs of locations that are not connected will have $a _ { i j } = 0$ , so that infeasible routes do not affect the sum. The portion of the expression in the parentheses is equivalent to the $( i , j )$ element of the weighted adjacency matrix to the power $K$ , i.e.: 

$$
\tau_ {i j} ^ {- \theta} = \sum_ {K = 0} ^ {\infty} \mathbf {A} _ {i j} ^ {K},
$$

where $\mathbf { A } ^ { K } = \left\lfloor A _ { i j } ^ { K } \right\rfloor$ , i.e. $A _ { i j } ^ { K }$ is the $( i , j )$ element of the matrix A to the matrix power $K$ . As we note in the paper, for a matrix A with spectral radius less than one, the geometric sum can be expressed as: 

$$
\sum_ {K = 0} ^ {\infty} \mathbf {A} ^ {K} = (\mathbf {I} - \mathbf {A}) ^ {- 1} \equiv \mathbf {B},
$$

where $\mathbf { B } = [ b _ { i j } ]$ is the Leontief inverse of the weighted adjacency matrix, so the transportation cost from $_ i$ to $j$ can be expressed as a function of the infrastructure matrix: 

$$
\tau_ {i j} ^ {- \theta} = b _ {i j},
$$

as in equation (21). 

# A.2 Section 3.2: Traffic Flows

Beginning with equation (22) we have: 

$$
\pi_ {i j} ^ {k l} = \sum_ {r \in \Re_ {i j}} \frac {\pi_ {i j , r}}{\sum_ {r ^ {\prime} \in \Re_ {i j}} \pi_ {i j , r ^ {\prime}}} n _ {r} ^ {k l} \iff
$$

$$
\pi_ {i j} ^ {k l} = \sum_ {r \in \Re_ {i j}} \frac {\left(\prod_ {l = 1} ^ {K} t _ {r _ {l - 1} , r _ {l}} ^ {- \theta}\right)}{\sum_ {r \in \Re_ {i j}} \left(\prod_ {l = 1} ^ {K} t _ {r _ {l - 1} , r _ {l}} ^ {- \theta}\right)} n _ {r} ^ {k l} \iff
$$

$$
\pi_ {i j} ^ {k l} = \tau_ {i j} ^ {\theta} \sum_ {r \in \Re_ {i j}} \prod_ {l = 1} ^ {K} t _ {r _ {l - 1}, r _ {l}} ^ {- \theta} n _ {r} ^ {k l},
$$

where the second line used either equation (2) (for the economic geography model) or equation (12) (for the urban model), and the third line used the definition of $\tau _ { i j }$ from equation (4). 

For each route in multiplied by the num $\boldsymbol { r } \in \Re _ { i j }$ , the value  imes the ro $\begin{array} { r } { \prod _ { l = 1 } ^ { K } t _ { r _ { l - 1 } , r _ { l } } ^ { - \theta } n _ { r } ^ { k l } } \end{array}$ k ansportation costs incurred along the route. To calculate this, we proceed by summing $\{ k , l \}$ across all possible traverses that occur on all routes from $_ i$ to $j$ . To do so, note for any $r \in \Re _ { i j }$ of length $K$ (the set of which we denote as $\Re _ { i j , K }$ ), a traverse is possible at any point $B \in [ 1 , 2 , . . . , K - 1 ]$ in the route. 

Defining $\mathbf { A } \equiv [ a _ { k l } ] = \left[ t _ { k l } ^ { - \theta } \right]$ and ${ \bf B } \equiv [ b _ { i j } ] = \left[ \tau _ { i j } ^ { - \theta } \right]$ as above, we can write: 

$$
\pi_ {i j} ^ {k l} = \frac {1}{b _ {i j}} \sum_ {K = 0} ^ {\infty} \sum_ {B = 0} ^ {K - 1} \left(\sum_ {r \in \Re_ {i k, B}} \prod_ {n = 1} ^ {B} a _ {r _ {n - 1}, r _ {n}}\right) \times a _ {k l} \times \left(\sum_ {r \in \Re_ {k j, K - B - 1}} \prod_ {n = 1} ^ {K - B - 1} a _ {r _ {n - 1}, r _ {n}}\right)
$$

This can in turn allows us to explicitly enumerate all possible paths from $_ i$ to $k$ of length $B$ and all possible paths from $\it l$ to $j$ of length $K - B - 1$ : 

$$
\pi_ {i j} ^ {k l} = \frac {1}{b _ {i j}} \sum_ {K = 0} ^ {\infty} \sum_ {B = 0} ^ {K - 1} \left(\sum_ {n _ {1} = 1} ^ {N} \dots \sum_ {n _ {B - 1} = 1} ^ {N} a _ {i, n _ {1}} \times \ldots \times a _ {n _ {B - 1}, k}\right) \times a _ {k l} \times \left(\sum_ {n _ {1} = 1} ^ {N} \dots \sum_ {n _ {K - B - 1} = 1} ^ {N} a _ {l, n _ {1}} \times \ldots \times a _ {n _ {K - B - 1}, j}\right),
$$

which can be expressed more succinctly as elements of matrix powers of A : 

$$
\pi_ {i j} ^ {k l} = \frac {1}{b _ {i j}} \sum_ {K = 0} ^ {\infty} \sum_ {B = 0} ^ {K - 1} A _ {i k} ^ {B} \times a _ {k l} \times A _ {l j} ^ {K - B - 1}.
$$

A result from matrix calculus (see e.g. Weber and Arfken (2003)) is for any $N \times N$ matrix $\mathbf { C }$ we have: 

$$
\sum_ {K = 0} ^ {\infty} \sum_ {B = 0} ^ {K - 1} \mathbf {A} ^ {B} \mathbf {C} \mathbf {A} ^ {K - B - 1} = (\mathbf {I} - \mathbf {A}) ^ {- 1} \mathbf {C} (\mathbf {I} - \mathbf {A}) ^ {- 1}. \tag {A.1}
$$

Define $\mathbf { C }$ to be an $N \times N$ matrix that takes the value of $a _ { k l }$ at row $k$ and column $\it l$ and zeros everywhere else. Using equation (A.1) we obtain our result: 

$$
\pi_ {i j} ^ {k l} = \frac {b _ {i k} a _ {k l} b _ {l j}}{b _ {i j}} \Longleftrightarrow
$$

$$
\pi_ {i j} ^ {k l} = \frac {\tau_ {i k} ^ {- \theta} t _ {k l} ^ {- \theta} \tau_ {l j} ^ {- \theta}}{\tau_ {i j} ^ {- \theta}}, \tag {A.2}
$$

as in equation (23). 

We now derive gravity equations for traffic over a link for both economic geography and commuting models. For trade, we sum over all trade between all origins and destinations, and all routes taken by that trade, to get: 

$$
\begin{array}{l} \Xi_ {k l} = \sum_ {i \in \mathcal {N}} \sum_ {j \in \mathcal {N}} \sum_ {r \in \Re_ {i j}} \pi_ {i j, r} n _ {r} ^ {k l} E _ {j} \Longleftrightarrow \\ \Xi_ {k l} = \sum_ {i \in \mathcal {N}} \sum_ {j \in \mathcal {N}} \pi_ {i j} ^ {k l} X _ {i j} \iff \\ \Xi_ {k l} = \sum_ {i \in \mathcal {N}} \sum_ {j \in \mathcal {N}} \frac {\tau_ {i k} ^ {- \theta} t _ {k l} ^ {- \theta} \tau_ {l j} ^ {- \theta}}{\tau_ {i j} ^ {- \theta}} \times \tau_ {i j} ^ {- \theta} \frac {Y _ {i}}{\Pi_ {i} ^ {- \theta}} \frac {E _ {j}}{P _ {j} ^ {- \theta}} \Longleftrightarrow \\ \Xi_ {k l} = t _ {k l} ^ {- \theta} \sum_ {i \in \mathcal {N}} \tau_ {i k} ^ {- \theta} \frac {Y _ {i}}{\Pi_ {i} ^ {- \theta}} \sum_ {j \in \mathcal {N}} \tau_ {l j} ^ {- \theta} \frac {E _ {j}}{P _ {j} ^ {- \theta}}, \\ \end{array}
$$

where the second line used equations (2) and (22), the third lines used equation (23), and the fourth lined rearranged. Recalling our definition of the consumer and producer market access terms (equations (8) and (7)) in the text, this becomes: 

$$
\Xi_ {k l} = t _ {k l} ^ {- \theta} \times P _ {k} ^ {- \theta} \times \Pi_ {l} ^ {- \theta},
$$

as in equation (24). 

Turning to the commuting model, we proceed similarly, summing over all commuting flow pairs and the routes they take: 

$$
\begin{array}{l} \Xi_ {k l} = \sum_ {i \in \mathcal {N}} \sum_ {j \in \mathcal {N}} \sum_ {r \in \Re_ {i j}} \pi_ {i j, r} n _ {r} ^ {k l} \bar {L} \iff \\ \Xi_ {k l} = \sum_ {i \in \mathcal {N}} \sum_ {j \in \mathcal {N}} \pi_ {i j} ^ {k l} L _ {i j} \iff \\ \Xi_ {k l} = \sum_ {i \in \mathcal {N}} \sum_ {j \in \mathcal {N}} \frac {\tau_ {i k} ^ {- \theta} t _ {k l} ^ {- \theta} \tau_ {l j} ^ {- \theta}}{\tau_ {i j} ^ {- \theta}} \times \tau_ {i j} ^ {- \theta} \frac {L _ {i} ^ {R}}{\Pi_ {i} ^ {- \theta}} \frac {L _ {j} ^ {F}}{P _ {j} ^ {- \theta}} \iff \\ \Xi_ {k l} = \frac {\bar {L}}{W ^ {\theta}} \times t _ {k l} ^ {- \theta} \times \left(\sum_ {i \in \mathcal {N}} \tau_ {i k} ^ {- \theta} \frac {L _ {i} ^ {R}}{\Pi_ {i} ^ {- \theta}}\right) \times \left(\sum_ {j \in \mathcal {N}} \tau_ {l j} ^ {- \theta} \frac {L _ {j} ^ {F}}{P _ {j} ^ {- \theta}}\right), \\ \end{array}
$$

where the second line used equations (12) and and (22), the third lines used equation (23), and the fourth lined rearranged. We substitute in the consumer and producer market access defined in (17) and (18) to generate traffic gravity for the commuting framework: 

$$
\Xi_ {k l} = t _ {k l} ^ {- \theta} \times (P _ {k}) ^ {- \theta} \times (\Pi_ {l}) ^ {- \theta},
$$

again as in equation (24). 

# A.3 Section 4.1: Equilibrium

Trade Model In this Appendix section, we derive the equilibrium conditions for the economic geography and commuting frameworks. 

For the trade equilibrium conditions, we start with equation (10) from the paper. Note that $\tau _ { i j } ^ { - \theta } =$ $[ \mathbf { I } - \mathbf { A } ] _ { i j } ^ { - 1 }$ , where $\mathbf { A } = [ a _ { i j } ] = \left[ t _ { i j } ^ { - \theta } \right]$ is the adjacency matrix, so with a change of notation, we can rewrite the summation term as a matrix product: 

$$
\begin{array}{l} \left[ \bar {A} _ {i} ^ {- \theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {- \theta (1 + \alpha)} \right] = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \times \left[ \tau_ {i j} ^ {- \theta} \right] \times \left[ \bar {u} _ {j} ^ {\theta} y _ {j} ^ {1 + \theta} l _ {j} ^ {\theta (\beta - 1)} \right] \Longleftrightarrow \\ \left[ \bar {A} _ {i} ^ {- \theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {- \theta (1 + \alpha)} \right] = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \times [ \mathbf {I} - \mathbf {A} ] ^ {- 1} \times \left[ \bar {u} _ {j} ^ {\theta} y _ {j} ^ {1 + \theta} l _ {j} ^ {\theta (\beta - 1)} \right] \\ \end{array}
$$

where $\left[ \bar { A } _ { i } ^ { - \theta } y _ { i } ^ { 1 + \theta } l _ { i } ^ { - \theta ( 1 + \alpha ) } \right]$ and $\left[ \bar { u } _ { j } ^ { \theta } { y _ { j } } ^ { 1 + \theta } l _ { j } ^ { \theta ( \beta - 1 ) } \right]$ are column vectors. Taking a matrix inversion and converting back to summation notation: 

$$
\begin{array}{l} [ \mathbf {I} - \mathbf {A} ] \times \left[ \bar {A} _ {i} ^ {- \theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {- \theta (1 + \alpha)} \right] = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \times \left[ \bar {u} _ {i} ^ {\theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {\theta (\beta - 1)} \right] \Longleftrightarrow \\ \left[ \bar {A} _ {i} ^ {- \theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {- \theta (1 + \alpha)} \right] - \mathbf {A} \times \left[ \bar {A} _ {j} ^ {- \theta} y _ {j} ^ {1 + \theta} l _ {j} ^ {- \theta (1 + \alpha)} \right] = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \times \left[ \bar {u} _ {i} ^ {\theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {\theta (\beta - 1)} \right] \Longleftrightarrow \\ \bar {A} _ {i} ^ {- \theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {- \theta (1 + \alpha)} - \sum_ {j} a _ {i j} \bar {A} _ {j} ^ {- \theta} y _ {j} ^ {1 + \theta} l _ {j} ^ {- \theta (1 + \alpha)} = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {u} _ {i} ^ {\theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {\theta (\beta - 1)} \\ \end{array}
$$

The second equilibrium condition, equation (11), can also be written as a matrix multiplication, where $\left[ \bar { u } _ { i } ^ { - \theta } y _ { i } ^ { - \theta } l _ { i } ^ { \theta ( 1 - \beta ) } \right]$ −θ and $\left[ \bar { A } _ { j } ^ { \theta } y _ { j } ^ { - \theta } l _ { j } ^ { \theta ( \alpha + 1 ) } \right]$ are row vectors. Applying the same matrix inversion we did to equiu¯i librium equation we did to the first equilibrium condition, we get: 

$$
\bar {u} _ {i} ^ {- \theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (1 - \beta)} = \frac {\bar {L} ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \sum_ {j = 1} ^ {N} \tau_ {j i} ^ {- \theta} \bar {A} _ {j} ^ {\theta} y _ {j} ^ {- \theta} l _ {j} ^ {\theta (\alpha + 1)} \iff
$$

$$
\left[ \bar {u} _ {i} ^ {- \theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (1 - \beta)} \right] = \frac {\bar {L} ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \times \left[ \bar {A} _ {j} ^ {\theta} y _ {j} ^ {- \theta} l _ {j} ^ {\theta (\alpha + 1)} \right] \times \left[ \tau_ {j i} ^ {- \theta} \right] \Longleftrightarrow
$$

$$
\left[ \bar {u} _ {i} ^ {- \theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (1 - \beta)} \right] = \frac {\bar {L} ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \times \left[ \bar {A} _ {j} ^ {\theta} y _ {j} ^ {- \theta} l _ {j} ^ {\theta (\alpha + 1)} \right] \times \left[ \mathbf {I} - \mathbf {A} ^ {T} \right] ^ {- 1} \Longleftrightarrow
$$

$$
\left[ \bar {u} _ {i} ^ {- \theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (1 - \beta)} \right] \times \left[ \mathbf {I} - \mathbf {A} ^ {T} \right] = \frac {\bar {L} ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \times \left[ \bar {A} _ {i} ^ {\theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (\alpha + 1)} \right] \Longleftrightarrow
$$

$$
\left[ \bar {u} _ {i} ^ {- \theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (1 - \beta)} \right] - \left[ \bar {u} _ {j} ^ {- \theta} y _ {j} ^ {- \theta} l _ {j} ^ {\theta (1 - \beta)} \right] \times \mathbf {A} ^ {T} = \frac {\bar {L} ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \times \left[ \bar {A} _ {i} ^ {\theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (\alpha + 1)} \right] \Longleftrightarrow
$$

$$
\bar {u} _ {i} ^ {- \theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (1 - \beta)} - \sum_ {j} a _ {j i} \bar {u} _ {j} ^ {- \theta} y _ {j} ^ {- \theta} l _ {j} ^ {\theta (1 - \beta)} = \frac {\bar {L} ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {A} _ {i} ^ {\theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (\alpha + 1)}
$$

Recalling that $a _ { i j } \equiv t _ { i j } ^ { - \theta }$ , we have for our two equilibrium conditions (before incorporating traffic congestion): 

$$
\bar {A} _ {i} ^ {- \theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {- \theta (1 + \alpha)} = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {u} _ {i} ^ {\theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {\theta (\beta - 1)} + \sum_ {j} t _ {i j} ^ {- \theta} \bar {A} _ {j} ^ {- \theta} y _ {j} ^ {1 + \theta} l _ {j} ^ {- \theta (1 + \alpha)}
$$

$$
\bar {u} _ {i} ^ {- \theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (1 - \beta)} = \frac {\bar {L} ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {A} _ {i} ^ {\theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (\alpha + 1)} + \sum_ {j} t _ {j i} ^ {- \theta} \bar {u} _ {j} ^ {- \theta} y _ {j} ^ {- \theta} l _ {j} ^ {\theta (1 - \beta)}
$$

To incorporate congestion, we combine these two equations with the expression (26), converting from market access terms to equilibrium $\{ y _ { i } \}$ and $\{ l _ { i } \}$ (as in Appendix C.1). Starting with the first equilibrium condition: 

$$
\begin{array}{l} \bar {A} _ {i} ^ {- \theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {- \theta (1 + \alpha)} = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {u} _ {i} ^ {\theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {\theta (\beta - 1)} + \\ \sum_ {j} \left(\left(\bar {t} _ {i j} \bar {L} ^ {\lambda}\right) ^ {\frac {1}{1 + \theta \lambda}} \left(\frac {\bar {L} ^ {- (\alpha + \beta) \theta}}{W ^ {- \theta}}\right) ^ {\frac {\lambda}{1 + \theta \lambda}} \bar {A} _ {j} ^ {- \frac {\theta \lambda}{1 + \theta \lambda}} \bar {u} _ {i} ^ {- \frac {\theta \lambda}{1 + \theta \lambda}} l _ {i} ^ {- \frac {\theta \lambda (\beta - 1)}{1 + \theta \lambda}} l _ {j} ^ {- \frac {\theta \lambda (1 + \alpha)}{1 + \theta \lambda}} y _ {i} ^ {- \frac {\theta \lambda}{1 + \theta \lambda}} y _ {j} ^ {\frac {\lambda (1 + \theta)}{1 + \theta \lambda}}\right) ^ {- \theta} \times \\ \bar {A} _ {j} ^ {- \theta} y _ {j} ^ {1 + \theta} l _ {j} ^ {- \theta (1 + \alpha)} \Longleftrightarrow \\ \end{array}
$$

$$
\begin{array}{l} \bar {A} _ {i} ^ {- \theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {- \theta (1 + \alpha)} = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {u} _ {i} ^ {\theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {\theta (\beta - 1)} + \\ \sum_ {j} \left(\bar {t} _ {i j} \bar {L} ^ {\lambda}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \left(\frac {\bar {L} ^ {(\alpha + \beta) \theta}}{W ^ {\theta}}\right) ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \bar {A} _ {j} ^ {- \frac {\theta}{1 + \theta \lambda}} \bar {u} _ {i} ^ {\theta} \frac {\frac {\theta \lambda}{1 + \theta \lambda}}{l _ {i} ^ {\theta (\beta - 1)}} \frac {\frac {\theta \lambda}{1 + \theta \lambda}}{y _ {i} ^ {\theta} \frac {\frac {\theta \lambda}{1 + \theta \lambda}}{l _ {j} ^ {1 + \theta \lambda}}} y _ {j} ^ {\frac {1 + \theta}{1 + \theta \lambda}} l _ {j} ^ {- \frac {\theta (1 + \alpha)}{1 + \theta \lambda}} \iff \\ \end{array}
$$

$$
\begin{array}{l} y _ {i} ^ {\frac {1 + \theta \lambda + \theta}{1 + \theta \lambda}} l _ {i} ^ {\frac {- \theta (1 + \alpha + \theta \lambda (\alpha + \beta))}{1 + \theta \lambda}} = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {A} _ {i} ^ {\theta} \bar {u} _ {i} ^ {\theta} y _ {i} ^ {\frac {1 + \theta \lambda + \theta}{1 + \theta \lambda}} l _ {i} ^ {\frac {\theta (\beta - 1)}{1 + \theta \lambda}} + \\ \left(\frac {\bar {L} ^ {(\alpha + \beta) \theta}}{W ^ {\theta}}\right) ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j} \left(\bar {t} _ {i j} \bar {L} ^ {\lambda}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \times \bar {A} _ {i} ^ {\theta} \bar {u} _ {i} ^ {\theta} \frac {\theta \lambda}{1 + \theta \lambda} \bar {A} _ {j} ^ {- \frac {\theta}{1 + \theta \lambda}} y _ {j} ^ {\frac {1 + \theta}{1 + \theta \lambda}} l _ {j} ^ {- \frac {\theta (1 + \alpha)}{1 + \theta \lambda}} \iff \\ \end{array}
$$

$$
y _ {i} ^ {\frac {1 + \theta \lambda + \theta}{1 + \theta \lambda}} l _ {i} ^ {\frac {- \theta (1 + \alpha + \theta \lambda (\alpha + \beta))}{1 + \theta \lambda}} = \chi \bar {A} _ {i} ^ {\theta} \bar {u} _ {i} ^ {\theta} y _ {i} ^ {\frac {1 + \theta \lambda + \theta}{1 + \theta \lambda}} l _ {i} ^ {\frac {\theta (\beta - 1)}{1 + \theta \lambda}} + \chi^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j} \left(\bar {t} _ {i j} \bar {L} ^ {\lambda}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \bar {A} _ {i} ^ {\theta} \bar {u} _ {i} ^ {\theta} \frac {\frac {\theta \lambda}{1 + \theta \lambda}}{\frac {1}{1 + \theta \lambda}} \bar {A} _ {j} ^ {- \frac {\theta}{1 + \theta \lambda}} y _ {j} ^ {\frac {1 + \theta}{1 + \theta \lambda}} l _ {j} ^ {- \frac {\theta (1 + \alpha)}{1 + \theta \lambda}},
$$

where $\begin{array} { r } { \chi = \left( \frac { L ^ { \alpha + \beta } } { W } \right) ^ { \theta } } \end{array}$ , as in equation (28). For the second equilibrium condition, we proceed similarly: 

$$
\begin{array}{l} \bar {u} _ {i} ^ {- \theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (1 - \beta)} = \frac {\bar {L} ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {A} _ {i} ^ {\theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (\alpha + 1)} + \\ \sum_ {j} \left(\left(\bar {t} _ {j i} \bar {L} ^ {\lambda}\right) ^ {\frac {1}{1 + \theta \lambda}} \left(\frac {\bar {L} ^ {- (\alpha + \beta) \theta}}{W ^ {- \theta}}\right) ^ {\frac {\lambda}{1 + \theta \lambda}} \bar {A} _ {i} ^ {- \frac {\theta \lambda}{1 + \theta \lambda}} \bar {u} _ {j} ^ {- \frac {\theta \lambda}{1 + \theta \lambda}} l _ {j} ^ {- \frac {\theta \lambda (\beta - 1)}{1 + \theta \lambda}} l _ {i} ^ {- \frac {\theta \lambda (1 + \alpha)}{1 + \theta \lambda}} y _ {j} ^ {- \frac {\theta \lambda}{1 + \theta \lambda}} y _ {i} ^ {\frac {\lambda (1 + \theta)}{1 + \theta \lambda}}\right) ^ {- \theta} \times \\ \bar {u} _ {j} ^ {- \theta} y _ {j} ^ {- \theta} l _ {j} ^ {\theta (1 - \beta)} \Longleftrightarrow \\ \end{array}
$$

$$
\begin{array}{l} \bar {u} _ {i} ^ {- \theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (1 - \beta)} = \frac {\bar {L} ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {A} _ {i} ^ {\theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (\alpha + 1)} + \\ \sum_ {j} \left(\bar {t} _ {j i} \bar {L} ^ {\lambda}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \times \left(\frac {\bar {L} ^ {(\alpha + \beta) \theta}}{W ^ {\theta}}\right) ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \bar {A} _ {i} ^ {\theta \frac {\theta \lambda}{1 + \theta \lambda}} \bar {u} _ {j} ^ {- \frac {\theta}{1 + \theta \lambda}} l _ {j} ^ {\frac {\theta (1 - \beta)}{1 + \theta \lambda}} l _ {i} ^ {\theta (1 + \alpha) \frac {\theta \lambda}{1 + \theta \lambda}} y _ {j} ^ {- \frac {\theta}{1 + \theta \lambda}} y _ {i} ^ {- \frac {\theta \lambda (1 + \theta)}{1 + \theta \lambda}} \iff \\ \end{array}
$$

$$
\begin{array}{l} y _ {i} ^ {- \frac {\theta (1 - \lambda)}{1 + \theta \lambda}} l _ {i} ^ {\frac {\theta (1 - \beta - \theta \lambda (\alpha + \beta))}{1 + \theta \lambda}} = \frac {\bar {L} ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {A} _ {i} ^ {\theta} \bar {u} _ {i} ^ {\theta} y _ {i} ^ {- \frac {\theta (1 - \lambda)}{1 + \theta \lambda}} l _ {i} ^ {\frac {\theta (\alpha + 1)}{1 + \theta \lambda}} + \\ \left(\frac {\bar {L} ^ {(\alpha + \beta) \theta}}{W ^ {\theta}}\right) ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j} \left(\bar {t} _ {j i} \bar {L} ^ {\lambda}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \bar {A} _ {i} ^ {\theta} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \bar {u} _ {i} ^ {\theta} \bar {u} _ {j} ^ {- \frac {\theta}{1 + \theta \lambda}} y _ {j} ^ {- \frac {\theta}{1 + \theta \lambda}} l _ {j} ^ {\frac {\theta (1 - \beta)}{1 + \theta \lambda}} \iff \\ \end{array}
$$

$$
y_{i}^{-\frac{\theta(1 - \lambda)}{1 + \theta\lambda}}l_{i}^{\frac{\theta(1 - \beta - \theta\lambda(\alpha + \beta))}{1 + \theta\lambda}} = \chi \bar{A}_{i}^{\theta}\bar{u}_{i}^{\theta}y_{i}^{-\frac{\theta(1 - \lambda)}{1 + \theta\lambda}}l_{i}^{\frac{\theta(\alpha + 1)}{1 + \theta\lambda}}\\ +\chi^{\frac{\theta\lambda}{1 + \theta\lambda}}\sum_{j}\left(\bar{t}_{ji}\bar{L}^{\lambda}\right)^{-\frac{\theta}{1 + \theta\lambda}}\bar{A}_{i}^{\theta}\frac{\theta\lambda}{1 + \theta\lambda}\bar{u}_{i}^{\theta}\bar{u}_{j}^{-\frac{\theta}{1 + \theta\lambda}}y_{j}^{-\frac{\theta}{1 + \theta\lambda}}l_{j}^{\frac{\theta(1 - \beta)}{1 + \theta\lambda}}
$$

where again χ =  Lα+βW  $\begin{array} { r } { \chi = \left( \frac { L ^ { \alpha + \beta } } { W } \right) ^ { \theta } } \end{array}$ , as in equation (29). 

Commuting Model The derivations for the commuting model follow a very similar process to that of the economic geography model. We rewrite the first equilibrium condition, equation (19), as a matrix product and invert: 

$$
\left[ \bar {u} _ {i} ^ {- \theta} \left(l _ {i} ^ {R}\right) ^ {- \theta \beta + 1} \right] = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \times \left[ \tau_ {i j} ^ {- \theta} \right] \times \left[ \bar {A} _ {j} ^ {\theta} \left(l _ {j} ^ {F}\right) ^ {\theta \alpha} \right] \Longleftrightarrow
$$

$$
\left[ \bar {u} _ {i} ^ {- \theta} \left(l _ {i} ^ {R}\right) ^ {- \theta \beta + 1} \right] = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \times [ \mathbf {I} - \mathbf {A} ] ^ {- 1} \times \left[ \bar {A} _ {j} ^ {\theta} \left(l _ {j} ^ {F}\right) ^ {\theta \alpha} \right]
$$

$$
[ \mathbf {I} - \mathbf {A} ] \times \left[ \bar {u} _ {i} ^ {- \theta} \left(l _ {i} ^ {R}\right) ^ {- \theta \beta + 1} \right] = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \times \left[ \bar {A} _ {j} ^ {\theta} \left(l _ {j} ^ {F}\right) ^ {\theta \alpha} \right] \Longleftrightarrow
$$

$$
\left[ \bar {u} _ {i} ^ {- \theta} \left(l _ {i} ^ {R}\right) ^ {- \theta \beta + 1} \right] - \mathbf {A} \times \left[ \bar {u} _ {j} ^ {- \theta} \left(l _ {j} ^ {R}\right) ^ {- \theta \beta + 1} \right] = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \times \left[ \bar {A} _ {i} ^ {\theta} \left(l _ {i} ^ {F}\right) ^ {\theta \alpha} \right] \Longleftrightarrow
$$

$$
\bar {u} _ {i} ^ {- \theta} \left(l _ {i} ^ {R}\right) ^ {- \theta \beta + 1} - \sum_ {j} a _ {i j} \bar {u} _ {j} ^ {- \theta} \left(l _ {j} ^ {R}\right) ^ {- \theta \beta + 1} = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {A} _ {i} ^ {\theta} \left(l _ {i} ^ {F}\right) ^ {\theta \alpha}
$$

where $\left[ \bar { u } _ { i } ^ { \theta } \left( l _ { i } ^ { R } \right) ^ { - \theta \beta + 1 } \right]$ and $\left[ T _ { j } ^ { \theta } \left( l _ { j } ^ { F } \right) ^ { \theta \alpha } \right]$ are column vectors. 

Applying the same steps to equilibrium equation (20), where $\left[ \bar { A } _ { i } ^ { - \theta } \left( l _ { i } ^ { F } \right) ^ { - \theta \alpha + 1 } \right]$ and $\left[ \hat { u } _ { j } ^ { \theta } \left( l _ { j } ^ { R } \right) ^ { \theta \beta } \right]$ are row vectors: 

$$
\left(l _ {i} ^ {F}\right) ^ {- \theta \alpha + 1} = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \sum_ {j} \bar {u} _ {j} ^ {\theta} \tau_ {j i} ^ {- \theta} \bar {A} _ {i} ^ {\theta} \left(l _ {j} ^ {R}\right) ^ {\theta \beta} \Longleftrightarrow
$$

$$
\left[ \bar {A} _ {i} ^ {- \theta} \left(l _ {i} ^ {F}\right) ^ {- \theta \alpha + 1} \right] = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \times \left[ \bar {u} _ {j} ^ {\theta} \left(l _ {j} ^ {R}\right) ^ {\theta \beta} \right] \times \left[ \tau_ {j i} ^ {- \theta} \right] \Longleftrightarrow
$$

$$
\left[ \bar {A} _ {i} ^ {- \theta} \left(l _ {i} ^ {F}\right) ^ {- \theta \alpha + 1} \right] = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \times \left[ \bar {u} _ {j} ^ {\theta} \left(l _ {j} ^ {R}\right) ^ {\theta \beta} \right] \times \left[ \mathbf {I} - \mathbf {A} ^ {T} \right] ^ {- 1} \Longleftrightarrow
$$

$$
\left[ \bar {A} _ {i} ^ {- \theta} \left(l _ {i} ^ {F}\right) ^ {- \theta \alpha + 1} \right] \times \left[ \mathbf {I} - \mathbf {A} ^ {T} \right] = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \times \left[ \bar {u} _ {i} ^ {\theta} \left(l _ {i} ^ {R}\right) ^ {\theta \beta} \right] \Longleftrightarrow
$$

$$
\begin{array}{l} \left[ \bar {A} _ {i} ^ {- \theta} \left(l _ {i} ^ {F}\right) ^ {- \theta \alpha + 1} \right] - \left[ \bar {A} _ {j} ^ {- \theta} \left(l _ {j} ^ {F}\right) ^ {- \theta \alpha + 1} \right] \times \mathbf {A} ^ {T} = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \times \left[ \bar {u} _ {i} ^ {\theta} \left(l _ {i} ^ {R}\right) ^ {\theta \beta} \right] \Longleftrightarrow \\ T _ {i} ^ {- \theta} \left(l _ {i} ^ {F}\right) ^ {- \theta \alpha + 1} - \sum_ {j} a _ {j i} \bar {A} _ {j} ^ {- \theta} \left(l _ {j} ^ {F}\right) ^ {- \theta \alpha + 1} = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {u} _ {i} ^ {\theta} \left(l _ {i} ^ {R}\right) ^ {\theta \beta} \\ \end{array}
$$

Recalling $a _ { i j } \equiv t _ { i j } ^ { - \theta }$ ,we have two equilibrium conditions for our commuting model: 

$$
\bar {u} _ {i} ^ {- \theta} \left(l _ {i} ^ {R}\right) ^ {- \theta \beta + 1} = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {A} _ {i} ^ {\theta} \left(l _ {i} ^ {F}\right) ^ {\theta \alpha} + \sum_ {j} t _ {i j} ^ {- \theta} \bar {u} _ {j} ^ {- \theta} \left(l _ {j} ^ {R}\right) ^ {- \theta \beta + 1}
$$

$$
\bar {A} _ {i} ^ {- \theta} \left(l _ {i} ^ {F}\right) ^ {- \theta \alpha + 1} = \frac {L ^ {(\alpha + \beta) \theta}}{\bar {W} ^ {\theta}} \bar {u} _ {i} ^ {\theta} \left(l _ {i} ^ {R}\right) ^ {\theta \beta} + \sum_ {j} t _ {j i} ^ {- \theta} \bar {A} _ {j} ^ {- \theta} \left(l _ {j} ^ {F}\right) ^ {- \theta \alpha + 1}
$$

As above, substituting in our expression for the iceberg transportation costs along a link using equation (26), and converting from market access terms to equilibrium $\left\{ l _ { i } ^ { F } \right\}$ and $\left\{ l _ { i } ^ { R } \right\}$ (as in Appendix C.1), incorporates endogenous traffic congestion. For the first equilibrium condition (30), we have: 

$$
\begin{array}{l} \bar {u} _ {i} ^ {- \theta} \left(l _ {i} ^ {R}\right) ^ {- \theta \beta + 1} = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {A} _ {i} ^ {\theta} \left(l _ {i} ^ {F}\right) ^ {\theta \alpha} + \\ \sum_ {j} \left(\left(\bar {t} _ {i j} \bar {L} ^ {\lambda}\right) ^ {\frac {1}{1 + \theta \lambda}} \bar {A} _ {i} ^ {\frac {- \theta \lambda}{1 + \theta \lambda}} \left(l _ {i} ^ {F}\right) ^ {\frac {(1 - \alpha \theta) \lambda}{1 + \theta \lambda}} \bar {u} _ {j} ^ {\frac {- \theta \lambda}{1 + \theta \lambda}} \left(l _ {j} ^ {R}\right) ^ {\frac {(1 - \beta \theta) \lambda}{1 + \theta \lambda}} \times W ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \bar {L} ^ {\frac {- \theta \lambda (\alpha + \beta)}{1 + \theta \lambda}}\right) ^ {- \theta} \bar {u} _ {j} ^ {- \theta} \left(l _ {j} ^ {R}\right) ^ {- \theta \beta + 1} \Longleftrightarrow \\ \end{array}
$$

$$
\begin{array}{l} \bar {u} _ {i} ^ {- \theta} \left(l _ {i} ^ {R}\right) ^ {- \theta \beta + 1} = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {A} _ {i} ^ {\theta} \left(l _ {i} ^ {F}\right) ^ {\theta \alpha} + \\ \sum_ {j} \left(\bar {t} _ {i j} \bar {L} ^ {\lambda}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \bar {u} _ {j} ^ {\theta \frac {\theta \lambda}{1 + \theta \lambda} - \theta} \left(l _ {j} ^ {R}\right) ^ {- \frac {\theta \lambda (1 - \beta \theta)}{1 + \theta \lambda} + (1 - \beta \theta)} \bar {A} _ {i} ^ {\theta \frac {\theta \lambda}{1 + \theta \lambda}} \left(l _ {i} ^ {F}\right) ^ {- \frac {\theta \lambda (1 - \alpha \theta)}{1 + \theta \lambda}} W ^ {- \theta \frac {\theta \lambda}{1 + \theta \lambda}} \bar {L} ^ {\frac {\theta \lambda (\alpha \theta + \beta \theta)}{1 + \theta \lambda}} \iff \\ \end{array}
$$

$$
\begin{array}{l} \left(l _ {i} ^ {R}\right) ^ {- \theta \beta + 1} \left(l _ {i} ^ {F}\right) ^ {\frac {\theta \lambda (1 - \alpha \theta)}{1 + \theta \lambda}} = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {A} _ {i} ^ {\theta} \bar {u} _ {i} ^ {\theta} \left(l _ {i} ^ {F}\right) ^ {\theta \alpha + \frac {\theta \lambda (1 - \alpha \theta)}{1 + \theta \lambda}} + \\ \sum_ {j} \left(\bar {t} _ {i j} \bar {L} ^ {\lambda}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \bar {u} _ {j} ^ {\theta} \frac {\theta \lambda}{1 + \theta \lambda} ^ {- \theta} \bar {u} _ {i} ^ {\theta} \left(l _ {j} ^ {R}\right) ^ {- \frac {\theta \lambda (1 - \beta \theta)}{1 + \theta \lambda} + (1 - \beta \theta)} \bar {A} _ {i} ^ {\theta} \frac {\theta \lambda}{1 + \theta \lambda} W ^ {- \theta} \frac {\theta \lambda}{1 + \theta \lambda} \bar {L} ^ {\frac {\theta \lambda (\alpha \theta + \beta \theta)}{1 + \theta \lambda}} \Longleftrightarrow \\ \end{array}
$$

$$
\begin{array}{l} \left(l _ {i} ^ {R}\right) ^ {- \theta \beta + 1} \left(l _ {i} ^ {F}\right) ^ {\frac {\theta \lambda (1 - \alpha \theta)}{1 + \theta \lambda}} = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {A} _ {i} ^ {\theta} \bar {u} _ {i} ^ {\theta} \left(l _ {i} ^ {F}\right) ^ {\theta \alpha + \frac {\theta \lambda (1 - \alpha \theta)}{1 + \theta \lambda}} + \\ \sum_ {j} \left(\bar {t} _ {i j} \bar {L} ^ {\lambda}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \bar {u} _ {j} ^ {\theta} \frac {\theta \lambda}{1 + \theta \lambda} ^ {- \theta} \bar {u} _ {i} ^ {\theta} \left(l _ {j} ^ {R}\right) ^ {- \frac {\theta \lambda (1 - \beta \theta)}{1 + \theta \lambda} + (1 - \beta \theta)} \bar {A} _ {i} ^ {\theta} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} W ^ {- \theta} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \bar {L} ^ {\frac {\theta \lambda (\alpha \theta + \beta \theta)}{1 + \theta \lambda}} \Longleftrightarrow \\ \end{array}
$$

$$
\left(l _ {i} ^ {R}\right) ^ {- \theta \beta + 1} \left(l _ {i} ^ {F}\right) ^ {\frac {\theta \lambda (1 - \alpha \theta)}{1 + \theta \lambda}} = \chi \bar {A} _ {i} ^ {\theta} \bar {u} _ {i} ^ {\theta} \left(l _ {i} ^ {F}\right) ^ {\frac {\theta (\alpha + \lambda)}{1 + \theta \lambda}} + \chi^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j} \left(\bar {t} _ {i j} \bar {L} ^ {\lambda}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \bar {A} _ {i} ^ {\theta} \frac {\theta \lambda}{1 + \theta \lambda} \bar {u} _ {j} ^ {- \frac {\theta}{1 + \theta \lambda}} \bar {u} _ {i} ^ {\theta} \left(l _ {j} ^ {R}\right) ^ {\frac {1 - \beta \theta}{1 + \theta \lambda}}
$$

where $\begin{array} { r } { \chi = \left( \frac { L ^ { \alpha + \beta } } { W } \right) ^ { \theta } } \end{array}$  Lα+β  θ , as claimed. For the second equilibrium condition (31): 

$$
\begin{array}{l} \bar {A} _ {i} ^ {- \theta} \left(l _ {i} ^ {F}\right) ^ {- \theta \alpha + 1} = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {u} _ {i} ^ {\theta} \left(l _ {i} ^ {R}\right) ^ {\theta \beta} + \\ \sum_ {j} \left(\left(\bar {t} _ {j i} \bar {L} ^ {\lambda}\right) ^ {\frac {1}{1 + \theta \lambda}} \bar {u} _ {i} ^ {- \frac {\theta \lambda}{1 + \theta \lambda}} \left(l _ {i} ^ {R}\right) ^ {\frac {(1 - \beta \theta) \lambda}{1 + \theta \lambda}} \bar {A} _ {j} ^ {- \frac {\theta \lambda}{1 + \theta \lambda}} \left(l _ {j} ^ {F}\right) ^ {\frac {(1 - \alpha \theta) \lambda}{1 + \theta \lambda}} W ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \bar {L} ^ {- \frac {\theta \lambda (\alpha + \beta)}{1 + \theta \lambda}}\right) ^ {- \theta} \bar {A} _ {j} ^ {- \theta} \left(l _ {j} ^ {F}\right) ^ {- \theta \alpha + 1} \Longleftrightarrow \\ \end{array}
$$

$$
\begin{array}{l} \bar {A} _ {i} ^ {- \theta} \left(l _ {i} ^ {F}\right) ^ {- \theta \alpha + 1} = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {u} _ {i} ^ {\theta} \left(l _ {i} ^ {R}\right) ^ {\theta \beta} + \\ \sum_ {j} \left(\bar {t} _ {j i} \bar {L} ^ {\lambda}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \bar {u} _ {i} ^ {\theta} \frac {\theta \lambda}{1 + \theta \lambda} \left(l _ {i} ^ {R}\right) ^ {- \frac {\theta \lambda (1 - \beta \theta)}{1 + \theta \lambda}} \bar {A} _ {j} ^ {\theta} \frac {\theta \lambda}{1 + \theta \lambda} ^ {- \theta} \left(l _ {j} ^ {F}\right) ^ {- \frac {\theta \lambda (1 - \alpha \theta)}{1 + \theta \lambda} + (1 - \alpha \theta)} W ^ {- \theta} \frac {\theta \lambda}{1 + \theta \lambda} \bar {L} ^ {\frac {\theta \lambda (\alpha \theta + \beta \theta)}{1 + \theta \lambda}} \iff \\ \end{array}
$$

$$
\begin{array}{l} \left(l _ {i} ^ {F}\right) ^ {- \theta \alpha + 1} \left(l _ {i} ^ {R}\right) ^ {\frac {\theta \lambda (1 - \beta \theta)}{1 + \theta \lambda}} = \frac {L ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \bar {A} _ {i} ^ {\theta} \bar {u} _ {i} ^ {\theta} \left(l _ {i} ^ {R}\right) ^ {\theta \beta + \frac {\theta \lambda (1 - \beta \theta)}{1 + \theta \lambda}} + \\ \sum_ {j} \left(\bar {t} _ {j i} \bar {L} ^ {\lambda}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \bar {u} _ {i} ^ {\theta} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \bar {A} _ {i} ^ {\theta} \bar {A} _ {j} ^ {\theta} ^ {\frac {\theta \lambda}{1 + \theta \lambda} - \theta} \left(l _ {j} ^ {F}\right) ^ {- \frac {\theta \lambda (1 - \alpha \theta)}{1 + \theta \lambda} + (1 - \alpha \theta)} W ^ {- \theta} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \bar {L} ^ {\frac {\theta \lambda (\alpha \theta + \beta \theta)}{1 + \theta \lambda}} \Longleftrightarrow \\ \end{array}
$$

$$
\left(l _ {i} ^ {F}\right) ^ {- \theta \alpha + 1} \left(l _ {i} ^ {R}\right) ^ {\frac {\theta \lambda (1 - \beta \theta)}{1 + \theta \lambda}} = \chi \bar {A} _ {i} ^ {\theta} \bar {u} _ {i} ^ {\theta} \left(l _ {i} ^ {R}\right) ^ {\frac {\theta (\beta + \lambda)}{1 + \theta \lambda}} + \chi^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j} \left(\bar {t} _ {j i} \bar {L} ^ {\lambda}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \bar {A} _ {i} ^ {\theta} \bar {A} _ {j} ^ {- \frac {\theta}{1 + \theta \lambda}} \bar {u} _ {i} ^ {\theta} \frac {\frac {\theta \lambda}{1 + \theta \lambda}}{1 + \theta \lambda} \left(l _ {j} ^ {F}\right) ^ {\frac {1 - \alpha \theta}{1 + \theta \lambda}},
$$

as claimed. 

# Part

# Online Appendix

# Table of Contents

# B Additional Derivations 66

B.1 Section 2.1: An Economic Geography Model with Optimal Routing . . . . 66 

B.2 Section 2.2: An Urban Model with Optimal Routing . 67 

B.3 Section 4.3: Scale Dependence . . 68 

B.4 Section 5.1: Traffic, Trade, and Commuting flows . . 69 

# C Proofs 73

C.1 Proof of Proposition 1 74 

C.2 Proof of Proposition 2 86 

# D Extensions 93

D.1 Additive Transportation Costs 93 

D.2 Nested Route Choice . 94 

D.3 An Armington trade model 95 

D.4 Traffic congestion in quantities of labor 96 

D.5 Traffic congestion in quantities of goods 99 

# E Algorithm for Conducting Counterfactuals 100

# F Data Construction 101

F.1 The U.S. Highway System Network . . 102 

F.2 The Seattle Road Network . 104 

F.3 Trade, Commuting, and Distance . 107 

# G Alternative Parameter Constellations 107

# B Additional Derivations

This section presents the full derivations of the results mentioned in the text. It is organized sequentially by section of the main text. 

# B.1 Section 2.1: An Economic Geography Model with Optimal Routing

Below, we derive the equilibrium conditions for the economic geography model described in equations (10) and (11) in the paper. We start with the first market clearing condition defined in equation (5) and combine it with the gravity equation described in equation (6): 

$$
Y _ {i} = \sum_ {j = 1} ^ {N} X _ {i j} \iff
$$

$$
Y _ {i} = \sum_ {j = 1} ^ {N} \tau_ {i j} ^ {- \theta} w _ {i} ^ {- \theta} A _ {i} ^ {\theta} E _ {j} P _ {j} ^ {\theta} \Longleftrightarrow
$$

$$
\frac {Y _ {i}}{A _ {i} ^ {\theta}} w _ {i} ^ {\theta} = \sum_ {j = 1} ^ {N} \tau_ {i j} ^ {- \theta} E _ {j} P _ {j} ^ {\theta} \iff
$$

$$
\bar {A} _ {i} ^ {- \theta} L _ {i} ^ {1 - \alpha \theta} w _ {i} ^ {\theta + 1} = \sum_ {j = 1} ^ {N} \tau_ {i j} ^ {- \theta} w _ {j} L _ {j} P _ {j} ^ {\theta}.
$$

Assuming welfare equalization, the above becomes: 

$$
\bar {A} _ {i} ^ {- \theta} L _ {i} ^ {1 - \alpha \theta} w _ {i} ^ {\theta + 1} = \sum_ {j = 1} ^ {N} \tau_ {i j} ^ {- \theta} w _ {j} L _ {j} w _ {j} ^ {\theta} \bar {u} _ {j} ^ {\theta} L _ {j} ^ {\beta \theta} W ^ {- \theta} \Longleftrightarrow
$$

$$
\bar {A} _ {i} ^ {- \theta} L _ {i} ^ {1 - \alpha \theta} w _ {i} ^ {\theta + 1} = \sum_ {j = 1} ^ {N} \tau_ {i j} ^ {- \theta} w _ {j} ^ {\theta + 1} \bar {u} _ {j} ^ {\theta} L _ {j} ^ {\beta \theta + 1} W ^ {- \theta}.
$$

Now, defining $\begin{array} { r } { y _ { i } = \frac { Y _ { i } } { Y ^ { W } } = \frac { w _ { i } L _ { i } } { Y ^ { W } } } \end{array}$ and $\begin{array} { r } { l _ { i } = \frac { L _ { i } } { \bar { L } } } \end{array}$ 

$$
\bar {A} _ {i} ^ {- \theta} L _ {i} ^ {1 - \alpha \theta} w _ {i} ^ {\theta + 1} = \sum_ {j = 1} ^ {N} \tau_ {i j} ^ {- \theta} w _ {j} ^ {\theta + 1} \bar {u} _ {j} ^ {\theta} L _ {j} ^ {\beta \theta + 1} W ^ {- \theta} \iff
$$

$$
\bar {A} _ {i} ^ {- \theta} l _ {i} ^ {1 - \alpha \theta} \bar {L} ^ {1 - \alpha \theta} \left(\frac {y _ {i} Y ^ {W}}{l _ {i} \bar {L}}\right) ^ {\theta + 1} = \sum_ {j = 1} ^ {N} \tau_ {i j} ^ {- \theta} \left(\frac {y _ {i} Y ^ {W}}{l _ {i} \bar {L}}\right) ^ {\theta + 1} \bar {u} _ {j} ^ {\theta} L _ {j} ^ {\beta \theta + 1} W ^ {- \theta} \iff
$$

$$
\bar {A} _ {i} ^ {- \theta} y _ {i} ^ {\theta + 1} l _ {i} ^ {- \theta (1 + \alpha)} \bar {L} ^ {\theta (1 + \alpha)} \left(Y ^ {W}\right) ^ {\theta + 1} = \left(Y ^ {W}\right) ^ {\theta + 1} \bar {L} ^ {\theta (\beta - 1)} W ^ {- \theta} \sum_ {j = 1} ^ {N} \tau_ {i j} ^ {- \theta} \bar {u} _ {j} ^ {\theta} y _ {j} ^ {1 + \theta} l _ {j} ^ {\theta (\beta - 1)} \iff
$$

$$
\bar {A} _ {i} ^ {- \theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {- \theta (1 + \alpha)} = \frac {\bar {L} ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \sum_ {j = 1} ^ {N} \tau_ {i j} ^ {- \theta} \bar {u} _ {j} ^ {\theta} y _ {j} ^ {1 + \theta} l _ {j} ^ {\theta (\beta - 1)}.
$$

This concludes our derivation of the equilibrium condition we describe in equation (10). 

Moving on to the second market clearing condition defined in equation (5) and combining with the gravity equation described in equation (6): 

$$
E _ {i} = \sum_ {j = 1} ^ {N} X _ {j i} \iff
$$

$$
E _ {i} = \sum_ {j = 1} ^ {N} \tau_ {j i} ^ {- \theta} w _ {j} ^ {- \theta} A _ {j} ^ {\theta} E _ {i} P _ {i} ^ {\theta} \iff
$$

$$
P _ {i} ^ {- \theta} = \sum_ {j = 1} ^ {N} \tau_ {j i} ^ {- \theta} w _ {j} ^ {- \theta} \bar {A} _ {j} ^ {\theta} L _ {j} ^ {\alpha \theta}.
$$

Assuming welfare equalization, this definition of the consumer price index becomes: 

$$
W ^ {\theta} w _ {i} ^ {- \theta} \bar {u} _ {i} ^ {- \theta} L _ {i} ^ {- \beta \theta} = \sum_ {j = 1} ^ {N} \tau_ {j i} ^ {- \theta} w _ {j} ^ {- \theta} \bar {A} _ {j} ^ {\theta} L _ {j} ^ {\alpha \theta}.
$$

Defining the same income and labor shares, $\begin{array} { r } { y _ { i } = \frac { Y _ { i } } { Y ^ { W } } = \frac { w _ { i } L _ { i } } { Y ^ { W } } } \end{array}$ and $\begin{array} { r } { l _ { i } = \frac { L _ { i } } { L } } \end{array}$ , we used for the first equilibrium condition, we get the following: 

$$
W ^ {\theta} w _ {i} ^ {- \theta} \bar {u} _ {i} ^ {- \theta} L _ {i} ^ {- \beta \theta} = \sum_ {j = 1} ^ {N} \tau_ {j i} ^ {- \theta} w _ {j} ^ {- \theta} \bar {A} _ {j} ^ {\theta} L _ {j} ^ {\alpha \theta} \iff
$$

$$
W ^ {\theta} \left(\frac {y _ {i} Y ^ {W}}{l _ {i} \bar {L}}\right) ^ {- \theta} \bar {u} _ {i} ^ {- \theta} l _ {i} ^ {- \beta \theta} \bar {L} ^ {- \beta \theta} = \sum_ {j = 1} ^ {N} \tau_ {j i} ^ {- \theta} \left(\frac {y _ {j} Y ^ {W}}{l _ {j} \bar {L}}\right) ^ {- \theta} \bar {A} _ {j} ^ {\theta} l _ {j} ^ {\alpha \theta} \bar {L} ^ {\alpha \theta} \iff
$$

$$
W ^ {\theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (1 - \beta)} \bar {u} _ {i} ^ {- \theta} \left(Y ^ {W}\right) ^ {- \theta} \bar {L} ^ {- \beta \theta} = \left(Y ^ {W}\right) ^ {- \theta} \bar {L} ^ {\theta (\alpha + 1)} \sum_ {j = 1} ^ {N} \tau_ {j i} ^ {- \theta} y _ {j} ^ {- \theta} l _ {j} ^ {\theta (\alpha + 1)} \bar {A} _ {j} ^ {\theta} \iff
$$

$$
\bar {u} _ {i} ^ {- \theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (1 - \beta)} = \frac {\bar {L} ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \sum_ {j = 1} ^ {N} \tau_ {j i} ^ {- \theta} \bar {A} _ {j} ^ {\theta} y _ {j} ^ {- \theta} l _ {j} ^ {\theta (\alpha + 1)}.
$$

At this point, we have reached the second equilibrium condition described in equation (11). 

# B.2 Section 2.2: An Urban Model with Optimal Routing

In this section, we derive the equilibrium conditions for the urban model described in equations (19) and (20). We start by combining commuting gravity with the adding-up constraint on the residential population: 

$$
\begin{array}{l} L _ {i} ^ {R} = \sum_ {j} L _ {i j} \iff \\ = \sum_ {j} \tau_ {i j} ^ {- \theta} u _ {i} ^ {\theta} A _ {j} ^ {\theta} \left(\frac {\bar {L}}{W ^ {\theta}}\right). \\ \end{array}
$$

We substitute in for the spillovers as defined in equation (14): 

$$
L _ {i} ^ {R} = \sum_ {j} \tau_ {i j} ^ {- \theta} \bar {u} _ {i} ^ {\theta} \left(L _ {i} ^ {R}\right) ^ {\beta \theta} \bar {A} _ {j} ^ {\theta} \left(L _ {j} ^ {F}\right) ^ {\alpha \theta} \left(\frac {\bar {L}}{W ^ {\theta}}\right).
$$

Next, we define residential labor shares and commercial labor shares as $l _ { i } ^ { R } = L _ { i } ^ { R } / L , l _ { i } ^ { F } = L _ { i } ^ { F } / L$ , and putting the above equation in terms of the shares, we get equilibrium equation (19): 

$$
L _ {i} ^ {R} = \sum_ {j} \tau_ {i j} ^ {- \theta} \bar {u} _ {i} ^ {\theta} \left(L _ {i} ^ {R}\right) ^ {\beta \theta} \bar {A} _ {j} ^ {\theta} \left(L _ {j} ^ {F}\right) ^ {\alpha \theta} \left(\frac {\bar {L}}{W ^ {\theta}}\right) \Longleftrightarrow
$$

$$
l _ {i} ^ {R} \bar {L} = \sum_ {j} \tau_ {i j} ^ {- \theta} \bar {u} _ {i} ^ {\theta} \left(l _ {i} ^ {R} \bar {L}\right) ^ {\beta \theta} \bar {A} _ {j} ^ {\theta} \left(l _ {j} ^ {F} \bar {L}\right) ^ {\alpha \theta} \left(\frac {\bar {L}}{W ^ {\theta}}\right) \Longleftrightarrow
$$

$$
\left(l _ {i} ^ {R}\right) ^ {1 - \beta \theta} = \frac {\bar {L} ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \sum_ {j} \tau_ {i j} ^ {- \theta} \bar {u} _ {i} ^ {\theta} \bar {A} _ {j} ^ {\theta} \left(l _ {j} ^ {F}\right) ^ {\alpha \theta}.
$$

Moving to the derivation of the second equilibrium condition, we start by combining commuting gravity with the adding up constraint on commercial population: 

$$
\begin{array}{l} L _ {i} ^ {F} = \sum_ {j} L _ {j i} \iff \\ = \sum_ {j} \tau_ {j i} ^ {- \theta} u _ {j} ^ {\theta} A _ {i} ^ {\theta} \left(\frac {\bar {L}}{W ^ {\theta}}\right). \\ \end{array}
$$

We substitute in for the spillovers to arrive at the following characterization of the commercial labor force: 

$$
L _ {i} ^ {F} = \sum_ {j} \tau_ {j i} ^ {- \theta} \bar {u} _ {j} ^ {\theta} \left(L _ {j} ^ {R}\right) ^ {\beta \theta} \bar {A} _ {i} ^ {\theta} \left(L _ {j} ^ {F}\right) ^ {\alpha \theta} \left(\frac {\bar {L}}{W ^ {\theta}}\right).
$$

Finally, we define the above expression in terms of residential and commercial labor shares: 

$$
L _ {i} ^ {F} = \sum_ {j} \tau_ {j i} ^ {- \theta} \bar {u} _ {j} ^ {\theta} \left(L _ {j} ^ {R}\right) ^ {\beta \theta} \bar {A} _ {i} ^ {\theta} \left(L _ {j} ^ {F}\right) ^ {\alpha \theta} \left(\frac {\bar {L}}{W ^ {\theta}}\right) \Longleftrightarrow
$$

$$
l _ {i} ^ {F} \bar {L} = \sum_ {j} \tau_ {j i} ^ {- \theta} \bar {u} _ {j} ^ {\theta} \left(l _ {j} ^ {R} \bar {L}\right) ^ {\beta \theta} \bar {A} _ {i} ^ {\theta} \left(l _ {j} ^ {F} \bar {L}\right) ^ {\alpha \theta} \left(\frac {\bar {L}}{W ^ {\theta}}\right) \Longleftrightarrow
$$

$$
\left(l _ {i} ^ {F}\right) ^ {1 - \alpha \theta} = \frac {\bar {L} ^ {(\alpha + \beta) \theta}}{W ^ {\theta}} \sum_ {j} \tau_ {j i} ^ {- \theta} \bar {A} _ {i} ^ {\theta} \bar {u} _ {j} ^ {\theta} \left(l _ {j} ^ {R}\right) ^ {\beta \theta},
$$

which is the second equilibrium condition defined in equation (20). 

# B.3 Section 4.3: Scale Dependence

In this section, we present a derivation of the partial derivative of trade costs about $c = 1$ . Let’s define the matrix of bilateral trade costs as ${ \bf B } \equiv ( { \bf I } - ( \exp \left( - \theta \ln c \right) { \bf A } ) ) ^ { - 1 }$ . We then have from, by matrix calculus, that: 

$$
\frac {\partial \mathbf {B}}{\partial \ln c} = - (\mathbf {B}) \frac {\partial (\mathbf {I} - (\exp (- \theta \ln c) \mathbf {A}))}{\partial \ln c} (\mathbf {B}) \Longleftrightarrow
$$

$$
\frac {\partial \mathbf {B}}{\partial \ln c} = - \theta (\mathbf {B}) c ^ {- \theta} \mathbf {A} (\mathbf {B})
$$

so that: 

$$
\left[ \frac {\partial \mathbf {B}}{\partial \ln c} \right] _ {i j} = - \theta \sum_ {k} \sum_ {l} \mathbf {B} _ {i k} \bar {t} _ {k l} ^ {- \theta} \mathbf {B} _ {l j}
$$

or in our notation: 

$$
\frac {\partial \tau_ {i j} ^ {- \theta} (c)}{\partial \ln c} | _ {c = 1} = - \theta \sum_ {k = 1} ^ {N} \sum_ {l = 1} ^ {N} \tau_ {i k} ^ {- \theta} \bar {t} _ {k l} ^ {- \theta} \tau_ {l j} ^ {- \theta} \iff
$$

$$
\frac {\partial \ln \tau_ {i j} (c)}{\partial \ln c} | _ {c = 1} = \sum_ {k = 1} ^ {N} \sum_ {l = 1} ^ {N} \frac {\tau_ {i k} ^ {- \theta} \bar {t} _ {k l} ^ {- \theta} \tau_ {l j} ^ {- \theta}}{\tau_ {i j} ^ {- \theta}}
$$

Recall: 

$$
\pi_ {i j} ^ {k l} = \left(\frac {\tau_ {i j}}{\tau_ {i k} t _ {k l} \tau_ {l j}}\right) ^ {\theta}
$$

so that: 

$$
\frac {\partial \ln \tau_ {i j} (c)}{\partial \ln c} | _ {c = 1} = \sum_ {k = 1} ^ {N} \sum_ {l = 1} ^ {N} \pi_ {i j} ^ {k l}
$$

# B.4 Section 5.1: Traffic, Trade, and Commuting flows

In this section, we derive an analytical mapping between traffic and gravity flows of trade and commuting. We begin by writing trade, commuting and traffic gravity from equations (3), (13), and (24), respectively in matrix form: 

$$
\mathbf {X} = \left(\frac {Y}{\Pi}\right) (\mathbf {I} - \mathbf {A}) ^ {- 1} \left(\frac {E}{P}\right)
$$

$$
\mathbf {L} = \left(\frac {L ^ {R}}{\Pi}\right) (\mathbf {I} - \mathbf {A}) ^ {- 1} \left(\frac {L ^ {F}}{P}\right)
$$

$$
\boldsymbol {\Xi} = \mathbf {P A I I},
$$

where for the gravity equations we used equation (21) and where 

$$
\mathbf {X} = [ X _ {i j} ]
$$

$$
\Xi = \left[ \Xi_ {i j} \right]
$$

$$
\mathbf {P} = \operatorname {d i a g} \left(P _ {i} ^ {- \theta}\right)
$$

$$
Y = \operatorname {d i a g} \left(Y _ {i}\right)
$$

$$
L ^ {R} = \operatorname {d i a g} \left(L _ {i} ^ {R}\right)
$$

$$
\frac {Y}{\Pi} = \operatorname {d i a g} \left(\frac {Y _ {i}}{\Pi_ {i} ^ {- \theta}}\right)
$$

$$
\frac {L ^ {R}}{\Pi} = \operatorname {d i a g} \left(\frac {L _ {i} ^ {R}}{\Pi_ {i} ^ {- \theta}}\right)
$$

$$
\mathbf {L} = \left[ L _ {i j} \right]
$$

$$
\mathbf {A} = \left[ t _ {i j} ^ {- \theta} \right]
$$

$$
\boldsymbol {\Pi} = \operatorname {d i a g} \left(\Pi_ {i} ^ {- \theta}\right)
$$

$$
E = \operatorname {d i a g} \left(E _ {i}\right)
$$

$$
L ^ {F} = \operatorname {d i a g} \left(L _ {i} ^ {F}\right)
$$

$$
\frac {E}{P} = \operatorname {d i a g} \left(\frac {E _ {i}}{P _ {i} ^ {- \theta}}\right)
$$

$$
\frac {L ^ {F}}{P} = \operatorname {d i a g} \left(\frac {L _ {i} ^ {F}}{P _ {i} ^ {- \theta}}\right)
$$

are each $N \times N$ matrices. 

Solving for the adjacency matrix A, we have: 

$$
\Xi = \mathrm {P A} \Pi \Longleftrightarrow
$$

$$
\mathbf {A} = \mathbf {P} ^ {- 1} \boldsymbol {\Xi} \boldsymbol {\Pi} ^ {- 1}
$$

which we can substitute into our trade gravity equation: 

$$
\mathbf {X} = \left(\frac {Y}{\Pi}\right) (\mathbf {I} - \mathbf {A}) ^ {- 1} \left(\frac {E}{\mathbf {P}}\right) \Longleftrightarrow
$$

$$
\mathbf {X} = \left(\frac {Y}{\Pi}\right) \left(\mathbf {I} - \mathbf {P} ^ {- 1} \boldsymbol {\Xi} \boldsymbol {\Pi} ^ {- 1}\right) ^ {- 1} \left(\frac {E}{\mathbf {P}}\right) \Longleftrightarrow
$$

$$
\mathbf {X} ^ {- 1} = \left(\frac {E}{\mathbf {P}}\right) ^ {- 1} \left(\mathbf {I} - \mathbf {P} ^ {- 1} \boldsymbol {\Xi} \boldsymbol {\Pi} ^ {- 1}\right) \left(\frac {Y}{\Pi}\right) ^ {- 1} \Longleftrightarrow
$$

$$
\mathbf {X} ^ {- 1} = \left(\frac {E}{\mathbf {P}}\right) ^ {- 1} \left(\frac {Y}{\Pi}\right) ^ {- 1} - (E) ^ {- 1} \Xi (Y) ^ {- 1} \Longleftrightarrow
$$

$$
\mathbf {X} ^ {- 1} = (E) ^ {- 1} \left(\mathbf {P} \boldsymbol {\Pi} - \boldsymbol {\Xi}\right) (Y) ^ {- 1} \Longleftrightarrow
$$

$$
\mathbf {X} = (Y) (\mathbf {P} \boldsymbol {\Pi} - \boldsymbol {\Xi}) ^ {- 1} (E)
$$

and our commuting gravity equation: 

$$
\mathbf {L} = \left(\frac {L ^ {R}}{\Pi}\right) (\mathbf {I} - \mathbf {A}) ^ {- 1} \left(\frac {L ^ {F}}{P}\right) \Longleftrightarrow
$$

$$
\mathbf {L} = \left(\frac {L ^ {R}}{\Pi}\right) \left(\mathbf {I} - \mathbf {P} ^ {- 1} \boldsymbol {\Xi} \boldsymbol {\Pi} ^ {- 1}\right) ^ {- 1} \left(\frac {L ^ {F}}{P}\right) \Longleftrightarrow
$$

$$
\mathbf {L} ^ {- 1} = \left(\frac {L ^ {F}}{P}\right) ^ {- 1} \left(\mathbf {I} - \mathbf {P} ^ {- 1} \boldsymbol {\Xi} \boldsymbol {\Pi} ^ {- 1}\right) \left(\frac {L ^ {R}}{\Pi}\right) ^ {- 1} \Longleftrightarrow
$$

$$
\mathbf {L} ^ {- 1} = \left(\left(\frac {L ^ {F}}{P}\right) ^ {- 1} \left(\frac {L ^ {R}}{\Pi}\right) ^ {- 1} - \left(L ^ {f}\right) ^ {- 1} \boldsymbol {\Xi} \left(L ^ {R}\right) ^ {- 1}\right) \Longleftrightarrow
$$

$$
\mathbf {L} ^ {- 1} = \left(L ^ {F}\right) ^ {- 1} \left(\mathbf {P I I} - \boldsymbol {\Xi}\right) \left(L ^ {R}\right) ^ {- 1} \Longleftrightarrow
$$

$$
\mathbf {L} = \left(L ^ {R}\right) \left(\mathbf {P} \boldsymbol {\Pi} - \boldsymbol {\Xi}\right) ^ {- 1} \left(L ^ {F}\right)
$$

Now all that remains is to define diagonal matrix $\mathbf { P I I }$ in terms of traffic $\Xi$ and other observables. For the trade model, we have the following, where $P$ and $\frac { Y } { \Pi }$ are column vectors: 

$$
P _ {i} ^ {- \theta} = \sum_ {j} \tau_ {j i} ^ {- \theta} \frac {Y _ {j}}{\Pi_ {j} ^ {- \theta}} \Longleftrightarrow
$$

$$
P = \left(\mathbf {I} - \mathbf {A} ^ {T}\right) ^ {- 1} \left( \begin{array}{c} Y \\ \overline {{\Pi}} \end{array} \right) \Longleftrightarrow
$$

$$
\left(\mathbf {I} - \Pi^ {- 1} \Xi^ {T} P ^ {- 1}\right) P = \frac {Y}{\Pi} \iff
$$

$$
P - \Pi^ {- 1} \Xi^ {T} 1 = \frac {Y}{\Pi}
$$

$$
\Pi P = Y + \Xi^ {T} 1
$$

and: 

$$
\Pi_ {i} ^ {- \theta} = \sum_ {j} \tau_ {i j} ^ {- \theta} \frac {E _ {j}}{P _ {j} ^ {- \theta}} \Longleftrightarrow
$$

$$
\Pi = (\mathbf {I} - \mathbf {A}) ^ {- 1} \left(\frac {E}{P}\right) \Longleftrightarrow
$$

$$
\left(\mathbf {I} - P ^ {- 1} \Xi \Pi^ {- 1}\right) \Pi = \left( \begin{array}{c} E \\ \overline {{P}} \end{array} \right)
$$

$$
\Pi - P ^ {- 1} \Xi \mathbf {1} = \frac {E}{P}
$$

$$
\Pi P = E + \Xi \mathbf {1}
$$

where $\mathbf { I }$ is the $N \times N$ identity matrix, 1 is an $N \times 1$ column vector of ones, and $\mathrm { I I }$ and $\textstyle { \frac { E } { P } }$ are column vectors. 

Since we have two definitions of column vector $1 1 P$ , we average them: 

$$
\Pi P = \frac {1}{2} (E + \Xi 1) + \frac {1}{2} (Y + \Xi^ {T} 1) = \frac {1}{2} (E + Y) + \frac {1}{2} (\Xi 1 + \Xi^ {T} 1)
$$

and plug that definition into our matrix product characterization of trade flows, where the diagonal matrix $\mathbf { P } \mathbf { I I } = \mathrm { d i a g } ( \Pi P )$ : 

$$
\mathbf {X} = (Y) \left(\operatorname {d i a g} \left(\frac {1}{2} (E + Y) + \frac {1}{2} (\Xi 1 + \Xi^ {T} 1)\right) - \Xi\right) ^ {- 1} (E) \Longleftrightarrow
$$

$$
X _ {i j} = \left[ \mathbf {D} ^ {X} - \boldsymbol {\Xi} \right] _ {i j} ^ {- 1} \times Y _ {i} \times E _ {j},
$$

where ${ \bf D } ^ { X } \equiv \mathrm { d i a g } \left( \frac { _ 1 } { ^ 2 } \left( E + Y \right) + \frac { _ 1 } { ^ 2 } \left( \Xi 1 + \Xi ^ { T } 1 \right) \right)$ , as in equation (34). 

For the commuting model, we derive the following for vector $P \mathrm { I I }$ 

$$
P _ {i} ^ {- \theta} = \sum_ {j} \tau_ {j i} ^ {- \theta} \frac {L _ {j} ^ {R}}{\Pi_ {j} ^ {- \theta}} \iff
$$

$$
P = \left(I - A ^ {T}\right) ^ {- 1} \left(\frac {L ^ {R}}{\Pi}\right) \Longleftrightarrow
$$

$$
\left(I - \Pi^ {- 1} \Xi^ {T} P ^ {- 1}\right) P = \frac {L ^ {R}}{\Pi} \Longleftrightarrow
$$

$$
P - \Pi^ {- 1} \Xi^ {T} 1 = \frac {L ^ {R}}{\Pi} \Longleftrightarrow
$$

$$
\Pi P = \frac {L ^ {R}}{\Pi} + \Xi^ {T} 1
$$

where $P$ and L R $\frac { L ^ { R } } { \Pi }$ are column vectors. We also derive another definition: 

$$
\Pi_ {i} ^ {- \theta} = \sum_ {j} \tau_ {i j} ^ {- \theta} \frac {L _ {j} ^ {F}}{P _ {j} ^ {- \theta}} \iff
$$

$$
\Pi = (I - A) ^ {- 1} \left(\frac {L ^ {F}}{P}\right) \Longleftrightarrow
$$

$$
\left(I - P ^ {- 1} \Xi \Pi^ {- 1}\right) \Pi = \left(\frac {L ^ {F}}{P}\right)
$$

$$
\Pi - P ^ {- 1} \Xi 1 = \frac {L ^ {F}}{P}
$$

$$
\Pi P = L ^ {F} + \Xi 1
$$

where $\mathrm { I I }$ and L F $\frac { L ^ { F } } { P }$ are column vectors. 

Like with the trade case, we average over the two different definitions for traffic to define the vector $1 1 P$ 

$$
\Pi P = \frac {1}{2} \left(L ^ {F} + \Xi 1\right) + \frac {1}{2} \left(L ^ {R} + \Xi^ {T} 1\right) = \frac {1}{2} \left(L ^ {F} + L ^ {R}\right) + \frac {1}{2} \left(\Xi 1 + \Xi^ {T} 1\right)
$$

and combining this average definition with our matrix product characterization of commuting flows, with the equality $\mathbf { P } \mathbf { I I } = \mathrm { d i a g } ( \Pi P ) \mathrm { w e }$ get: 

$$
\mathbf {L} = \left(L ^ {R}\right) \left(\operatorname {d i a g} \left(\frac {1}{2} \left(L ^ {F} + L ^ {R}\right) + \frac {1}{2} \left(\Xi 1 + \Xi^ {T} 1\right)\right) - \Xi\right) ^ {- 1} \left(L ^ {F}\right)
$$

$$
L _ {i j} = \left[ \mathbf {D} ^ {L} - \Xi \right] _ {i j} ^ {- 1} \times L _ {i} ^ {R} \times L _ {j} ^ {F},
$$

as in equation (35). 

# Recovering the level of trade costs from observed traffic flows

As mentioned in footnote 23, we can use a similar methodology to recover the observed level of (symmetric) travel costs on each link $\left\{ t _ { i j } ^ { - \theta } \right\}$ from observed traffic flows (and measures of economic activity in each location). 

We being by noting that in both the economic geography and urban variants of our model, we can write the equilibrium inward and outward market access in each location as satisfying the following two equations: 

$$
\Pi_ {i} ^ {- \theta} = \sum_ {j = 1} ^ {N} \tau_ {i j} ^ {- \theta} X _ {j} ^ {i n} P _ {j} ^ {\theta} \tag {B.1}
$$

$$
P _ {j} ^ {- \theta} = \sum_ {i = 1} ^ {N} \tau_ {i j} ^ {- \theta} X _ {i} ^ {\text {o u t}} \Pi_ {i} ^ {\theta} \tag {B.2}
$$

where $X _ { j } ^ { i n } \equiv E _ { j }$ in an economic geography model or $X _ { j } ^ { i n } \equiv L _ { j } ^ { F }$ in an urban model an, similarly, $X _ { i } ^ { o u t } \equiv Y _ { i }$ in a economic geography model or $X _ { i } ^ { o u t } \equiv L _ { i } ^ { R }$ in an urban model. 

Using equation (21) that expresses the equilibrium bilateral transportation costs as a function of the transportation network, we can rewrite these market access equilibrium conditions as: 

$$
\Pi_ {i} ^ {- \theta} = X _ {i} ^ {\text {i n}} P _ {i} ^ {\theta} + \sum_ {j = 1} ^ {N} t _ {i j} ^ {- \theta} \Pi_ {j} ^ {- \theta} \tag {B.3}
$$

$$
P _ {i} ^ {- \theta} = X _ {i} ^ {\text {o u t}} \Pi_ {i} ^ {\theta} + \sum_ {j = 1} ^ {N} t _ {j i} ^ {- \theta} P _ {j} ^ {- \theta} \tag {B.4}
$$

Then from the traffic equation (24), we can express the product of the equilibrium market access terms solely as a function of observed traffic flows – i.e. $\{ \Xi _ { i j } \}$ – and the measures of economic activity in each location, i.e. $\{ X _ { i } ^ { i n } , X _ { i } ^ { o u t } \}$ : 

$$
\Pi_ {i} ^ {- \theta} = X _ {i} ^ {i n} P _ {i} ^ {\theta} + \sum_ {j = 1} ^ {N} t _ {i j} ^ {- \theta} \Pi_ {j} ^ {- \theta} \iff
$$

$$
\Pi_ {i} ^ {- \theta} P _ {i} ^ {- \theta} = X _ {i} ^ {i n} + \sum_ {j = 1} ^ {N} \Xi_ {i j} \tag {B.5}
$$

and: 

$$
P _ {i} ^ {- \theta} = X _ {i} ^ {o u t} \Pi_ {i} ^ {\theta} + \sum_ {j = 1} ^ {N} t _ {j i} ^ {- \theta} P _ {j} ^ {- \theta} \iff
$$

$$
P _ {i} ^ {- \theta} \Pi_ {i} ^ {- \theta} = X _ {i} ^ {\text {o u t}} + \sum_ {j = 1} ^ {N} \Xi_ {j i}, \tag {B.6}
$$

or, by combining the two expressions: 

$$
P _ {i} ^ {- \theta} \Pi_ {i} ^ {- \theta} = \frac {1}{2} \left(X _ {i} ^ {\text {o u t}} + X _ {i} ^ {\text {i n}}\right) + \frac {1}{2} \left(\sum_ {j = 1} ^ {N} \Xi_ {i j} + \sum_ {j = 1} ^ {N} \Xi_ {j i}\right) \tag {B.7}
$$

Substituting this expression into traffic equation (24) and multiplying the bilateral link costs in the two directions results in: 

$$
t _ {i j} ^ {- \theta} \times t _ {j i} ^ {- \theta} = \frac {\Xi_ {i j}}{P _ {i} ^ {- \theta} \times \Pi_ {j} ^ {- \theta}} \times \frac {\Xi_ {j i}}{P _ {j} ^ {- \theta} \times \Pi_ {i} ^ {- \theta}} \Longleftrightarrow
$$

$$
t _ {i j} ^ {- \theta} \times t _ {j i} ^ {- \theta} = \frac {\Xi_ {i j} \times \Xi_ {j i}}{\left(\frac {1}{2} \left(X _ {i} ^ {o u t} + X _ {i} ^ {i n}\right) + \frac {1}{2} \left(\sum_ {k = 1} ^ {N} \Xi_ {i k} + \sum_ {k = 1} ^ {N} \Xi_ {k i}\right)\right) \times \left(\frac {1}{2} \left(X _ {j} ^ {o u t} + X _ {j} ^ {i n}\right) + \frac {1}{2} \left(\sum_ {k = 1} ^ {N} \Xi_ {j k} + \sum_ {k = 1} ^ {N} \Xi_ {k j}\right)\right)}
$$

Finally, imposing symmetry – i.e. $t _ { i j } = t _ { j i }$ – immediately implies: 

$$
t _ {i j} ^ {- \theta} = \sqrt {\frac {\Xi_ {i j} \times \Xi_ {j i}}{\left(\frac {1}{2} \left(X _ {i} ^ {o u t} + X _ {i} ^ {i n}\right) + \frac {1}{2} \left(\sum_ {k = 1} ^ {N} \Xi_ {i k} + \sum_ {k = 1} ^ {N} \Xi_ {k i}\right)\right) \times \left(\frac {1}{2} \left(X _ {j} ^ {o u t} + X _ {j} ^ {i n}\right) + \frac {1}{2} \left(\sum_ {k = 1} ^ {N} \Xi_ {j k} + \sum_ {k = 1} ^ {N} \Xi_ {k j}\right)\right)}},
$$

i.e. we can express the cost of traversing a link solely as a function of observed traffic flows and the observed economic activity in each location. 

Finally, we note that equation (B.8) can be combined with equation (23) to calculate the link intensity $\pi _ { i j } ^ { k l }$ , i.e. the expected number of times a route from any origin $i$ to any destination $j$ traverses link $k l \mathrm { ~ - ~ }$ something that we show in Figure 4. Note that this expression does not require any assumptions on model parameters. 

# C Proofs

This section presents the proofs of Propositions 1 and 2. 

# C.1 Proof of Proposition 1

# C.1.1 Preliminaries

In this subsection, we show how the economic geography and urban models both are special cases of the following mathematical system: 

$$
x _ {i, 1} = D _ {i, 1} \frac {C _ {i , 1} x _ {i , 1} ^ {b _ {1 1}} x _ {i , 2} ^ {b _ {1 2}}}{\sum_ {j \in \mathcal {N}} C _ {j , 1} x _ {j , 1} ^ {b _ {1 1}} x _ {j , 2} ^ {b _ {1 2}}} x _ {i, 1} ^ {\frac {a}{a + 1}} x _ {i, 2} ^ {- \frac {1}{a + 1}} + \sum_ {j \in \mathcal {N}} K _ {i j} x _ {j, 1} ^ {\frac {1}{a + 1}} x _ {j, 2} ^ {- \frac {a}{a + 1}} \forall i \in \mathcal {N} \tag {C.1}
$$

$$
x _ {i, 2} = D _ {i, 2} \frac {C _ {i , 2} x _ {i , 1} ^ {b _ {2 1}} x _ {i , 2} ^ {b _ {2 2}}}{\sum_ {j \in \mathcal {N}} C _ {j , 2} x _ {j , 1} ^ {b _ {2 1}} x _ {j , 2} ^ {b _ {2 2}}} x _ {i, 1} ^ {- \frac {1}{a + 1}} x _ {i, 2} ^ {\frac {a}{a + 1}} + \sum_ {j \in \mathcal {N}} K _ {j i} x _ {j, 1} ^ {- \frac {a}{a + 1}} x _ {j, 2} ^ {\frac {1}{a + 1}} \forall i \in \mathcal {N} \tag {C.2}
$$

In the following sections, we then prove the existence and provide conditions for the uniqueness of any equilibrium characterized by equations (C.1) and (C.2). 

To begin, note that gravity equation in both frameworks can be written as follows: 

$$
F _ {i j} = \tau_ {i j} ^ {- \theta} \times \frac {\gamma_ {i}}{\Pi_ {i} ^ {- \theta}} \times \frac {\delta_ {j}}{P _ {j} ^ {- \theta}}, \tag {C.3}
$$

where in the economic geography model $F _ { i j } \equiv X _ { i j }$ , $\gamma _ { i } \equiv Y _ { i }$ , and $\delta _ { j } \equiv E _ { j }$ and in the urban model $F _ { i j } \equiv L _ { i j }$ , $\gamma _ { i } \equiv L _ { i } ^ { R }$ , and $\delta _ { j } \equiv L _ { i } ^ { F }$ . Written like this, the market market clearing conditions in both frameworks can be expressed identically as follows: 

$$
\gamma_ {i} = \sum_ {j} F _ {i j} \tag {C.4}
$$

$$
\delta_ {j} = \sum_ {i} F _ {i j} \tag {C.5}
$$

Substituting the gravity equation (C.3) into the two market clearing conditions yields: 

$$
\Pi_ {i} ^ {- \theta} = \sum_ {j} \tau_ {i j} ^ {- \theta} \times \frac {\delta_ {j}}{P _ {j} ^ {- \theta}}
$$

$$
P _ {j} ^ {- \theta} = \sum_ {i} \tau_ {i j} ^ {- \theta} \times \frac {\gamma_ {i}}{\Pi_ {i} ^ {- \theta}}
$$

Substituting in the expression for the endogenous trade costs $\tau _ { i j } ^ { - \theta }$ as a function of the transportation network from equation (21) and inverting each linear equation, the system becomes: 

$$
P _ {i} ^ {- \theta} \Pi_ {i} ^ {- \theta} = \delta_ {i} + \sum_ {j} t _ {i j} ^ {- \theta} P _ {i} ^ {- \theta} \Pi_ {j} ^ {- \theta}
$$

$$
P _ {i} ^ {- \theta} \Pi_ {i} ^ {- \theta} = \gamma_ {i} + \sum_ {j} t _ {j i} ^ {- \theta} P _ {j} ^ {- \theta} \Pi_ {i} ^ {- \theta}
$$

Finally, expressing the transportation network as a function of equilibrium traffic and the infrastructure network from equation (26): 

$$
\left(P _ {i} ^ {- \theta}\right) ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \Pi_ {i} ^ {- \theta} = \delta_ {i} \left(P _ {i} ^ {- \theta}\right) ^ {- \frac {1}{1 + \theta \lambda}} + \sum_ {j} \left(\bar {t} _ {i j} ^ {\frac {1 - \theta}{1 + \theta \lambda}}\right) \left(\Pi_ {j} ^ {- \theta}\right) ^ {\frac {1}{1 + \theta \lambda}}
$$

$$
P _ {i} ^ {- \theta} \left(\Pi_ {i} ^ {- \theta}\right) ^ {\frac {\theta \lambda}{1 + \theta \lambda}} = \gamma_ {i} \left(\Pi_ {i} ^ {- \theta}\right) ^ {- \frac {1}{1 + \theta \lambda}} + \sum_ {j} \bar {t} _ {j i} ^ {- \frac {\theta}{1 + \theta \lambda}} \left(P _ {j} ^ {- \theta}\right) ^ {\frac {1}{1 + \theta \lambda}},
$$

or, by defining $p _ { i } \equiv P _ { i } ^ { - \theta }$ , $\pi _ { i } \equiv \Pi _ { i } ^ { - \theta }$ , $\begin{array} { r } { a \equiv \frac { \theta \lambda } { 1 + \theta \lambda } } \end{array}$ , and $K _ { i j } \equiv \bar { t } _ { i j } ^ { - \frac { \theta } { 1 + \theta \lambda } }$ we can write this more succinctly as: 

$$
p _ {i} ^ {a} \pi_ {i} = \delta_ {i} p _ {i} ^ {- (1 - a)} + \sum_ {j} K _ {i j} \pi_ {j} ^ {(1 - a)} \tag {C.6}
$$

$$
p _ {i} \pi_ {i} ^ {a} = \gamma_ {i} \pi_ {i} ^ {- (1 - a)} + \sum_ {j} K _ {j i} p _ {j} ^ {(1 - a)}. \tag {C.7}
$$

We proceed by one final change of variables to get equations (C.6) and (C.7) into a form more amenable to define an operator to establish existence and uniqueness. Define $x _ { i , 1 } \equiv p _ { i } ^ { a } \pi _ { i }$ and $x _ { i , 2 } \equiv p _ { i } \pi _ { i } ^ { a }$ . Note this in turn implies pi = xi,1 $p _ { i } = x _ { i , 1 } ^ { \frac { a } { a ^ { 2 } - 1 } } x _ { i , 2 } ^ { - \frac { 1 } { a ^ { 2 } - 1 } }$ aa2−1 x i, 2 and πi = x ai,1 $\pi _ { i } = x _ { i , 1 } ^ { - { \frac { 1 } { a ^ { 2 } - 1 } } } x _ { i , 2 } ^ { \frac { a } { a ^ { 2 } - 1 } }$ x i, 2 , so that equations (C.6) and (C.7) become: 

$$
x _ {i, 1} = \delta_ {i} x _ {i, 1} ^ {\frac {a}{a + 1}} x _ {i, 2} ^ {- \frac {1}{a + 1}} + \sum_ {j} K _ {i j} x _ {j, 1} ^ {\frac {1}{a + 1}} x _ {j, 2} ^ {- \frac {a}{a + 1}} \tag {C.8}
$$

$$
x _ {i, 2} = \gamma_ {i} x _ {i, 1} ^ {- \frac {1}{a + 1}} x _ {i, 2} ^ {\frac {a}{a + 1}} + \sum_ {j} K _ {j i} x _ {j, 1} ^ {- \frac {a}{a + 1}} x _ {j, 2} ^ {\frac {1}{a + 1}}. \tag {C.9}
$$

The final step is to write $\{ \gamma _ { i } , \delta _ { i } \}$ as functions of $\{ x _ { i , 1 } , x _ { i , 2 } \}$ . As mentioned in the text, this mapping between the endogenous economic activity and the market access variables differs depending on the model considered, so we do it separately for each. 

The Economic Geography Model In the economic geography model in equilibrium, we have $\gamma _ { i } = \delta _ { i } = Y _ { i }$ . From welfare equalization equation (9) we have: 

$$
P _ {i} = \frac {1}{\bar {W}} \bar {u} _ {i} L _ {i} ^ {\beta - 1} Y _ {i}
$$

and from the definition of the producer price index equation (7) we have: 

$$
\Pi_ {i} = \bar {A} _ {i} L _ {i} ^ {1 + \alpha} Y _ {i} ^ {- \frac {\theta + 1}{\theta}}.
$$

Combining these two equations yields a log-linear system, which can be inverted to write $Y _ { i }$ as a function of $p _ { i }$ and $\pi _ { i }$ , as follows: 

$$
\left( \begin{array}{c} \ln p _ {i} \\ \ln \pi_ {i} \end{array} \right) = \left( \begin{array}{c} \theta \ln \bar {W} - \theta \ln \bar {u} _ {i} \\ - \theta \ln \bar {A} \end{array} \right) + \left( \begin{array}{c c} \theta   (1 - \beta) & - \theta \\ - \theta   (1 + \alpha) & (1 + \theta) \end{array} \right) \left( \begin{array}{c} \ln L _ {i} \\ \ln Y _ {i} \end{array} \right) \iff
$$

$$
\left( \begin{array}{c} \ln L _ {i} \\ \ln Y _ {i} \end{array} \right) = \left( \begin{array}{c c} \theta   (1 - \beta) & - \theta \\ - \theta   (1 + \alpha) & (1 + \theta) \end{array} \right) ^ {- 1} \left( \begin{array}{c} \ln p _ {i} - \theta \ln \bar {W} + \theta \ln \bar {u} _ {i} \\ \ln \pi_ {i} + \theta \ln \bar {A} \end{array} \right) \Longleftrightarrow
$$

$$
\left( \begin{array}{c} \ln L _ {i} \\ \ln Y _ {i} \end{array} \right) = \left( \begin{array}{c c} \frac {\left(\frac {1 + \theta}{\theta}\right)}{1 - \beta - \theta (\alpha + \beta)} & \frac {1}{1 - \beta - \theta (\alpha + \beta)} \\ \frac {(1 + \alpha)}{1 - \beta - \theta (\alpha + \beta)} & \frac {(1 - \beta)}{1 - \beta - \theta (\alpha + \beta)} \end{array} \right) \left( \begin{array}{c} \ln p _ {i} - \theta \ln \bar {W} + \theta \ln \bar {u} _ {i} \\ \ln \pi_ {i} + \theta \ln \bar {A} \end{array} \right),
$$

so that: 

$$
Y _ {i} = p _ {i} ^ {\frac {(1 + \alpha)}{1 - \beta - \theta (\alpha + \beta)}} \pi_ {i} ^ {\frac {(1 - \beta)}{1 - \beta - \theta (\alpha + \beta)}} \left(\frac {\bar {u} _ {i}}{W}\right) ^ {\theta^ {\frac {(1 + \alpha)}{1 - \beta - \theta (\alpha + \beta)}}} \bar {A} _ {i} ^ {\theta^ {\frac {(1 - \beta)}{1 - \beta - \theta (\alpha + \beta)}}}.
$$

We can then express $Y _ { i }$ as a function of $\{ x _ { 1 , i } , x _ { 2 , i } \}$ as follows: 

$$
Y _ {i} = x _ {i, 1} ^ {\frac {(a (1 + \alpha) - (1 - \beta))}{(a ^ {2} - 1) (1 - \beta - \theta (\alpha + \beta))}} x _ {i, 2} ^ {\frac {a (1 - \beta) - (1 + \alpha)}{(a ^ {2} - 1) (1 - \beta - \theta (\alpha + \beta))}} \left(\frac {\bar {u} _ {i}}{W}\right) ^ {\theta^ {\frac {(1 + \alpha)}{1 - \beta - \theta (\alpha + \beta)}}} \bar {A} _ {i} ^ {\theta^ {\frac {(1 - \beta)}{1 - \beta - \theta (\alpha + \beta)}}},
$$

or, more succinctly, as: 

$$
Y _ {i} = W ^ {\rho} C _ {i} x _ {i, 1} ^ {b _ {1}} x _ {i, 2} ^ {b _ {2}}, \tag {C.10}
$$

where $\begin{array} { r } { b _ { 1 } \equiv \frac { a ( 1 + \alpha ) - ( 1 - \beta ) } { ( a ^ { 2 } - 1 ) ( 1 - \beta - \theta ( \alpha + \beta ) ) } } \end{array}$ , $\begin{array} { r } { b _ { 2 } \equiv \frac { a ( 1 - \beta ) - ( 1 + \alpha ) } { ( a ^ { 2 } - 1 ) ( 1 - \beta - \theta ( \alpha + \beta ) ) } } \end{array}$ , $C _ { i } \equiv \bar { u } _ { i } ^ { \theta { \frac { \left( 1 + \alpha \right) } { 1 - \beta - \theta \left( \alpha + \beta \right) } } } \bar { A } _ { i } ^ { \theta { \frac { \left( 1 - \beta \right) } { 1 - \beta - \theta \left( \alpha + \beta \right) } } }$ u¯θ 1−β−θ(α+β) (1+α) Ai ¯θ (1−β)1−β−θ(α+β) , and ρ ≡ −θ 1−β−θ(α+β) . $\begin{array} { r } { \mathrm { a n d } \rho \equiv - \theta \frac { ( 1 + \alpha ) } { 1 - \beta - \theta ( \alpha + \beta ) } . } \end{array}$ (1+α) Substituting equation (C.10) into the equilibrium system (C.8) and (C.9) yields: 

$$
x _ {i, 1} = W ^ {\rho} C _ {i} x _ {i, 1} ^ {\frac {a}{a + 1} + b _ {1}} x _ {i, 2} ^ {- \frac {1}{a + 1} + b _ {2}} + \sum_ {j} K _ {i j} x _ {j, 1} ^ {\frac {1}{a + 1}} x _ {j, 2} ^ {- \frac {a}{a + 1}} \tag {C.11}
$$

$$
x _ {i, 2} = W ^ {\rho} C _ {i} x _ {i, 1} ^ {- \frac {1}{a + 1} + b _ {1}} x _ {i, 2} ^ {\frac {a}{a + 1} + b _ {2}} + \sum_ {j} K _ {j i} x _ {j, 1} ^ {- \frac {a}{a + 1}} x _ {j, 2} ^ {\frac {1}{a + 1}}. \tag {C.12}
$$

Equations (C.11) and (C.12) will form the basis of the uniqueness analysis below (where we additionally impose a symmetric infrastructure matrix). For the analysis of existence, we proceed one step further and impose the normalization that $\textstyle \sum _ { i } Y _ { i } = L$ , which then allows us to write: 

$$
Y _ {i} = \bar {L} \frac {C _ {i} x _ {i , 1} ^ {b _ {1}} x _ {i , 2} ^ {b _ {2}}}{\sum_ {j} C _ {j} x _ {j , 1} ^ {b _ {1}} x _ {j , 2} ^ {b _ {2}}}
$$

so that the system of equations becomes: 

$$
x _ {i, 1} = \bar {L} \frac {C _ {i} x _ {i , 1} ^ {b _ {1}} x _ {i , 2} ^ {b _ {2}}}{\sum_ {j} C _ {j} x _ {j , 1} ^ {b _ {1}} x _ {j , 2} ^ {b _ {2}}} x _ {i, 1} ^ {\frac {a}{a + 1}} x _ {i, 2} ^ {- \frac {1}{a + 1}} + \sum_ {j} K _ {i j} x _ {j, 1} ^ {\frac {1}{a + 1}} x _ {j, 2} ^ {- \frac {a}{a + 1}} \tag {C.13}
$$

$$
x _ {i, 2} = \bar {L} \frac {C _ {i} x _ {i , 1} ^ {b _ {1}} x _ {i , 2} ^ {b _ {2}}}{\sum_ {j} C _ {j} x _ {j , 1} ^ {b _ {1}} x _ {j , 2} ^ {b _ {2}}} x _ {i, 1} ^ {- \frac {1}{a + 1}} x _ {i, 2} ^ {\frac {a}{a + 1}} + \sum_ {j} K _ {j i} x _ {j, 1} ^ {- \frac {a}{a + 1}} x _ {j, 2} ^ {\frac {1}{a + 1}}, \tag {C.14}
$$

which is a special case of equations (C.1) and (C.2), as claimed. 

The Urban Model In the urban model, we have $\gamma _ { i } \equiv L _ { i } ^ { R }$ and $\delta _ { i } \equiv L _ { i } ^ { F }$ and, from equations (17) and (18) we have: 

$$
P _ {i} ^ {- \theta} = \bar {A} _ {i} ^ {- \theta} \left(\frac {\bar {L}}{\bar {W} ^ {\theta}}\right) ^ {\frac {1}{2}} \left(L _ {i} ^ {F}\right) ^ {1 - \theta \alpha}
$$

$$
\Pi_ {i} ^ {- \theta} = \bar {u} _ {i} ^ {- \theta} \left(\frac {\bar {L}}{\bar {W} ^ {\theta}}\right) ^ {\frac {1}{2}} \left(L _ {i} ^ {R}\right) ^ {1 - \theta \beta},
$$

so that $L _ { i } ^ { F } = p _ { i } ^ { \frac { 1 } { 1 - \theta \alpha } } \bar { A } _ { i } ^ { \frac { \theta } { 1 - \theta \alpha } } \bigg ( \frac { \bar { W } ^ { \theta } } { \bar { L } } \bigg ) ^ { \frac { 1 } { 2 ( 1 - \theta \alpha ) } }$ A¯ 1i  W¯ θ  and $L _ { i } ^ { R } = \pi _ { i } ^ { \frac { 1 } { 1 - \theta \beta } } \bar { u } _ { i } ^ { \frac { \theta } { 1 - \theta \beta } } \bigg ( \frac { \bar { W } ^ { \theta } } { \bar { L } } \bigg ) ^ { \frac { 1 } { 2 ( 1 - \theta \beta ) } }$ = π 1i u¯ 1i We can then express $L _ { i } ^ { F ^ { \prime } }$ and $L _ { i } ^ { R }$ as a function of $\{ x _ { 1 , i } , x _ { 2 , i } \}$ as follows: 

$$
L _ {i} ^ {F} = x _ {i, 1} ^ {\frac {a}{a ^ {2} - 1} \frac {1}{1 - \theta \alpha}} x _ {i, 2} ^ {- \frac {1}{a ^ {2} - 1} \frac {1}{1 - \theta \alpha}} \bar {A} _ {i} ^ {\frac {\theta}{1 - \theta \alpha}} \left(\frac {\bar {W} ^ {\theta}}{\bar {L}}\right) ^ {\frac {1}{2 (1 - \theta \alpha)}},
$$

$$
L _ {i} ^ {R} = x _ {i, 1} ^ {- \frac {1}{a ^ {2} - 1} \frac {1}{1 - \theta \beta}} x _ {i, 2} ^ {\frac {a}{a ^ {2} - 1} \frac {1}{1 - \theta \beta}} \bar {u} _ {i} ^ {\frac {\theta}{1 - \theta \beta}} \left(\frac {\bar {W} ^ {\theta}}{\bar {L}}\right) ^ {\frac {1}{2 (1 - \theta \beta)}}
$$

or more succinctly as: 

$$
L _ {i} ^ {F} = \bar {W} ^ {\rho_ {1}} C _ {i, 1} x _ {i, 1} ^ {b _ {1 1}} x _ {i, 2} ^ {b _ {1 2}} \tag {C.15}
$$

$$
L _ {i} ^ {R} = \bar {W} ^ {\rho_ {2}} C _ {i, 2} x _ {i, 1} ^ {b _ {2 1}} x _ {i, 2} ^ {b _ {2 2}} \tag {C.16}
$$

where aa2−1 1−θβ , ρ1 ≡ $\frac { a } { a ^ { 2 } - 1 } \frac { 1 } { 1 - \theta \beta }$ 1 $\begin{array} { r l r } { C _ { i , 1 } } & { { } \equiv } & { \frac { \bar { A } _ { i } ^ { \frac { \theta } { 1 - \theta \alpha } } } { \bar { L } ^ { \frac { 1 } { 2 } } \frac { 1 } { 1 - \theta \alpha } } , C _ { i , 2 } \equiv \frac { \bar { u } _ { i } ^ { \frac { \theta } { 1 - \theta \beta } } } { \bar { L } ^ { \frac { 1 } { 2 } } \frac { 1 } { 1 - \theta \beta } } } \end{array}$ $\begin{array} { r } { \rho _ { 1 } \equiv \frac { \theta } { 2 ( 1 - \theta \alpha ) } } \end{array}$ A¯ θ1−θαi L¯ 12 11 − θ α ,Ci,2 and $\begin{array} { r } { \rho _ { 2 } \equiv \frac { \theta } { 2 ( 1 - \theta \beta ) } } \end{array}$ u¯ 1 θ−θβ b11 $\begin{array} { r } { b _ { 1 1 } \equiv \frac { a } { a ^ { 2 } - 1 } \frac { 1 } { 1 - \theta \alpha } } \end{array}$ . Substituting equations (C.15) and (C.16) into equations into the , $\begin{array} { r } { b _ { 1 2 } ~ \equiv ~ - \frac { 1 } { a ^ { 2 } - 1 } \frac { 1 } { 1 - \theta \alpha } } \end{array}$ , $\begin{array} { r } { b _ { 2 1 } \equiv - \frac { 1 } { a ^ { 2 } - 1 } \frac { 1 } { 1 - \theta \beta } } \end{array}$ − a2−1 11−θβ , b22 ≡ 1 1 $b _ { 2 2 } \equiv$ 

equilibrium system (C.8) and (C.9) yields: 

$$
x _ {i, 1} = W ^ {\rho_ {1}} C _ {i, 1} x _ {i, 1} ^ {\frac {a}{a + 1} + b _ {1 1}} x _ {i, 2} ^ {- \frac {1}{a + 1} + b _ {1 2}} + \sum_ {j} K _ {i j} x _ {j, 1} ^ {\frac {1}{a + 1}} x _ {j, 2} ^ {- \frac {a}{a + 1}} \tag {C.17}
$$

$$
x _ {i, 2} = W ^ {\rho_ {2}} C _ {i, 2} x _ {i, 1} ^ {- \frac {1}{a + 1} + b _ {2 1}} x _ {i, 2} ^ {\frac {a}{a + 1} + b _ {2 2}} + \sum_ {j} K _ {j i} x _ {j, 1} ^ {- \frac {a}{a + 1}} x _ {j, 2} ^ {\frac {1}{a + 1}}. \tag {C.18}
$$

Equations (C.17) and (C.18) will form the basis of the uniqueness analysis below. For the analysis of existence, we proceed one step further and impose the normalization that $\begin{array} { r } { \sum _ { i } L _ { i } ^ { F } = \sum _ { i } L _ { i } ^ { R } = L } \end{array}$ ,  which then allows us to write: 

$$
L _ {i} ^ {F} = \frac {C _ {i , 1} x _ {i , 1} ^ {b _ {1 1}} x _ {i , 2} ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} x _ {j , 1} ^ {b _ {1 1}} x _ {j , 2} ^ {b _ {1 2}}} \bar {L}
$$

$$
L _ {i} ^ {R} = \frac {C _ {i , 2} x _ {i , 1} ^ {b _ {2 1}} x _ {i , 2} ^ {b _ {2 2}}}{\sum_ {j} C _ {j , 2} x _ {j , 1} ^ {b _ {2 1}} x _ {j , 2} ^ {b _ {2 2}}} \bar {L}
$$

so that the system of equations becomes: 

$$
x _ {i, 1} = \bar {L} \frac {C _ {i , 1} x _ {i , 1} ^ {b _ {1 1}} x _ {i , 2} ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} x _ {j , 1} ^ {b _ {1 1}} x _ {j , 2} ^ {b _ {1 2}}} x _ {i, 1} ^ {\frac {a}{a + 1}} x _ {i, 2} ^ {- \frac {1}{a + 1}} + \sum_ {j} K _ {i j} x _ {j, 1} ^ {\frac {1}{a + 1}} x _ {j, 2} ^ {- \frac {a}{a + 1}} \tag {C.19}
$$

$$
x _ {i, 2} = \bar {L} \frac {C _ {i , 2} x _ {i , 1} ^ {b _ {2 1}} x _ {i , 2} ^ {b _ {2 2}}}{\sum_ {j} C _ {j , 2} x _ {j , 1} ^ {b _ {2 1}} x _ {j , 2} ^ {b _ {2 2}}} x _ {i, 1} ^ {- \frac {1}{a + 1}} x _ {i, 2} ^ {\frac {a}{a + 1}} + \sum_ {j} K _ {j i} x _ {j, 1} ^ {- \frac {a}{a + 1}} x _ {j, 2} ^ {\frac {1}{a + 1}}. \tag {C.20}
$$

which is a special case of equations (C.1) and (C.2), as claimed. 

# C.1.2 Part 1 (Existence)

We first turn to the existence of the system. To do so, we rely on Brouwer’s fixed point theorem. A natural way to proceed would be to define the operator $T ( x ) : \mathbb { R } _ { + + } ^ { 2 N }  \mathbb { R } _ { + + } ^ { 2 N }$ such that: 

$$
T \left(x\right) = \left( \begin{array}{c} \left(T _ {1} \left(x\right)\right) _ {i \in \mathcal {N}} \\ \left(T _ {2} \left(x\right)\right) _ {i \in \mathcal {N}} \end{array} \right) \equiv \left( \begin{array}{c} \left(D _ {i, 1} \frac {C _ {i , 1} x _ {i , 1} ^ {b _ {1 1}} x _ {i , 2} ^ {b _ {1 2}}}{\sum_ {j \in \mathcal {N}} C _ {j , 1} x _ {j , 1} ^ {b _ {1 1}} x _ {j , 2} ^ {b _ {1 2}}} x _ {i, 1} ^ {\frac {a}{a + 1}} x _ {i, 2} ^ {- \frac {1}{a + 1}} + \sum_ {j \in \mathcal {N}} K _ {i j} x _ {j, 1} ^ {\frac {1}{a + 1}} x _ {j, 2} ^ {- \frac {a}{a + 1}}\right) _ {i} \\ \left(D _ {i, 2} \frac {C _ {i , 2} x _ {i , 1} ^ {b _ {2 1}} x _ {i , 2} ^ {b _ {2 2}}}{\sum_ {j \in \mathcal {N}} C _ {j , 2} x _ {j , 1} ^ {b _ {2 1}} x _ {j , 2} ^ {b _ {2 2}}} x _ {i, 1} ^ {- \frac {1}{a + 1}} x _ {i, 2} ^ {\frac {a}{a + 1}} + \sum_ {j \in \mathcal {N}} K _ {j i} x _ {j, 1} ^ {- \frac {a}{a + 1}} x _ {j, 2} ^ {\frac {1}{a + 1}}\right) _ {i} \end{array} \right)
$$

While it is immediately obvious that $T$ is continuous, unfortunately, $T$ does not operate on a compact space, and hence Brouwer’s theorem cannot be directly applied to it. Instead, following Allen, Arkolakis, and Takahashi (2020), we consider an alternative “scaled” system, the equilibrium of which will turn out to be a scaled version of the equilibrium of our system. 

Define the operator: $\tilde { T } ( x ) \equiv ( \tilde { T } _ { 1 } ( x ) _ { i } ) : \mathbb { R } _ { + + } ^ { 2 N }  \mathbb { R } _ { + + } ^ { 2 N }$ T˜1 (x)i T˜2 (x)i 1 such that: 

$$
\begin{array}{l} \tilde {T} _ {1} (x) _ {i} \equiv \frac {D _ {i , 1} \frac {C _ {i , 1} x _ {i , 1} ^ {b _ {1 1}} x _ {i , 2} ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} x _ {j , 1} ^ {b _ {1 1}} x _ {j , 2} ^ {b _ {1 2}}} x _ {i , 1} ^ {\frac {a}{a + 1}} x _ {i , 2} ^ {- \frac {1}{a + 1}}}{\left(\sum_ {i} \left(D _ {i , 1} \frac {C _ {i , 1} x _ {i , 1} ^ {b _ {1 1}} x _ {i , 2} ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} x _ {j , 1} ^ {b _ {1 1}} x _ {j , 2} ^ {b _ {1 1 2}}} x _ {i , 1} ^ {\frac {a}{a + 1}} x _ {i , 2} ^ {- \frac {1}{a + 1}} + \sum_ {j} K _ {i j} x _ {j , 1} ^ {\frac {1}{a + 1}} x _ {j , 2} ^ {- \frac {a}{a + 1}}\right) + 1\right) ^ {c _ {1 1}}} \\ + \frac {\sum_ {j} K _ {i j} x _ {j , 1} ^ {\frac {1}{a + 1}} x _ {j , 2} ^ {- \frac {a}{a + 1}}}{\left(\sum_ {i} \left(D _ {i , 1} \frac {C _ {i , 1} x _ {i , 1} ^ {b _ {1 1}} x _ {i , 2} ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} x _ {j , 1} ^ {b _ {1 1}} x _ {j , 2} ^ {b _ {1 2}}} x _ {i , 1} ^ {\frac {a}{a + 1}} x _ {i , 2} ^ {- \frac {1}{a + 1}} + \sum_ {j} K _ {i j} x _ {j , 1} ^ {\frac {1}{a + 1}} x _ {j , 2} ^ {- \frac {a}{a + 1}}\right) + 1\right) ^ {c _ {1 2}}} \\ \end{array}
$$

$$
\begin{array}{l} \tilde {T} _ {2} (x) _ {i} \equiv \frac {D _ {i , 2} \frac {C _ {i , 2} x _ {i , 1} ^ {b _ {2 1}} x _ {i , 2} ^ {b _ {2 2}}}{\sum_ {j} C _ {j , 2} x _ {j , 1} ^ {b _ {2 1}} x _ {j , 2} ^ {b _ {2 2}}} x _ {i , 1} ^ {- \frac {1}{a + 1}} x _ {i , 2} ^ {\frac {a}{a + 1}}}{\left(\sum_ {i} \left(D _ {i , 2} \frac {C _ {i , 2} x _ {i , 1} ^ {b _ {2 1}} x _ {i , 2} ^ {b _ {2 2}}}{\sum_ {j} C _ {j , 2} x _ {j , 1} ^ {b _ {2 1}} x _ {j , 2} ^ {b _ {2 1}}} x _ {i , 1} ^ {- \frac {1}{a + 1}} x _ {i , 2} ^ {\frac {a}{a + 1}} + \sum_ {j} K _ {j i} x _ {j , 1} ^ {- \frac {a}{a + 1}} x _ {j , 2} ^ {\frac {1}{a + 1}}\right) + 1\right) ^ {c _ {2 1}}} \\ + \frac {\sum_ {j} K _ {j i} x _ {j , 1} ^ {- \frac {a}{a + 1}} x _ {j , 2} ^ {\frac {1}{a + 1}}}\left(\sum_ {i} \left(D _ {i , 2} \frac {C _ {i , 2} x _ {i , 1} ^ {b _ {2 1}} x _ {i , 2} ^ {b _ {2 2}}}{\sum_ {j} C _ {j , 2} x _ {j , 1} ^ {b _ {2 1}} x _ {j , 2} ^ {b _ {2 2}}} x _ {i , 1} ^ {- \frac {1}{a + 1}} x _ {i , 2} ^ {\frac {a}{a + 1}} + \sum_ {j} K _ {j i} x _ {j , 1} ^ {- \frac {a}{a + 1}} x _ {j , 2} ^ {\frac {1}{a + 1}}\right) + 1\right) ^ {c _ {2 2}}, \\ \end{array}
$$

where $x \equiv \binom { \left( x _ { i , 1 } \right) _ { i } } { \left( x _ { i , 2 } \right) _ { i } } \in \mathbb { R } _ { + + } ^ { 2 N }$ (xi,1)i and ${ \binom { c _ { 1 1 } } { c _ { 2 1 } } } c _ { 1 2 } ) \equiv { \binom { \frac { 1 } { a } } { \frac { 1 } { a } } } 1 )$ are positive scalars weakly greater than one. Like with $T$ , it is immediately obvious that $\ddot { T }$ is continuous. We show now that it also maps a compact space to itself, thereby allowing us to establish the existence of a fixed point of $\tilde { T }$ . 

Define the compact space $\mathcal { M } \subset \mathbb { R } _ { + + } ^ { 2 N } \equiv \left\{ x \in \mathbb { R } _ { + + } ^ { 2 N } | x _ { i l } \in [ 0 , M _ { l } ] \ \forall l \in \left\{ 1 , 2 \right\} , i \in \mathcal { N } \right\}$ , where: 

$$
M _ {1} \equiv \max \left\{\max _ {i} \frac {D _ {i , 1}}{\sum_ {l} D _ {l , 1}}, \max _ {i} \sum_ {j} \frac {K _ {i j}}{\sum_ {l} K _ {l j}} \right\}
$$

$$
M _ {2} \equiv \max  \left\{\max  _ {i} \frac {D _ {i , 2}}{\sum_ {l} D _ {l , 2}}, \max  _ {i} \sum_ {j} \frac {K _ {j i}}{\sum_ {l} K _ {j l}} \right\}
$$

We claim that for all $x \in \mathcal { M }$ , $\tilde { T } \left( x \right) \in \mathcal { M }$ ,i.e. $\ddot { T }$ operates from a compact space $\mathcal { M }$ to itself. It is immediately evident that for all $x \in \mathcal { M }$ , $\tilde { T } \left( x \right) \geq 0$ , so it only remains to show that $\ddot { T } _ { 1 } \left( x \right) \leq M _ { 1 }$ and $\tilde { T } _ { 1 } \left( x \right) \leq M _ { 2 }$ . We show the former; the argument for the latter proceeds similarly. Consider the first term of the $\tilde { T } _ { 1 } \left( x \right)$ , 

$$
\tilde {T} _ {1} ^ {A} (x) \equiv \frac {D _ {i , 1} \frac {C _ {i , 1} x _ {i , 1} ^ {b _ {1 1}} x _ {i , 2} ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} x _ {j , 1} ^ {b _ {1 1}} x _ {j , 2} ^ {b _ {1 2}}} x _ {i , 1} ^ {\frac {a}{a + 1}} x _ {i , 2} ^ {- \frac {1}{a + 1}}}{\left(\sum_ {i} \left(D _ {i , 1} \frac {C _ {i , 1} x _ {i , 1} ^ {b _ {1 1}} x _ {i , 2} ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} x _ {j , 1} ^ {b _ {1 1}} x _ {j , 2} ^ {b _ {1 1}}} x _ {i , 1} ^ {\frac {a}{a + 1}} x _ {i , 2} ^ {- \frac {1}{a + 1}} + \sum_ {j} K _ {i j} x _ {j , 1} ^ {\frac {1}{a + 1}} x _ {j , 2} ^ {- \frac {a}{a + 1}}\right) + 1\right) ^ {c _ {1 1}}}, \text {w h i c h c a n b e b o u n d e d a s f o l l o w s :}
$$

$$
\begin{array}{l} \tilde {T} _ {1} ^ {A} (x) = \frac {D _ {i , 1}}{\sum_ {l} D _ {l , 1}} \times \frac {\sum_ {l} D _ {l , 1}}{D _ {i , 1}} \times \frac {D _ {i , 1} \frac {C _ {i , 1} x _ {i , 1} ^ {b _ {1 1}} x _ {i , 2} ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} x _ {j , 1} ^ {b _ {1 1}} x _ {j , 2} ^ {b _ {1 2}}} x _ {i , 1} ^ {\frac {a}{a + 1}} x _ {i , 2} ^ {- \frac {1}{a + 1}}}{\left(\sum_ {i} \left(D _ {i , 1} \frac {C _ {i , 1} x _ {i , 1} ^ {b _ {1 1}} x _ {i , 2} ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} x _ {j , 1} ^ {b _ {1 1}} x _ {j , 2} ^ {b _ {1 1}}}\right) x _ {i , 1} ^ {\frac {a}{a + 1}} x _ {i , 2} ^ {- \frac {1}{a + 1}} + \sum_ {j} K _ {i j} x _ {j , 1} ^ {\frac {1}{a + 1}} x _ {j , 2} ^ {- \frac {a}{a + 1}}\right) + 1\left. \right) ^ {c _ {1 1}}} \iff \\ = \frac {D _ {i , 1}}{\sum_ {l} D _ {l , 1}} \times \frac {\sum_ {l} D _ {l , 1} \frac {C _ {i , 1} x _ {i , 1} ^ {b _ {1 1}} x _ {i , 2} ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} x _ {j , 1} ^ {b _ {1 1}} x _ {j , 2} ^ {b _ {1 2}}} x _ {i , 1} ^ {\frac {a}{a + 1}} x _ {i , 2} ^ {- \frac {1}{a + 1}}}{\left(\sum_ {i} \left(D _ {i , 1} \frac {C _ {i , 1} x _ {i , 1} ^ {b _ {1 1}} x _ {i , 2} ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} x _ {j , 1} ^ {b _ {1 1}} x _ {j , 2} ^ {b _ {1 2}}} x ^ {\frac {a}{a + 1}} x ^ {- \frac {1}{a + 1}} + \sum_ {j} K _ {i j} x _ {j , 1} ^ {\frac {1}{a + 1}} x _ {j , 2} ^ {- \frac {a}{a + 1}}\right) + 1\right) ^ {c _ {1 1}}} \Rightarrow \\ \leq \left(\max  _ {i} \frac {D _ {i , 1}}{\sum_ {l} D _ {l , 1}}\right) \times \frac {\sum_ {l} D _ {l , 1} \frac {C _ {i , 1} x _ {i , 1} ^ {b _ {1 1}} x _ {i , 2} ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} x _ {j , 1} ^ {b _ {1 1}} x _ {j , 2} ^ {b _ {1 2}}} x _ {i , 1} ^ {\frac {a}{a + 1}} x _ {i , 2} ^ {- \frac {1}{a + 1}}}{\left(\sum_ {i} \left(D _ {i , 1} \frac {C _ {i , 1} x _ {i , 1} ^ {b _ {1 1}} x _ {i , 2} ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} x _ {j , 1} ^ {b _ {1 1}} x _ {j , 2} ^ {b _ {1 2}}} x ^ {\frac {a}{a + 1}} x _ {i , 1} ^ {- \frac {1}{a + 1}} + \sum_ {j} K _ {i j} x _ {j , 1} ^ {\frac {1}{a + 1}} x _ {j , 2} ^ {- \frac {a}{a + 1}}\right) + 1\right) ^ {c _ {1 1}}} \Rightarrow \\ \leq \left(\max _ {i} \frac {D _ {i , 1}}{\sum_ {l} D _ {l , 1}}\right), \\ \end{array}
$$

where the last line follows from the fact that the denominator of the second term is strictly larger than the numerator (since $c _ { 1 1 } > 1$ and the denominator value is weakly greater than one). Similarly, consider the 

$$
\text {s e c o n d} \tilde {T} _ {1} (x), \tilde {T} _ {1} ^ {B} (x) \equiv \frac {\sum_ {j} K _ {i j} x _ {j , 1} ^ {\frac {1}{a + 1}} x _ {j , 2} ^ {- \frac {a}{a + 1}}}{\left(\sum_ {i} \left(D _ {i , 1} \frac {C _ {i , 1} x _ {i , 1} ^ {b _ {1 1}} x _ {i , 2} ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} x _ {j , 1} ^ {b _ {1 1}} x _ {j , 2} ^ {b _ {1 2}}} x _ {i , 1} ^ {\frac {a}{a + 1}} x _ {i , 2} ^ {- \frac {1}{a + 1}} + \sum_ {j} K _ {i j} x _ {j , 1} ^ {\frac {1}{a + 1}} x _ {j , 2} ^ {- \frac {a}{a + 1}}\right) + 1\right) ^ {c _ {1 2}}}, \text {w h i c h c a n b e}
$$

bounded as follows:: 

$$
\begin{array}{l} \tilde {T} _ {1} ^ {B} (x) = \sum_ {j} \frac {K _ {i j}}{\sum_ {j} K _ {l j}} \times \frac {\sum_ {l} K _ {l j}}{K _ {i j}} \times \frac {K _ {i j} x _ {j , 1} ^ {\frac {1}{a + 1}} x _ {j , 2} ^ {- \frac {a}{a + 1}}}{\left(\sum_ {i} \left(D _ {i , 1} \frac {C _ {i , 1} x _ {i , 1} ^ {b _ {1 1}} x _ {i , 2} ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} x _ {j , 1} ^ {b _ {1 1}} x _ {j , 2} ^ {b _ {1 2}}} x _ {i , 1} ^ {\frac {a}{a + 1}} x _ {i , 2} ^ {- \frac {1}{a + 1}} + \sum_ {j} K _ {i j} x _ {j , 1} ^ {\frac {1}{a + 1}} x _ {j , 2} ^ {- \frac {a}{a + 1}}\right) + 1\right) ^ {c _ {1 2}}} \iff \\ = \sum_ {j} \frac {K _ {i j}}{\sum_ {l} K _ {l j}} \times \frac {\sum_ {l} \sum_ {j} K _ {l j} x _ {j , 1} ^ {\frac {1}{a + 1}} x _ {j , 2} ^ {- \frac {a}{a + 1}}}{\left(\sum_ {i} \left(D _ {i , 1} \frac {C _ {i , 1} x _ {i , 1} ^ {b _ {1 1}} x _ {i , 2} ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} x _ {j , 1} ^ {b _ {1 1}} x _ {j , 2} ^ {b _ {1 2}}} x _ {i , 1} ^ {\frac {a}{a + 1}} x _ {i , 2} ^ {- \frac {1}{a + 1}} + \sum_ {j} K _ {i j} x _ {j , 1} ^ {\frac {1}{a + 1}} x _ {j , 2} ^ {- \frac {a}{a + 1}}\right) + 1\right) ^ {c _ {1 2}}} \Rightarrow \\ \leq \max  _ {i} \sum_ {j} \frac {K _ {i j}}{\sum_ {l} K _ {l j}} \times \frac {\left(\sum_ {l} \sum_ {j} K _ {l j} x _ {j , 1} ^ {\frac {1}{a + 1}} x _ {j , 2} ^ {- \frac {a}{a + 1}}\right)}{\left(\sum_ {i} \left(D _ {i , 1} \frac {C _ {i , 1} x _ {i , 1} ^ {b _ {1 1}} x _ {i , 2} ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} x _ {j , 1} ^ {b _ {1 1}} x _ {j , 2} ^ {b _ {1 2}}} x _ {i , 1} ^ {\frac {a}{a + 1}} x _ {i , 2} ^ {- \frac {1}{a + 1}} + \sum_ {j} K _ {i j} x _ {j , 1} ^ {\frac {1}{a + 1}} x _ {j , 2} ^ {- \frac {a}{a + 1}}\right) + 1\right) ^ {c _ {1 2}}} \Rightarrow \\ \leq \max  _ {i} \sum_ {j} \frac {K _ {i j}}{\sum_ {l} K _ {l j}}, \\ \end{array}
$$

where again the last line follows from the fact that the denominator of the second term is strictly larger than its numerator. Together we have: 

$$
\tilde {T} _ {1} (x) _ {i} \leq \max  \left\{\max  _ {i} \frac {D _ {i , 1}}{\sum_ {l} D _ {l , 1}}, \max  _ {i} \sum_ {j} \frac {K _ {i j}}{\sum_ {l} K _ {l j}} \right\},
$$

as claimed. 

Since $\tilde { T }$ is a continuous operator mapping a compact set $\mathcal { M }$ to itself, by Brouwer’s fixed point, there exists a fixed point $\tilde { x } ^ { * } \in \mathcal { M }$ such that $\tilde { T } \left( \tilde { x } ^ { * } \right) = \tilde { x } ^ { * }$ . 

The next step of the existence proof is to show that there exists an equilibrium of the system characterized by equations (C.1) and (C.2), or, equivalently, there exists a fixed point $x ^ { * }$ of the (un-scaled) operator $T$ , i.e. $x ^ { * } = T \left( x ^ { * } \right)$ . To do so, we will show that there exists a scalar $t > 0$ such that $x ^ { * } = t \tilde { x } ^ { * }$ . Before we do this, however, we first need the following result that says the denominators (inside the exponents $\{ c _ { k l } \}$ ) of the $\ddot { T } _ { 1 } \left( \tilde { x } ^ { * } \right) _ { i }$ and $\ddot { T } _ { 2 } \left( \tilde { x } ^ { * } \right) _ { i }$ are equal, i.e.: 

$$
\sum_ {i} \left(D _ {i, 1} \frac {C _ {i , 1} \left(\tilde {x} _ {i , 1} ^ {*}\right) ^ {b _ {1 1}} \left(\tilde {x} _ {i , 2} ^ {*}\right) ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} \left(\tilde {x} _ {j , 1} ^ {*}\right) ^ {b _ {1 1}} \left(\tilde {x} _ {j , 2} ^ {*}\right) ^ {b _ {1 2}}} \left(\tilde {x} _ {i, 1} ^ {*}\right) ^ {\frac {a}{a + 1}} \left(\tilde {x} _ {i, 2} ^ {*}\right) ^ {- \frac {1}{a + 1}} + \sum_ {j} K _ {i j} \left(\tilde {x} _ {j, 1} ^ {*}\right) ^ {\frac {1}{a + 1}} \left(\tilde {x} _ {j, 2} ^ {*}\right) ^ {- \frac {a}{a + 1}}\right) =
$$

$$
\sum_ {i} \left(D _ {i, 2} \frac {C _ {i , 2} \left(\tilde {x} _ {i , 1} ^ {*}\right) ^ {b _ {2 1}} \left(\tilde {x} _ {i , 2} ^ {*}\right) ^ {b _ {2 2}}}{\sum_ {j} C _ {j , 2} \left(\tilde {x} _ {i , 1} ^ {*}\right) ^ {b _ {2 1}} \left(\tilde {x} _ {i , 2} ^ {*}\right) ^ {b _ {2 2}}} \left(\tilde {x} _ {i, 1} ^ {*}\right) ^ {- \frac {1}{a + 1}} \left(\tilde {x} _ {i, 2} ^ {*}\right) ^ {\frac {a}{a + 1}} + \sum_ {j} K _ {j i} \left(\tilde {x} _ {j, 1} ^ {*}\right) ^ {- \frac {a}{a + 1}} \left(\tilde {x} _ {j, 2} ^ {*}\right) ^ {\frac {1}{a + 1}}\right) \tag {C.21}
$$

The easiest way to see this is to re-write the scaled equilibrium conditions as functions of the endogenous variables $\{ p _ { i } , \pi _ { i } \}$ (rather than $\{ x _ { i , 1 } , x _ { i , 2 } \} )$ , which becomes: 

$$
p _ {i} ^ {a} \pi_ {i} = \frac {\delta_ {i} p _ {i} ^ {- (1 - a)} + \sum_ {j} K _ {i j} \pi_ {j} ^ {(1 - a)}}{\sum_ {i} \left(\delta_ {i} p _ {i} ^ {- (1 - a)} + \sum_ {j} K _ {i j} \pi_ {j} ^ {(1 - a)}\right)}
$$

$$
p _ {i} \pi_ {i} ^ {a} = \frac {\gamma_ {i} \pi_ {i} ^ {- (1 - a)} + \sum_ {j} K _ {j i} p _ {j} ^ {(1 - a)}}{\sum_ {i} \left(\gamma_ {i} \pi_ {i} ^ {- (1 - a)} + \sum_ {j} K _ {j i} p _ {j} ^ {(1 - a)}\right)}
$$

Multiplying the first condition by $p _ { i } ^ { 1 - a }$ and summing over $i$ yields: 

$$
\sum_ {i} p _ {i} \pi_ {i} = \frac {\sum_ {i} \delta_ {i} + \sum_ {i} \sum_ {j} K _ {i j} \pi_ {j} ^ {(1 - a)} p _ {i} ^ {(1 - a)}}{\sum_ {i} \left(\delta_ {i} p _ {i} ^ {- (1 - a)} + \sum_ {j} K _ {i j} \pi_ {j} ^ {(1 - a)}\right)}
$$

Similarly, multiplying the second condition by $\pi _ { i } ^ { 1 - a }$ and summing over $i$ yields: 

$$
\sum_ {i} p _ {i} \pi_ {i} = \frac {\sum_ {i} \gamma_ {i} + \sum_ {i} \sum_ {j} K _ {j i} p _ {j} ^ {(1 - a)} \pi_ {i} ^ {(1 - a)}}{\sum_ {i} \left(\gamma_ {i} \pi_ {i} ^ {- (1 - a)} + \sum_ {j} K _ {j i} p _ {j} ^ {(1 - a)}\right)}
$$

Since the left hand side of both conditions are the same and the numerators on the right hand side are both the same (as recall $\begin{array} { r } { \sum _ { i } \delta _ { i } = \sum _ { i } \gamma _ { i } = L } \end{array}$ in both the economic geography and urban models), the denominators must also be the same. Hence we have: 

$$
\sum_ {i} \left(\delta_ {i} p _ {i} ^ {- (1 - a)} + \sum_ {j} K _ {i j} \pi_ {j} ^ {(1 - a)}\right) = \sum_ {i} \left(\gamma_ {i} \pi_ {i} ^ {- (1 - a)} + \sum_ {j} K _ {j i} p _ {j} ^ {(1 - a)}\right),
$$

or in $\{ x _ { i , 1 } , x _ { i , 2 } \}$ space, equation (C.21) holds, as claimed. 

Armed with this result, we are now prepared to construct a solution of equations (C.1) and (C.2). We posit that for all $i \in \mathcal N$ and $l \in \{ 1 , 2 \}$ we have: 

$$
x _ {i, l} ^ {*} = t \tilde {x} _ {i, l} ^ {*} \tag {C.22}
$$

where: 

$$
t \equiv \left(\sum_ {i} \left(D _ {i, 1} \frac {C _ {i , 1} \left(\tilde {x} _ {i , 1} ^ {*}\right) ^ {b _ {1 1}} \left(\tilde {x} _ {i , 2} ^ {*}\right) ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} \left(\tilde {x} _ {j , 1} ^ {*}\right) ^ {b _ {1 1}} \left(\tilde {x} _ {j , 2} ^ {*}\right) ^ {b _ {1 2}}} \left(\tilde {x} _ {i, 1} ^ {*}\right) ^ {\frac {a}{a + 1}} \left(\tilde {x} _ {i, 2} ^ {*}\right) ^ {- \frac {1}{a + 1}} + \sum_ {j} K _ {i j} \left(\tilde {x} _ {j, 1} ^ {*}\right) ^ {\frac {1}{a + 1}} \left(\tilde {x} _ {j, 2} ^ {*}\right) ^ {- \frac {a}{a + 1}}\right) + 1\right) ^ {\frac {a + 1}{2 a}} \Longleftrightarrow \tag {C.23}
$$

$$
t = \left(\sum_ {i} \left(D _ {i, 2} \frac {C _ {i , 2} \left(\tilde {x} _ {i , 1} ^ {*}\right) ^ {b _ {2 1}} \left(\tilde {x} _ {i , 2} ^ {*}\right) ^ {b _ {2 2}}}{\sum_ {j} C _ {j , 2} \left(\tilde {x} _ {i , 1} ^ {*}\right) ^ {b _ {2 1}} \left(\tilde {x} _ {i , 2} ^ {*}\right) ^ {b _ {2 2}}} \left(\tilde {x} _ {i, 1} ^ {*}\right) ^ {- \frac {1}{a + 1}} \left(\tilde {x} _ {i, 2} ^ {*}\right) ^ {\frac {a}{a + 1}} + \sum_ {j} K _ {j i} \left(\tilde {x} _ {j, 1} ^ {*}\right) ^ {- \frac {a}{a + 1}} \left(\tilde {x} _ {j, 2} ^ {*}\right) ^ {\frac {1}{a + 1}}\right) + 1\right) ^ {\frac {a + 1}{2 a}}, \tag {C.24}
$$

where the second line immediately follows from equation (C.21). 

Since $\tilde { x } ^ { * } = \ddot { T } \left( \tilde { x } ^ { * } \right)$ , we first consider $\tilde { x } _ { 1 } ^ { * } = \tilde { T } _ { 1 } \left( \tilde { x } ^ { * } \right)$ . Imposing (C.22) and (C.23) yields: 

$$
\frac {1}{t} x _ {i, 1} ^ {*} = \frac {D _ {i , 1} \frac {C _ {i , 1} \left(\frac {1}{t} x _ {i , 1} ^ {*}\right) ^ {b _ {1 1}} \left(\frac {1}{t} x _ {i , 2} ^ {*}\right) ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} \left(\frac {1}{t} x _ {j , 1} ^ {*}\right) ^ {b _ {1 1}} \left(\frac {1}{t} x _ {j , 2} ^ {*}\right) ^ {b _ {1 2}}} \left(\frac {1}{t} x _ {i , 1} ^ {*}\right) ^ {\frac {a}{a + 1}} \left(\frac {1}{t} x _ {i , 2} ^ {*}\right) ^ {- \frac {1}{a + 1}}}\left(\sum_ {i} \left(D _ {i , 1} \frac {C _ {i , 1} \left(\tilde {x} _ {i , 1} ^ {*}\right) ^ {b _ {1 1}} \left(\tilde {x} _ {i , 2} ^ {*}\right) ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} \left(\tilde {x} _ {j , 1} ^ {*}\right) ^ {b _ {1 1}} \left(\tilde {x} _ {j , 2} ^ {*}\right) ^ {b _ {1 2}}} \left(\tilde {x} _ {i , 1} ^ {*}\right) ^ {\frac {a}{a + 1}} \left(\tilde {x} _ {i , 2} ^ {*}\right) ^ {- \frac {1}{a + 1}} + \sum_ {j} K _ {i j} \left(\tilde {x} _ {j , 1} ^ {*}\right) ^ {\frac {1}{a + 1}} \left(\tilde {x} _ {j , 2} ^ {*}\right) ^ {- \frac {a}{a + 1}}\right) + 1\right) ^ {c _ {1 1}} +
$$

$$
\frac {\sum_ {j} K _ {i j} \left(\frac {1}{t} x _ {j , 1} ^ {*}\right) ^ {\frac {1}{a + 1}} \left(\frac {1}{t} x _ {j , 2} ^ {*}\right) ^ {- \frac {a}{a + 1}}}{\left(\sum_ {i} \left(D _ {i , 1} \frac {C _ {i , 1} \left(\tilde {x} _ {i , 1} ^ {*}\right) ^ {b _ {1 1}} \left(\tilde {x} _ {i , 2} ^ {*}\right) ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} \left(\tilde {x} _ {j , 1} ^ {*}\right) ^ {b _ {1 1}} \left(\tilde {x} _ {j , 2} ^ {*}\right) ^ {b _ {1 2}}} \left(\tilde {x} _ {i , 1} ^ {*}\right) ^ {\frac {a}{a + 1}} \left(\tilde {x} _ {i , 2} ^ {*}\right) ^ {- \frac {1}{a + 1}} + \sum_ {j} K _ {i j} \left(\tilde {x} _ {j , 1} ^ {*}\right) ^ {\frac {1}{a + 1}} \left(\tilde {x} _ {j , 2} ^ {*}\right) ^ {- \frac {a}{a + 1}}\right) + 1\right) ^ {c _ {1 2}}} \Longleftrightarrow
$$

$$
\begin{array}{l} x _ {i, 1} ^ {*} = t ^ {\frac {2}{a + 1} - \frac {2 a}{a + 1} c _ {1 1}} D _ {i, 1} \frac {C _ {i , 1} \left(x _ {i , 1} ^ {*}\right) ^ {b _ {1 1}} \left(x _ {i , 2} ^ {*}\right) ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} \left(x _ {j , 1} ^ {*}\right) ^ {b _ {1 1}} \left(x _ {j , 2} ^ {*}\right) ^ {b _ {1 2}}} \left(x _ {i, 1} ^ {*}\right) ^ {\frac {a}{a + 1}} \left(x _ {i, 2} ^ {*}\right) ^ {- \frac {1}{a + 1}} \\ + t ^ {\frac {2 a}{a + 1} - \frac {2 a}{a + 1} c _ {1 2}} \sum_ {j} K _ {i j} \left(x _ {j, 1} ^ {*}\right) ^ {\frac {1}{a + 1}} \left(x _ {j, 2} ^ {*}\right) ^ {- \frac {a}{a + 1}} \iff \\ \end{array}
$$

$$
x _ {i, 1} ^ {*} = D _ {i, 1} \frac {C _ {i , 1} \left(x _ {i , 1} ^ {*}\right) ^ {b _ {1 1}} \left(x _ {i , 2} ^ {*}\right) ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} \left(x _ {j , 1} ^ {*}\right) ^ {b _ {1 1}} \left(x _ {j , 2} ^ {*}\right) ^ {b _ {1 2}}} \left(x _ {i, 1} ^ {*}\right) ^ {\frac {a}{a + 1}} \left(x _ {i, 2} ^ {*}\right) ^ {- \frac {1}{a + 1}} + \sum_ {j} K _ {i j} \left(x _ {j, 1} ^ {*}\right) ^ {\frac {1}{a + 1}} \left(x _ {j, 2} ^ {*}\right) ^ {- \frac {a}{a + 1}},
$$

i.e. $x _ { 1 } ^ { * } = \ddot { T } _ { 1 } \left( x ^ { * } \right)$ holds as claimed. 

Similarly, considering $\tilde { x } _ { 2 } ^ { * } = \tilde { T } _ { 2 } \left( \tilde { x } ^ { * } \right)$ and imposing (C.22) and (C.24) yields: 

$$
\begin{array}{l} \frac {1}{t} x _ {i, 2} ^ {*} = \frac {D _ {i , 2} \frac {C _ {i , 2} \left(\frac {1}{t} x _ {i , 1} ^ {*}\right) ^ {b _ {2 1}} \left(\frac {1}{t} x _ {i , 2} ^ {*}\right) ^ {b _ {2 2}}}{\sum_ {j} C _ {j , 2} \left(\frac {1}{t} x _ {j , 1} ^ {*}\right) ^ {b _ {2 1}} \left(\frac {1}{t} x _ {j , 2} ^ {*}\right) ^ {b _ {2 2}}} \left(\frac {1}{t} x _ {i , 1} ^ {*}\right) ^ {- \frac {1}{a + 1}} \left(\frac {1}{t} x _ {i , 2} ^ {*}\right) ^ {\frac {a}{a + 1}}}\left(\sum_ {i} \left(D _ {i , 2} \frac {C _ {i , 2} \left(\tilde {x} _ {i , 1} ^ {*}\right) ^ {b _ {2 1}} \left(\tilde {x} _ {i , 2} ^ {*}\right) ^ {b _ {2 2}}}{\sum_ {j} C _ {j , 2} \left(\tilde {x} _ {i , 1} ^ {*}\right) ^ {b _ {2 1}} \left(\tilde {x} _ {i , 2} ^ {*}\right) ^ {b _ {2 2}}} \left(\tilde {x} _ {i , 1} ^ {*}\right) ^ {- \frac {1}{a + 1}} \left(\tilde {x} _ {i , 2} ^ {*}\right) ^ {\frac {a}{a + 1}} + \sum_ {j} K _ {j i} \left(\tilde {x} _ {j , 1} ^ {*}\right) ^ {- \frac {a}{a + 1}} \left(\tilde {x} _ {j , 2} ^ {*}\right) ^ {\frac {1}{a + 1}}\right) + 1\right) ^ {c _ {2 1}} + \\ \frac {\sum_ {j} K _ {j i} \left(\frac {1}{t} x _ {j , 1} ^ {*}\right) ^ {- \frac {a}{a + 1}} \left(\frac {1}{t} x _ {j , 2} ^ {*}\right) ^ {\frac {1}{a + 1}}}{\left(\sum_ {i} \left(D _ {i , 2} \frac {C _ {i , 2} \left(\tilde {x} _ {i , 1} ^ {*}\right) ^ {b _ {2 1}} \left(\tilde {x} _ {i , 2} ^ {*}\right) ^ {b _ {2 2}}}{\sum_ {j} C _ {j , 2} \left(\tilde {x} _ {i , 1} ^ {*}\right) ^ {b _ {2 1}} \left(\tilde {x} _ {i , 2} ^ {*}\right) ^ {b _ {2 2}}} \left(\tilde {x} _ {i , 1} ^ {*}\right) ^ {- \frac {1}{a + 1}} \left(\tilde {x} _ {i , 2} ^ {*}\right) ^ {\frac {a}{a + 1}} + \sum_ {j} K _ {j i} \left(\tilde {x} _ {j , 1} ^ {*}\right) ^ {- \frac {a}{a + 1}} \left(\tilde {x} _ {j , 2} ^ {*}\right) ^ {\frac {1}{a + 1}}\right) + 1\right) ^ {c _ {2 2}}} \Longleftrightarrow \\ \end{array}
$$

$$
\begin{array}{l} x _ {i, 2} ^ {*} = t ^ {\frac {2}{a + 1} - \frac {2 a}{a + 1} c _ {2 1}} D _ {i, 2} \frac {C _ {i , 2} \left(x _ {i , 1} ^ {*}\right) ^ {b _ {2 1}} \left(x _ {i , 2} ^ {*}\right) ^ {b _ {2 2}}}{\sum_ {j} C _ {j , 2} \left(x _ {j , 1} ^ {*}\right) ^ {b _ {2 1}} \left(x _ {j , 2} ^ {*}\right) ^ {b _ {2 2}}} \left(x _ {i, 1} ^ {*}\right) ^ {- \frac {1}{a + 1}} \left(x _ {i, 2} ^ {*}\right) ^ {\frac {a}{a + 1}} + \\ t ^ {\frac {2 a}{a + 1} - \frac {2 a}{a + 1} c _ {2 2}} \frac {\sum_ {j} K _ {j i} \left(x _ {j , 1} ^ {*}\right) ^ {- \frac {a}{a + 1}} \left(x _ {j , 2} ^ {*}\right) ^ {\frac {1}{a + 1}}}{\sum_ {i} \left(D _ {i , 2} \frac {C _ {i , 2} \left(\tilde {x} _ {i , 1} ^ {*}\right) ^ {b _ {2 1}} \left(\tilde {x} _ {i , 2} ^ {*}\right) ^ {b _ {2 2}}}{\sum_ {j} C _ {j , 2} \left(\tilde {x} _ {i , 1} ^ {*}\right) ^ {b _ {2 1}} \left(\tilde {x} _ {i , 2} ^ {*}\right) ^ {b _ {2 2}}} \left(\tilde {x} _ {i , 1} ^ {*}\right) ^ {- \frac {1}{a + 1}} \left(\tilde {x} _ {i , 2} ^ {*}\right) ^ {\frac {a}{a + 1}} + \sum_ {j} K _ {j i} \left(\tilde {x} _ {j , 1} ^ {*}\right) ^ {- \frac {a}{a + 1}} \left(\tilde {x} _ {j , 2} ^ {*}\right) ^ {\frac {1}{a + 1}}\right)} \Longleftrightarrow \\ x _ {i, 2} ^ {*} = D _ {i, 2} \frac {C _ {i , 2} \left(x _ {i , 1} ^ {*}\right) ^ {b _ {2 1}} \left(x _ {i , 2} ^ {*}\right) ^ {b _ {2 2}}}{\sum_ {j} C _ {j , 2} \left(x _ {j , 1} ^ {*}\right) ^ {b _ {2 1}} \left(x _ {j , 2} ^ {*}\right) ^ {b _ {2 2}}} \left(x _ {i, 1} ^ {*}\right) ^ {- \frac {1}{a + 1}} \left(x _ {i, 2} ^ {*}\right) ^ {\frac {a}{a + 1}} + \sum_ {j} K _ {j i} \left(x _ {j, 1} ^ {*}\right) ^ {- \frac {a}{a + 1}} \left(x _ {j, 2} ^ {*}\right) ^ {\frac {1}{a + 1}}, \\ \end{array}
$$

i.e. $x _ { 2 } ^ { * } = \tilde { T } _ { 2 } \left( x ^ { * } \right)$ as well. Together, this means we have successfully constructed a fixed point $x ^ { * }$ of the (un-scaled) operator $T$ , i.e. $x ^ { * } = T \left( x ^ { * } \right)$ , i.e. there exists a solution to the system characterized by equations (C.1) and (C.2). 

Finally, we show that the equilibrium we have found is strictly positive, i.e. $x _ { i , 1 } > 0$ for all $i \in \mathcal N$ and $l \in \{ 1 , 2 \}$ . (Note that all $x _ { i , l } ^ { * }$ are trivially finite since $t \in ( 0 , \infty )$ and $\tilde { x } ^ { \ast } \in \mathcal { M }$ , so that $x _ { i , l } ^ { * } \leq t M _ { l } < \infty$ ). We proceed by contradiction: suppose not. Then there exists an $i \in \mathcal N$ and $l \in \{ 1 , 2 \}$ such that $x _ { i , l } ^ { * } = 0$ . Suppose that ${ \mathit { l } } = 1$ (the case of $l = 2$ proceeds analogously). We then have: 

$$
x _ {i, 1} ^ {*} = D _ {i, 1} \frac {C _ {i , 1} \left(x _ {i , 1} ^ {*}\right) ^ {b _ {1 1}} \left(x _ {i , 2} ^ {*}\right) ^ {b _ {1 2}}}{\sum_ {j} C _ {j , 1} \left(x _ {j , 1} ^ {*}\right) ^ {b _ {1 1}} \left(x _ {j , 2} ^ {*}\right) ^ {b _ {1 2}}} \left(x _ {i, 1} ^ {*}\right) ^ {\frac {a}{a + 1}} \left(x _ {i, 2} ^ {*}\right) ^ {- \frac {1}{a + 1}} + \sum_ {j} K _ {i j} \left(x _ {j, 1} ^ {*}\right) ^ {\frac {1}{a + 1}} \left(x _ {j, 2} ^ {*}\right) ^ {- \frac {a}{a + 1}} = 0
$$

Since $\{ x _ { j , 1 } ^ { * } \}$ are weakly positive and finite, since $K _ { i j } \geq 0$ , and $a \in [ 0 , 1 ] ,$ for this to be true, this means that for all $j$ connected to $i$ (i.e. the set of $j$ such that $K _ { i j } > 0$ ), it must be the case that $x _ { j , 1 } ^ { * } = 0$ as well. The same argument implies that all $j ^ { \prime }$ connected to any of these $j$ also have $x _ { j ^ { \prime } , 1 } ^ { * } = 0$ . Since the matrix $\mathbf { K } = [ K _ { i j } ]$ is connected, repeating this argument iteratively then implies that $x _ { j , 1 } ^ { * } = 0$ for all $j$ . But if $x _ { j , 1 } ^ { * } = 0$ for all $j$ , then from equation (C.2) $\boldsymbol { x } _ { j , 2 } ^ { * }$ would be infinite for all $j$ , which is a contradiction, as $x _ { i , 2 } ^ { * } \leq t M _ { 2 }$ . Hence, the equilibrium is strictly positive. 

# C.1.3 Part 2 (Uniqueness)

We now proceed to study the uniqueness of the system. We consider first an economic geography with a symmetric infrastructure matrix; we then consider an urban model. 

Part 2(a): Uniqueness in an Economic Geography Model with Symmetric Infrastructure Matrix Consider an economic geography model with a symmetric infrastructure matrix, i.e. $t _ { k l } = t _ { l k }$ for all $l , k \in \mathcal { N } \times \mathcal { N }$ . It is well known (see e.g. Anderson and Van Wincoop (2003); Allen, Arkolakis, and Takahashi (2020)) that with symmetric transportation costs, the market access terms are equal up to scale). It turns out that this is also true about a symmetric infrastructure matrix in the presence of endogenous traffic congestion. To see this, note a symmetric infrastructure matrix implies $K _ { i j } = K _ { j i }$ for all $i$ and $j$ , so that equations (C.6) and (C.7) can be written in the economic geography case as: 

$$
p _ {i} ^ {a} \pi_ {i} = Y _ {i} p _ {i} ^ {- (1 - a)} + \sum_ {j} K _ {i j} \pi_ {j} ^ {(1 - a)} \tag {C.25}
$$

$$
p _ {i} \pi_ {i} ^ {a} = Y _ {i} \pi_ {i} ^ {- (1 - a)} + \sum_ {j} K _ {i j} p _ {j} ^ {(1 - a)}. \tag {C.26}
$$

It is obvious by inspection of equations (C.25) and (C.26) that a solution to this system is $\pi _ { i } = \kappa p _ { i } $ for some $\kappa > 0$ . (Indeed, one can show using the same tools applied to the urban model in the next section that this is the unique solution). As a result, it is sufficient to focus on the solution of the single equation case, which can be written as follows: 

$$
\kappa^ {a} p _ {i} ^ {1 + a} = \kappa^ {a - 1} Y _ {i} p _ {i} ^ {- (1 - a)} + \sum_ {j} K _ {i j} p _ {j} ^ {(1 - a)}
$$

Recall from Section C.1.1 that we can write: 

$$
Y _ {i} = p _ {i} ^ {\frac {(1 + \alpha)}{1 - \beta - \theta (\alpha + \beta)}} \pi_ {i} ^ {\frac {(1 - \beta)}{1 - \beta - \theta (\alpha + \beta)}} \left(\frac {\bar {u} _ {i}}{W}\right) ^ {\theta \frac {(1 + \alpha)}{1 - \beta - \theta (\alpha + \beta)}} \bar {A} _ {i} ^ {\theta \frac {(1 - \beta)}{1 - \beta - \theta (\alpha + \beta)}}.
$$

so that the full system becomes: 

$$
\kappa^ {a} p _ {i} ^ {1 + a} = \kappa^ {(a - 1)} \bar {W} ^ {- \theta \frac {(1 + \alpha)}{1 - \beta - \theta (\alpha + \beta)}} \bar {u} _ {i} ^ {\theta \frac {(1 + \alpha)}{1 - \beta - \theta (\alpha + \beta)}} \bar {A} _ {i} ^ {\theta \frac {(1 - \beta)}{1 - \beta - \theta (\alpha + \beta)}} p _ {i} ^ {\frac {(1 + \alpha) + (1 - \beta)}{1 - \beta - \theta (\alpha + \beta)} - (1 - a)} p _ {i} ^ {- (1 - a)} + \sum_ {j} K _ {i j} p _ {j} ^ {(1 - a)}.
$$

Define s 1 $x _ { i } \equiv \kappa ^ { a } p _ { i } ^ { 1 + a }$ o that pi = κ− a1+a x 1+i $p _ { i } = \kappa ^ { - \frac { a } { 1 + a } } x _ { i } ^ { \frac { 1 } { 1 + a } }$ a so that we can write: 

$$
x _ {i} = \kappa^ {- \frac {a}{1 + a} \frac {(1 + \alpha) + (1 - \beta)}{1 - \beta - \theta (\alpha + \beta)} - \left(\frac {1 - a}{1 + a}\right)} \bar {W} ^ {- \theta} \frac {(1 + \alpha)}{1 - \beta - \theta (\alpha + \beta)} \bar {u} _ {i} ^ {\theta} \frac {(1 + \alpha)}{1 - \beta - \theta (\alpha + \beta)} \bar {A} _ {i} ^ {\theta} \frac {(1 - \beta)}{1 - \beta - \theta (\alpha + \beta)} x _ {i} ^ {\frac {1}{1 + a} \frac {(1 + \alpha) + (1 - \beta)}{1 - \beta - \theta (\alpha + \beta)}} x _ {i} ^ {- \left(\frac {1 - a}{1 + a}\right)} + \kappa^ {- \frac {a (1 - a)}{1 + a}} \sum_ {j} K _ {i j} x _ {j} ^ {\frac {1 - a}{1 + a}},
$$

or, written slightly more succinctly: 

$$
x _ {i} = \kappa^ {\rho_ {1}} \bar {W} ^ {\rho_ {2}} C _ {i} x _ {i} ^ {b - \left(\frac {1 - a}{1 + a}\right)} + \kappa^ {\rho_ {3}} \sum_ {j} K _ {i j} x _ {j} ^ {\frac {1 - a}{1 + a}}, \tag {C.27}
$$

where $\begin{array} { r } { \rho _ { 1 } \equiv - \frac { a } { 1 + a } \frac { ( 1 + \alpha ) + ( 1 - \beta ) } { 1 - \beta - \theta ( \alpha + \beta ) } - \left( \frac { 1 - a } { 1 + a } \right) } \end{array}$ , $\begin{array} { r } { \rho _ { 2 } \equiv - \theta \frac { ( 1 + \alpha ) } { 1 - \beta - \theta ( \alpha + \beta ) } } \end{array}$ θ 1−β−θ(α+β) , ρ3 ≡ (1+α) $\begin{array} { r } { \rho _ { 3 } \equiv - \frac { a ( 1 - a ) } { 1 + a } , b \equiv \frac { 1 } { 1 + a } \frac { ( 1 + \alpha ) + ( 1 - \beta ) } { 1 - \beta - \theta ( \alpha + \beta ) } } \end{array}$ , and $C _ { i } \equivq$ $\bar { u } _ { i } ^ { \theta \frac { ( 1 + \alpha ) } { 1 - \beta - \theta ( \alpha + \beta ) } } \bar { A } _ { i } ^ { \theta \frac { ( 1 - \beta ) } { 1 - \beta - \theta ( \alpha + \beta ) } }$ u¯θ 1−β−θ(α+β) i (1+α) A A¯θ 1−β−θ(α+β) . (1−β) ， 

We proceed by contradiction. Suppose there are two solutions $\left\{ W _ { y } , \kappa _ { y } , \left\{ y _ { i } \right\} _ { i \in \mathcal { N } } \right\}$ and $\left\{ W _ { x } , \kappa _ { x } , \left\{ x _ { i } \right\} _ { i \in \mathcal { N } } \right\}$ that both solve the system of equations (C.27) and that the two equations are distinct, i.e. they are not equal up to scale. To proceed, we take ratios of the two assumed solutions: 

$$
\frac {y _ {i}}{x _ {i}} = \frac {\kappa_ {y} ^ {\rho_ {1}} \bar {W} _ {y} ^ {\rho_ {2}} C _ {i} y _ {i} ^ {b - \left(\frac {1 - a}{1 + a}\right)} + \kappa_ {y} ^ {\rho_ {3}} \sum_ {j} K _ {i j} y _ {j} ^ {\frac {1 - a}{1 + a}}}{\kappa_ {x} ^ {\rho_ {1}} \bar {W} _ {x} ^ {\rho_ {2}} C _ {i} x _ {i} ^ {b - \left(\frac {1 - a}{1 + a}\right)} + \kappa_ {x} ^ {\rho_ {3}} \sum_ {j} K _ {i j} x _ {j} ^ {\frac {1 - a}{1 + a}}} \iff
$$

$$
\hat {x} _ {i} = D _ {i} \hat {\kappa} ^ {\rho_ {1}} \hat {W} ^ {\rho_ {2}} \hat {x} _ {i} ^ {b - \left(\frac {1 - a}{1 + a}\right)} + \sum_ {j} F _ {i j} \hat {\kappa} ^ {\rho_ {3}} \hat {x} _ {j} ^ {\frac {1 - a}{1 + a}}, \tag {C.28}
$$

where $\begin{array} { r } { \hat { x } _ { i } \equiv \frac { y _ { i } } { x _ { i } } } \end{array}$ $\begin{array} { r } { \hat { \kappa } \equiv \frac { \kappa _ { y } } { \kappa _ { x } } } \end{array}$ $\begin{array} { r } { \hat { W } \equiv \frac { W _ { y } } { W _ { x } } } \end{array}$ Wy $\begin{array} { r } { D _ { i } \equiv \frac { \kappa _ { x } ^ { \rho _ { 1 } } \bar { W _ { x } } ^ { \rho _ { 2 } } C _ { i } x _ { i } ^ { b - \left( \frac { 1 - a } { 1 + a } \right) } } { \kappa _ { x } ^ { \rho _ { 1 } } \bar { W _ { x } } ^ { \rho _ { 2 } } C _ { i } x _ { i } ^ { b - \left( \frac { 1 - a } { 1 + a } \right) } + \kappa _ { x } ^ { \rho _ { 3 } } \sum _ { j } K _ { i j } x _ { j } ^ { \frac { 1 - a } { 1 + a } } } } \end{array}$ −( 1−a1+a )+κρ3x Pj Kij x 11j , and $\begin{array} { r } { F _ { i j } \equiv \frac { \kappa _ { x } ^ { \rho _ { 3 } } K _ { i j } x _ { j } ^ { \frac { 1 - a } { 1 + a } } } { \kappa _ { x } ^ { \rho _ { 1 } } \bar { W _ { x } } ^ { \rho _ { 2 } } C _ { i } x _ { i } ^ { b - \left( \frac { 1 - a } { 1 + a } \right) } + \kappa _ { x } ^ { \rho _ { 3 } } \sum _ { j } K _ { i j } x _ { j } ^ { \frac { 1 - a } { 1 + a } } } . } \end{array}$ ( 1−a1+a )+κρ3x Pj Kij x 1−1+j 

We now construct a maximum and minimum bound to equation (C.28). Define $M \equiv \operatorname* { m a x } _ { i } { \hat { x } } _ { i }$ , m ≡ $\operatorname* { m i n } _ { i } { \hat { x } } _ { i }$ , and $\mu \equiv M / m$ . Note that because $a \in [ 0 , 1 ] ,$ we have: 

$$
\hat {x} _ {i} \leq D _ {i} \hat {\kappa} ^ {\rho_ {1}} \hat {W} ^ {\rho_ {2}} \hat {x} _ {i} ^ {b - \left(\frac {1 - a}{1 + a}\right)} + \sum_ {j} F _ {i j} \hat {\kappa} ^ {\rho_ {3}} \max  _ {j} x _ {j} ^ {\frac {1 - a}{1 + a}}. \tag {C.29}
$$

Indeed (and importantly), the inequality is strict. How do you see this? Well, suppose the inequality is not strict. Then there exists an $i \in S$ such that: 

$$
\hat {x} _ {i, 1} = D _ {i} \hat {\kappa} ^ {\rho_ {1}} \hat {W} ^ {\rho_ {2}} \hat {x} _ {i} ^ {b - \left(\frac {1 - a}{1 + a}\right)} + \sum_ {j} F _ {i j} \hat {\kappa} ^ {\rho_ {3}} \max _ {j} x _ {j} ^ {\frac {1 - a}{1 + a}}.
$$

For this to be true, it must be the case that for all $j$ such that $F _ { i j } > 0$ that $\hat { x } _ { j } ^ { \frac { 1 } { a + 1 } } = \operatorname* { m a x } _ { i } \hat { x } _ { j } ^ { \frac { 1 } { a + 1 } }$ 1+1 = maxi xˆ 1a+1j , ony of these j as ou r equivalently, ${ \hat { x } } _ { j } = \operatorname* { m a x } _ { j } { \hat { x } } _ { j } = { \hat { x } } _ { i }$ . Since all the locations are connected, we can choose a $j$ r new $i$ , they must also satisfy the equation with equality, and the argument can continue to the point that we have $\hat { x } _ { i } = \hat { x } _ { j }$ for all $_ i$ and $j$ . But this is a contradiction, since we assumed that there are two distinct solutions. Hence the inequality is strict. As a result, equation (C.29) implies: 

$$
M <   D _ {i} \hat {\kappa} ^ {\rho_ {1}} \hat {W} ^ {\rho_ {2}} M ^ {1 ^ {+} (b - (\frac {1 - a}{1 + a}))} m ^ {1 ^ {-} (b - (\frac {1 - a}{1 + a}))} + \sum_ {j} F _ {i j} \hat {\kappa} ^ {\rho_ {3}} M ^ {\frac {1 - a}{a + 1}}, \tag {C.30}
$$

where the (apologetically cumbersome) notation $1 ^ { + } \left( x \right) \equiv { \left\{ \begin{array} { l l } { x } & { { \mathrm { i f ~ } } x > 0 } \\ { 0 } & { { \mathrm { i f ~ } } x \leq 0 } \end{array} \right. }$ and $1 ^ { - } \left( x \right) \equiv { \left\{ \begin{array} { l l } { x } & { { \mathrm { i f ~ } } x < 0 } \\ { 0 } & { { \mathrm { i f ~ } } x \geq 0 } \end{array} \right. }$ , is necessary to consider multiple cases of the signs of the exponents at once. A similar argument can be made to establish the following lower bound, resulting in: 

$$
m > D _ {i} \hat {\kappa} ^ {\rho_ {1}} \hat {W} ^ {\rho_ {2}} m ^ {1 + (b - (\frac {1 - a}{1 + a}))} M ^ {1 - (b - (\frac {1 - a}{1 + a}))} + \sum_ {j} F _ {i j} \hat {\kappa} ^ {\rho_ {2}} m ^ {\frac {1 - a}{a + 1}}. \tag {C.31}
$$

Taking ratios of the upper to lower bounds achieves a strict upper bound on $\mu$ 

$$
\mu <   \frac {D _ {i} \hat {\kappa} ^ {\rho_ {1}} \hat {W} ^ {\rho_ {2}} M ^ {1 + (b - (\frac {1 - a}{1 + a}))} m ^ {1 - (b - (\frac {1 - a}{1 + a}))} + \sum_ {j} F _ {i j} \hat {\kappa} ^ {\rho_ {3}} M ^ {\frac {1 - a}{a + 1}}}{D _ {i} \hat {\kappa} ^ {\rho_ {1}} \hat {W} ^ {\rho_ {2}} m ^ {1 + (b - (\frac {1 - a}{1 + a}))} M ^ {1 - (b - (\frac {1 - a}{1 + a}))} + \sum_ {j} F _ {i j} \hat {\kappa} ^ {\rho_ {2}} m ^ {\frac {1 - a}{a + 1}}} \iff
$$

$$
\mu <   \beta \mu^ {\left| b - \left(\frac {1 - a}{1 + a}\right) \right|} + (1 - \beta) \mu^ {\frac {1 - a}{a + 1}}, \tag {C.32}
$$

where β ≡ ${ \cal A } \ = \ - \ \frac { D _ { i } \hat { \kappa } ^ { \rho _ { 1 } } \hat { W } ^ { \rho _ { 2 } } m ^ { 1 + } \big ( b - \big ( \frac { 1 - a } { 1 + a } \big ) \big ) _ { M ^ { 1 } } - \big ( b - \big ( \frac { 1 - a } { 1 + a } \big ) \big ) } { \big . \big . }$ Equation (C.32) says that is bounded $\begin{array} { r } { \rho = \frac { 1 } { D _ { i } \hat { \kappa } ^ { \rho _ { 1 } } \hat { W } ^ { \rho _ { 2 } } m ^ { 1 + \left( b - \left( \frac { 1 - a } { 1 + a } \right) \right) } M ^ { 1 - \left( b - \left( \frac { 1 - a } { 1 + a } \right) \right) } + \sum _ { j } F _ { i j } \hat { \kappa } ^ { \rho _ { 2 } } m ^ { \frac { 1 - a } { a + 1 } } } } \end{array}$ Diκˆρ1 Wˆ ρ2 m1+(b−( $\mu$ above by a weighted average of two terms. For this to be true, it must either be the case that $\mu$ is smaller than the first of the two terms or $\mu$ is smaller than the second of the two terms (or both). Since $\mu \geq 1$ , we then have a contradiction if both 

$$
\left| \frac {1 - a}{a + 1} \right| \leq 1,
$$

which is assured since $a \in [ 0 , 1 ]$ and: 

$$
\left| b - \left(\frac {1 - a}{1 + a}\right) \right| \leq 1 \tag {C.33}
$$

Hence a contradiction arises (and therefore uniqueness is assured) if equation (C.33) holds. Recall from above that a ≡ $\begin{array} { r } { a \equiv \frac { \theta \lambda } { 1 + \theta \lambda } } \end{array}$ and $\begin{array} { r } { b \equiv \frac { 1 } { 1 + a } \frac { ( 1 + \alpha ) + ( 1 - \beta ) } { 1 - \beta - \theta ( \alpha + \beta ) } } \end{array}$ . It is straightforward (but tedious) to verify that for all $\alpha \in [ - 1 , 1 ]$ , $\beta \in [ - 1 , 1 ]$ , $\lambda \geq 0$ , and $\theta \geq 0$ the following condition ensures equation (C.33) is satisfied an uniqueness is assured: 

$$
\alpha + \beta \leq 0,
$$

as claimed. 

Part 2(b): Uniqueness in an Urban Model From equations (C.17) and (C.18), we have that the equilibrium of the urban model satisfies the following system of equations: 

$$
x _ {i, 1} = W ^ {\rho_ {1}} C _ {i, 1} x _ {i, 1} ^ {\frac {a}{a + 1} + b _ {1 1}} x _ {i, 2} ^ {- \frac {1}{a + 1} + b _ {1 2}} + \sum_ {j} K _ {i j} x _ {j, 1} ^ {\frac {1}{a + 1}} x _ {j, 2} ^ {- \frac {a}{a + 1}} \tag {C.34}
$$

$$
x _ {i, 2} = W ^ {\rho_ {2}} C _ {i, 2} x _ {i, 1} ^ {- \frac {1}{a + 1} + b _ {2 1}} x _ {i, 2} ^ {\frac {a}{a + 1} + b _ {2 2}} + \sum_ {j} K _ {j i} x _ {j, 1} ^ {- \frac {a}{a + 1}} x _ {j, 2} ^ {\frac {1}{a + 1}}. \tag {C.35}
$$

As in the economic geography model, we prove uniqueness by contradiction. Suppose there exists a $\left\{ W _ { y } , \{ y _ { i , 1 } \} _ { i } , \{ y _ { i , 2 } \} _ { i } \right\}$ and $\left\{ W _ { x } , \{ x _ { i , 1 } \} _ { i } , \{ x _ { i , 2 } \} _ { i } \right\}$ that both solve equations (C.34) and (C.35), where solutions $\{ y _ { i , 1 } , y _ { i , 2 } \}$ and $\{ x _ { i , 1 } , x _ { i , 2 } \}$ are not equal up to scale. We will derive conditions under which this implies a contradiction. 

The first step is to re-write the system in terms of ratios of the two proposed solutions. Beginning with equation (C.34) , we have: 

$$
\frac {y _ {i , 1}}{x _ {i , 1}} = \frac {W _ {y} ^ {\rho_ {1}} C _ {i , 1} y _ {i , 1} ^ {\frac {a}{a + 1} + b _ {1 1}} y _ {i , 2} ^ {- \frac {1}{a + 1} + b _ {1 2}} + \sum_ {j} K _ {i j} y _ {j , 1} ^ {\frac {1}{a + 1}} y _ {j , 2} ^ {- \frac {a}{a + 1}}}{W _ {x} ^ {\rho_ {1}} C _ {i , 1} x _ {i , 1} ^ {\frac {a}{a + 1} + b _ {1 1}} x _ {i , 2} ^ {- \frac {1}{a + 1} + b _ {1 2}} + \sum_ {j} K _ {i j} x _ {j , 1} ^ {\frac {1}{a + 1}} x _ {j , 2} ^ {- \frac {a}{a + 1}}} \iff
$$

$$
\hat {x} _ {i, 1} = D _ {i, 1} \hat {W} ^ {\rho_ {1}} \left(\hat {x} _ {i, 1}\right) ^ {\frac {a}{a + 1} + b _ {1 1}} \left(\hat {x} _ {i, 2}\right) ^ {- \frac {1}{a + 1} + b _ {1 2}} + \sum_ {j} F _ {i j, 1} \left(\hat {x} _ {j, 1}\right) ^ {\frac {1}{a + 1}} \left(\hat {x} _ {j, 2}\right) ^ {- \frac {a}{a + 1}}, \tag {C.36}
$$

where ˆxi,1 ≡ $\begin{array} { r } { \hat { x } _ { i , 1 } \equiv \frac { y _ { i , 1 } } { x _ { i , 1 } } } \end{array}$ yi,1 , $\begin{array} { r } { \hat { W } \equiv \frac { W _ { y } } { W _ { x } } } \end{array}$ A Wy 

$$
D _ {i, 1} \equiv \frac {W _ {x} ^ {\rho_ {1}} C _ {i , 1} x _ {i , 1} ^ {\frac {a}{a + 1} + b _ {1 1}} x _ {i , 2} ^ {- \frac {1}{a + 1} + b _ {1 2}}}{W _ {x} ^ {\rho_ {1}} C _ {i , 1} x _ {i , 1} ^ {\frac {a}{a + 1} + b _ {1 1}} x _ {i , 2} ^ {- \frac {1}{a + 1} + b _ {1 2}} + \sum_ {j} K _ {i j} x _ {j , 1} ^ {\frac {1}{a + 1}} x _ {j , 2} ^ {- \frac {a}{a + 1}}},
$$

$$
F _ {i j, 1} \equiv \frac {K _ {i j}}{W _ {x} ^ {\rho_ {1}} C _ {i , 1} x _ {i , 1} ^ {\frac {a}{a + 1} + b _ {1 1}} x _ {i , 2} ^ {- \frac {1}{a + 1} + b _ {1 2}} + \sum_ {j} K _ {i j} x _ {j , 1} ^ {\frac {1}{a + 1}} x _ {j , 2} ^ {- \frac {a}{a + 1}}}.
$$

Note that $\begin{array} { r } { D _ { i , 1 } + \sum _ { j } F _ { i j , 1 } = 1 } \end{array}$ for all $i \in \mathcal N$ . Similarly, equation (C.35) can be written in ratios as follows: 

$$
\hat {x} _ {i, 2} = D _ {i, 2} \hat {W} ^ {\rho_ {2}} \left(\hat {x} _ {i, 1}\right) ^ {- \frac {1}{a + 1} + b _ {2 1}} \left(\hat {x} _ {i, 2}\right) ^ {\frac {a}{a + 1} + b _ {2 2}} + \sum_ {j} F _ {i j, 2} \left(\hat {x} _ {j, 1}\right) ^ {- \frac {a}{a + 1}} \left(\hat {x} _ {j, 2}\right) ^ {\frac {1}{a + 1}}, \tag {C.37}
$$

where $\begin{array} { r } { \hat { x } _ { i , 2 } \equiv \frac { y _ { i , 2 } } { x _ { i , 2 } } } \end{array}$ , 

$$
D _ {i, 2} \equiv \frac {W ^ {\rho_ {2}} C _ {i , 2} x _ {i , 1} ^ {- \frac {1}{a + 1} + b _ {2 1}} x _ {i , 2} ^ {\frac {a}{a + 1} + b _ {2 2}}}{W ^ {\rho_ {2}} C _ {i , 2} x _ {i , 1} ^ {- \frac {1}{a + 1} + b _ {2 1}} x _ {i , 2} ^ {\frac {a}{a + 1} + b _ {2 2}} + \sum_ {j} K _ {j i} x _ {j , 1} ^ {- \frac {a}{a + 1}} x _ {j , 2} ^ {\frac {1}{a + 1}}},
$$

$$
F _ {i j, 2} \equiv \frac {K _ {j i}}{W ^ {\rho_ {2}} C _ {i , 2} x _ {i , 1} ^ {- \frac {1}{a + 1} + b _ {2 1}} x _ {i , 2} ^ {\frac {a}{a + 1} + b _ {2 2}}} + \sum_ {j} K _ {j i} x _ {j, 1} ^ {- \frac {a}{a + 1}} x _ {j, 2} ^ {\frac {1}{a + 1}}.
$$

Note that $\begin{array} { r } { D _ { i , 2 } + \sum _ { j } F _ { i j , 2 } = 1 } \end{array}$ for all $i \in \mathcal N$ . 

In what follows, we focus on the system of equations (C.36) and (C.37) written in ratios. Because we are assuming the solutions $\{ y _ { i , 1 } , y _ { i , 2 } \}$ and $\{ x _ { i , 1 } , x _ { i , 2 } \}$ are not equal up to scale, it must be the case that there is at least one $\hat { x } _ { i , l } \neq 1$ . 

To proceed, we start bounding the ratios of the solutions. Define $M _ { l } \equiv \operatorname* { m a x } _ { i } { \hat { x } } _ { i , l }$ , $m _ { l } \equiv \operatorname* { m i n } _ { i } { \hat { x } } _ { i , l }$ , and $\begin{array} { r } { \mu _ { l } \equiv \frac { M _ { l } } { m _ { l } } } \end{array}$ . From equation (C.36) we have for all $i \in \mathcal N$ : 

$$
\hat {x} _ {i, 1} \leq D _ {i, 1} \hat {W} ^ {\rho_ {1}} \left(\hat {x} _ {i, 1}\right) ^ {\frac {a}{a + 1} + b _ {1 1}} \left(\hat {x} _ {i, 2}\right) ^ {- \frac {1}{a + 1} + b _ {1 2}} + \sum_ {j} F _ {i j, 1} \max  _ {i} \left(\left(\hat {x} _ {j, 1}\right) ^ {\frac {1}{a + 1}}\right) \max  _ {j} \left(\left(\hat {x} _ {j, 2}\right) ^ {- \frac {a}{a + 1}}\right). \tag {C.38}
$$

Indeed (and importantly), the inequality is strict. How do you see this? Well, suppose the inequality is not strict. Then there exists an $i \in S$ such that: 

$$
\hat {x} _ {i, 1} = D _ {i, 1} \hat {W} ^ {\rho_ {1}} (\hat {x} _ {i, 1}) ^ {\frac {a}{a + 1} + b _ {1 1}} (\hat {x} _ {i, 2}) ^ {- \frac {1}{a + 1} + b _ {1 2}} + \sum_ {j} F _ {i j, 1} \max _ {i} \left((\hat {x} _ {j, 1}) ^ {\frac {1}{a + 1}}\right) \max _ {j} \left((\hat {x} _ {j, 2}) ^ {- \frac {a}{a + 1}}\right).
$$

For this to be true, it must be the case that for all $j$ such that $F _ { i j , 1 } > 0$ that $( \hat { x } _ { j , 1 } ) ^ { \frac { 1 } { a + 1 } } = \operatorname* { m a x } _ { i } \left( ( \hat { x } _ { j , 1 } ) ^ { \frac { 1 } { a + 1 } } \right)$ , or equivalently, $\begin{array} { r } { \hat { x } _ { j , 1 } = \operatorname* { m a x } _ { j } \hat { x } _ { j , 1 } = \hat { x } _ { i , 1 } } \end{array}$ . Similarly, it must be the case that $\hat { x } _ { j , 2 } = \operatorname* { m i n } _ { j } \hat { x } _ { j , 2 }$ . Since all the locations are connected, we can choose any of these $j$ as our new $_ i$ , they must also satisfy the equation with equality, and the argument can continue to the point that we have $\hat { x } _ { i , 1 } = \hat { x } _ { j , 1 }$ and $\hat { x } _ { i , 2 } = \hat { x } _ { j , 2 }$ for all $i$ and $j$ . But this is a contradiction, since we assumed that there are two distinct solutions. Hence the inequality is strict. As a result, we can write equation (C.38) as implying: 

$$
M _ {1} <   D _ {i, 1} \hat {W} ^ {\rho_ {1}} M _ {1} ^ {1 ^ {+} \left(\frac {a}{a + 1} + b _ {1 1}\right)} m _ {1} ^ {1 ^ {-} \left(\frac {a}{a + 1} + b _ {1 1}\right)} M _ {2} ^ {1 ^ {+} \left(- \frac {1}{a + 1} + b _ {1 2}\right)} m _ {2} ^ {1 ^ {+} \left(- \frac {1}{a + 1} + b _ {1 2}\right)} + \sum_ {j} F _ {i j, 1} M _ {1} ^ {\frac {1}{a + 1}} m _ {2} ^ {- \frac {a}{a + 1}}, \tag {C.39}
$$

where the (apologetically cumbersome) notation $1 ^ { + } \left( x \right) \equiv { \left\{ \begin{array} { l l } { x } & { { \mathrm { i f ~ } } x > 0 } \\ { 0 } & { { \mathrm { i f ~ } } x \leq 0 } \end{array} \right. }$ and $1 ^ { - } \left( x \right) \equiv { \left\{ \begin{array} { l l } { x } & { { \mathrm { i f ~ } } x < 0 } \\ { 0 } & { { \mathrm { i f ~ } } x \geq 0 } \end{array} \right. }$ , is necessary to consider multiple cases of the signs of the exponents at once. 

We can proceed similarly for the minimum bound of equation (C.36), yielding: 

$$
m _ {1} > D _ {i, 1} \hat {W} ^ {\rho_ {1}} m _ {1} ^ {1 + \left(\frac {a}{a + 1} + b _ {1 1}\right)} M _ {1} ^ {1 - \left(\frac {a}{a + 1} + b _ {1 1}\right)} m _ {2} ^ {1 + \left(- \frac {1}{a + 1} + b _ {1 2}\right)} M _ {2} ^ {1 + \left(- \frac {1}{a + 1} + b _ {1 2}\right)} + \sum_ {j} F _ {i j, 1} m _ {1} ^ {\frac {1}{a + 1}} M _ {2} ^ {- \frac {a}{a + 1}}. \tag {C.40}
$$

Taking ratios of equations (C.39) and (C.40) yields: 

$$
\frac {M _ {1}}{m _ {1}} <   \frac {D _ {i , 1} \hat {W} ^ {\rho_ {1}} M _ {1} ^ {1 + \left(\frac {a}{a + 1} + b _ {1 1}\right)} m _ {1} ^ {1 - \left(\frac {a}{a + 1} + b _ {1 1}\right)} M _ {2} ^ {1 + \left(- \frac {1}{a + 1} + b _ {1 2}\right)} m _ {2} ^ {1 + \left(- \frac {1}{a + 1} + b _ {1 2}\right)} + \sum_ {j} F _ {i j , 1} M _ {1} ^ {\frac {1}{a + 1}} m _ {2} ^ {- \frac {a}{a + 1}}}{D _ {i , 1} \hat {W} ^ {\rho_ {1}} m _ {1} ^ {1 + \left(\frac {a}{a + 1} + b _ {1 1}\right)} M _ {1} ^ {1 - \left(\frac {a}{a + 1} + b _ {1 1}\right)} m _ {2} ^ {1 + \left(- \frac {1}{a + 1} + b _ {1 2}\right)} M _ {2} ^ {1 + \left(- \frac {1}{a + 1} + b _ {1 2}\right)} + \sum_ {j} F _ {i j , 1} m _ {1} ^ {\frac {1}{a + 1}} M _ {2} ^ {- \frac {a}{a + 1}}} \iff
$$

$$
\mu_ {1} <   \beta_ {1} \mu_ {1} ^ {\left| \frac {a}{a + 1} + b _ {1 1} \right|} \mu_ {2} ^ {\left| - \frac {1}{a + 1} + b _ {1 2} \right|} + (1 - \beta_ {1}) \mu_ {1} ^ {\frac {1}{a + 1}} \mu_ {2} ^ {\frac {a}{a + 1}}, \tag {C.41}
$$

where $\begin{array} { r } { \beta _ { 1 } \equiv \frac { D _ { i , 1 } \hat { W } ^ { \rho _ { 1 } } m _ { 1 } ^ { 1 + \left( \frac { a } { a + 1 } + b _ { 1 1 } \right) } M _ { 1 } ^ { 1 - \left( \frac { a } { a + 1 } + b _ { 1 1 } \right) } m _ { 2 } ^ { 1 + \left( - \frac { 1 } { a + 1 } + b _ { 1 2 } \right) } M _ { 2 } ^ { 1 + \left( - \frac { 1 } { a + 1 } + b _ { 1 2 } \right) } } { D _ { i , 1 } \hat { W } ^ { \rho _ { 1 } } m _ { 1 } ^ { 1 + \left( \frac { a } { a + 1 } + b _ { 1 1 } \right) } M _ { 1 } ^ { 1 - \left( \frac { a } { a + 1 } + b _ { 1 1 } \right) } m _ { 2 } ^ { 1 + \left( - \frac { 1 } { a + 1 } + b _ { 1 2 } \right) } M _ { 2 } ^ { 1 + \left( - \frac { 1 } { a + 1 } + b _ { 1 2 } \right) } + \sum _ { j } F _ { i j , 1 } m _ { 1 } ^ { \frac { 1 } { a + 1 } } M _ { 2 } ^ { - \frac { a } { a + 1 } } } . } \end{array}$ m 2 M 2 1+( m 2 

Proceeding identically for equation (C.37) yields the corresponding bound: 

$$
\mu_ {2} <   \beta_ {2} \mu_ {1} ^ {\left| - \frac {1}{a + 1} + b _ {2 1} \right|} \mu_ {2} ^ {\left| \frac {a}{a + 1} + b _ {2 2} \right|} + (1 - \beta_ {2}) \mu_ {1} ^ {\frac {a}{a + 1}} \mu_ {2} ^ {\frac {1}{a + 1}}, \tag {C.42}
$$

where $\delta _ { 2 } \equiv { \frac { D _ { i , 2 } \hat { W } ^ { \rho 2 } m _ { 1 } ^ { 1 } \left( - { \frac { 1 } { a + 1 } } + b _ { 2 1 } \right) _ { M _ { 1 } ^ { 1 } } \left( - { \frac { 1 } { a + 1 } } + b _ { 2 1 } \right) _ { m _ { 2 } ^ { 1 } } + \left( { \frac { a } { a + 1 } } + b _ { 2 2 } \right) _ { M _ { 2 } ^ { 1 } } \left( { \frac { a } { a + 1 } } + { \frac { b _ { 2 2 } } { a + 2 } } \right) } { \int _ { - \infty } ^ { \infty } { \frac { 1 } { a } } \left( - { \frac { 1 } { a + 1 } } + b _ { 2 2 } \right) _ { M _ { 1 } ^ { 1 } } \left( - { \frac { 1 } { a + 1 } } + b _ { 2 1 } \right) _ { M _ { 2 } ^ { 1 } } } } .$ $\begin{array} { r } { \beta _ { 2 } \equiv \frac { D _ { i , 2 } W ^ { \prime 2 } m _ { 1 } } { D _ { i , 2 } \hat { W } ^ { \prime 2 } m _ { 1 } ^ { 1 + } \left( - \frac { 1 } { a + 1 } + b _ { 2 1 } \right) _ { M _ { 1 } ^ { 1 } } ^ { - } \left( - \frac { 1 } { a + 1 } + b _ { 2 1 } \right) _ { m _ { 2 } ^ { 1 + } } \left( \frac { a } { a + 1 } + b _ { 2 2 } \right) _ { M _ { 2 } ^ { 1 } } ^ { - } \left( \frac { a } { a + 1 } + b _ { 2 2 } \right) _ { + \sum _ { j } F _ { i , 2 } M _ { 1 } ^ { - } } ^ { - } \frac { a } { a + 1 } _ { m _ { 2 } ^ { \frac { 1 } { a + 1 } } } ^ { - } } . } \end{array}$ 1− m 2 1+( M M1 m 2 M 2 m a2 

Both equations (C.41) and (C.42) bound $\mu _ { l }$ above by a weighted average of two terms. For each equation to be true, it must then be the case that either $\mu _ { l }$ is bounded above by the first term or it is bounded above by the second term (or both). Considering all possible combinations, for there not to be a contradiction, we require at least one of the following conditions to be true: 

$$
\left( \begin{array}{c} \ln \mu_ {1} \\ \ln \mu_ {2} \end{array} \right) <   \left( \begin{array}{c c} \left| \frac {a}{a + 1} + b _ {1 1} \right| & \left| - \frac {1}{a + 1} + b _ {1 2} \right| \\ \left| - \frac {1}{a + 1} + b _ {2 1} \right| & \left| \frac {a}{a + 1} + b _ {2 2} \right| \end{array} \right) \left( \begin{array}{c} \ln \mu_ {1} \\ \ln \mu_ {2} \end{array} \right)
$$

$$
\left( \begin{array}{c} \ln \mu_ {1} \\ \ln \mu_ {2} \end{array} \right) <   \left( \begin{array}{c c} \frac {1}{a + 1} & \frac {a}{a + 1} \\ \frac {a}{a + 1} & \frac {1}{a + 1} \end{array} \right) \left( \begin{array}{c} \ln \mu_ {1} \\ \ln \mu_ {2} \end{array} \right)
$$

$$
\left( \begin{array}{c} \ln \mu_ {1} \\ \ln \mu_ {2} \end{array} \right) <   \left( \begin{array}{c c} \left| \frac {a}{a + 1} + b _ {1 1} \right| & \left| - \frac {1}{a + 1} + b _ {1 2} \right| \\ \frac {a}{a + 1} & \frac {1}{a + 1} \end{array} \right) \left( \begin{array}{c} \ln \mu_ {1} \\ \ln \mu_ {2} \end{array} \right)
$$

$$
\left( \begin{array}{c} \ln \mu_ {1} \\ \ln \mu_ {2} \end{array} \right) <   \left( \begin{array}{c c} \frac {1}{a + 1} & \frac {a}{a + 1} \\ \left| - \frac {1}{a + 1} + b _ {2 1} \right| & \left| \frac {a}{a + 1} + b _ {2 2} \right| \end{array} \right) \left( \begin{array}{c} \ln \mu_ {1} \\ \ln \mu_ {2} \end{array} \right)
$$

By the Collatz-Wielandt formula, note the each inequality can hold only if its matrix has a spectral radius greater than one. By the Gershgorin circle theorem, if all row sums of a matrix is no greater than one, its spectral radius is also no greater than one. As a result, none of the four inequalities will hold if: 

$$
\left| \frac {a}{a + 1} + b _ {1 1} \right| + \left| - \frac {1}{a + 1} + b _ {1 2} \right| \leq 1 \tag {C.43}
$$

and: 

$$
\left| - \frac {1}{a + 1} + b _ {2 1} \right| + \left| \frac {a}{a + 1} + b _ {2 2} \right| \leq 1 \tag {C.44}
$$

Hence if both inequalities (C.43) and (C.44) are satisfied, then we have a contradiction, thereby establishing uniqueness. Recall that $\begin{array} { r } { a \ = \ \frac { \theta \lambda } { 1 + \theta \lambda } } \end{array}$ , $\begin{array} { r } { b _ { 1 1 } \equiv \frac { a } { a ^ { 2 } - 1 } \frac { 1 } { 1 - \theta \alpha } } \end{array}$ , $\begin{array} { r } { b _ { 1 2 } \equiv - \frac { 1 } { a ^ { 2 } - 1 } \frac { 1 } { 1 - \theta \alpha } } \end{array}$ , $\begin{array} { r } { b _ { 2 1 } \equiv - \frac { 1 } { a ^ { 2 } - 1 } \frac { 1 } { 1 - \theta \beta } } \end{array}$ , and $b _ { 2 2 } \equiv$ $\frac { a } { a ^ { 2 } - 1 } \frac { 1 } { 1 - \theta \beta }$ . Although the calculations are tedious, for all $\alpha \in [ 0 , 1 ]$ , $\beta \in [ 0 , 1 ]$ , $\theta \geq 0$ ,and $\lambda \geq 0$ , one can verify that inequality (C.43) holds if: 

$$
\alpha \leq \frac {1}{2} \left(\frac {1}{\theta} - \lambda\right).
$$

Similarly, one can verify that inequality (C.44) holds if: 

$$
\beta \leq \frac {1}{2} \left(\frac {1}{\theta} - \lambda\right),
$$

as required. 

# C.2 Proof of Proposition 2

In this section, we describe how to conduct counterfactuals (and relate the resulting equilibrium system to that of Proposition 1). We begin by summarizing the results that follow. 

Suppose the observed economy has infrastructure network $\mathbf { T } \equiv [ \bar { t } _ { k l } ]$ and is in equilibrium. To determine how the distribution of economic activity and the welfare will change under an alternative infrastructure 

network $\bar { \bf T } ^ { \prime } \equiv [ \bar { t } _ { k l } ^ { \prime } ]$ , we follow the “exact hat algebra” approach pioneered by Dekle, Eaton, and Kortum (2008), where we denote with hats the change in variables, e.g. $\hat { \bar { t } } _ { k l } \equiv \vec { t } _ { k l } / \bar { t } _ { k l }$ , $\begin{array} { r } { \hat { \gamma } _ { i } \equiv \frac { \gamma _ { i } ^ { \prime } } { \gamma _ { i } } } \end{array}$ . For the economic geography model, one can write the equilibrium system of equations (10) and (11) in changes as: 

$$
\hat {y} _ {i} ^ {\frac {1 + \theta \lambda + \theta}{1 + \theta \lambda}} \hat {l} _ {i} ^ {- \frac {\theta (1 + \alpha + \theta \lambda (\beta + \alpha))}{1 + \theta \lambda}} = \hat {\chi} \left(\frac {E _ {i}}{E _ {i} + \sum_ {k} \Xi_ {i k}}\right) \hat {y} _ {i} ^ {\frac {1 + \theta \lambda + \theta}{1 + \theta \lambda}} \hat {l} _ {i} ^ {\frac {\theta (\beta - 1)}{1 + \theta \lambda}} + \hat {\chi} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j} \left(\frac {\Xi_ {i j}}{E _ {i} + \sum_ {k} \Xi_ {i k}}\right) \hat {\tilde {t}} _ {i j} ^ {- \frac {\theta}{1 + \theta \lambda}} \hat {y} _ {j} ^ {\frac {1 + \theta}{1 + \theta \lambda}} \hat {l} _ {j} ^ {- \frac {\theta (1 + \alpha)}{1 + \theta \lambda}}
$$

$$
\hat {y} _ {i} ^ {\frac {- \theta (1 - \lambda)}{1 + \theta \lambda}} \hat {l} _ {i} ^ {\frac {\theta (1 - \beta - \theta \lambda (\alpha + \beta))}{1 + \theta \lambda}} = \hat {\chi} \left(\frac {Y _ {i}}{Y _ {i} + \sum_ {k} \Xi_ {k i}}\right) \hat {y} _ {i} ^ {\frac {- \theta (1 - \lambda)}{1 + \theta \lambda}} \hat {l} _ {i} ^ {\frac {\theta (\alpha + 1)}{1 + \theta \lambda}} + \hat {\chi} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j = 1} ^ {N} \left(\frac {\Xi_ {j i}}{Y _ {i} + \sum_ {k} \Xi_ {k i}}\right) \hat {\tilde {t}} _ {j i} ^ {- \frac {\theta}{1 + \theta \lambda}} \hat {y} _ {j} ^ {- \frac {\theta}{1 + \theta \lambda}} \hat {l} _ {j} ^ {\frac {\theta (1 - \beta)}{1 + \theta \lambda}}
$$

Similarly, for the urban model, the equilibrium system defined by equations (19) and (20) can be written in changes as: 

$$
\left(\hat {l} _ {i} ^ {R}\right) ^ {1 - \theta \beta} \left(\hat {l} _ {i} ^ {F}\right) ^ {\frac {\theta \lambda (1 - \alpha \theta)}{1 + \theta \lambda}} = \hat {\chi} \left(\frac {L _ {i} ^ {F}}{L _ {i} ^ {F} + \sum_ {k} \Xi_ {i k}}\right) \left(\hat {l} _ {i} ^ {F}\right) ^ {\frac {\theta (\alpha + \lambda)}{1 + \theta \lambda}} + \hat {\chi} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j = 1} ^ {N} \left(\frac {\Xi_ {i j}}{L _ {i} ^ {F} + \sum_ {k} \Xi_ {i k}}\right) \hat {\bar {t}} _ {i j} ^ {- \frac {\theta}{1 + \theta \lambda}} \left(\hat {l} _ {j} ^ {R}\right) ^ {\frac {1 - \theta \beta}{1 + \theta \lambda}}
$$

$$
\left(\hat {l} _ {i} ^ {R}\right) ^ {\frac {\theta \lambda (1 - \beta \theta)}{1 + \theta \lambda}} \left(\hat {l} _ {i} ^ {F}\right) ^ {1 - \theta \alpha} = \hat {\chi} \left(\frac {L _ {i} ^ {R}}{L _ {i} ^ {R} + \sum_ {k} \Xi_ {k i}}\right) \left(\hat {l} _ {i} ^ {R}\right) ^ {\frac {\theta (\beta + \lambda)}{1 + \theta \lambda}} + \hat {\chi} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j = 1} ^ {N} \left(\frac {\Xi_ {j i}}{L _ {i} ^ {R} + \sum_ {k} \Xi_ {k i}}\right) \hat {\tilde {t}} _ {j i} ^ {- \frac {\theta}{1 + \theta \lambda}} \left(\hat {l} _ {j} ^ {F}\right) ^ {\frac {1 - \theta \alpha}{1 + \theta \lambda}}
$$

As noted in the main text, both systems of equilibrium changes bear a close resemblance to their level variants above; the key distinction, however, is that the local geography and the infrastructure matrix are replaced with shares that depend only on the observed traffic flows and observed distribution of economic activity. As a result, it should come as no surprise that they inherit the mathematical properties of their level variants. 

We now turn to deriving both expressions. 

# C.2.1 Preliminaries

To simplify the derivation of counterfactual expressions for the trade and commuting models, we present a generalized version of models in term of market access terms and convert that to the “exact hat” notation developed in Dekle, Eaton, and Kortum (2008). Substituting in for the specific definitions of market access of each model returns the relevant counterfactual equilibrium conditions for both models. 

For both models, flows $X _ { i j }$ can be described in gravity form as: 

$$
X _ {i j} = \tau_ {i j} ^ {- \theta} \times \frac {\gamma_ {i}}{\Pi_ {i} ^ {- \theta}} \times \frac {\delta_ {j}}{P _ {j} ^ {- \theta}}
$$

where $\gamma _ { i }$ and $\delta _ { j }$ are cumulative flows out of an origin and into a destination, respectively, and $\Pi _ { i }$ and $P _ { j }$ are origin and destination market access terms. Trade costs can be represented as: 

$$
\left[ \tau_ {i j} ^ {- \theta} \right] = \left(I - \left[ t _ {i j} \right] ^ {- \theta}\right) ^ {- 1}
$$

And the iceberg cost of traveling along a link can be described as: 

$$
t _ {i j} = \bar {t} _ {i j} [ \Xi_ {i j} ] ^ {- \lambda} \Longleftrightarrow
$$

$$
t _ {i j} = \bar {t} _ {i j} \left[ t _ {i j} ^ {- \theta} \times P _ {i} ^ {- \theta} \times \Pi_ {j} ^ {- \theta} \right] ^ {- \lambda} \Longleftrightarrow
$$

$$
t _ {i j} = \bar {t} _ {i j} ^ {\frac {1}{1 + \theta \lambda}} \times P _ {i} ^ {- \frac {\theta \lambda}{1 + \theta \lambda}} \times \Pi_ {j} ^ {- \frac {\theta \lambda}{1 + \theta \lambda}}
$$

For both models, we have the equilibrium conditions: 

$$
\gamma_ {i} = \sum_ {j} X _ {i j}
$$

$$
\delta_ {i} = \sum_ {j} X _ {j i}
$$

Starting with the first equilibrium conditions, we substitute in for our gravity model of flows and solve for market access term $\Pi _ { i }$ : 

$$
\gamma_ {i} = \sum_ {j} X _ {i j} \iff
$$

$$
\gamma_ {i} = \sum_ {j} \tau_ {i j} ^ {- \theta} \times \frac {\gamma_ {i}}{\Pi_ {i} ^ {- \theta}} \times \frac {\delta_ {j}}{P _ {j} ^ {- \theta}} \Longleftrightarrow
$$

$$
\gamma_ {i} = \sum_ {j} \left(I - \left[ t _ {i j} \right] ^ {- \theta}\right) ^ {- 1} \times \frac {\gamma_ {i}}{\Pi_ {i} ^ {- \theta}} \times \frac {\delta_ {j}}{P _ {j} ^ {- \theta}} \Longleftrightarrow
$$

$$
\Pi_ {i} ^ {- \theta} = \sum_ {j} \left(I - \left[ t _ {i j} \right] ^ {- \theta}\right) ^ {- 1} \times \frac {\delta_ {j}}{P _ {j} ^ {- \theta}} \Longleftrightarrow
$$

$$
\left(I - \left[ t _ {i j} \right] ^ {- \theta}\right) \Pi_ {i} ^ {- \theta} = \frac {\delta_ {i}}{P _ {i} ^ {- \theta}} \iff
$$

$$
\Pi_ {i} ^ {- \theta} = \frac {\delta_ {i}}{P _ {i} ^ {- \theta}} + \sum_ {j} t _ {i j} ^ {- \theta} \Pi_ {j} ^ {- \theta}
$$

For the second equilibrium condition, we do the same, but solving for market access term $P _ { i }$ : 

$$
\delta_ {i} = \sum_ {j} X _ {j i} \iff
$$

$$
\delta_ {j} = \sum_ {j} \tau_ {j i} ^ {- \theta} \times \frac {\gamma_ {j}}{\Pi_ {j} ^ {- \theta}} \times \frac {\delta_ {i}}{P _ {i} ^ {- \theta}} \Longleftrightarrow
$$

$$
P _ {i} ^ {- \theta} = \sum_ {j} \tau_ {j i} ^ {- \theta} \times \frac {\gamma_ {j}}{\Pi_ {j} ^ {- \theta}} \Longleftrightarrow
$$

$$
P _ {i} ^ {- \theta} = \sum_ {j} \left(I - \left[ t _ {j i} \right] ^ {- \theta}\right) ^ {- 1} \times \frac {\gamma_ {j}}{\Pi_ {j} ^ {- \theta}} \Longleftrightarrow
$$

$$
\left(I - \left[ t _ {j i} \right] ^ {- \theta}\right) ^ {- 1} P _ {i} ^ {- \theta} = \frac {\gamma_ {i}}{\Pi_ {i} ^ {- \theta}} \iff
$$

$$
P _ {i} ^ {- \theta} = \frac {\gamma_ {i}}{\Pi_ {i} ^ {- \theta}} + \sum_ {j} t _ {j i} ^ {- \theta} P _ {j} ^ {- \theta}
$$

Expressed in changes, these two equilibrium conditions become 

$$
\hat {\Pi} _ {i} ^ {- \theta} = \left(\frac {\frac {\delta_ {i}}{P _ {i} ^ {- \theta}}}{\frac {\delta_ {i}}{P _ {i} ^ {- \theta}} + \sum_ {j} t _ {i j} ^ {- \theta} \Pi_ {j} ^ {- \theta}}\right) \frac {\hat {\delta} _ {i}}{\hat {P} _ {i} ^ {- \theta}} + \sum_ {j} \left(\frac {t _ {i j} ^ {- \theta} \Pi_ {j} ^ {- \theta}}{\frac {\delta_ {i}}{P _ {i} ^ {- \theta}} + \sum_ {j} t _ {i j} ^ {- \theta} \Pi_ {j} ^ {- \theta}}\right) \hat {t} _ {i j} ^ {- \theta} \hat {\Pi} _ {j} ^ {- \theta}
$$

and: 

$$
\hat {P} _ {i} ^ {- \theta} = \left(\frac {\frac {\gamma_ {i}}{\Pi_ {i} ^ {- \theta}}}{\frac {\gamma_ {i}}{\Pi_ {i} ^ {- \theta}} + \sum_ {j} t _ {j i} ^ {- \theta} P _ {j} ^ {- \theta}}\right) \frac {\hat {\gamma} _ {i}}{\hat {\Pi} _ {i} ^ {- \theta}} + \sum_ {j} \left(\frac {t _ {j i} ^ {- \theta} P _ {j} ^ {- \theta}}{\frac {\gamma_ {i}}{\Pi_ {i} ^ {- \theta}} + \sum_ {j} t _ {j i} ^ {- \theta} P _ {j} ^ {- \theta}}\right) \hat {t} _ {j i} ^ {- \theta} \hat {P} _ {j} ^ {- \theta}
$$

We multiply both the numerator and denominator by their appropriate market access term so that we 

can substitute in our expression for traffic, resulting in: 

$$
\hat {\Pi} _ {i} ^ {- \theta} = \left(\frac {\delta_ {i}}{\delta_ {i} + \sum_ {j} t _ {i j} ^ {- \theta} P _ {i} ^ {- \theta} \Pi_ {j} ^ {- \theta}}\right) \frac {\hat {\delta} _ {i}}{\hat {P} _ {i} ^ {- \theta}} + \sum_ {j} \left(\frac {t _ {i j} ^ {- \theta} P _ {i} ^ {- \theta} \Pi_ {j} ^ {- \theta}}{\delta_ {i} + \sum_ {j} t _ {i j} ^ {- \theta} P _ {i} ^ {- \theta} \Pi_ {j} ^ {- \theta}}\right) \hat {t} _ {i j} ^ {- \theta} \hat {\Pi} _ {j} ^ {- \theta} \iff
$$

$$
\hat {\Pi} _ {i} ^ {- \theta} = \left(\frac {\delta_ {i}}{\delta_ {i} + \sum_ {j} \Xi_ {i j}}\right) \frac {\hat {\delta} _ {i}}{\hat {P} _ {i} ^ {- \theta}} + \sum_ {j} \left(\frac {\Xi_ {i j}}{\delta_ {i} + \sum_ {j} \Xi_ {i j}}\right) \hat {t} _ {i j} ^ {- \theta} \hat {\Pi} _ {j} ^ {- \theta}
$$

and: 

$$
\hat {P} _ {i} ^ {- \theta} = \left(\frac {\gamma_ {i}}{\gamma_ {i} + \sum_ {j} t _ {j i} ^ {- \theta} \Pi_ {i} ^ {- \theta} P _ {j} ^ {- \theta}}\right) \frac {\hat {\gamma} _ {i}}{\hat {\Pi} _ {i} ^ {- \theta}} + \sum_ {j} \left(\frac {t _ {j i} ^ {- \theta} \Pi_ {i} ^ {- \theta} P _ {j} ^ {- \theta}}{\gamma_ {i} + \sum_ {j} t _ {j i} ^ {- \theta} \Pi_ {i} ^ {- \theta} P _ {j} ^ {- \theta}}\right) \hat {t} _ {j i} ^ {- \theta} \hat {P} _ {j} ^ {- \theta} \iff
$$

$$
\hat {P} _ {i} ^ {- \theta} = \left(\frac {\gamma_ {i}}{\gamma_ {i} + \sum_ {j} \Xi_ {j i}}\right) \frac {\hat {\gamma} _ {i}}{\hat {\Pi} _ {i} ^ {- \theta}} + \sum_ {j} \left(\frac {\Xi_ {j i}}{\gamma_ {i} + \sum_ {j} \Xi_ {j i}}\right) \hat {t} _ {j i} ^ {- \theta} \hat {P} _ {j} ^ {- \theta}
$$

Finally, substituting in our expression for iceberg trade costs along a link, 

$$
\hat {t} _ {i j} = \hat {\tilde {t}} _ {i j} ^ {\frac {1}{1 + \theta \lambda}} \times \hat {P} _ {i} ^ {- \frac {\theta \lambda}{1 + \theta \lambda}} \times \hat {\Pi} _ {j} ^ {- \frac {\theta \lambda}{1 + \theta \lambda}}
$$

and multiplying each equation by the other market access term, we obtain the following for the two equilibrium conditions: 

$$
\hat {\Pi} _ {i} ^ {- \theta} \hat {P} _ {i} ^ {- \theta} = \left(\frac {\delta_ {i}}{\delta_ {i} + \sum_ {j} \Xi_ {i j}}\right) \hat {\delta} _ {i} + \sum_ {j} \left(\frac {\Xi_ {i j}}{\delta_ {i} + \sum_ {j} \Xi_ {i j}}\right) \hat {\tilde {t}} _ {i j} ^ {- \frac {\theta}{1 + \theta \lambda}} \hat {P} _ {i} ^ {- \frac {\theta}{1 + \theta \lambda}} \hat {\Pi} _ {j} ^ {- \frac {\theta}{1 + \theta \lambda}}
$$

$$
\hat {\Pi} _ {i} ^ {- \theta} \hat {P} _ {i} ^ {- \theta} = \left(\frac {\gamma_ {i}}{\gamma_ {i} + \sum_ {j} \Xi_ {j i}}\right) \hat {\gamma} _ {i} + \sum_ {j} \left(\frac {\Xi_ {j i}}{\gamma_ {i} + \sum_ {j} \Xi_ {j i}}\right) \hat {t} _ {j i} ^ {- \frac {\theta}{1 + \theta \lambda}} \hat {P} _ {j} ^ {- \frac {\theta}{1 + \theta \lambda}} \hat {\Pi} _ {i} ^ {- \frac {\theta}{1 + \theta \lambda}}
$$

Now that we have defined the counterfactual equations generally, we turn to the specific cases of the trade and traffic models. 

# C.2.2 The Economic Geography Model

For the economic geography model, we have the following definitions for the fixed effects: 

$$
\delta_ {i} = E _ {i}
$$

$$
\gamma_ {i} = Y _ {i}
$$

We also derive the following for the price indices: 

$$
P _ {i} = \frac {w _ {i} u _ {i}}{W} \iff
$$

$$
P _ {i} = Y _ {i} \bar {u} _ {i} L _ {i} ^ {\beta - 1} W ^ {- 1} \Longrightarrow
$$

$$
\hat {P} _ {i} = \hat {y} _ {i} \hat {l} _ {i} ^ {\beta - 1} \hat {W} ^ {- 1}
$$

and 

$$
\Pi_ {i} = A _ {i} L _ {i} Y _ {i} ^ {- \frac {\theta + 1}{\theta}} \iff
$$

$$
\Pi_ {i} = \bar {A} _ {i} L _ {i} ^ {\alpha + 1} Y _ {i} ^ {- \frac {\theta + 1}{\theta}} \Rightarrow
$$

$$
\hat {\Pi} _ {i} = \hat {l} _ {i} ^ {\alpha + 1} \hat {y} _ {i} ^ {- \frac {\theta + 1}{\theta}}
$$

Substituting into the equilibrium conditions, we get: 

$$
\begin{array}{l} \left(\hat {l} _ {i} ^ {\alpha + 1} \hat {y} _ {i} ^ {- \frac {\theta + 1}{\theta}}\right) ^ {- \theta} \left(\hat {y} _ {i} \hat {l} _ {i} ^ {\beta - 1} \hat {W} ^ {- 1}\right) ^ {- \theta} = \left(\frac {E _ {i}}{E _ {i} + \sum_ {j} \Xi_ {i j}}\right) \hat {y} _ {i} + \\ \sum_ {j} \left(\frac {\Xi_ {i j}}{E _ {i} + \sum_ {j} \Xi_ {i j}}\right) \hat {\tilde {t}} _ {i j} ^ {- \frac {\theta}{1 + \theta \lambda}} \left(\hat {y} _ {i} \hat {l} _ {i} ^ {\beta - 1} \hat {W} ^ {- 1}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \left(\hat {l} _ {j} ^ {\alpha + 1} \hat {y} _ {j} ^ {- \frac {\theta + 1}{\theta}}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \iff \\ \hat {y} _ {i} \hat {l} _ {i} ^ {- \theta (\alpha + \beta)} \left(\hat {y} _ {i} \hat {l} _ {i} ^ {\beta - 1}\right) ^ {\frac {\theta}{1 + \theta \lambda}} = \left(\frac {E _ {i}}{E _ {i} + \sum_ {j} \Xi_ {i j}}\right) \hat {W} ^ {- \theta} \hat {y} _ {i} \left(\hat {y} _ {i} \hat {l} _ {i} ^ {\beta - 1}\right) ^ {\frac {\theta}{1 + \theta \lambda}} + \\ \hat {W} ^ {\frac {\theta}{1 + \theta \lambda} - \theta} \sum_ {j} \left(\frac {\Xi_ {i j}}{E _ {i} + \sum_ {j} \Xi_ {i j}}\right) \hat {y} _ {i} \hat {\tilde {t}} _ {i j} ^ {- \frac {\theta}{1 + \theta \lambda}} \left(\hat {l} _ {j} ^ {\alpha + 1} \hat {y} _ {j} ^ {- \frac {\theta + 1}{\theta}}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \iff \\ \hat {y} _ {i} ^ {\frac {1 + \theta \lambda + \theta}{1 + \theta \lambda}} \hat {l} _ {i} ^ {- \frac {\theta (1 + \alpha + \theta \lambda (\beta + \alpha))}{1 + \theta \lambda}} = \hat {\chi} \left(\frac {E _ {i}}{E _ {i} + \sum_ {k} \Xi_ {i k}}\right) \hat {y} _ {i} ^ {\frac {1 + \theta \lambda + \theta}{1 + \theta \lambda}} \hat {l} _ {i} ^ {\frac {\theta (\beta - 1)}{1 + \theta \lambda}} + \hat {\chi} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j} \left(\frac {\Xi_ {i j}}{E _ {i} + \sum_ {k} \Xi_ {i k}}\right) \hat {\tilde {t}} _ {i j} ^ {- \frac {\theta}{1 + \theta \lambda}} \hat {y} _ {j} ^ {\frac {1 + \theta}{1 + \theta \lambda}} \hat {l} _ {j} ^ {- \frac {\theta (1 + \alpha)}{1 + \theta \lambda}}, \\ \end{array}
$$

and: 

$$
\begin{array}{l} \left(\hat {l} _ {i} ^ {\alpha + 1} \hat {y} _ {i} ^ {- \frac {\theta + 1}{\theta}}\right) ^ {- \theta} \left(\hat {y} _ {i} \hat {l} _ {i} ^ {\beta - 1} \hat {W} ^ {- 1}\right) ^ {- \theta} = \left(\frac {Y _ {i}}{Y _ {i} + \sum_ {j} \Xi_ {j i}}\right) \hat {y} _ {i} + \\ \sum_ {j} \left(\frac {\Xi_ {j i}}{Y _ {i} + \sum_ {j} \Xi_ {j i}}\right) \hat {t} _ {j i} ^ {- \frac {\theta}{1 + \theta \lambda}} \left(\hat {y} _ {j} \hat {l} _ {j} ^ {\beta - 1} \hat {W} ^ {- 1}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \left(\hat {l} _ {i} ^ {\alpha + 1} \hat {y} _ {i} ^ {- \frac {\theta + 1}{\theta}}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \iff \\ \hat {y} _ {i} \hat {l} _ {i} ^ {- \theta (\alpha + \beta)} \left(\hat {l} _ {i} ^ {\alpha + 1} \hat {y} _ {i} ^ {- \frac {\theta + 1}{\theta}}\right) ^ {\frac {\theta}{1 + \theta \lambda}}. = \left(\frac {Y _ {i}}{Y _ {i} + \sum_ {j} \Xi_ {j i}}\right) \hat {W} ^ {\theta} \hat {y} _ {i} \left(\hat {l} _ {i} ^ {\alpha + 1} \hat {y} _ {i} ^ {- \frac {\theta + 1}{\theta}}\right) ^ {\frac {\theta}{1 + \theta \lambda}} + \\ \hat {W} ^ {\frac {\theta}{1 + \theta \lambda} - \theta} \sum_ {j} \left(\frac {\Xi_ {j i}}{Y _ {i} + \sum_ {j} \Xi_ {j i}}\right) \hat {\tilde {t}} _ {j i} ^ {- \frac {\theta}{1 + \theta \lambda}} \left(\hat {y} _ {j} \hat {l} _ {j} ^ {\beta - 1} \hat {W} ^ {- 1}\right) ^ {- \frac {\theta}{1 + \theta \lambda}} \iff \\ \hat {y} _ {i} ^ {\frac {- \theta (1 - \lambda)}{1 + \theta^ {\lambda}}} \hat {l} _ {i} ^ {\frac {\theta (1 - \beta - \theta \lambda (\alpha + \beta))}{1 + \theta^ {\lambda}}} = \hat {\chi} \left(\frac {Y _ {i}}{Y _ {i} + \sum_ {k} \Xi_ {k i}}\right) \hat {y} _ {i} ^ {\frac {- \theta (1 - \lambda)}{1 + \theta^ {\lambda}}} \hat {l} _ {i} ^ {\frac {\theta (\alpha + 1)}{1 + \theta^ {\lambda}}} + \hat {\chi} ^ {\frac {\theta \lambda}{1 + \theta^ {\lambda}}} \sum_ {j = 1} ^ {N} \left(\frac {\Xi_ {j i}}{Y _ {i} + \sum_ {k} \Xi_ {k i}}\right) \hat {y} _ {j} ^ {- \frac {\theta}{1 + \theta^ {\lambda}}} \hat {l} _ {j} ^ {\frac {\theta (1 - \beta)}{1 + \theta^ {\lambda}}}. \\ \end{array}
$$

To summarize the economic geography model, one can write the equilibrium system of equations (10) and (11) in changes as: 

$$
\begin{array}{l} \hat {y} _ {i} ^ {\frac {1 + \theta \lambda + \theta}{1 + \theta \lambda}} \hat {l} _ {i} ^ {- \frac {\theta (1 + \alpha + \theta \lambda (\beta + \alpha))}{1 + \theta \lambda}} = \hat {\chi} \left(\frac {E _ {i}}{E _ {i} + \sum_ {k} \Xi_ {i k}}\right) \hat {y} _ {i} ^ {\frac {1 + \theta \lambda + \theta}{1 + \theta \lambda}} \hat {l} _ {i} ^ {\frac {\theta (\beta - 1)}{1 + \theta \lambda}} + \hat {\chi} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j} \left(\frac {\Xi_ {i j}}{E _ {i} + \sum_ {k} \Xi_ {i k}}\right) \hat {\tilde {t}} _ {i j} ^ {- \frac {\theta}{1 + \theta \lambda}} \hat {y} _ {j} ^ {\frac {1 + \theta}{1 + \theta \lambda}} \hat {l} _ {j} ^ {- \frac {\theta (1 + \alpha)}{1 + \theta \lambda}} (C.45) \\ \hat {y} _ {i} ^ {\frac {- \theta (1 - \lambda)}{1 + \theta \lambda}} \hat {l} _ {i} ^ {\frac {\theta (1 - \beta - \theta \lambda (\alpha + \beta))}{1 + \theta \lambda}} = \hat {\chi} \left(\frac {Y _ {i}}{Y _ {i} + \sum_ {k} \Xi_ {k i}}\right) \hat {y} _ {i} ^ {\frac {- \theta (1 - \lambda)}{1 + \theta \lambda}} \hat {l} _ {i} ^ {\frac {\theta (\alpha + 1)}{1 + \theta \lambda}} + \hat {\chi} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j = 1} ^ {N} \left(\frac {\Xi_ {j i}}{Y _ {i} + \sum_ {k} \Xi_ {k i}}\right) \hat {\tilde {t}} _ {j i} ^ {- \frac {\theta}{1 + \theta \lambda}} \hat {y} _ {j} ^ {- \frac {\theta}{1 + \theta \lambda}} \hat {l} _ {j} ^ {\frac {\theta (1 - \beta)}{1 + \theta \lambda}} (C.46) \\ \end{array}
$$

It is immediately evident that the mathematical structure is identical to the one considered in Proposition 

1 (see equations (C.11) and (C.12)). As a result, the results of that Proposition apply. 

# C.2.3 The Urban Model

For the commuting model, we define the following for fixed effects: 

$$
\gamma_ {i} = L _ {i} ^ {R}
$$

$$
\delta_ {i} = L _ {i} ^ {F}
$$

We also derive the following for the price indices: 

$$
P _ {i} = A _ {i} \left(L _ {i} ^ {F}\right) ^ {- \frac {1}{\theta}} W ^ {- 1} \Longleftrightarrow
$$

$$
P _ {i} = \bar {A} _ {i} \left(L _ {i} ^ {F}\right) ^ {\alpha - \frac {1}{\theta}} W ^ {- 1} \Longrightarrow
$$

$$
\hat {P} _ {i} = \left(\hat {l} _ {i} ^ {F}\right) ^ {\alpha - \frac {1}{\theta}} \hat {W} ^ {- 1}
$$

and: 

$$
\Pi_ {i} = u _ {i} \left(L _ {i} ^ {R}\right) ^ {- \frac {1}{\theta}} \iff
$$

$$
\Pi_ {i} = \bar {u} _ {i} \left(L _ {i} ^ {R}\right) ^ {\beta - \frac {1}{\theta}} \Longrightarrow
$$

$$
\hat {\Pi} _ {i} = \left(\hat {l} _ {i} ^ {R}\right) ^ {\beta - \frac {1}{\theta}}
$$

Substituting these into our generalized equilibrium conditions, we get: 

$$
\begin{array}{l} \left(\hat {l} _ {i} ^ {R}\right) ^ {1 - \beta \theta} \left(\hat {l} _ {i} ^ {F}\right) ^ {1 - \alpha \theta} \hat {W} ^ {\theta} = \left(\frac {L _ {i} ^ {F}}{L _ {i} ^ {F} + \sum_ {j} \Xi_ {i j}}\right) \hat {l} _ {i} ^ {F} + \sum_ {j} \left(\frac {\Xi_ {i j}}{L _ {i} ^ {F} + \sum_ {j} \Xi_ {i j}}\right) \hat {\bar {t}} _ {i j} ^ {- \frac {\theta}{1 + \theta \lambda}} \left(\hat {l} _ {i} ^ {F}\right) ^ {\frac {1 - \alpha \theta}{1 + \theta \lambda}} \hat {W} ^ {\frac {\theta}{1 + \theta \lambda}} \left(\hat {l} _ {j} ^ {R}\right) ^ {\frac {1 - \beta \theta}{1 + \theta \lambda}} \iff \\ \left(\hat{l}_{i}^{R}\right)^{1 - \beta \theta}\left(\hat{l}_{i}^{F}\right)^{1 - \alpha \theta -\frac{1 - \alpha\theta}{1 + \theta \lambda}} = \left(\frac{L_{i}^{F}}{L_{i}^{F} + \sum_{j}\Xi_{ij}}\right)\hat{W}^{-\theta}\left(\hat{l}_{i}^{F}\right)^{1 - \frac{1 - \alpha\theta}{1 + \theta \lambda}} + \hat{W}^{\frac{\theta}{1 + \theta \lambda} -\theta}\sum_{j}\left(\frac{\Xi_{ij}}{L_{i}^{F} + \sum_{j}\Xi_{ij}}\right)\hat{\tilde{t}}_{ij}^{-\frac{\theta}{1 + \theta \lambda}}\left(\hat{l}_{j}^{R}\right)^{\frac{1 - \beta\theta}{1 + \theta \lambda}}\iff \\ \left(\hat {l} _ {i} ^ {R}\right) ^ {1 - \beta \theta} \left(\hat {l} _ {i} ^ {F}\right) ^ {\frac {\theta \lambda (1 - \alpha \theta)}{1 + \theta \lambda}} = \hat {\chi} \left(\frac {L _ {i} ^ {F}}{L _ {i} ^ {F} + \sum_ {j} \Xi_ {i j}}\right) \left(\hat {l} _ {i} ^ {F}\right) ^ {\frac {\theta (\alpha + \lambda)}{1 + \theta \lambda}} + \hat {\chi} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j} \left(\frac {\Xi_ {i j}}{L _ {i} ^ {F} + \sum_ {j} \Xi_ {i j}}\right) \hat {\bar {t}} _ {i j} ^ {- \frac {\theta}{1 + \theta \lambda}} \left(\hat {l} _ {j} ^ {R}\right) ^ {\frac {1 - \beta \theta}{1 + \theta \lambda}}, \\ \end{array}
$$

and 

$$
\begin{array}{l} \left(\hat {l} _ {i} ^ {R}\right) ^ {1 - \beta \theta} \left(\hat {l} _ {i} ^ {F}\right) ^ {1 - \alpha \theta} \hat {W} ^ {\theta} = \left(\frac {L _ {i} ^ {R}}{L _ {i} ^ {R} + \sum_ {j} \Xi_ {j i}}\right) \hat {l} _ {i} ^ {R} + \sum_ {j} \left(\frac {\Xi_ {j i}}{L _ {i} ^ {R} + \sum_ {j} \Xi_ {j i}}\right) \hat {\tilde {t}} _ {j i} ^ {- \frac {\theta}{1 + \theta \lambda}} \left(\hat {l} _ {j} ^ {F}\right) ^ {\frac {1 - \alpha \theta}{1 + \theta \lambda}} \hat {W} ^ {\frac {\theta}{1 + \theta \lambda}} \left(\hat {l} _ {i} ^ {R}\right) ^ {\frac {1 - \beta \theta}{1 + \theta \lambda}} \iff \\ \left(\hat {l} _ {i} ^ {R}\right) ^ {1 - \beta \theta - \frac {1 - \beta \theta}{1 + \theta \lambda}} \left(\hat {l} _ {i} ^ {F}\right) ^ {1 - \alpha \theta} = \left(\frac {L _ {i} ^ {R}}{L _ {i} ^ {R} + \sum_ {j} \Xi_ {j i}}\right) \hat {W} ^ {- \theta} \left(\hat {l} _ {i} ^ {R}\right) ^ {1 - \frac {1 - \beta \theta}{1 + \theta \lambda}} + \hat {W} ^ {\frac {\theta}{1 + \theta \lambda} - \theta} \sum_ {j} \left(\frac {\Xi_ {j i}}{L _ {i} ^ {R} + \sum_ {j} \Xi_ {j i}}\right) \hat {\tilde {t}} _ {j i} ^ {- \frac {\theta}{1 + \theta \lambda}} \left(\hat {l} _ {j} ^ {F}\right) ^ {\frac {1 - \alpha \theta}{1 + \theta \lambda}} \iff \\ \left(\hat {l} _ {i} ^ {R}\right) ^ {\frac {\theta \lambda (1 - \beta \theta)}{1 + \theta \lambda}} \left(\hat {l} _ {i} ^ {F}\right) ^ {1 - \alpha \theta} = \hat {\chi} ^ {\theta} \left(\frac {L _ {i} ^ {R}}{L _ {i} ^ {R} + \sum_ {j} \Xi_ {j i}}\right) \left(\hat {l} _ {i} ^ {R}\right) ^ {\frac {\theta (\beta + \lambda)}{1 + \theta \lambda}} + \hat {\chi} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j} \left(\frac {\Xi_ {j i}}{L _ {i} ^ {R} + \sum_ {j} \Xi_ {j i}}\right) \hat {\tilde {t}} _ {j i} ^ {- \frac {\theta}{1 + \theta \lambda}} \left(\hat {l} _ {j} ^ {F}\right) ^ {\frac {1 - \alpha \theta}{1 + \theta \lambda}}. \\ \end{array}
$$

To summarize, for the urban model, the equilibrium system defined by equations (19) and (20) can be written 

in changes as: 

$$
\left(\hat {l} _ {i} ^ {R}\right) ^ {1 - \theta \beta} \left(\hat {l} _ {i} ^ {F}\right) ^ {\frac {\theta \lambda (1 - \alpha \theta)}{1 + \theta \lambda}} = \hat {\chi} \left(\frac {L _ {i} ^ {F}}{L _ {i} ^ {F} + \sum_ {k} \Xi_ {i k}}\right) \left(\hat {l} _ {i} ^ {F}\right) ^ {\frac {\theta (\alpha + \lambda)}{1 + \theta \lambda}} + \hat {\chi} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j = 1} ^ {N} \left(\frac {\Xi_ {i j}}{L _ {i} ^ {F} + \sum_ {k} \Xi_ {i k}}\right) \hat {t} _ {i j} ^ {- \frac {\theta}{1 + \theta \lambda}} \left(\hat {l} _ {j} ^ {R}\right) ^ {\frac {1 - \theta \beta}{1 + \theta \lambda}} \tag {C.47}
$$

$$
\left(\hat {l} _ {i} ^ {R}\right) ^ {\frac {\theta \lambda (1 - \beta \theta)}{1 + \theta \lambda}} \left(\hat {l} _ {i} ^ {F}\right) ^ {1 - \theta \alpha} = \hat {\chi} \left(\frac {L _ {i} ^ {R}}{L _ {i} ^ {R} + \sum_ {k} \Xi_ {k i}}\right) \left(\hat {l} _ {i} ^ {R}\right) ^ {\frac {\theta (\beta + \lambda)}{1 + \theta \lambda}} + \hat {\chi} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j = 1} ^ {N} \left(\frac {\Xi_ {j i}}{L _ {i} ^ {R} + \sum_ {k} \Xi_ {k i}}\right) \hat {t} _ {j i} ^ {- \frac {\theta}{1 + \theta \lambda}} \left(\hat {l} _ {j} ^ {F}\right) ^ {\frac {1 - \theta \alpha}{1 + \theta \lambda}}. \tag {C.48}
$$

Again, it is immediately evident that the mathematical structure is identical to the one considered in Proposition 1 (see equations (C.17) and (C.18)). As a result, the results of that Proposition apply. 

# D Extensions

In this section, we present two model extensions mentioned in the main text. 

# D.1 Additive Transportation Costs

In this subsection, we describe an alternative framework where transportation costs are additive across segments (rather than multiplicative, as assumed in the main text). Let $\mathit { L } _ { i j }$ be the worker hours employed in the prothe good, and shipment of goods from  and hours the workers spen $i$ $j$ , which is splitping the good the hours workers spend producing: $L _ { i j } ^ { p r o d u c e }$ $L _ { i j } ^ { s h i p }$ 

$$
L _ {i j} = L _ {i j} ^ {p r o d u c e} + L _ {i j} ^ {s h i p}
$$

The total number of goods being sent from $i$ to $j$ is: 

$$
Q _ {i j} = A _ {i} L _ {i j} ^ {p r o d u c e}
$$

Let $t _ { i j }$ be the expected travel time from $i$ to $j$ (see below). Suppose that each unit of good requires a separate truck so that the total amount of labor hours used in shipping goods is: 

$$
L _ {i j} ^ {s h i p} = Q _ {i j} t _ {i j} = A _ {i} L _ {i j} ^ {p r o d u c e} t _ {i j}.
$$

Hence, to produce each unit of a good to send from $i$ to $j$ requires $\frac { ( 1 + t _ { i j } ) } { A _ { i } }$ units of labor and costs $p _ { i j } =$ (1 + tij ) wiAi . $\begin{array} { r } { ( 1 + t _ { i j } ) \frac { w _ { i } } { A _ { i } } } \end{array}$ W. If we define $\tau _ { i j } \equiv 1 + t _ { i j }$ i, the model is identical to the economic geography presented in the main text. 

Now we calculate the expected cost. Let $\mu _ { i j }$ the direct travel time between $i$ and $j$ . The total aggregate travel time from $i$ to $j$ on path $p$ of length $K$ is: 

$$
\tilde {t} _ {i j} (p) = \sum_ {k = 1} ^ {K} \mu_ {p _ {k - 1}, p _ {k}}
$$

Suppose each worker $\nu \in [ 1 , L _ { i } ]$ is heterogeneous in her preferences of routes so that she chooses the path $p$ to minimize: 

$$
\hat {t} _ {i j} \left(\nu\right) = \min  _ {K \geq 0, p \in P _ {i j} ^ {K}} \tilde {t} _ {i j} \left(p\right) + \varepsilon_ {i j} \left(p, \nu\right),
$$

where $\varepsilon _ { i j } \left( p , \nu \right)$ is Gumbel distributed with shape parameter $\theta$ . 

The expected trade cost can then be written as: 

$$
\begin{array}{l} t _ {i j} = - \frac {1}{\theta} \ln \sum_ {K \geq 0, p \in P _ {i j} ^ {K}} \exp \left(\tilde {t} _ {i j} (p)\right) ^ {- \theta} \Longleftrightarrow \\ t _ {i j} = - \frac {1}{\theta} \ln \sum_ {K \geq 0, p \in P _ {i j} ^ {K}} \exp \left(\sum_ {k = 1} ^ {K} \mu_ {p _ {k - 1}, p _ {k}}\right) ^ {- \theta} \Longleftrightarrow \\ t _ {i j} = - \frac {1}{\theta} \ln \sum_ {K = 0} ^ {\infty} \sum_ {p \in P _ {i j} ^ {K}} \left(\prod_ {k = 1} ^ {K} \exp \left(- \theta \mu_ {p _ {k - 1}, p _ {k}}\right)\right) \Longleftrightarrow \\ t _ {i j} = - \frac {1}{\theta} \ln \sum_ {K = 0} ^ {\infty} \left(\sum_ {k _ {1} = 1} ^ {N} \sum_ {k _ {2} = 1} ^ {N} \dots \sum_ {k _ {K - 1} = 1} ^ {N} \left(a _ {i, k _ {1}} \times a _ {k _ {1}, k _ {2}} \times \dots \times a _ {k _ {K - 2}, k _ {K - 1}} \times a _ {k _ {K - 1}, j}\right)\right) \iff , \\ t _ {i j} = - \frac {1}{\theta} \ln \sum_ {K = 0} ^ {\infty} A _ {i j} ^ {K} \\ \end{array}
$$

where $a _ { i j } \equiv \exp \left( - \theta \mu _ { i j } \right)$ , $\mathbf { A } \equiv [ a _ { i j } ]$ and $\mathbf { A } ^ { K } \equiv \left\lfloor A _ { i j } ^ { K } \right\rfloor$ . Define $\textstyle \mathbf { B } \equiv \left( \mathbf { I } - \mathbf { A } \right) ^ { - 1 } = \sum _ { K = 0 } ^ { \infty } { A ^ { K } }$ . We then have: 

$$
t _ {i j} = - \frac {1}{\theta} \ln b _ {i j}.
$$

Hence, the iceberg transportation cost with additive transportation costs can be written as: 

$$
\tau_ {i j} ^ {\prime} = 1 + \ln \left(b _ {i j} ^ {- \frac {1}{\theta}}\right) \approx b _ {i j} ^ {- \frac {1}{\theta}},
$$

i.e. it is approximately equal to the iceberg transportation cost defined in equation (21). 

# D.2 Nested Route Choice

In this subsection, we present an alternative economic geography model where agents first choose from which location to source the good and then choose on which route to ship the good. Suppose that each location $i \in \mathcal N$ is endowed with a constant returns to scale technology for producing and shipping each good $\nu \in [ 0 , 1 ]$ to each destination $j \in \mathcal N$ along each route $r \in \Re _ { i j }$ , which under perfect competition yields the following price of good $\nu \in [ 0 , 1 ]$ in destination $j \in \mathcal N$ from origin $i \in \mathcal N$ along route $r \in \Re _ { i j }$ : 

$$
p _ {i j, r} (\nu) = \frac {w _ {i}}{\varepsilon_ {i} (\nu)} \times \frac {\prod_ {k = 1} ^ {K} t _ {r _ {k - 1} , r _ {k}}}{\epsilon_ {i j , r} (\nu)},
$$

where $\varepsilon _ { i } \left( \nu \right)$ is an independently and identically Frechet distributed across goods distributed with scale parameter $1 / A _ { i }$ and shape parameter $\theta _ { g }$ and $\epsilon _ { i j , r } \left( \nu \right)$ is independently and identically Frechet distributed across routes with a scale parameter equal to $\Gamma \left( \frac { \theta _ { r } - 1 } { \theta _ { r } } \right) ^ { - 1 }$ (which is done without loss of generality and for convenience alone) and shape parameter $\theta _ { r }$ . The timing is as follows: first, individuals observe $\varepsilon _ { i } \left( \nu \right)$ and choose a source to purchase the good; second, individuals observe $\epsilon _ { i j , r } \left( \nu \right)$ and choose the route to ship the good. (For simplicity, individuals are not permitted to alter their decision of where to source the good once $\epsilon _ { i j , r } \left( \nu \right)$ is revealed). 

To solve the model, we proceed by backwards induction. Conditional on choosing to source good $\nu$ from location $i$ , a consumer in location $j$ will choose a route from $_ i$ to $j$ to minimize the shipping cost incurred, so that the probability she selects a route $r \in \Re _ { i j }$ is: 

$$
\pi_ {i j, r | i} = \frac {\left(\prod_ {k = 1} ^ {K} t _ {r _ {k - 1} , r _ {k}}\right) ^ {\theta_ {r}}}{\sum_ {r ^ {\prime} \in \Re_ {i j}} \left(\prod_ {k = 1} ^ {K} t _ {r _ {k - 1} ^ {\prime} , r _ {k} ^ {\prime}}\right) ^ {\theta_ {r}}}
$$

and the expected cost she incurs in shipping a good from $i$ to $j$ is: 

$$
\tau_ {i j} = \left(\sum_ {r ^ {\prime} \in \Re_ {i j}} \left(\prod_ {k = 1} ^ {K} t _ {r _ {k - 1} ^ {\prime}, r _ {k} ^ {\prime}}\right) ^ {- \theta_ {r}}\right) ^ {- \frac {1}{\theta_ {r}}} \tag {D.1}
$$

Apart from the different $\theta$ , equation (D.1) is equivalent to equation (4) and so a similar be expressed equivalently as in equation (21): 

$$
\tau_ {i j} = b _ {i j} ^ {- \frac {1}{\theta_ {r}}}, \tag {D.2}
$$

where $b _ { i j }$ is the $( i , j )$ element of the matrix $\mathbf { B } = \left( \mathbf { I } - \mathbf { A } \right) ^ { - 1 }$ and $\mathbf { A } \equiv [ a _ { i j } ] = \left[ t _ { i j } ^ { - \theta _ { r } } \right]$ 

Knowing that this is the expected cost she will incur, a consumer in location $j$ will first choose the location $i$ to source the good $\nu$ from in order to minimize the expected total cost, so that the probability she 

sources from $i \in \mathcal N$ is: 

$$
\pi_ {i j} = \frac {\tau_ {i j} ^ {- \theta_ {g}} \left(w _ {i} / A _ {i}\right) ^ {- \theta_ {g}}}{\sum_ {k \in \mathcal {N}} \tau_ {k j} ^ {- \theta_ {g}} \left(w _ {k} / A _ {k}\right) ^ {- \theta_ {g}}} \tag {D.3}
$$

and the total value of trade from $i$ to $j$ can be written as: 

$$
X _ {i j} = \frac {\tau_ {i j} ^ {- \theta_ {g}} \left(w _ {i} / A _ {i}\right) ^ {- \theta_ {g}}}{\sum_ {k \in \mathcal {N}} \tau_ {k j} ^ {- \theta_ {g}} \left(w _ {k} / A _ {k}\right) ^ {- \theta_ {g}}} E _ {j}. \tag {D.4}
$$

It is immediately evident that when $\theta _ { r } = \theta _ { g } = \theta$ , equation (D.2) is identical to equation (21) and equation (D.4) is identical to equation (3), i.e. the model here becomes isomorphic to the one presented in the main text. More generally, (and to get a sense of where the tractability is lost when ${ \theta _ { r } } \neq \theta _ { g }$ ), combining equations (D.2) and (D.4) yields: 

$$
X _ {i j} = \frac {b _ {i j} ^ {\left(\frac {\theta_ {g}}{\theta_ {r}}\right)} (w _ {i} / A _ {i}) ^ {- \theta_ {g}}}{\sum_ {k \in \mathcal {N}} b _ {k j} ^ {\left(\frac {\theta_ {g}}{\theta_ {r}}\right)} (w _ {k} / A _ {k}) ^ {- \theta_ {g}}} E _ {j},
$$

so the bilateral trade flows are functions of elements of the Leontief inverse raised to a power (rather than the elements of the Leontief inverse themselves). 

# D.3 An Armington trade model

In the economic geography framework in the main paper, consumers have idiosyncratic preferences over the source and route of a particular good. Here we consider an alternative framework where (a) traders choose the least cost routes; (b) each location produces a unique differentiated variety; and (c) trade is the result of consumers having a constant elasticity of substitution preferences across these differentiated varieties. We show that our routing framework can deliver identical predictions as in the main paper in this Armington framework. 

Suppose that for each pair of locations $( i , j ) \in \mathcal { N } \times \mathcal { N }$ , there are measure of perfectly competitive traders $\nu \in [ 0 , 1 ]$ who choose the route $\boldsymbol { r } \in \Re _ { i j }$ from $i$ to $j$ to minimize their trade costs: 

$$
\tau_ {i j} \left(\nu\right) = \prod_ {k = 1} ^ {K} t _ {r _ {k - 1}, r _ {k}} \varepsilon_ {i j, r} \left(\nu\right),
$$

where $\varepsilon _ { i j , r } \left( \nu \right)$ is an idiosyncratic preference that is i.i.d. across traders and Frechet distributed with shape parameter $\theta$ . Suppose too that each unit value of consumption by consumers is transported by a randomly selected trader; an equivalent derivation to that in the text implies that a consumer’s expected trade cost from $_ i$ to $j$ is as in equation (21), i.e.: 

$$
\tau_ {i j} = \left(\left[ \mathbf {I} - \mathbf {A} \right] _ {i j} ^ {- 1}\right) ^ {\frac {1}{\theta}},
$$

where $\mathbf { A } \equiv \left[ t _ { k l } ^ { - \theta } \right]$ is the $N \times N$ adjacency matrix. Similarly, the expected intensity with which a trade from $i$ to $j$ uses a link $k l$ can be written as in equation (23), i.e.: 

$$
\pi_ {i j} ^ {k l} = \left(\frac {\tau_ {i j}}{\tau_ {i k} t _ {k l} \tau_ {l j}}\right) ^ {\theta}. \tag {D.5}
$$

That is, both the expected trade costs and the link intensity are identical to those derived in the main paper. 

The expression for trade flows, however, is distinct, as it depends on the elasticity of substitution $\sigma$ rather than the Frechet elasticity $\theta$ . Given an elasticity of substitution $\sigma$ , the value of equilibrium trade flows from $i$ to $j$ , $X _ { i j }$ , can be written as: 

$$
X _ {i j} = \tau_ {i j} ^ {1 - \sigma} \times w _ {i} ^ {1 - \sigma} A _ {i} ^ {\sigma - 1} \times P _ {j} ^ {\sigma - 1} E _ {j}, \tag {D.6}
$$

or, equivalently, as: 

$$
X _ {i j} = \tau_ {i j} ^ {1 - \sigma} \times \frac {Y _ {i}}{\Pi_ {i} ^ {1 - \sigma}} \times \frac {E _ {j}}{P _ {j} ^ {1 - \sigma}}, \tag {D.7}
$$

where $\Pi _ { i } ^ { 1 - \sigma }$ and $P _ { i } ^ { 1 - \sigma }$ are the outward market access and inward market access, respectively, and are implicitly defined as the solution to the following system of equations: 

$$
\Pi_ {i} ^ {1 - \sigma} = \sum_ {j \in \mathcal {N}} \tau_ {i j} ^ {1 - \sigma} \times \frac {E _ {j}}{P _ {j} ^ {1 - \sigma}}
$$

$$
P _ {j} ^ {1 - \sigma} = \sum_ {i \in \mathcal {N}} \tau_ {i j} ^ {1 - \sigma} \times \frac {Y _ {i}}{\Pi_ {i} ^ {1 - \sigma}}.
$$

It is immediately clear by inspection of (D.6) and (D.7) that $\Pi _ { i } ^ { 1 - \sigma } = Y _ { i } / w _ { i } ^ { 1 - \sigma } A _ { i } ^ { \sigma - 1 }$ and $P _ { j } ^ { 1 - \sigma } = E _ { j } / P _ { j } ^ { \sigma - 1 } E _ { j }$ . Notice also that given the functional form of $\tau _ { i j } = \left( \left[ \mathbf { I } - \mathbf { A } \right] _ { i j } ^ { - 1 } \right) ^ { \frac { 1 } { \theta } }$ assuming that $1 - \sigma = \theta$ makes the two equilibrium conditions linear functions of the infrastructure cost matrix, $t _ { k l }$ . This is a crucial step in deriving our equilibrium conditions via an inversion. 

Furthermore, let us now characterize the value of trade traversing over each link, i.e. the traffic flows. As in the main paper, equilibrium traffic on link $( k , l )$ can be calculated by summing the product of the link intensity and the total trade flows across all origins and destinations: 

$$
\Xi_ {k l} = \sum_ {i \in \mathcal {N}} \sum_ {j \in \mathcal {N}} \pi_ {i j} ^ {k l} X _ {i j} \tag {D.8}
$$

Substituting equations (D.5) and (D.7) into (D.8) yields: 

$$
\Xi_ {k l} = t _ {k l} ^ {- \theta} \times P _ {k} ^ {1 - \sigma} \times \Pi_ {l} ^ {1 - \sigma} \times \sum_ {i \in \mathcal {N}} \frac {\left(\tau_ {i k} ^ {(\sigma - 1) - \theta}\right) X _ {i k}}{Y _ {k}} \times \sum_ {j \in \mathcal {N}} \frac {\left(\tau_ {l j} ^ {(\sigma - 1) - \theta}\right) X _ {l j}}{Y _ {l}} \times \left(\tau_ {i j} ^ {\theta - (\sigma - 1)}\right)
$$

In the special case where $\sigma - 1 = \theta - \mathrm { i . e }$ . the route and trade elasticities are equal to each other – this simplifies to: 

$$
\Xi_ {k l} = t _ {k l} ^ {- \theta} \times P _ {k} ^ {1 - \sigma} \times \Pi_ {l} ^ {1 - \sigma}, \tag {D.9}
$$

i.e. the equilibrium traffic flows in the Armington framework are identical to those in the main paper. Intuitively, the restriction that $\sigma - 1 = \theta$ here is equivalent to the assumption in the main paper that individuals preferences are distributed over routes-source pairs. 

# D.4 Traffic congestion in quantities of labor

In the economic geography framework presented in the main paper, it is assumed that the cost of traversing a link is increasing in the traffic congestion along that link, where traffic is measured in terms of the value of goods flowing over the link. In this section, we show how to derive the traffic (and traffic congestion) in terms of the quantity of labor flowing over a link. In Online Appendix D.5, we show how the framework can be altered so that traffic (and traffic congestion) are functions of the quantity of goods flowing over a link in an Armington variant of the model. 

The setup is identical to the economic geography framework in the paper. We note that the quantity of labor required in $i$ to produce the value of goods exported to $j$ , $\boldsymbol { L } _ { i j }$ , can be written as: 

$$
L _ {i j} = \frac {X _ {i j}}{w _ {i}} = \tau_ {i j} ^ {- \theta} \times w _ {i} ^ {- \theta - 1} A _ {i} ^ {\theta} \times P _ {j} ^ {\theta} w _ {j} L _ {j}, \tag {D.10}
$$

or, equivalently, as: 

$$
L _ {i j} = \tau_ {i j} ^ {- \theta} \times \frac {\sum_ {j \in \mathcal {N}} L _ {i j}}{M A _ {i} ^ {L , o u t}} \times \frac {\sum_ {i \in \mathcal {N}} L _ {i j}}{M A _ {j} ^ {L , i n}}, \tag {D.11}
$$

where $M A _ { i } ^ { L , o u t }$ and $M A _ { j } ^ { L , i n }$ are outward and inward labor-based market access measures, respectively, and are implicitly defined as the solution to the following system of equations: 

$$
M A _ {i} ^ {L, o u t} = \sum_ {j} \tau_ {i j} ^ {- \theta} \times \frac {\sum_ {i \in \mathcal {N}} L _ {i j}}{M A _ {j} ^ {L , i n}} \tag {D.12}
$$

$$
M A _ {j} ^ {L, i n} = \sum_ {i \in \mathcal {N}} \tau_ {i j} ^ {- \theta} \times \frac {\sum_ {j \in \mathcal {N}} L _ {i j}}{M A _ {i} ^ {L , o u t}}. \tag {D.13}
$$

Comparison of equations (D.10) and (D.11) immediately imply $\begin{array} { r } { M A _ { i } ^ { L , o u t } = \sum _ { j \in \mathcal { N } } { L _ { i j } } / { w _ { i } ^ { - \theta - 1 } } A _ { i } ^ { \theta } } \end{array}$ and $M A _ { j } ^ { L , i n } =$ $\sum _ { i \in \mathcal { N } } { L _ { i j } } / P _ { j } ^ { \theta } w _ { j } L _ { j }$ . 

Unlike in the paper where traffic is measured in the value of goods flowing over a link, we now measure traffic in terms of the quantity of labor used to produce the value of goods exported from $_ i$ to $j$ . That is, we suppose that trade is accomplished – and traffic through the network is generated – by traveling salesmen (or management consultants) residing in $i$ who travel to $j$ to provide their time to consumers in $j$ . Then through a similar derivation as that in Section 3.2 of the paper, we can calculate the traffic flows in quantities of labor as: 

$$
\Xi_ {k l} ^ {L} = \sum_ {i \in \mathcal {N}} \sum_ {j \in \mathcal {N}} \pi_ {i j} ^ {k l} L _ {i j}. \tag {D.14}
$$

Substituting the link intensity equation (23) from the main paper and the labor gravity equation D.11 into equation D.14 yields: 

$$
\Xi_ {k l} ^ {L} = t _ {k l} ^ {- \theta} \times M A _ {k} ^ {L, i n} \times M A _ {l} ^ {L, o u t} \tag {D.15}
$$

i.e. the equilibrium traffic flows in quantities of labor follow a gravity equation just as in equation (24) in the main paper. Moreover, the intuition for equation (D.15) is the same as in the main paper: traffic is declining in the cost of traversing the link, the greater the inward market access, the more traffic that flows into a link $k$ , and the greater the outward market access $\left( \Pi _ { l } ^ { - \theta } \right)$ , the more traffic that flows out of link $\it { \Delta } l$ . Indeed, the only difference between the traffic gravity equation (D.15) here and the traffic gravity equation (24) in the paper is that the measures of inward and outward market access are in terms of quantities of labor rather than values of trade. 

Given the similarity between the traffic gravity equation here and in the paper, it is also straightforward to incorporate traffic congestion (in quantities). To do so, suppose that the cost of traversing a link is increasing in the quantity of goods flowing over the link, i.e.: 

$$
t _ {k l} = \bar {t} _ {k l} \times \left(\Xi_ {k l} ^ {L}\right) ^ {\lambda}. \tag {D.16}
$$

Combining equations (D.16) and (D.15) yield: 

$$
\Xi_ {k l} ^ {L} = \bar {t} _ {k l} ^ {- \frac {\theta}{1 + \lambda \theta}} \times \left(M A _ {k} ^ {L, i n}\right) ^ {\frac {1}{1 + \lambda \theta}} \times \left(M A _ {l} ^ {L, o u t}\right) ^ {\frac {1}{1 + \lambda \theta}} \tag {D.17}
$$

$$
t _ {k l} = \bar {t} _ {k l} ^ {- \frac {1}{1 + \lambda \theta}} \times \left(M A _ {k} ^ {L, i n}\right) ^ {\frac {\lambda}{1 + \lambda \theta}} \times \left(M A _ {l} ^ {L, o u t}\right) ^ {\frac {\lambda}{1 + \lambda \theta}}, \tag {D.18}
$$

i.e. just as in the framework presented in the paper, we can express the equilibrium traffic flows (and trade costs) as analytical functions of the underlying infrastructure matrix and the (endogenous) market access. Hence, as the discussion here hopefully makes clear, it is straightforward to adjust the framework to incorporate traffic congestion in quantities of labor, even without any alteration to the underlying framework (unlike in Online Appendix D.5, which is based on an Armington variant of the model). 

As in Online Appendix D.5, however, assuming that traffic congestion arises from traffic in quantities does increase the complexity of calculating the model equilibrium. We now discuss this additional complexity in detail. 

The additional complexity arises from the fact that it is no longer possible to express the labor market access terms solely as functions of the endogenous variables $\left\{ l _ { i } , y _ { i } , \chi \right\}$ and local geography. As a result, we require an additional equilibrium condition – that for the inward flow of the quantity of labor used elsewhere 

in the goods consumed locally – to characterize the equilibrium. Define $L _ { i } ^ { i n }$ to be the total labor used throughout the economy for the goods consumed in location $i$ : 

$$
L _ {i} ^ {i n} \equiv \sum_ {j \in \mathcal {N}} L _ {j i} \iff
$$

$$
l _ {i} ^ {i n} l _ {i} ^ {\theta (1 - \beta)} y _ {i} ^ {- (1 + \theta)} = \chi \sum_ {j \in \mathcal {N}} \tau_ {j i} ^ {- \theta} y _ {j} ^ {- \theta - 1} l _ {j} ^ {1 + (1 + \alpha) \theta} \bar {A} _ {j} ^ {\theta} \bar {u} _ {i} ^ {\theta}, \tag {D.19}
$$

where $l _ { i } ^ { i n } \equiv L _ { i } ^ { i n } / \bar { L }$ is the share of inward flows of labor and recall $\begin{array} { r } { \chi \equiv \left( \frac { \hat { L } ^ { ( \alpha + \beta ) } } { W } \right) ^ { \theta } } \end{array}$  L¯(α+β) θ Equation (D.19) can be inverted after substituting for $\tau _ { j i }$ (using equation (21)) to yield an equilibrium relationship in terms of the infrastructure network rather than the transportation costs: 

$$
l _ {i} ^ {i n} l _ {i} ^ {\theta (1 - \beta)} y _ {i} ^ {- (1 + \theta)} \bar {u} _ {i} ^ {- \theta} = \chi y _ {i} ^ {- \theta - 1} l _ {i} ^ {1 + (1 + \alpha) \theta} \bar {u} _ {i} ^ {\theta} \bar {A} _ {i} ^ {\theta} + \sum_ {j \in \mathcal {N}} t _ {j i} ^ {- \theta} l _ {j} ^ {i n} l _ {j} ^ {\theta (1 - \beta)} y _ {j} ^ {- (1 + \theta)} \bar {u} _ {j} ^ {- \theta}. \tag {D.20}
$$

The existing equilibrium equations (10) and (11) after substituting for $\tau _ { j i }$ and inverting (again, using equation (21)) remain unchanged 

$$
\bar {A} _ {i} ^ {- \theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {- \theta (1 + \alpha)} = \chi \bar {u} _ {i} ^ {\theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {\theta (\beta - 1)} + \sum_ {j = 1} ^ {N} t _ {i j} ^ {- \theta} \bar {A} _ {j} ^ {- \theta} y _ {j} ^ {1 + \theta} l _ {j} ^ {- \theta (1 + \alpha)} \tag {D.21}
$$

$$
\bar {u} _ {i} ^ {- \theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (1 - \beta)} = \chi \bar {A} _ {i} ^ {\theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (\alpha + 1)} + \sum_ {j = 1} ^ {N} t _ {j i} ^ {- \theta} \bar {u} _ {j} ^ {- \theta} y _ {j} ^ {- \theta} l _ {j} ^ {\theta (1 - \beta)}. \tag {D.22}
$$

comprise our equilibrium system, albeit in terms of the endogenous infrastructure network $[ t _ { k l } ]$ . To write the equilibrium system in terms of the fundamental infrastructure network $[ t _ { k l } ]$ , we take the definitions of the labor market access terms above combined with equation (D.18) to write the endogenous infrastructure network $[ t _ { k l } ]$ solely as a function of the endogenous variables $\left\{ l _ { i } ^ { i n } , l _ { i } , y _ { i } , \chi \right\}$ and model fundamentals. With some algebra, this becomes: 

$$
t _ {i j} = \chi^ {- \frac {\lambda}{1 + \lambda \theta}} \times \left(\bar {t} _ {i j} \bar {L} ^ {\lambda}\right) ^ {- \frac {1}{1 + \lambda \theta}} \times \left(l _ {i} ^ {i n} y _ {i} ^ {- (1 + \theta)} l _ {i} ^ {\theta (1 - \beta)} \bar {u} _ {i} ^ {- \theta}\right) ^ {\frac {\lambda}{1 + \lambda \theta}} \times \left(l _ {j} ^ {- \theta (1 + \alpha)} y _ {j} ^ {1 + \theta} \bar {A} _ {j} ^ {- \theta}\right) ^ {\frac {\lambda}{1 + \lambda \theta}}. \tag {D.23}
$$

Substituting equation (D.23) into equations (D.21), (D.22), and the new equilibrium equation (D.20), we arrive at our full equilibrium system written solely as a function of fundamental geographic variables: 

$$
\begin{array}{l} \bar {A} _ {i} ^ {- \theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {- \theta (1 + \alpha)} = \chi \bar {u} _ {i} ^ {\theta} y _ {i} ^ {1 + \theta} l _ {i} ^ {\theta (\beta - 1)} + \\ \chi^ {\frac {\lambda \theta}{1 + \lambda \theta}} \sum_ {j = 1} ^ {N} \left(\left(\bar {t} _ {i j} \bar {L} ^ {\lambda}\right) \left(l _ {i} ^ {i n} y _ {i} ^ {- (1 + \theta)} l _ {i} ^ {\theta (1 - \beta)} \bar {u} _ {i} ^ {- \theta} l _ {j} ^ {- \theta (1 + \alpha)} y _ {j} ^ {1 + \theta} \bar {A} _ {j} ^ {- \theta}\right) ^ {- \lambda}\right) ^ {\frac {\theta}{1 + \lambda \theta}} \bar {A} _ {j} ^ {- \theta} y _ {j} ^ {1 + \theta} l _ {j} ^ {- \theta (1 + \alpha)} \\ \end{array}
$$

$$
\begin{array}{l} \bar {u} _ {i} ^ {- \theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (1 - \beta)} = \chi \bar {A} _ {i} ^ {\theta} y _ {i} ^ {- \theta} l _ {i} ^ {\theta (\alpha + 1)} + \\ \chi^ {\frac {\lambda \theta}{1 + \lambda \theta}} \sum_ {j = 1} ^ {N} \left(\left(\bar {t} _ {j i} \bar {L} ^ {\lambda}\right) \left(l _ {j} ^ {i n} y _ {j} ^ {- (1 + \theta)} l _ {j} ^ {\theta (1 - \beta)} \bar {u} _ {j} ^ {- \theta} l _ {i} ^ {- \theta (1 + \alpha)} y _ {i} ^ {1 + \theta} \bar {A} _ {i} ^ {- \theta}\right) ^ {- \lambda}\right) ^ {\frac {\theta}{1 + \lambda \theta}} \bar {u} _ {j} ^ {- \theta} y _ {j} ^ {- \theta} l _ {j} ^ {\theta (1 - \beta)} \\ \end{array}
$$

$$
\begin{array}{l} l _ {i} ^ {i n} l _ {i} ^ {\theta (1 - \beta)} y _ {i} ^ {- (1 + \theta)} \bar {u} _ {i} ^ {- \theta} = \chi y _ {i} ^ {- \theta - 1} l _ {i} ^ {1 + (1 + \alpha) \theta} \bar {u} _ {i} ^ {\theta} \bar {A} _ {i} ^ {\theta} + \\ \chi^ {\frac {\lambda \theta}{1 + \lambda \theta}} \sum_ {j \in \mathcal {N}} \left(\left(\bar {t} _ {j i} \bar {L} ^ {\lambda}\right) \left(l _ {j} ^ {i n} y _ {j} ^ {- (1 + \theta)} l _ {j} ^ {\theta (1 - \beta)} \bar {u} _ {j} ^ {- \theta} l _ {i} ^ {- \theta (1 + \alpha)} y _ {i} ^ {1 + \theta} \bar {A} _ {i} ^ {- \theta}\right) ^ {- \lambda}\right) ^ {\frac {\theta}{1 + \lambda \theta}} l _ {j} ^ {i n} l _ {j} ^ {\theta (1 - \beta)} y _ {j} ^ {- (1 + \theta)} \bar {u} _ {j} ^ {- \theta}. \\ \end{array}
$$

These $3 N$ equations (along with imposing that the shares sum to one, i.e. $\begin{array} { r } { \sum _ { i \in \mathcal { N } } l _ { i } = \sum _ { i \in \mathcal { N } } l _ { i } ^ { N } = \sum _ { i \in \mathcal { N } } y _ { i } = } \end{array}$ 1) can be solved simultaneously for the $3 N + 1$ equilibrium variables $\left\{ l _ { i } ^ { i n } , l _ { i } , y _ { i } , \chi \right\}$ . 

# D.5 Traffic congestion in quantities of goods

In the economic geography framework presented in the main paper, it is assumed that the cost of traversing a link is increasing in the traffic congestion along that link, where traffic is measured in terms of the value of goods flowing over the link. In this section, we show how the framework can be altered so that traffic (and traffic congestion) are functions of the quantity of goods flowing over a link. 

We begin by noting that in the framework presented in the main paper, the value of bilateral trade flows between origins and destinations are aggregated across a continuum of different varieties. As the probabilistic nature of that setup precludes an analytical expression for the quantity traded of each distinct variety, here we instead focus on an Armington trade framework (discussed in Online Appendix D.3) where each location produces a single differentiated variety. As the focus here is on quantities traded, we now suppose too that each unit of quantity consumed (rather each unit of value) is transported by a randomly selected trader. Then the quantity shipped from $i$ to $j$ can be written simply as: 

$$
Q _ {i j} = \frac {X _ {i j}}{p _ {i j}} = \tau_ {i j} ^ {- \sigma} \times w _ {i} ^ {- \sigma} A _ {i} ^ {\sigma} \times P _ {j} ^ {\sigma - 1} E _ {j}. \tag {D.24}
$$

As with values, we can equivalently express the bilateral quantity traded as follows: 

$$
Q _ {i j} = \tau_ {i j} ^ {- \sigma} \times \frac {\sum_ {j \in \mathcal {N}} Q _ {i j}}{M A _ {i} ^ {Q , o u t}} \times \frac {\sum_ {i \in \mathcal {N}} Q _ {i j}}{M A _ {j} ^ {Q , i n}}, \tag {D.25}
$$

where $M A _ { i } ^ { Q , o u t }$ and $M A _ { j } ^ { Q , i n }$ are outward and inward quantity-based market access measures, respectively, and are implicitly defined as the solution to the following system of equations: 

$$
M A _ {i} ^ {Q, o u t} = \sum_ {j} \tau_ {i j} ^ {- \sigma} \times \frac {\sum_ {i \in \mathcal {N}} Q _ {i j}}{M A _ {j} ^ {Q , i n}} \tag {D.26}
$$

$$
M A _ {j} ^ {Q, i n} = \sum_ {i \in \mathcal {N}} \tau_ {i j} ^ {- \sigma} \times \frac {\sum_ {j \in \mathcal {N}} Q _ {i j}}{M A _ {i} ^ {Q , o u t}}. \tag {D.27}
$$

It is immediately clear by inspection of (D.24) and (D.25) that $\begin{array} { r } { M A _ { i } ^ { Q , o u t } = \sum _ { j \in \mathcal { N } } Q _ { i j } / w _ { i } ^ { - \sigma } A _ { i } ^ { \sigma } } \end{array}$ and $M A _ { j } ^ { Q , i n } =$ $\sum _ { i \in \mathcal { N } } Q _ { i j } / P _ { j } ^ { \sigma - 1 } E _ { j }$ . We remark that these quantity-based market access measures are distinct from the (typical) value-based market access measures discussed in the paper (and in Online Appendix D.3). 

Now we consider the quantity of goods that traverse a particular link, i.e. the traffic. As traders are matched randomly with each unit of good transported, the equilibrium traffic on link $( k , l )$ can be calculated by summing the product of the link intensity and the quantity of trade flows across all origins and destinations: 

$$
\Xi_ {k l} ^ {Q} = \sum_ {i \in \mathcal {N}} \sum_ {j \in \mathcal {N}} \pi_ {i j} ^ {k l} Q _ {i j}. \tag {D.28}
$$

Substituting the link intensity equation (D.5) (unchanged from Online Appendix D.3, as the traders route 

choice problem remains the same) and the quantity gravity equation (D.25) into equation (D.28) yields: 

$$
\Xi_ {k l} ^ {Q} = t _ {k l} ^ {- \theta} \times M A _ {k} ^ {Q, i n} \times M A _ {l} ^ {Q, o u t} \times \sum_ {i \in \mathcal {N}} \left(\frac {\left(\tau_ {i k} ^ {\sigma - \theta}\right) Q _ {i k}}{\sum_ {i \in \mathcal {N}} Q _ {i k}}\right) \sum_ {j \in \mathcal {N}} \left(\frac {\left(\tau_ {l j} ^ {\sigma - \theta}\right) Q _ {l j}}{\sum_ {j \in \mathcal {N}} Q _ {l j}}\right) \left(\tau_ {i j} ^ {\theta - \sigma}\right)
$$

In the special case where $\sigma = \theta - \mathrm { i . e }$ . the route elasticity and the elasticity of substitution are equal to each other – this simplifies to: 

$$
\Xi_ {k l} ^ {Q} = t _ {k l} ^ {- \theta} \times M A _ {k} ^ {Q, i n} \times M A _ {l} ^ {Q, o u t}, \tag {D.29}
$$

i.e. the equilibrium traffic flows in quantities follow a gravity equation just as in equation (24) in the main paper. Moreover, the intuition for equation (D.29) is the same as in the main paper: traffic is declining in the cost of traversing the link, the greater the inward market access, the more traffic that flows into a link $k$ , and the greater the outward market access, the more traffic that flows out of link $\it l$ . Indeed, the only difference between the traffic gravity equation (D.29) here and the traffic gravity equation (24) in the paper is that the measures of inward and outward market access are quantity rather than value based. 

Given the similarity between the traffic gravity equation (D.29) here and the traffic gravity equation (24) in the paper, it is also straightforward to incorporate traffic congestion (in quantities). To do so, suppose that the cost of traversing a link is increasing in the quantity of goods flowing over the link, i.e.: 

$$
t _ {k l} = \bar {t} _ {k l} \times \left(\Xi_ {k l} ^ {Q}\right) ^ {\lambda}. \tag {D.30}
$$

Combining equations (D.30) and (D.29) yield: 

$$
\Xi_ {k l} ^ {Q} = \bar {t} _ {k l} ^ {- \frac {\theta}{1 + \lambda \theta}} \times \left(M A _ {k} ^ {Q, i n}\right) ^ {\frac {1}{1 + \lambda \theta}} \times \left(M A _ {l} ^ {Q, o u t}\right) ^ {\frac {1}{1 + \lambda \theta}} \tag {D.31}
$$

$$
t _ {k l} = \bar {t} _ {k l} ^ {\frac {1}{1 + \lambda \theta}} \times \left(M A _ {k} ^ {Q, i n}\right) ^ {\frac {\lambda}{1 + \lambda \theta}} \times \left(M A _ {l} ^ {Q, o u t}\right) ^ {\frac {\lambda}{1 + \lambda \theta}}, \tag {D.32}
$$

i.e. just as in the framework presented in the paper, we can express the equilibrium traffic flows (and trade costs) as analytical functions of the underlying infrastructure matrix and the (endogenous) market access. Hence, as the discussion here hopefully makes clear, it is straightforward to adjust the framework to incorporate traffic congestion in quantities. However, we should note that doing so does increase the complexity of the model, as the quantity based market access measures constitute additional endogenous variables in the equilibrium system of equations to solve. That is, in addition to solving the $2 N$ equilibrium conditions (10) and (11), one has to also solve the $2 N$ quantity market access equations (D.26) and (D.27). This additional complexity is not present in the paper where traffic congestion arises from the value of traffic, as the traditional market access measures (in values) can be written solely in terms of local variables (see e.g. (7)). See Online Appendix D.4 for a detailed discussion of this additional complexity in a closely related problem where traffic congestion is assumed to be depend on the quantity of labor used to produce the bilateral goods traded. 

# E Algorithm for Conducting Counterfactuals

The algorithm we use to find the equilibrium in our counterfactual simulations consists of an outer loop, where we guess a $\hat { \chi }$ , and an inner loop, where, given a $\hat { \chi }$ , we solve for vectors of $\{ \hat { y } _ { i } \}$ and $\left\{ \hat { l } _ { i } \right\}$ in the economic geography case and vectors of $\left\{ \hat { l } _ { i } ^ { R } \right\}$ and $\left\{ \hat { l } _ { i } ^ { F } \right\}$ in the commuting case for which the system is equal up to scale. We can see that for equilibrium equations (36) and (37) (the economic geography model) and equilibrium equations (38) and (39) (the commuting model), given a $\hat { \chi }$ , the system forms a system of non-linear equations in $\left( \left\{ \hat { y } _ { i } \right\} , \left\{ \hat { l } _ { i } \right\} \right)$ and $\left( \left\{ \hat { l } _ { i } ^ { R } \right\} , \left\{ \hat { l } _ { i } ^ { F } \right\} \right)$ , respectively. At the same time, the term on the left-hand side of each of the equilibrium conditions is a log-linear combination of the endogenous variables. Because of these similarities, we use a very similar algorithm to solve both models. For concision, we will 

describe the algorithm in detail in terms of the endogenous variables of the economic geography model. 

Let us start with a detailed description of the inner loop. Given an initial guess of $\hat { \chi } _ { ( 0 ) } = 1$ , we plug in an initial guess of the endogenous variables $\{ \hat { y } _ { i } \} _ { ( 0 ) } = 1$ and $\left\{ \hat { l } _ { i } \right\} _ { ( 0 ) } = 1$ (0) into our equilibrium conditions (38) and (39). To help us update our guess, we define the following: 

$$
\hat {x} _ {1, i} \equiv \hat {y} _ {i} ^ {\frac {1 + \theta \lambda + \theta}{1 + \theta \lambda}} \hat {l} _ {i} ^ {- \theta \left(\frac {1 + \alpha + \theta \lambda (\alpha + \beta)}{1 + \theta \lambda}\right)}
$$

$$
\hat {x} _ {2, i} \equiv \hat {y} _ {i} ^ {- \theta \left(\frac {1 - \lambda}{1 + \theta \lambda}\right)} \hat {l} _ {i} ^ {- \theta \left(\frac {- 1 + \beta + \theta \lambda (\alpha + \beta)}{1 + \theta \lambda}\right)}
$$

so that: 

$$
\left( \begin{array}{c} \ln \hat {x} _ {1, i} \\ \ln \hat {x} _ {2, i} \end{array} \right) = \left( \begin{array}{c c} \frac {1 + \theta \lambda + \theta}{1 + \theta \lambda} & - \theta \left(\frac {1 + \alpha + \theta \lambda (\alpha + \beta)}{1 + \theta \lambda}\right) \\ - \theta \left(\frac {1 - \lambda}{1 + \theta \lambda}\right) & - \theta \left(\frac {- 1 + \beta + \theta \lambda (\alpha + \beta)}{1 + \theta \lambda}\right) \end{array} \right) \left( \begin{array}{c} \ln \hat {y} _ {i} \\ \ln \hat {l} _ {i} \end{array} \right) \Longleftrightarrow
$$

$$
\left( \begin{array}{c} \ln \hat {y} _ {i} \\ \ln \hat {l} _ {i} \end{array} \right) = \left( \begin{array}{c c} \frac {1 + \theta \lambda + \theta}{1 + \theta \lambda} & - \theta \left(\frac {1 + \alpha + \theta \lambda (\alpha + \beta)}{1 + \theta \lambda}\right) \\ - \theta \left(\frac {1 - \lambda}{1 + \theta \lambda}\right) & - \theta \left(\frac {- 1 + \beta + \theta \lambda (\alpha + \beta)}{1 + \theta \lambda}\right) \end{array} \right) ^ {- 1} \left( \begin{array}{c} \ln \hat {x} _ {1, i} \\ \ln \hat {x} _ {2, i} \end{array} \right)
$$

By this definition, for a guess of $\left( \{ \hat { y } _ { i } \} _ { ( 0 ) } , \left\{ \hat { l } _ { i } \right\} _ { ( 0 ) } \right)$ , the equilibrium conditions yield a set of vectors $\left( \{ \hat { x } _ { 1 , i } \} _ { ( 1 ) } , \{ \hat { x } _ { 2 , i } \} _ { ( 1 ) } \right)$ , which by the log-linear transformation defined above, imply an updated guess of $\left( \left\{ \hat { y } _ { i } \right\} _ { ( 1 ) } , \left\{ \hat { l } _ { i } \right\} _ { ( 1 ) } \right)$ . On each iteration, we rescale the vectors $\left( \{ \hat { x } _ { 1 , i } \} _ { ( 1 ) } , \{ \hat { x } _ { 2 , i } \} _ { ( 1 ) } \right)$ such that the secondperiod income and labor distribution still sum to 1, and update our guess of $\left( \{ \hat { y } _ { i } \} _ { ( 0 ) } , \left\{ \hat { l } _ { i } \right\} _ { ( 0 ) } \right)$ towards $\left( \{ \hat { y } _ { i } \} _ { ( 1 ) } , \left\{ \hat { l } _ { i } \right\} _ { ( 1 ) } \right)$ , We iterate through procedure this until the equilibrium conditions are solved up to scale. Therefore, the inner loop returns a $\left( \{ \hat { y } _ { i } \} _ { ( 0 ) } , \left\{ \hat { l } _ { i } \right\} _ { ( 0 ) } \right)$ for which: 

$$
\lambda_ {1, i} \left(\hat {x} _ {1, i}\right) _ {(1)} = \hat {\chi} _ {(0)} \left(\frac {E _ {i}}{E _ {i} + \sum_ {k} \Xi_ {i k}}\right) \left(\hat {y} _ {i}\right) _ {(0)} ^ {\frac {1 + \theta \lambda + \theta}{1 + \theta \lambda}} \left(\hat {l} _ {i}\right) _ {(0)} ^ {\frac {\theta (\beta - 1)}{1 + \theta \lambda}} + \hat {\chi} _ {(0)} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j} \left(\frac {\Xi_ {i j}}{E _ {i} + \sum_ {k} \Xi_ {i k}}\right) \hat {\tilde {t}} _ {i j} ^ {- \frac {\theta}{1 + \theta \lambda}} \left(\hat {y} _ {j}\right) _ {(0)} ^ {\frac {1 + \theta}{1 + \theta \lambda}} \left(\hat {l} _ {j}\right) _ {(0)} ^ {- \frac {\theta (1 + \alpha)}{1 + \theta \lambda}}
$$

$$
\lambda_ {2, i} \left(\hat {x} _ {2, i}\right) _ {(1)} = \hat {\chi} _ {(0)} \left(\frac {Y _ {i}}{Y _ {i} + \sum_ {k} \Xi_ {k i}}\right) (\hat {y} _ {i}) _ {(0)} ^ {\frac {- \theta (1 - \lambda)}{1 + \theta \lambda}} (\hat {l} _ {i}) _ {(0)} ^ {\frac {\theta (\alpha + 1)}{1 + \theta \lambda}} + \hat {\chi} _ {(0)} ^ {\frac {\theta \lambda}{1 + \theta \lambda}} \sum_ {j = 1} ^ {N} \left(\frac {\Xi_ {j i}}{Y _ {i} + \sum_ {k} \Xi_ {k i}}\right) (\hat {y} _ {j}) _ {(0)} ^ {- \frac {\theta}{1 + \theta \lambda}} (\hat {l} _ {j}) _ {(0)} ^ {\frac {\theta (1 - \beta)}{1 + \theta \lambda}}
$$

To solve the equilibrium conditions, we would like to have a $\hat { \chi }$ for which $\{ \lambda _ { 1 , i } \} = 1$ and $\left\{ \lambda _ { 2 , i } \right\} = 1$ . Therefore, in the outer loop we implement a version of the fmincon function, a non-linear constrained minimization function built into MATLAB, which finds a $\hat { \chi } _ { ( 1 ) }$ which minimizes the matrix norm of $\left[ \left\{ \ln \lambda _ { 1 , i } \right\} \ \left\{ \ln \lambda _ { 2 , i } \right\} \right]$ We use the output of fmincon to update our guess of $\hat { \chi }$ , and we iterate on this outer loop until $\hat { \chi } _ { ( 0 ) }$ converges with ˆχ(1). $\hat { \chi } _ { ( 1 ) }$ 

# F Data Construction

In this section, we provide details on the construction of the data used in the paper. 

# F.1 The U.S. Highway System Network

Our version of the interstate highway system consists of road segments, which are parts of an interstate which pass through a county. Each segment has data on mileage, traffic flows, and through lanes over the entirety of the segment, as well as over subsections of the segment which fall into different road type classifications under the Federal Highway Administration’s (FHA) Highway Economic Requirements System (HERS). To generate this dataset, we combine GIS data on the interstate highway network released by the FHA, elevation data from NASA’s Shuttle Radar Topography Mission (CGIAR-CSI, 2017), population data on and geographies of urban areas from (Edwards, 2017), Commodity Flow Survey Areas, and Census-Based Statistical Areas (CBSA) released by the Census Bureau and sourced through IPUMS NH-GIS (MPC, 2011), and trade information from the CFS (CFS, 2012). 

# F.1.1 Interstate Highway System & Traffic Data

Since 2011, the FHA has released shapefiles of the national road network, which are linked to data collected from the Highway Performance Monitoring System (HPMS). Critically, this shapefile contains information on the state and county a road segment is in, the roads’ Department of Transportation functional system classification, a road’s name or route number, whether it is in an urban or rural area and which urban area it is in, the mileage of a segment, the average annual daily traffic (AADT) that goes over it, and how many lanes the segment has. 

To generate our version of the interstate highway system, we begin with the HPMS release from 2012, trim it to only those roads within the contiguous United States, and remove all roads that are not classified as Interstates, resulting a road network with 334,040 segments. Beyond this, we resolve several data quality issues within the 2012 HPMS release. The most visible of these issues is the road reports from West Virginia, which are missing large sections of several interstate highways. To resolve these, we replace the reported data on West Virginia interstates from 2012 with that from 2013. Additionally, because the national HPMS dataset is sourced from reports prepared by state departments of transportation, there are some discrepancies in how roads are labeled. For instance, for Interstate Route 10, one department of transportation might code it as 10, another as I10, another as I-10, and another might use an entirely different identification system altogether. In order to consistently code each road segment by its integer route number (i.e. Interstate 10 as “10”), it was necessary to recode interstate segments in Arizona, California, Nevada, and Rhode Island, with reference to Google Maps. Finally, the HPMS has some road segments which are classified as interstates but are not coded with a route number at all. These “zeroes” come from three states: Alabama, Maine, and New York. We delete those in Alabama, being exceptionally small in length (<0.1 miles), and those in Maine, which are short ramps, insignificant to the overall connectivity of the road network. The “zero” in New York is a transitional ramp between I-90 and I-87 near Albany and is recoded as part of I-90. The resulting dataset consists of 333,021 road segments. 

Finally, notice that the CFS data are intended as end-to-end transactions in the sense that they are shipments from one establishment to another. The census explicitly asks in their questionnaires“An outbound shipment in this survey is defined as a movement of commodities from the location specified in Item A to another single location. For shipments traveling to multiple locations, report one shipment for each location”.49 Notice that to avoid counting third-party shipment the census indicates: “Certain wholesalers (i.e., manufacturers’ sales offices, agents, and brokers) and certain importers (i.e., own brand importers and marketers) were excluded from the frame. These wholesalers do not maintain inventory at their office location but rather arrange for products to be shipped to a buyer from some other location”.50 

# F.1.2 Road Type Classifications

Using the cleaned version of the interstate highway system, we proceed to match each segment of the highway system to its underlying terrain. An elevation raster of the United States was sourced from DIVA-GIS, which gathers the underlying elevation data for the raster from the CGIAR’s Center for Spatial Information SRTM 


Table F.1: Interstate Highway System: Cost of Adding a Lane-Mile


<table><tr><td>Road Type</td><td>Cost ($m)</td></tr><tr><td>Rural - Flat</td><td>1.923</td></tr><tr><td>Rural - Rolling</td><td>2.085</td></tr><tr><td>Rural - Mountainous</td><td>6.492</td></tr><tr><td>Small Urban</td><td>3.061</td></tr><tr><td>Small Urbanized</td><td>3.345</td></tr><tr><td>Large Urbanized</td><td>5.598</td></tr><tr><td>Major Urbanized</td><td>11.197</td></tr></table>


Source: Federal Highway Administration (2015) 


data. Using this elevation raster and ArcGIS’ 3D Analyst toolkit, we extract the average grade of each segment of the interstate highway system. 

Each road segment is then matched with the population of the urban area it passes through based on its urban code in the HPMS data, which are the same codes that the Census Bureau uses to identify its urban areas and urbanized centers. Populations for urban areas are sourced from the 2010 Census Urban-Rural Classification, which was released in March 2012. The urban area codes in use in the HPMS data differ slightly from those in the census release, so a handful of urban areas in the HPMS data are recoded to match their updated codes in 2010 Census Urban-Rural Classification. Based on this terrain and population data, each segment was then classified into one of the seven HERS urban-rural road types below. 

Urban road segments are classified based on population, per standards outlined by the FHA in its HPMS field manual Federal Highway Administration (2016). Rural road segments, which are all segments which pass through areas of less than 5,000 in population, are classified based on the average grade. The FHA offers only general guidance on how to classify roads by terrain. Based upon the guidance that Level terrain “generally includes short grades of no more than 2 percent” Federal Highway Administration (2016), all roads of grade below 2% were classified as Level, and based upon the maximum grade for Interstate Highways going over rolling terrain with a speed limit of 60 mph being set at 4% American Association of State Highway and Transportation Officials (2016), all roads of grades between 2% and 4% were classified as Rolling. The remainder of roads (those over 4% in grade) were classified as Mountainous. 

For each section, a measure of its length is generated by subtracting the mileage marker of its endpoint from the mileage marker of its beginning point. Then, we generate a measure of vehicle miles traveled (VMT) by multiplying AADT by the length of the segment and a measure of lane-miles over the segment by multiplying through lanes (the total two-way lane width of a road) by the length of the segment. Each of these three measures—length, VMT, and lane-miles—is also interacted with the seven dummy variables that code road type. 

# F.1.3 Observed Network of the Interstate System

To simplify the geometry of the interstate network, we aggregate road segments based on their state, county, and route number, summing length, VMT, lane-miles, and all road type-interactions. This reduces the number of road segments from 333,021 to 1,761. Finally, we join segments within a radius of 3500 meters of each other together. This links together geometries which were either not precisely connected in the shapefile or were connected by shorter roads not coded as interstates in the HPMS dataset, and therefore removed in the initial data cleaning. This dataset forms the links of our interstate network. 

For our network analysis of the interstate system, we choose to place nodes at every intersection between two different interstates and endpoint in the interstate highway system. This results in a set of 616 nodes, each of which we geocode with its latitude, longitude, distance from and name of nearest CBSA, and CFS area. We also identify adjacent nodes as nodes which can be reached from each another without passing through another node on the network. 

# F.1.4 Simplified Network of the Interstate System

In ArcMap, we set the nodes as origins and destinations in a symmetric origin-destination matrix; identify length, VMT, lane-miles and all road type-interactions as accumulation attributes (that is, fields that the ArcMap Network Analyst toolbox should integrate over as it calculates the least-cost route); and set length as the impedance, which is the field that the Network Analyst toolbox minimizes to identify the least-cost route. Solving this symmetric origin-destination matrix between 616 origins and 616 destination leads to 379,456 bilateral connections. Each of these bilateral connections contains data on the length of the route, total VMT over the route, total lane-miles over the route, and the length, total VMT, and total lane-miles over sections of the route which fall into each of the seven road types. 

Using the CBSA information geocoded to the nodes, we code a node as being within its nearest CBSA if it is less than 3000 m away from the boundary of that CBSA. This is necessary to address the presence of nodes which just barely fall outside of a CBSA. The set of bilateral connections generated by solving the origin-destination matrix has several “redundant” connections because our approach to generating the nodes of the interstate network yields clusters of nodes near cities, with large ring roads or through which many interstate highways pass. To eliminate these redundancies, we consolidate all nodes coded to the same CBSA into one “CBSA node,” coded with the average latitude and longitude of all nodes within that CBSA and its relevant CBSA and CFS area. The bilateral connections from this “CBSA node” to all other nodes in the interstate network contain the average distance units (length, VMT, lane-miles, and road type interactions) for each unique connection from that CBSA to another node. This process of consolidating nodes within CBSAs yields a simplified adjacency matrix, where the clusters of nodes around major cities due to ring roads are absorbed into one node. The 616 nodes and 379,456 links of the original OD matrix are reduced down to 228 nodes and 51,984 links, of which 704 are between adjacent nodes. 

# F.1.5 Estimated Cost of Improvements and Congestion

For the simplified network, we calculate the cost of adding an additional lane-mile along an link by identifying the share of each link that goes over each road type and using those shares as weights in a weighted average of the different cost figures for adding a lane-mile in each terrain type estimated by the Federal Highway Administration (2015). To identify congestion measures, we divide the total VMT along an link by the total lane-miles along that link. This gives us a measure of traffic per lane-mile over the course of the road. 

# F.1.6 Consistent Measures of Node Population and Income

To generate a consistent measure of population and income at each node, we sum the population of and average the median income of all cities within 25 miles of a node, conditional on a city within that radius not being closer to another node. We name each node (for readability) after the city with the largest population in the aforementioned 25 mile radius. Population and median income data come from a purchased dataset from USCitiesList.org. Although consistent, this way of measuring population naturally tends to understate populations for less densely populated areas. 

# F.1.7 Trade Flows: Observed and Imputed

Using the CFS area coded to each node, we link links with trucking flows aggregated to the origin-destination level from the 2012 CFS. CFS areas are generally larger than CBSA’s, so to get more granular trade flows, we imputed commodity flows between CBSAs by assuming that, for CFS areas which consist of more than one CBSA, each CBSA receives and sends out a portion of flows proportional to its share of the CFS area’s total GDP. Both the observed CFS area-CFS area flows and imputed CBSA-CBSA flows are included in the data. 9,801 of the 51,984 links are linked to CFS flows; 9,651 of those links can be further disaggregated into CBSA-CBSA flows. 

# F.2 The Seattle Road Network

Our version of Seattle’s commuting network combines the road system reported in the Seattle-Tacoma-Bellevue CBSA in the FHA’s 2016 Highway Performance Monitoring System (HMPS, 2016) for the state 


Table F.2: Seattle Road Network: Cost of Adding a Lane-Mile


<table><tr><td>Road Type</td><td>Standard Cost ($m)</td><td>High Cost ($m)</td></tr><tr><td>Freeway/Interstate</td><td>11.197</td><td>46.691</td></tr><tr><td>Other Principal Arterial</td><td>8.252</td><td>31.988</td></tr><tr><td>Minor Arterial/Collector</td><td>5.614</td><td>31.988</td></tr><tr><td>Local*</td><td>5.614</td><td>N/A</td></tr></table>


*Local costs imputed as identical to Minor Arterial/Collector costs 



Source: Federal Highway Administration (2015) 


of Washington with commuting flow data from the 2017 LEHD LODES release from the Census Bureau (LODES, 2017). We also, similarly, to the our version of the interstate highway network rely on elevation data from the SRTM (CGIAR-CSI, 2017) and population data for urban areas and Census block groups released by the Census Bureau (MPC, 2011). 

# F.2.1 Local Road Data

We trim the data from the 2016 HPMS release for the state of Washington to cover all roads within the municipal boundaries of Seattle as specified by (Seattle, 2017), creating a dataset of 9,188 road segments. This dataset contains information on a road segment’s Department of Transportation functional system classification, which authority owns it, how many lanes it has, whether additional lanes could be easily added, and the AADT that flows over it, and whether it is a ramp or not. 

There are a handful of roads for which one or several of these datapoints are blank, so we impute those based upon the features of the surrounding roads. For functional system classification and ownership, which are both categorical variables with an associated hierarchy, we fill in these blanks with the “highest” level of the hierarchy that a road comes into contact with; for example, if a road touches an interstate highway and a minor arterial, it is classified as belonging to the interstate functional classification, and if a road touches a road that is owned by the state and another that is owned by the county, that road is classified as being owned by the state. Generally, for lane width, we fill in blanks with the maximum lane width among roads that a segment comes into contact with, and for traffic flows, we fill in blanks with the mean traffic flows among roads that a segment comes into contact with. For a subset of roads which happens to have blank lane widths and traffic flows because they represent the other lane of a dual lane road way with that data, like a large highway or boulevard, we simply impute traffic flows and lane width from the its parallel counterpart. The only other exception to the aforementioned general rule is surface streets, where we impute that each surface street has width of two lanes and an AADT of 120. Finally, we used geoprocessing tools in ArcGIS to fix connectivity issues in the road data. 

Taking the cleaned traffic data, we prepared the road network for network analysis by aggregating road segment-level data to what we call the “road section”-level—road sections being defined as a continuous segment of road not interrupted by an intersection with another road. This increases the number of road segments in our data to 17,261. We used ArcGIS tools to measure the length of each road section and code it to the Census urban area it belongs to. We further classified each road as being a high cost road to add a lane to or a not, based on whether the HPMS release says that additional lanes could be easily added to it. Using a road section’s functional system classification and its high-cost classification, we coded each road section with the cost of adding a lane-mile to it, based on the costs estimated by the FHA Federal Highway Administration (2015). 

Finally, for road sections which are missing posted speed limits, we filled in speed limits based on whether it is a ramp (ramps are given a speed limit of 30 mph) and who owns the road. Washington state has default statutory speed limits for city and town streets, county roads, and state highways Washington State Legislature (1965). Local roads are assigned a speed limit of 20 miles per hour, as per local Seattle traffic regulations. From these road type classifications, the length of the road, and the traffic and through-lane capacity data from the HPMS release, we created measures of VMT, lane-miles, improvement cost—generated 

by multiplying the cost of adding a lane-mile to a road by its length—and unimpeded travel time—generated by dividing the length of the road by its speed limit. 

# F.2.2 The Observed Road Network

We convert this local data into the observed road network by generating a set of nodes for network analysis, solving for the least-cost path between those nodes using ArcMap’s Network Analyst tool, and paring the resulting set of bilateral paths down to only those between adjacent nodes. 

To generate the nodes, we grid the city of Seattle into 224 1 km x 1 km parcels and set the center points of those parcels as nodes for network analysis. We restrict participation in the network analysis to only those nodes which are within a third of a kilometer of the road network; we find that this distance restriction does a good job of resolving the tradeoff between capturing the overall structure of Seattle’s roads and limiting the nodes to those within a reasonable distance of the road system. Overall, 217 nodes participate in the network analysis. 

Using the OD Matrix feature of the Network Analyst tool, we solve for the path which minimizes unimpeded travel time from each of the 217 nodes to all the other nodes. Along each path, we also sum over VMT, lane-miles, improvement cost, the number of intersections crossed, the number of turns taken, and the length traveled along arterial and local roads. While the network dataset has highly detailed data on all the roads in Seattle, this step means that we observe only those roads along which at least one least cost route between nodes travel. Solving this optimization problem yields a dataset with 47,089 bilateral connections between nodes. We define a node as adjacent to another node if it is in one of eight parcels which surrounds the other node’s parcel and is not separated from the other node by a body of water, without being connected by a bridge.51 There are 1,384 bilateral connections between adjacent nodes in our dataset. 

# F.2.3 Node Populations and Incomes

Since our nodes are not linked to any existing administrative dataset on populations and incomes, we need to generate population and income figure for each one. To do so, we use an ArcMap tool to identify the area of each intersection between a block group and a parcel and calculate the share of each parcel’s area that comes from a particular block group. We also calculate the population density of each block group. Then, assuming that the population of each block group is uniformly distributed within that block group, we estimate the population density of each parcel by finding the weighted average of the block groups it overlaps with, where the weights are provided by the share of each parcel’s area that comes from a block group. We calculate the per capita income of each parcel using a similar method. Further assuming that the residential and working population within each parcel is uniformly distributed, we calculate the total residential and working population of each parcel by multiplying the respective population density by its area. 

# F.2.4 Commuting Flows

We take commuting flow data from the LEHD LODES dataset (LODES, 2017). We narrow down the commuting flows, which originally are at the census block-to-census block level, to only those flows which begin and end in the Seattle Metro Area. Using the areas of the intersections between each block group and parcel, we calculate the share of each block group’s area that falls into a particular parcel. We distribute a block group’s residents and labor force among its intersections according to these shares, and then aggregate those intersections up to the parcel level to find the number of residents and the size of the labor force for each parcel. 

To distribute commuting flows, we estimate the commuting flows between block group-parcel intersections by multiplying total commuting flows between two block groups by the shares of the origin block group and the destination block group taken up by the origin intersection and the destination intersection. We then aggregate these commuting flows up to the parcel-to-parcel level. 


Figure F.1: Route Complexity: An Example


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/ddd7641805cbe1af22b9f0834b227721b3a7d813be7f4be332820bdda0e5634b.jpg)


Notes: This figure provides an example of how the instrument for traffic based on the route complexity is constructed for the Seattle road network. On this link, there are five turns and 19 intersections. 

# F.2.5 Instrument Construction

For our IV estimation of the congestion parameter in Seattle, we rely upon the number of turns along a route, conditional on the number of intersections traversed and origin and destination fixed effects, as an instrument. We define any deviation from the current bearing of the route by more than 30 degrees, in either direction, as a turn, and we use the Global Turn Delay within ArcMap’s Network Analyst to count the number of turns along the least-cost route between two nodes. We also use the Network Analyst to count the number of intersections traversed. The below figure presents an example of how this process works and what kind of data it results in. Between these two nodes, the least-cost path makes five turns and traverses 19 intersections. 

# F.3 Trade, Commuting, and Distance

# G Alternative Parameter Constellations

In the section, we compare the estimated welfare elasticities and returns on investment for each segment of the U.S. highway network and the Seattle road network to equivalent results under three different parameter constellations: (1) no externalities ( $\alpha = \beta = 0$ ); (2) lower trade elasticity ( $\theta = 4$ ); and (3) greater traffic congestion, which we calculate by estimating $\delta _ { 0 }$ from a gravity regression of either trade or commuting on travel times. We summarize the results in three figures, corresponding to each of the three parameter constellations. 


Figure F.2: Trade, Commuting, and Distance



(a) U.S. Trade Flows and Distance: Data and Calibrated Model


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/2e00a7ce2884b9a4c4c4592a29e94a61da755f7bfd44a182d185cccfaf788420.jpg)



(b) Commuting Flows and Distance in Seattle: Data and Calibrated Model


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/9b41e4208aa8856e10c9f94e7575ab76887b0e9a21eecf7376205729f8e0c9c7.jpg)


Notes: This figure compares the share of trade and commuting by pair of locations for different distance bins in our data to the one predicted by our model using our preferred parameter constellation. 


Figure G.1: Alternative parameter constellation: No externalities



(a) U.S. Highway Network


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/ffdff656f2f48ee027ecf9ae1f08bffd974d620958093abdc5bb763cce53c6ba.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/a8b14f6f69b34e2d0e064ca2bbacbfbda8b1a7e80a7a107fd5d3c52d0f5d06e4.jpg)



(b) Seattle Road Network


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/f9d68737b4f58079a331f45039105367315324f02050e23a1cacbf09ae8eecd9.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/4829b67761fd7b659760a970fe6b43d24dbb3b9dcce295894a52056844c6c465.jpg)


Notes: This figure compares the welfare elasticity (on the left) and return on investment (on the right) elasticity for each link in the U.S. highway network (panel a) and the Seattle road network (panel b) calculated using our preferred parameter constellation (on the x-axis) to an alternative parameter constellation where we assume no externalities, i.e. $\alpha = \beta = 0$ , (on the y-axis). 


Figure G.2: Alternative parameter constellation: Lower gravity elasticity



(a) U.S. Highway Network


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/2db4d8076543ef63527461c37860e5a9c1fc888118f41d1a40d31e4066e3ec43.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/ac1a5b02168a4f1c9251585f794225f0115913de5ccb8125a7fdbdf780e71bf1.jpg)



(b) Seattle Road Network


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/a1da485646fc3d38ca82ace2991253d39c34a3fe93f7bd3fd271d80107103a5a.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/325ce1be3d7aa34ab13fcbaec6ef0e1ff3ccf2a6b235097f20014a9c60dcc6b4.jpg)


Notes: This figure compares the welfare elasticity (on the left) and return on investment (on the right) elasticity for each link in the U.S. highway network (panel a) and the Seattle road network (panel b) calculated using our preferred parameter constellation (on the x-axis) to an alternative parameter constellation where we assume a smaller gravity elasticity, i.e. $\theta = 4$ , (on the y-axis). 


Figure G.3: Alternative parameter constellation: Greater traffic congestion



(a) U.S. Highway Network


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/06550b44b7a3d12cf89b374b1c9359a7565373b1fca62f4292d4e8f8f902debd.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/0ed9fa69f01062567af7a141eb6117a488976314afc983410d0e798857bf8b13.jpg)



(b) Seattle Road Network


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/bc14b47ed891ee0c640f36139de7ca970102d5fe19f635aebdf673c345e13b36.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-03-08/cf2f931b-5b7c-412d-a282-458c1e2a5d04/281be33923af754f3860575db49a4359b7cd4d556b8243ff1bbd17d6f34086b5.jpg)


Notes: This figure compares the welfare elasticity (on the left) and return on investment (on the right) elasticity for each link in the U.S. highway network (panel a) and the Seattle road network (panel b) calculated using our preferred parameter constellation (on the x-axis) to an alternative parameter constellation where we assume greater traffic congestion, i.e. higher $\lambda$ , (on the y-axis). The $\lambda$ is calculated here by estimating $\delta _ { 0 } \theta$ based on a gravity regression of trade (panel a) or commuting (panel b) flows on travel times, rather than setting $\delta _ { 0 } \theta = 1$ , as above. 