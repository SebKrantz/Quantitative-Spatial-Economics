# Quantitative Regional Economics*

Treb Allen and Costas Arkolakis 

January 2025 

# Abstract

This handbook chapter presents the major advances made in the field of economic geography over the past decade. It starts by documenting a number of motivating empirical facts. It then shows how a quantitative regional model that combines the insights from two seminal models from an earlier generation can explain these facts. It then presents a unified quantitative framework that incorporates this and many other economic geography models. This unified framework is sufficiently tractable to characterize its equilibrium properties while flexible enough to be combined with detailed spatial economic data to estimate the model parameters, conduct counterfactuals, and perform welfare analysis. The chapter concludes by discussing many extensions of the framework, some of which have already been explored and others which have not. 

# 1 Introduction

Space matters. Life in rural New Hampshire is different than life in New York City: it is cleaner, quieter, and less crowded. But the wages are lower and—unless you are an avid outdoors person—the amenities are worse. Why do some people choose to live in New Hampshire while others choose New York City? Why are wages higher in New York City? And how much of what we observe about the spatial distribution of people and economic activity today is due to innate geographical differences (e.g. New Hampshire has mountains, whereas New York City is a natural harbor) versus geographic location (e.g. New Hampshire is far from most major cities, New York City is not) versus other economic forces (e.g. agglomeration)? 

These questions form the bedrock of the study of regional economics. And while these topics have been studied for hundreds of years (see Smith (1776), Von Thunen (1826), and Marshall ¨ 

(1890)), it took until the 1980s for the field to develop formal mathematical models for characterizing the equilibrium distribution of economic activity across space with the work of Roback (1982) and a few years after with the work of Krugman (1991). And it was not until the past decade that these distinct mathematical frameworks have become unified into a single “quantitative” framework where the geography being considered is sufficiently realistic that it can be combined with real world data to provide precise answers to these questions. But the power of this new quantitative framework goes beyond explaining why the world looks the way it does: it also allows us to predict how changes in the world (say, to the cost of interacting across space) would change the equilibrium distribution of economic activity. In doing so, it opens up the possibility of using the framework to design better spatial policies. 

The purpose of this handbook chapter is to provide an overview of this new quantitative regional framework. While we situate the framework relative to earlier economic models and try to highlight several exciting paths toward which the framework may evolve, the purpose of the chapter is not a literature review. Instead, our goal is to offer the reader a stand alone toolkit for using the quantitative regional framework in their own research. Toward that end, the focus of the chapter is two-fold: first, we try to offer the reader sufficient technical detail in order to see how (and why) the framework works the way it does without overwhelming with technicalities; second, we try to discuss how the framework can be combined with data in order to perform counterfactual analysis. 

The chapter is organized as follows. We first motivate the key ingredients of the framework with a set of real world facts (Section 2). We then briefly summarize two seminal models developed in the 1980s and 1990s upon which the modern quantitative regional economic framework was built (Section 3). We next present the simplest version of the quantitative regional economy model and discuss a number of its properties (Section 4) before showing how this simple model has been extended in a number of different directions (Section 5). Section 6 presents a workhorse quantitative regional framework, for which all earlier models are special cases. Here we present the framework, characterize its equilibrium properties and derive expressions for its counterfactuals. Section 7 then provides a detailed description of how to combine the framework with spatial data to both estimate key model parameters and perform those counterfactuals. Finally, Section 9 discusses a number of exciting recent extensions to the framework before Section 10 identifies a number of interesting avenues for future research. 

Before proceeding, we also offer some advice on using this handbook chapter. Sections 2 and 3 are helpful for motivation and to see the genesis of the contemporary quantitative regional economic framework. Sections 4 and 5 offer particular micro-foundations for the workhorse framework. Section 6 is meant to offer a stand-alone description of the modern quantitative regional framework. This section presents the framework, characterizes its equilibrium properties, and derives expressions for its counterfactuals. Section 7 is then its empirical counterpart about how to take the framework to the data. Section 8 summarizes some of the lessons learned 

from the literature employing quantitative regional models, while sections 9 and 10 offer our opinions of ideas of how the field is evolving and may evolve in the future. We have separated more technical portions into appropriately labeled subsections, which can be skipped if the reader is interested in the results rather than the derivations themselves. 

Finally, we want to offer a brief note on the terminology that will be useful for this and other chapters of the book as well as for navigating the rapid progress in the field. With the term quantitative, we refer to frameworks that incorporate rich interactions between individuals across locations, which in turn enables a tight connection between economic theory and spatial data. We refer to the focus of this chapter—which summarizes the recent quantitative advances on the network of connections between cities and regions—as quantitative regional economics. We refer to the theme of its partner chapter Redding (2025)—which summarizes the recent quantitative advances on the network of connections within a city—as quantitative urban economics. We refer to the combination of regional economics and urban economics as economic geography. Finally, spatial economics comprises both economic geography and international economics, entailing the analysis of economic activity at any level and under any specification for the interactions across geographical space. 

# 2 Economic geography facts

We begin by offering three facts that motivate the need for a quantitative economic geography framework. 

# 2.1 Fact #1: Economic activity is distributed highly unequally across space

Economic activity across space is markedly uneven. Figure (1) depicts the spatial distribution of economic activity across counties in the United States (U.S.). The left panel depicts the population density, i.e. the number of workers per square kilometer. The population density differs by several orders of magnitude, with population densities below one worker per square kilometer in much of the mountain west and the population densities exceeding 1,000 workers per square kilometer along the eastern seaboard. 

The right panel depicts the distribution of economic output per worker, i.e. the average annualized payroll per worker. Here too there is substantial spatial variation. Many counties around city centers, such as San Francisco, Manhattan, Seattle, and others, have an annual payroll of more than $\$ 100,000$ while in counties a few dozens kilometers away the annual payroll quickly falls below $\$ 50,000$ . Because more highly populated locations have greater economic output per worker, the right panel demonstrates that total economic activity is even more spatially concentrated than the population density. 


Figure 1: The spatial distribution of economic activity



(a) Workers per square kilometer


![image](https://cdn-mineru.openxlab.org.cn/result/2026-02-25/4e4c25dd-e494-4502-8cd1-6a5c0f2e691d/ddf37f18010c30c9605a193f26945326f8f2620767fa46677f58062e7b178537.jpg)



(b) Annual payroll per worker


![image](https://cdn-mineru.openxlab.org.cn/result/2026-02-25/4e4c25dd-e494-4502-8cd1-6a5c0f2e691d/d10aff540c1fbb89ca278697e9f9beef12cb6a79d3d7336eaa7dd8f6bd6bdaeb.jpg)


Notes: The left panel depicts the number of workers per square kilometer; the right panel depicts the annual payroll per worker. The source for both panels is the 2017 U.S. County Business Patterns from the U.S. Census Bureau. 

One possible explanation for the substantial spatial variation in total economic activity is that we are observing an economy in transition: for example, perhaps the population density is the largest along the eastern seaboard is simply due to those areas being older. Figure 2 shows the distribution of the fraction of total workers (left panel) and total output (right panel) across U.S. counties in both 1986 and 2017. As is evident, the spatial distribution of both population and output has remained remarkably consistent over the 30 year period, suggesting that the substantial spatial inequality we observe is not a temporary phenomenon. 

This fact motivates the need for a framework that can explain the persistent presence of substantial spatial variation in both where people live and how much they earn. 

# 2.2 Fact #2: Regional trade is important

Our second stylized fact is that regions within the U.S. are strikingly economically integrated. To show this, we use the 2017 Commodity Flow Survey (CFS) Public Use File to calculate the fraction of observed spending by a Metropolitcan Statistical Area (MSA) on shipments that originated in the same MSA, i.e. its “own expenditure share.” There are 70 unique MSAs in our data. An own expenditure share of 1 indicates that an MSA is entirely in autarky, whereas a value of 0 indicates that an MSA only buys goods from elsewhere. Figure 3 depicts the spatial distribution of own expenditure shares across MSAs. As is evident, nearly all MSAs spend a majority on shipments from other MSAs, i.e. across-MSA trade comprises a majority of economic activity. That being said, own expenditure shares vary considerably. Major MSAs, such as New York or Miami have much larger own shares, close to a $4 0 \%$ or more. On the contrary, many small MSAs in the Midwest have own shares closer to 0.25. Geographic location appears 


Figure 2: The spatial distribution of economic activity over time



(a) Workers


![image](https://cdn-mineru.openxlab.org.cn/result/2026-02-25/4e4c25dd-e494-4502-8cd1-6a5c0f2e691d/03f623d504cf44c38ad03999e8fd753638c31834c22b4a6cfe7d9e8df4567eed.jpg)



(b) Payroll


![image](https://cdn-mineru.openxlab.org.cn/result/2026-02-25/4e4c25dd-e494-4502-8cd1-6a5c0f2e691d/562db540d4720d317da70689c7e90a7eeed32bd6ca9f77cbe825f5c18f0116f6.jpg)


Notes: This figure depicts the distribution across counties in employment (panel a) and output (panel b) for both 1986 and 2017. We calculate the fraction of total employment / output in each county, divided by overall employment and output, respectively, and then plot the resulting distribution across counties. The source for both panels are the 1986 and 2017 U.S. County Business Patterns from the U.S. Census Bureau. 

to matter too: One of the MSAs with the highest own share is Honolulu, in Hawaii, nearly four thousand kilometers away from an MSA in the mainland USA. 

This fact motivates the need for a framework where different regions interact economically through the flow of goods. 

# 2.3 Fact #3: Gravity holds

Our third and final stylized fact is that where a region is located and how large it is plays an important role in determining to extent to which it interacts with other regions, i.e. gravity holds. Figure 4 reveals that the “gravity” trade pattern very familiar to international trade economists (see e.g. Anderson (1979), Eaton and Kortum (2002), and Anderson and Van Wincoop (2003)) also holds within countries: trade systematically declines with distance. We regress the common logarithm of the trade value between all the pairs of MSAs in our data to routed distance bins of the pair, 0-250 kilometers (km), 250-500, 500-1000, 1000-2000, 2000-4000, or more than 4000 and origin and destination fixed effects. We provide the raw observations on the left panel and the coefficients of the regression in the right panel. The left panel shows that trade (measured as the ratio of exports to all shipments to the destination MSA) is greatest within-MSA (green), followed by trade between different MSAs that belong to the same state (blue), and then unrelated MSAs. Trade flows decline substantially with distance; for example, locations that are 250-500 kilometers away have about $1 0 \%$ $( 1 0 ^ { - 1 } )$ the value of trade as locations that are less than 250 kilometers away, while locations that are away more than 4000 kilometers away only about $1 \%$ . 


Figure 3: Own expenditure shares



(a) Spatial variation in own expenditure shares


![image](https://cdn-mineru.openxlab.org.cn/result/2026-02-25/4e4c25dd-e494-4502-8cd1-6a5c0f2e691d/782224f2f8a486be532bd45a01cfec24e4a0d7c8ea8e722756b6b67d7e6016c2.jpg)



(b) Distribution of own expenditure shares


![image](https://cdn-mineru.openxlab.org.cn/result/2026-02-25/4e4c25dd-e494-4502-8cd1-6a5c0f2e691d/45f2cf3f3819d373149f7923715806e776ae38b66d33109260fa9954a9351575.jpg)


Notes: Own expenditure share is calculated as the ratio of the value of within-MSA transactions to the total value of MSA transactions using the 2017 Commodity Flow Survey Public Use File from the U.S. Census Bureau. 


Figure 4: Distance matters



(a) Raw data


![image](https://cdn-mineru.openxlab.org.cn/result/2026-02-25/4e4c25dd-e494-4502-8cd1-6a5c0f2e691d/c7d58921c7064ffdf4d45e6c6bd1833156f084f7481c50fef41e0061673828ba.jpg)



(b) Gravity regression


![image](https://cdn-mineru.openxlab.org.cn/result/2026-02-25/4e4c25dd-e494-4502-8cd1-6a5c0f2e691d/3d905db174a3a5c30052afd730d25781f37e39d2d91c72c3566b7a93f50b75da.jpg)



Distance Bins


Notes: This figure uses the 2017 Commodity Flow Survey Dataset Public Use File to show that trade flows within the U.S. decline with distance. In the left panel, we present the scatter plot of the ratio of the value of bilateral trade flows to total MSA spending and the average distance between them. In the right panel, we report the results of a gravity regression of (log) trade flows on distance bins with origin and destination fixed effects. We construct bilateral trade flows and kilometric routed distance using the almost 6 million individual transactions provided by the Survey and aggregate them to bilateral MSA transactions using the weights provided by the CFS. 


Figure 5: Size matters


![image](https://cdn-mineru.openxlab.org.cn/result/2026-02-25/4e4c25dd-e494-4502-8cd1-6a5c0f2e691d/e41979903e7335532075244883e1c8d5eb4af1e0569e55ef7830b4c37eeeb25b.jpg)


Notes: This figure uses the 2017 Commodity Flow Survey Dataset Public Use File to show that trade flows within the U.S. increase with partner size. Each point represents the total trade between the New York MSA and another MSA; blue dots are exports from New York, red dots are imports into New York. We construct both average kilometric routed distance and bilateral trade value using the weights for each transaction provided by the CFS. GDP data at the MSA level are from the Bureau of Economic Analysis. 

Size matters too. Following the spirit of the analysis in Head and Mayer (2014), Figure 5 plots the total trade value of trade (of both imports and exports) between the New York MSA and the Gross Domestic Product of its partner MSAs. Both incoming and outgoing trade increases systematically with the size of the partner of New York MSA. 

This fact motivates the need for a framework where the degree to which different regions interact through trade depends on their respective sizes and the distances between them. 

# 3 Seminal models of the previous generation

We now present the two seminal models of economic geography from the 1980s and 1990s: the Rosen (1979)-Roback (1982) model of location choice and the Core-Periphery model developed by Krugman (1991). We first sketch out each model, then discuss their key insights, before 

noting the limitations of each framework, particularly in regard to explaining the stylized facts above. 

# 3.1 The Rosen (1979)-Roback (1982) model of location choice

We consider first the celebrated model of location choice developed by Rosen (1979) and Roback (1982). Our characterization of the model emphasizes the roles agglomeration and dispersion forces play in determining the equilibrium system of locations in a way reminiscent of Henderson (1974). This framework has found wide applicability in the urban and labor literature. For a more comprehensive review of the framework see Glaeser and Gottlieb (2009) and Moretti (2011). 

# 3.1.1 Model sketch

Consider a world of $N$ locations, which we index with $i$ and $j ,$ and let $\mathcal { N } \equiv \{ 1 , . . . , N \}$ be the set of locations, each of which produces an identical and costlessly traded numeraire good. Suppose there is a measure $\bar { L }$ of identical and perfectly mobile agents. Agents choose in which location to live (i.e. produce and consume the numeraire good). 

Agents who choose to live in location $j \in \mathcal N$ receive welfare: 

$$
W _ {j} = w _ {j} \times \bar {u} _ {j} \times L _ {j} ^ {\beta}, \tag {1}
$$

where $w _ { j }$ is their wage, $\bar { u } _ { j }$ is the innate amenity value of residing in location $j , L _ { j }$ is the measure of agents residing in location $j ,$ and $\beta < 0$ reflects a negative congestion externality (i.e. the amenity value of residing in location $j$ is lower the more people residing there). We note that treating the congestion force as a non-pecuniary externality is not an innocuous assumption; for example, if congestion were to instead arise through market forces like a land market, then the distributional consequences of the model would in general depend on the assumed ownership structure of the land. 

Free mobility implies that welfare in all locations is equalized, i.e. there exists a scalar $W > 0$ such that $W _ { j } = W$ for all $j \in \mathcal N$ . Combining welfare equalization with the welfare expression (1) allows us to derive the following (inverse) labor supply curve: 

$$
L _ {j} = w _ {j} ^ {- 1 / \beta} \times \bar {u} _ {j} ^ {- 1 / \beta} \times W ^ {1 / \beta}. \tag {2}
$$

Because $\beta < 0$ , equation (2) shows that labor supply is upward sloping in the wage of a location with constant elasticity $- 1 / \beta ,$ with the amenity of a location acting as a $( \log )$ shifter of the supply curve. This is intuitive: the greater the innate amenity of a location, the more the labor supply curve is shifted outward, as residents need not be compensated with as high of a wage 

in order to choose to live there. 

To produce the numeraire good, the Rosen (1979)-Roback (1982) model assumes that labor and an additional fixed factor (e.g. capital) are combined in production function $Q _ { i } \ =$ ${ \textstyle \frac { 1 } { 1 - \alpha } } \bar { A } _ { i } \bar { K } _ { i } ^ { \alpha } L _ { i } ^ { 1 - \alpha } .$ , where ${ \bar { A } } _ { i }$ is the innate productivity of a region $i , { \bar { K } } _ { i }$ is the capital endowment, and $\alpha$ is the capital share. Inverting the first order condition of the firm with respect to labor allows us to derive the following (inverse) labor demand curve: 

$$
L _ {i} = w _ {i} ^ {- 1 / \alpha} \times \bar {A} _ {i} ^ {1 / \alpha} \bar {K} _ {i}. \tag {3}
$$

Because $\alpha ~ > ~ 0 ,$ , equation (3) shows that labor demand is downward sloping with constant elasticity $- 1 / \alpha ,$ with the productivity and capital endowment of a location as a (log) shifter of the demand curve. This too is intuitive: the greater the productivity and/or capital endowment of a location, the greater the marginal product of labor in that location, increasing labor demand. 

The equilibrium distribution of economic activity—i.e. the equilibrium population and wages in each location—can be calculated by equating labor supply equation (2) with labor demand equation (3), where the equilibrium level of aggregate welfare W is determined by the aggregate labor market clearing constraint $\bar { L } = \textstyle \sum _ { j } L _ { j }$ . Doing so yields: 

$$
L _ {i} = \bar {u} _ {i} ^ {\frac {1}{\alpha - \beta}} \times (\bar {A} _ {i} \bar {K} _ {i} ^ {\alpha}) ^ {\frac {1}{\alpha - \beta}} \times W ^ {\frac {1}{\beta - \alpha}} \tag {4}
$$

$$
w _ {i} = \bar {u} _ {i} ^ {- \frac {\alpha}{\alpha - \beta}} \times \left(\bar {A} _ {i} \bar {K} _ {i} ^ {\alpha}\right) ^ {- \frac {\beta}{\alpha - \beta}} \times W ^ {\frac {\alpha}{\alpha - \beta}} \tag {5}
$$

$$
W = \left(\sum_ {i \in \mathcal {N}} \bar {u} _ {i} ^ {\frac {1}{\alpha - \beta}} \times (\bar {A} _ {i} \bar {K} _ {i} ^ {\alpha}) ^ {\frac {1}{\alpha - \beta}} / \bar {L}\right) ^ {\alpha - \beta}. \tag {6}
$$

Equations (4), (5), (6) completely specify the equilibrium spatial distribution of economic activity in terms of model fundamentals. 

# 3.1.2 Insights

Despite its simplicity, the Rosen (1979)-Roback (1982) model yields a number of insights into the spatial distribution of economic activity. Equation (4) shows that the equilibrium population is relatively larger in locations with better innate amenities, productivities, or capital endowments. Equation (5) shows that equilibrium wages are relatively greater in locations with better innate productivities or capital endowments, but relatively lower in places with higher amenities (since people are more willing to live in nice places even if the wages are low). And equation (6) demonstrates that the aggregate welfare is greater with greater amenities, productivities, and capital endowments world wide. All three equations also highlight the important role that the parameters $\alpha$ and $\beta$ play in determining how the underlying geography of the world shapes the equilibrium distribution of economic activity. 

# 3.1.3 Limitations

In principle, the Rosen (1979)-Roback (1982) location choice model is able to explain the high degree of inequality in the spatial distribution of economic activity highlighted in Fact 1 above by ascribing large differences in innate amenities, productivities, and endowments across locations. But it is also immediately obvious that the Rosen (1979)-Roback (1982) model is unable to capture Facts 2 and 3 for the simple reason that it abstracts entirely from the linkages between locations through the flow of goods. Given the important role that trade plays in the regional economy (Fact 2), this is a glaring omission. Moreover, without trade, there is no concept of distance or heterogeneity in who trades with whom (Fact 3). Indeed, as equations (4) and (5) highlight, apart from the endogenous scalar W, the equilibrium economic activity in one location is unaffected by the geography and endowments of every other location. 

# 3.2 The Krugman (1991) Core-Periphery model

We now consider the Core-Periphery model developed by Krugman (1991) (and further elaborated by Fujita, Krugman, and Venables (1999)). Unlike the Rosen (1979)-Roback (1982) residential choice model, the core-periphery model explicitly considers the economic linkages between locations. Yet it too has several limitations, which we will discuss below. 

# 3.2.1 Model sketch

As above, we consider a world with $N$ locations. Unlike above, however, we assume that each location is inhabited by an (endogenous) measure $N _ { i }$ of infinitesimal firms, each of which produces its own distinct variety. Because consumers will have a love of variety, this will create an incentive for trade across locations. 

We start with the consumer. An agent $l \in \mathcal { L } \equiv [ 0 , \bar { L } ]$ chooses (a) where to live and (b) how much to consume of each location-specific differentiated variety. As above, agents are assumed to be freely mobile, which again implies that welfare is equalized across all inhabited locations. Consumption, however, is very different than above. Agents are assumed to have constant elasticity of substitution (CES) demand and solve the following utility maximization problem: 

$$
\max  _ {j \in \mathcal {N}, \left\{c _ {i j} \right\} _ {i \in \mathcal {N}} \geq 0} \left(\sum_ {i \in \mathcal {N}} \int_ {\Omega_ {i}} c _ {i j} (\omega) ^ {\frac {\sigma - 1}{\sigma}} d \omega\right) ^ {\frac {\sigma}{\sigma - 1}} \text {s . t .} \sum_ {i \in \mathcal {N}} \int_ {\Omega_ {i}} p _ {i j} (\omega) c _ {i j} (\omega) d \omega \leq w _ {j}, \tag {7}
$$

where $\sigma > 1$ is the elasticity of substitution between varieties from different locations, $\Omega _ { i }$ is the set of firms in location i (with $| \Omega _ { i } | = N _ { i }$ ), $p _ { i j }$ is the price of a good produced in $i$ and sold in $j ,$ and $w _ { j }$ is the wage of an agent residing in location $j$ . 

This optimization results in the following consumer demand for an individual good $\omega$ in 

location j: 

$$
y _ {i j} (\omega) = p _ {i j} (\omega) ^ {1 - \sigma} P _ {j} ^ {\sigma - 1} w _ {j}. \tag {8}
$$

It also implies that the indirect utility function of the consumer is given by the wage, $w _ { j } ,$ divided by the price index in region $j , P _ { j }$ : 

$$
W _ {j} = w _ {j} / P _ {j}, \tag {9}
$$

where $\begin{array} { r } { P _ { j } ^ { 1 - \sigma } \equiv \sum _ { i \in \mathcal { N } } \int _ { \Omega _ { i } } p _ { i j } \left( \omega \right) ^ { 1 - \sigma } d \omega . } \end{array}$ . Following Redding and Venables (2004), we refer to $P _ { j } ^ { 1 - \sigma }$ as the (inward) market access of location $j$ . Intuitively, locations with better inward market access are those who face lower prices for the goods they purchase because they are closer (have lower trade costs) to origins with lower cost producers. 

We now consider the firm. Production uses one factor, labor. Firms are homogeneous but, unlike above, they produce differentiated products under monopolistic competition. Shipments to different regions are subject to trade costs that take the “iceberg” form such that $\tau _ { i j } \geq 1$ units of the goods need to be produced in i for one unit to arrive in region j. Given the demand of the consumer, the firm sets an optimal markup over the marginal cost and the price of a good in $i$ sold in region $j ,$ yielding the following equilibrium bilateral price: 

$$
p _ {i j} = \frac {\sigma}{\sigma - 1} \frac {w _ {i}}{A _ {i}} \tau_ {i j}, \tag {10}
$$

where, as above, $w _ { i }$ is the wage and $A _ { i }$ is the productivity. 

Substituting equation (10) into equation (8) and summing up across destinations we obtain total firm demand, $y _ { i } \left( \omega \right)$ , 

$$
y _ {i} (\omega) = \sum_ {j} \frac {\left(\frac {\sigma}{\sigma - 1} \frac {w _ {i}}{A _ {i}} \tau_ {i j}\right) ^ {1 - \sigma}}{P _ {j} ^ {1 - \sigma}} E _ {j},
$$

where $E _ { j }$ is the total expenditure in location $j$ . 

It is straightforward to show that the associated firm profits with CES are $y _ { i } \left( \omega \right) / \sigma$ . To generate their unique variety firms need to hire a fixed number of workers, $f _ { i } ^ { e }$ , to enter the market. The equilibrium number of firms is then pinned down by a zero-profit condition requiring that the expected profits of an entrant equal its fixed cost: 

$$
w _ {i} f _ {i} ^ {e} = \frac {1}{\sigma} \sum_ {j} \frac {\left(\frac {\sigma}{\sigma - 1} \frac {w _ {i}}{A _ {i}} \tau_ {i j}\right) ^ {1 - \sigma}}{P _ {j} ^ {1 - \sigma}} E _ {j}. \tag {11}
$$

In equilibrium, the total labor income in a location equals the sum of sales from firms of that location to all locations: 

$$
w _ {i} = \frac {N _ {i}}{L _ {i}} \sum_ {j} \frac {\left(\frac {\sigma}{\sigma - 1} \frac {w _ {i}}{A _ {i}} \tau_ {i j}\right) ^ {1 - \sigma}}{P _ {j} ^ {1 - \sigma}} E _ {j}. \tag {12}
$$

A comparison of equations (12) and (11) immediately shows that the equilibrium number of firms in a location is proportional to its population: 

$$
N _ {i} = \frac {L _ {i}}{\sigma f _ {i} ^ {e}}. \tag {13}
$$

We now can solve for the equilibrium. Substituting the price equation (10) into the indirect utility function of the consumer in equation (9) (using equation (13) to solve for the equilibrium number of firms) and applying welfare equalization yields the following (inverse) labor supply curve: 

$$
w _ {i} = W \left[ \sum_ {j} \frac {L _ {j}}{\sigma f _ {j} ^ {e}} \left(\frac {\sigma}{\sigma - 1} \frac {w _ {j}}{A _ {j}} \tau_ {j i}\right) ^ {1 - \sigma} \right] ^ {1 / (1 - \sigma)}. \tag {14}
$$

Similarly, re-writing equation (11) yields the following (inverse) labor demand curve: 

$$
w _ {i} = A _ {i} ^ {\frac {\sigma - 1}{\sigma}} \left(\frac {\sigma}{\sigma - 1}\right) ^ {\frac {1 - \sigma}{\sigma}} \left(\frac {1}{\sigma f _ {i} ^ {e}}\right) ^ {\frac {1}{\sigma}} \times \left(\sum_ {j} \frac {\left(\tau_ {i j}\right) ^ {1 - \sigma}}{P _ {j} ^ {1 - \sigma}} E _ {j}\right) ^ {\frac {1}{\sigma}}. \tag {15}
$$

A quick examination of equations (14) and (15) reveals a surprising discovery: both the labor supply and demand equations in the Core-Periphery model are perfectly inelastic to the labor in a location. This is in stark contrast to the Rosen-Roback framework above, where labor supply is upward sloping (as higher wages induce greater number of workers to reside in a location) and labor demand is downward sloping (as more labor entry leads to wages falling due to competitive pressures). Why are these forces not present here? On the supply side, there is no congestion force that directly reduces the utility of an agent residing in a populated location. On the demand side, competitive wage pressures are perfectly offset by the agglomeration force through firm entry. 

Instead of local labor, what matters here is the global market access. In the case of labor supply, workers are attracted to locations with better (inward) market access, as it increases their purchasing power; accordingly, such locations can offer lower wages and still attract workers. Similarly, labor demand depends on a weighted sum of the demand across all destinations $\begin{array} { r } { \Pi _ { i } ^ { 1 - \sigma } \equiv \sum _ { j } \frac { \left( \tau _ { i j } \right) ^ { 1 - \sigma } } { P _ { j } ^ { 1 - \sigma } } E _ { j } } \end{array}$ (τij )1−σ . Analogous to above, we refer to this term as (outward) market access. Locations with better outward market access are those for whom producers face higher demand for the goods they sell because they are closer (have lower trade costs) to destinations with higher demand. 

If you guessed that having a model with perfectly inelastic demand and supply curves may result in some issues, you are correct. It turns out that, absent of any additional intervention, the equilibrium is one in which there is complete concentration of all the workers in one location. Krugman (1991) solves this issue by adding an additional mechanism: the presence of an 

additional “agriculture” sector which is homogenous and costlessly traded across sectors. He assumes that consumers allocate a constant fraction of their income to this sector and the sector employs an exogenously given equal number of workers in each location. This additional mechanism is neither elegant nor necessary and so it has largely been left on the chopping block in the next generation of economic geography models. In the following subsection, however, we retain his original assumption in the following subsection to highlight several insights that arise from the model. 

# 3.2.2 Insights

The Core-Periphery model was a significant addition to the economic geography literature. It offered two major innovations: first, it introduced the concept of market access arising from the linkages between locations through the trade of goods, which Fact #2 and #3 above show are a crucial part of the spatial economy. Second, the presence of firm entry introduced an endogenous agglomeration force. Both the market access and agglomeration forces remain a cornerstone of modern quantitative economic geography, so it is worthwhile discussing briefly how they interact to determine the equilibrium distribution of economic activity in the Core-Periphery model. 

Innate differences and economic agglomeration First, economic agglomerations can exacerbate innate difference across locations. Figure 6 plots the population allocations for two locations under different scenarios. Here we fix the productivity of location 2 to 1 and slightly vary the productivity of location 1 from 1, to 1.1, and 1.2. Of course, in the scenario where locations are symmetric, population is equally distributed between the two locations. A small change of $1 0 \%$ in the productivity of location 1 implies a disproportionate change in the share of population in location 1, more than $2 0 \%$ . In fact, if we increase the productivity to 1.2 the solution indicates that the welfare difference between the location is always positive. In other words, the agglomeration force in this scenario is so strong that with just a small change in parameters the model is lead to a corner solution where all the population moves to location 1. 

This insight suggests that endogenous agglomeration forces may be part of the reason for such substantial variation in the distribution of economic activity across space presented in Fact #1 above. 

Market access and economic outcomes Second, a key innovation of the Core-Periphery model is the presence of the market access term that appears because of trade costs. As trade costs increase, wage difference across locations are more pronounced. In Figure 7, left panel, we plot the share of population in location 1 for the two locations with different productivities and vary trade costs across a range of values. Naturally, if trade costs are too high the two locations act 


Figure 6: The distribution of spatial economic activity in the Core-Periphery Setup


![image](https://cdn-mineru.openxlab.org.cn/result/2026-02-25/4e4c25dd-e494-4502-8cd1-6a5c0f2e691d/0b08c5cf6fb14ba84de0982b688c9d417228546eb3c0f76581cb342879f245d1.jpg)



Source: Author’s calculations based on the core periphery setup. We chose calibrated parameters to closely resemble the simulations in Fujita, Krugman, and Venables (1999), bearing the assymetry in productivities.


completely independently, but in the limit there are a few more workers in location 1, where the productivity is $1 0 \%$ higher. However, as the trade costs reduce the difference in productivity lead to more and more significant differences in the share of population in location 1. As trade costs decline the producer market access increases, which means higher labor demand and higher wages. This attracts population through entry and that creates an even stronger agglomeration force. The two forces coalesce at low trade costs leading to extreme levels of economic agglomeration and at trade costs lower than 2 the population share of location 1 is upward of 0.75. 

This insight suggests that the strength of the economic linkages plays an important role in determining the equilibrium spatial distribution of economic activity. 

Multiplicity of equilibria The final, and most intriguing feature of the Core-Periphery model is the possibility that arises for multiple equilibria. In Figure 7, right panel, we plot the equilibrium of the model with two symmetric locations, with different trade costs. This is the famous “Tomahawk graph” (see e.g. Chapter 5 in Fujita, Krugman, and Venables (1999)) that depicts the level of trade costs on the x-axis from starting from free trade and the share of population in location 1 on the y-axis. For large trade costs the equilibrium is unique and is evenly divided across the two locations. But below a threshold of trade costs the number of equilibria jumps to three, with the symmetric equilibria being still one possibility. The other two correspond to cases where population is concentrated in location 1 or location 2, due to economic 


Figure 7: The distribution of spatial economic activity in the Core-Periphery Setup


![image](https://cdn-mineru.openxlab.org.cn/result/2026-02-25/4e4c25dd-e494-4502-8cd1-6a5c0f2e691d/66e2b65c2e396d07dc73423d9ba22d5bd675b1ea127fdcc2f5953663e5a1ac6b.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-02-25/4e4c25dd-e494-4502-8cd1-6a5c0f2e691d/83cf0a12f33e884048e3159ed219234d160967a4917411d2d0d6e5b6946340fc.jpg)



Source: Author’s calculations based on the core periphery setup. We chose calibrated parameters to closely resemble the simulations in Fujita, Krugman, and Venables (1999), barring the asymmetry in productivities.


agglomeration. The intuition for the emergence of the alternative equilibria with concentration in population is that as trade costs lower, economic agglomeration increases the effect of market access to the manufacturing sector. After a critical point this effect is so strong that multiple equilibria arise. After that threshold the lower the trade costs lead to even more concentration, and eventually to only one location to be populated with manufacturing workers in equilibrium. 

The Core-Periphery model is silent on whether multiplicity is merely an academic novelty or actually reflects a fundamental property of the real world spatial equilibria. If multiplicity is indeed an empirical phenomenon, then conducting counterfactual analyses must account for the possibility that small changes in the geography may result in large changes in the spatial equilibria. This highlights the need for a clear understanding of when the uniqueness of the spatial equilibrium can be assured and when multiplicity can arise. 

# 3.2.3 Limitations

While the core-periphery model was pioneering in its introduction of spatial economic linkages and agglomeration forces, it has several notable limitations. First, because the equilibrium supply and demand only depends on market access and not local labor, even when considering simple geographies, an equilibrium may not exist or may be unrealistic (e.g. featuring complete concentration in a single location even if location fundamentals are not very different). Second, by relying entirely on firm entry to generate agglomeration forces, the model features a strong and inflexible agglomeration force that can generate economic responses to changes in the underlying geography that appear unrealistic (e.g. featuring large changes in the share 

of population to small changes in trade costs). Both these limitations are exacerbated when confronting real world data featuring many locations and complex geographies. 

# 3.3 The current generation

While both the location choice model of Rosen (1979) and Roback (1982) and the Core-Periphery model of Krugman (1991) offer important new insights, neither alone offer an appropriate starting point to build a quantitative regional economic framework. But, one might wonder what could happen if we combine the richness of the economic linkages in the Core-Periphery model with the tractability of the location choice model of Rosen and Roback. Is that possible? And what could happen if we replace the extreme agglomeration forces induced by entry with alternative ways of modeling agglomeration, borrowing ideas from the literature in labor economics that uses the Rosen-Roback model (see for example the review of Moretti (2011)) or from international trade work that models the presence of Marshallian externalities (e.g. Grossman and Rossi-Hansberg (2010) and Kucheryavyy, Lyn, and Rodr´ıguez-Clare (2023)). If we could add these features in a seamless way, we could revisit some of the facts in our introduction and, for example, ask how much of the increasingly unequal distribution of economic activity is due to reductions in trade costs, increases in trade openness, etc. This is what we turn to next. 

# 4 A simple quantitative regional model

We now turn to presenting a simple quantitative economic geography framework. This framework, which builds off both the location choice and core-periphery models of the previous section, is based on Allen and Arkolakis (2014). As we will then discuss in detail in Section 5, it also turns out to be mathematically equivalent to many alternative frameworks. 

# 4.1 Setup

We first describe the setup of the framework. We consider a world comprising a compact set $\mathcal { N }$ of locations, which we index with $i , j \in \mathcal { N }$ . We assume the world has perfectly competitive markets for both products and labor. 

# 4.1.1 Geography

Each location $i \in \mathcal N$ endowed with (a) a technology for producing a unique variety (which we will also index by i); (b) an innate productivity ${ \bar { A } } _ { i }$ for producing that variety; and (c) an innate amenity ${ \bar { u } } _ { i }$ that affects the welfare of agents residing in i. The assumption that each location produces a distinct differentiated variety—proposed by Armington (1969) and applied to trade flows by Anderson (1979)—is clearly an extreme simplification. Nevertheless, the Armington 

assumption turns out to be both surprisingly tractable and mathematically equivalent to more complicated (but more realistic) assumptions, as we will show below. In addition, each pair of locations $i , j \in \mathcal { N }$ are endowed with a technology for trading goods, which following Samuelson (1954) we assume take the “iceberg” form, where $\tau _ { i j } \geq 1$ units of the product must be shipped from $i \in \mathcal N$ in order for one unit to arrive in $j \in \mathcal N$ . Together, the $\left\{ \bar { A } _ { i } \right\} _ { i \in \mathcal { N } } , \left\{ \bar { u } _ { i } \right\} _ { i \in \mathcal { N } } ,$ and $\left\{ \tau _ { i j } \right\} _ { i , j \in \mathcal { N } }$ comprise the geography of the world. 

Our interest is to determine how the geography of the world shapes the equilibrium distribution of economic activity. To do so, we begin by populating the world with a measure $\bar { L }$ of infinitesimal agents. These agents are assumed to be perfectly mobile across locations. Wherever they choose to reside, they produce goods to exchange with other agents and consume goods that they have purchased from other agents. The particular value $\bar { L }$ turns out to not be particularly important in this simple framework: it will not affect the fraction of agents who live in a location or the share of economic activity occurring in each location. 

# 4.1.2 Production

Labor is the only factor of production. Given the assumption of perfect competition, the price of a good produced in $i \in \mathcal N$ and consumed in $j \in \mathcal N$ is equal to its marginal cost, namely: 

$$
p _ {i j} = \tau_ {i j} w _ {i} / A _ {i}, \tag {16}
$$

where $w _ { i }$ is the wage per unit of labor and $A _ { i }$ is the total productivity of workers in location $i \in \mathcal N$ . Note that the (innate) productivity ${ \bar { A } } _ { i }$ and the (total) productivity $A _ { i }$ are distinct objects; we will discuss their relationship below. 

# 4.1.3 Consumption

Consider an agent $l \in \mathcal { L } \equiv [ 0 , \bar { L } ]$ that supplies her unit of labor inelastically and chooses (a) where to live and (b) how much to consume of each location-specific differentiated variety in order to maximize her utility subject to her budget constraint by solving: 

$$
\max  _ {j \in \mathcal {N}, \left\{c _ {i j} \right\} _ {i \in \mathcal {N}} \geq 0} \left(\sum_ {i \in \mathcal {N}} c _ {i j} ^ {\frac {\sigma - 1}{\sigma}}\right) ^ {\frac {\sigma}{\sigma - 1}} \times u _ {j} \text {s . t .} \sum_ {i \in \mathcal {N}} p _ {i j} c _ {i j} \leq w _ {j}, \tag {17}
$$

where $\sigma \geq 0$ is the elasticity of substitution between varieties from different locations, $u _ { j }$ is the total amenity of residing in $j , p _ { i j }$ is the price of a good produced in i and sold in $j ,$ and $w _ { j }$ is the wage of an agent residing in location j. This problem is a variation of maximization in equation (7) with an additional location specific amenity term. We introduce this term to add flexibility in the model to fit the distribution of observed population in the data. As with the productivity, 

note that the (innate) amenity ${ \bar { u } _ { j } }$ and the (total) amenity $u _ { j }$ are distinct objects; we will discuss their relationship below. 

Optimization of equation (17) yields two important results. First, it allows us to write the welfare of any individual residing in location $i \in \mathcal N$ as an indirect function of her wage, the amenity, and the price index as follows: 

$$
W _ {i} = w _ {i} u _ {i} / P _ {i}. \tag {18}
$$

Second, the optimization yields the following expression for the quantity consumed by agent $l \in \mathcal L$ residing in $j \in \mathcal N$ : 

$$
c _ {i j} = p _ {i j} ^ {- \sigma} P _ {j} ^ {\sigma - 1} w _ {j}, \tag {19}
$$

where $\begin{array} { r } { P _ { j } ^ { 1 - \sigma } = \sum _ { i \in \mathcal { N } } p _ { i j } ^ { 1 - \sigma } } \end{array}$ is the Dixit-Stiglitz price index. Multiplying the quantity consumed given in equation (19) by the price and summing across the measure $L _ { j }$ agents residing in $j \in \mathcal N$ yields the total value of goods produced in $i \in \mathcal N$ and sold in $j \in \mathcal { N } , X _ { i j }$ : 

$$
X _ {i j} = p _ {i j} ^ {1 - \sigma} P _ {j} ^ {\sigma - 1} w _ {j} L _ {j}. \tag {20}
$$

Substituting equation (16) for the prices into the bilateral demand equation (20) yields the following “gravity” expression: 

$$
X _ {i j} = \tau_ {i j} ^ {1 - \sigma} \left(w _ {i} / A _ {i}\right) ^ {1 - \sigma} P _ {j} ^ {\sigma - 1} w _ {j} L _ {j}. \tag {21}
$$

Equation (21) states that, if goods are substitutes (i.e. $\sigma > 1$ ), then the value of bilateral trade flows decreases in the trade costs between the locations $( \tau _ { i j } )$ and the cost of production in the origin $( w _ { i } / A _ { i } ) .$ , and increasing in the total income of the destination $( w _ { j } L _ { j } )$ . It increases in the trade costs and costs of production in all other origins, summarized by the Dixit-Stiglitz price index $( P _ { j } )$ . In doing so, it succinctly summarizes the rich patterns of trade flows observed empirically described in Section 2.3 above. 

# 4.1.4 Agglomeration Forces

We now formally incorporate the agglomeration forces that we discussed in Section 3.2 in a flexible way. To do so we link the innate productivities and amenities to the total productivities and amenities. Suppose that the total productivity of a location is a function of its (exogenous) innate productivity and its (endogenous) labor population with constant elasticity $\alpha \in \mathbb { R }$ : 

$$
A _ {i} = \bar {A} _ {i} L _ {i} ^ {\alpha}. \tag {22}
$$

A positive $\alpha$ means that the productivity of a location is increasing with its population density. This could occur for a variety of reasons, e.g. firm entry, knowledge spillovers, better access to labor or inputs, etc. Conversely, a negative $\alpha$ means that the productivity of a location is decreasing with its population density, for example due to a fixed factor of production. Of course, none of these mechanisms are explicit in this simple framework; instead we interpret $\alpha$ as a “reduced form” way of capturing the effects of such forces have on the prices consumers pay for tradable goods. Below, however, we will show formally that there are multiple microfoundations incorporating such forces that justify the particular functional form of equation (22) in Section 5 below. 

Similarly, suppose that the total amenity of a location is a function of its innate amenity and its labor population with constant elasticity $\beta \in \mathbb { R }$ : 

$$
u _ {i} = \bar {u} _ {i} L _ {i} ^ {\beta}. \tag {23}
$$

As with the productivity spillovers, $\beta$ may be positive or negative. A positive value—which indicates that amenities are increasing with population—could arise if e.g. a greater population density enables a location to sustain a greater variety of restaurants, better schools, or better parks. Conversely, a negative value could arise if a greater population density e.g. increases the cost of fixed factors of consumption like housing. As with the productivity spillovers, we interpret $\beta$ as parameter meant to capture (in a reduced form way) the many possible forces through which the local population affects the amenity of living in a location other than the direct cost of tradable goods. We provide several examples of particular micro-foundations in Section 5 below. 

As we mentioned in Section 3.2, an earlier generation of economic geography theory highlighted the important role that agglomeration forces played in determining how the innate geography affected the equilibrium distribution of economic activity. As will become evident below, this is still true in quantitative economic geography frameworks: the particular value and magnitude of $\alpha$ and $\beta$ play an important role in determining the relationship between the innate geography and the equilibrium distribution of economic activity. 

# 4.2 Equilibrium

Our goal is to determine how geography affects the equilibrium distribution of economic activity across locations. To do so, we first derive the equilibrium conditions that define this mapping. 

# 4.2.1 Definition

We can now discuss the equilibrium of this simple quantitative framework. Formally, for any geography (which recall is the set of innate productivities $\left\{ { \bar { A } } _ { i } \right\} _ { i \in { \mathcal { N } } } ,$ amenities $\{ \bar { u } _ { i } \} _ { i \in \mathcal { N } } ,$ and trade costs $\left\{ \tau _ { i j } \right\} _ { i \in \mathcal { N } } ) .$ , set of model parameters $\sigma \geq 0$ , $\alpha \ \in \ \mathbb { R }$ , $\beta \in \mathbb { R }$ and $\bar { L } > 0 \AA$ , a spatial equilibrium is the distribution of labor $\{ L _ { i } \} _ { i \in { \mathcal { N } } } ,$ wages $\{ w _ { i } \} _ { i \in \mathcal { N } } ,$ and price indices $\{ P _ { i } \} _ { i \in \mathcal { N } }$ such that the following four conditions hold: 

1. Goods markets clear. In particular, the total labor income in each location is equal to its total sales to all locations: 

$$
w _ {i} L _ {i} = \sum_ {j \in \mathcal {N}} X _ {i j} \forall i \in \mathcal {N}. \tag {24}
$$

2. Budget constraints are satisfied. In particular, the total expenditure in each location is equal to its total income: 

$$
w _ {i} L _ {i} = \sum_ {j \in \mathcal {N}} X _ {j i} \forall i \in \mathcal {N}. \tag {25}
$$

3. Welfare is equalized. In particular, there exists a scalar $W > 0$ such that: 

$$
W _ {i} \leq W \forall i \in \mathcal {N}, \tag {26}
$$

with equality if $L _ { i } > 0$ . 

4. The total population of the world is equal to its endowed amount of labor: 

$$
\sum_ {i \in \mathcal {N}} L _ {i} = \bar {L}. \tag {27}
$$

Note that equilibrium conditions (24) and (25) together imply that trade is balanced in every location. While this is an abstraction from reality (where trade deficits and surpluses are commonplace), the static nature of the framework makes it difficult to incorporate such imbalances, although some prior work has simply imposed exogenous transfers, see e.g. Dekle, Eaton, and Kortum (2008). In Section 9, we briefly discuss recent dynamic extensions to this framework. We also refer the interested reader to Desmet and Parro (2025) for an in depth discussion. 

In what follows, we focus on an interior equilibrium where every location is inhabited, so that the welfare equalization condition (26) implies that $W _ { i } = W$ for all $i \in \mathcal N$ . We do so for two reasons. First, an interior solution is the empirically relevant one for the study of the regional economy: in an analysis of cities, states, counties, etc., all locations are inhabited. Second, the log linear functional forms in equations (22) and (23) for the agglomeration forces are not particularly realistic when locations are uninhabited, implying for a positive (negative) $\alpha$ or $\beta$ that an uninhabited location would have zero (infinite) productivity or amenity. 

# 4.2.2 Deriving the equilibrium conditions

Substituting the gravity equation (21) and agglomeration equations (22) and (23) into equilibrium conditions (24) and (25) and collecting like terms implies for all $i \in \mathcal N$ : 

$$
w _ {i} ^ {\sigma} L _ {i} ^ {1 - \alpha (\sigma - 1)} = \sum_ {j \in \mathcal {N}} \left(\bar {A} _ {i} / \tau_ {i j}\right) ^ {\sigma - 1} P _ {j} ^ {\sigma - 1} w _ {j} L _ {j} \tag {28}
$$

$$
P _ {i} ^ {1 - \sigma} = \sum_ {j \in \mathcal {N}} \left(\bar {A} _ {j} / \tau_ {j i}\right) ^ {\sigma - 1} w _ {j} ^ {1 - \sigma} L _ {j} ^ {\alpha (\sigma - 1)} \tag {29}
$$

Similarly, substituting the indirect utility equation (18) and agglomeration equations (22) and (23) into equilibrium condition (26) and solving for the price index yields for all $i \in \mathcal N$ : 

$$
P _ {i} = w _ {i} \bar {u} _ {i} L _ {i} ^ {\beta} / W. \tag {30}
$$

Notice that this expression can be reformulated to be written as a labor supply curve, a point we will discuss in detail below. 

Finally, substituting the welfare equalization equation (30) into the market clearing equation (28) and budget constraint equation (29) yields the following system of equations for all $i \in \mathcal N$ : 

$$
W ^ {\sigma - 1} w _ {i} ^ {\sigma} L _ {i} ^ {1 - \alpha (\sigma - 1)} = \sum_ {j \in \mathcal {N}} \left(\bar {A} _ {i} \bar {u} _ {j} / \tau_ {i j}\right) ^ {\sigma - 1} w _ {j} ^ {\sigma} L _ {j} ^ {1 + \beta (\sigma - 1)} \tag {31}
$$

$$
W ^ {\sigma - 1} w _ {i} ^ {1 - \sigma} L _ {i} ^ {\beta (1 - \sigma)} = \sum_ {j \in \mathcal {N}} \left(\bar {u} _ {i} \bar {A} _ {j} / \tau_ {j i}\right) ^ {\sigma - 1} w _ {j} ^ {1 - \sigma} L _ {j} ^ {\alpha (\sigma - 1)} \tag {32}
$$

Equations (31) and (32) define the relationship between the world geography and model parameters on the one hand and the equilibrium distribution of economic activity on the other. They provide $2 N$ equations (two for each location), which along with the aggregate labor market clearing condition (27), is the same number as the $2 N + 1$ endogenous outcomes (wages $w _ { i }$ and labor $L _ { i }$ in each location, along with the overall welfare $W$ of the world). This provides hope that this system of equations is sufficient for determining the equilibrium distribution of economic activity. This hope, as we will show below, turns out to be well founded. 

The particular mathematical structure of equations (31) and (32) deserves some discussion. Both equations state that a certain log-linear combination of endogenous outcomes in a particular location $i ,$ scaled appropriately, is equal to a weighted average of a different log-linear combination of endogenous outcomes in all other locations. The log-linear combinations depend on the model parameters $\sigma , \alpha ,$ and $\beta ,$ , whereas the weights depend on the geography $\left\{ \bar { A } _ { i } \right\} _ { i \in \mathcal { N } } , \{ \bar { u } _ { i } \} _ { i \in \mathcal { N } }$ , and $\left\{ \tau _ { i j } \right\} _ { i , j \in \mathcal { N } }$ of the system. That the geographic weights in equation (32) are the transpose of (31) will be important below. It turns out that this particular structure is ubiquitous among quantitative economic geography models, and so in Section 6 we discuss the properties of such a system in detail. Prior to doing so, however, we first consider three special 

cases of the framework. 

# 4.3 Market access, supply, and demand

We proceed by offering a slightly more intuitive expression for bilateral trade flows and, in doing so, revisit the important concept of market access. The following derivations come from Anderson and Van Wincoop (2003) (who referred to market access as “multilateral resistance.”) Substituting gravity equation (21) into the goods market clearing condition (24) and re-arranging yields: 

$$
\left(w _ {i} / A _ {i}\right) ^ {1 - \sigma} = \frac {w _ {i} L _ {i}}{\sum_ {j \in \mathcal {N}} \tau_ {i j} ^ {1 - \sigma} P _ {j} ^ {\sigma - 1} w _ {j} L _ {j}}. \tag {33}
$$

We then proceed by substituting equation (33) back into the gravity equation, yielding: 

$$
X _ {i j} = \tau_ {i j} ^ {1 - \sigma} \times \frac {w _ {i} L _ {i}}{\Pi_ {i} ^ {1 - \sigma}} \times \frac {w _ {j} L _ {j}}{P _ {j} ^ {1 - \sigma}}, \tag {34}
$$

where $\Pi _ { i } ^ { 1 - \sigma } \equiv \sum _ { j \in \mathcal { N } } \tau _ { i j } ^ { 1 - \sigma } P _ { j } ^ { \sigma - 1 } w _ { j } L _ { j }$ is the (outward) market access and $P _ { j } ^ { 1 - \sigma }$ was defined in equation (33) above is the (inward) market access. 

The market access terms play a key role in quantitative economic geography frameworks, as they mediate how economic activity elsewhere influences economic activity in a given location. To see this, we can substitute the productivity spillover equation (22) into equation (33), take logs, and solve for wages to construct an (inverse) labor demand curve: 

$$
\ln w _ {i} = - \left(\frac {1}{\sigma} - \alpha \left(\frac {\sigma - 1}{\sigma}\right)\right) \ln L _ {i} + \frac {1}{\sigma} \ln \Pi_ {i} ^ {1 - \sigma} + \frac {\sigma - 1}{\sigma} \ln \bar {A} _ {i}. \tag {35}
$$

Equation (35) highlights that, conditional on outward market access, the demand curve for labor in location $i \in \mathcal N$ is a straightforward log-linear function. As long as $\alpha$ is not too large (i.e. $\alpha < 1 / \left( \sigma - 1 \right) )$ , then the labor demand curve is downward sloping, and a larger agglomeration force attenuates the downward slope. Greater innate productivity or market access increase labor demand by shifting the demand curve upward. 

Similarly, we can take logs of equation (30) and solve for wages to construct an (inverse) labor supply curve: 

$$
\ln w _ {i} = - \beta \ln L _ {i} + \frac {1}{1 - \sigma} \ln P _ {i} ^ {1 - \sigma} + \ln W - \ln \bar {u} _ {i}. \tag {36}
$$

Like with the labor demand curve, equation (36) shows that, conditional on inward market access, the labor supply curve is also simple log-linear function. As long as $\beta < 0$ (i.e. there are congestion externalities in consumption), then the labor supply curve is upward sloping, with the strength of the externalities determining the slope. Greater innate amenities or better 

inward market access shifts the labor supply curve downward, since either implies that workers can maintain the same level of welfare with lower nominal wages. 

Note that the only place that economic activity in the rest of the world shows up in the labor demand curve (35) and labor supply curve (36) are in the outward and inward market access terms, respectively. Put another way, the market access terms together fully summarize how locations interact with each other. That means that if we could treat the market access terms as exogenous characteristics of a location, then solving for the equilibrium economic activity would be simple: we would simply equate supply and demand in each location separately and solve for the equilibrium wages and labor. But, of course, it is not so simple, as both the inward and outward market access terms are endogenous variables that depend on both the exogenous geography and the endogenous distribution of economic activity in all locations. 

Before we start exploring the role of the endogenous market access terms, it is worth pointing out that we already have a small reason to celebrate. Recall from Section 3.2 that the equilibrium of the core-periphery framework was quite fragile due to the labor supply and demand curves not depending directly on local labor. By modeling the amenity and productivity spillovers, we have solved that problem, with both the labor demand curve (35) and labor supply curve (36) now depending on the local labor while still retaining the important role that market access plays in the core-periphery model. This will prove crucial for the model to remain tractable even in a world with a rich geography, making it an ideal framework for quantitative applications. Indeed, it will provide the basis for the workhorse quantitative economic geography framework we present in Section 6. 

# 4.4 Special case #1: No agglomeration or congestion externalities

Suppose that there are not agglomeration or congestion externalities in either production or consumption, i.e. $\alpha = \beta = 0$ . In this case, the equilibrium conditions (31) and (32) simplify to the following: 

$$
W ^ {\sigma - 1} w _ {i} ^ {\sigma} L _ {i} = \sum_ {j \in \mathcal {N}} \left(\bar {A} _ {i} \bar {u} _ {j} / \tau_ {i j}\right) ^ {\sigma - 1} w _ {j} ^ {\sigma} L _ {j} \tag {37}
$$

$$
W ^ {\sigma - 1} w _ {i} ^ {1 - \sigma} = \sum_ {j \in \mathcal {N}} \left(\bar {u} _ {i} \bar {A} _ {j} / \tau_ {j i}\right) ^ {\sigma - 1} w _ {j} ^ {1 - \sigma}. \tag {38}
$$

In the absence of agglomeration or congestion externalities, the log linear combination of endogenous outcomes is the same on both the left hand side and right hand side of each equation. This allows us to write equations (37) and (38) more succinctly in matrix notation as follows: 

$$
\lambda \mathbf {x} = \mathbf {T x} \tag {39}
$$

$$
\lambda \mathbf {y} = \mathbf {T} ^ {\prime} \mathbf {y}, \tag {40}
$$

where $\lambda \equiv W ^ { \sigma - 1 }$ is an endogenous scalar, $\mathbf { x }$ is an $N \times 1$ vector with $i ^ { t h }$ element $w _ { i } ^ { \sigma } L _ { i } ,$ y is an $N \times 1$ vector with $i ^ { t h }$ element $w _ { i } ^ { 1 - \sigma } ,$ , and $\mathbf { T }$ is an $N \times N$ matrix with $( i , j ) ^ { t h }$ element $\big ( \bar { A _ { i } } \bar { u } _ { j } / \tau _ { i j } \big ) ^ { \sigma - 1 }$ . As equations (39) and (40) make clear, the endogenous distribution of economic activity in the absence of agglomeration or congestion externalities is simply determined by the left and right eigenvectors of the matrix T, which summarizes the geography of the world. 

Given the familiar mathematical structure of the equilibrium, a number of interesting properties of the equilibrium follow readily from well known mathematical results. 

First, as long as T is positive in all its elements and because x and y must be non-negative (since wages and populations must be weakly positive), by the Perron Frobenius theorem (Perron, 1907; Frobenius, 1912), there exists unique (to scale) and strictly positive eigenvectors corresponding to the largest (in absolute value) eigenvalue of equations (39) and (40). Given these eigenvectors, the equilibrium wages and populations can be immediately recovered by solving the two system $x _ { i } \propto w _ { i } ^ { \sigma } L _ { i }$ and $y _ { i } \propto w _ { i } ^ { 1 - \sigma }$ for $w _ { i }$ and $L _ { i } ,$ where the proportions are pinned down by a choice of numeraire (for wages) and the aggregate labor market clearing condition equation (27) (for labor). 

Second, because it is straightforward to calculate the largest eigenvalue and associated eigenvectors of a matrix (using e.g. a power iteration algorithm, see e.g. Mises and Pollaczek-Geiringer (1929)), determining the unique distribution of economic activity is very easy computationally and can be accomplished in very little time even with a large number of locations and a rich geography. 

Third, because economic activity is greater in locations with larger eigenvectors, if we view the geography matrix T as a weighted graph of a network between locations, then the concept of eigenvector centrality (see e.g. Jackson et al. (2008)) tells us that economic activity will be more concentrated in locations that are more central to this network. The weights of the graph, in this case, are the elements of the matrix, $T _ { i j }$ . Since $T _ { i j } = \big ( \bar { A } _ { i } \bar { u } _ { j } / \tau _ { i j } \big ) ^ { \sigma - 1 } , \sigma > 1 ,$ more central locations will (loosely speaking) be those with better productivities and amenities and lower trade costs to other locations (and those connected to other locations with similarly advantaged geographies). 

Fourth, because world welfare is $( \log )$ proportional to the largest eigenvalue of the geography matrix T, then world welfare will tend to be larger the lower the trade costs and the greater the innate productivities and amenities of locations. We can formalize this result by considering what happens to the largest eigenvalue when we perturb the $( i , j ) ^ { t h }$ element of the geography matrix T by a small amount. From Vahrenkamp (1976) we have that the resulting change in the largest eigenvalue is proportional to the product of the left and right eigenvectors as follows: 

$$
\frac {\partial \lambda}{\partial T _ {i j}} = y _ {i} x _ {j} / \sum_ {k \in \mathcal {N}} y _ {k} x _ {k}. \tag {41}
$$

Given the relationships of $\lambda , T _ { i j } , y _ { j }$ and $x _ { i }$ to the underlying economic variables, it is straightforward to show that equation (41) implies that the elasticity of aggregate welfare to a reduction in bilateral trade costs between locations $i$ and $j$ is equal to the value of trade flows between those locations, i.e.: 

$$
- \frac {\partial \ln W}{\partial \ln \tau_ {i j}} = X _ {i j} / Y ^ {W}, \tag {42}
$$

where $Y ^ { W }$ is the total world income. This result is a key building block to the “social savings” approach proposed by Fogel (1964) for advocating the welfare effects of changes in transportation costs. Intuitively, in the absence of agglomeration or congestion forces, the competitive equilibrium is efficient, and so the first order effects of technology improvements like the reduction of trade costs on aggregate welfare are proportional to the share of aggregate economic activity on those links, as in Hulten (1978). 

To summarize, in the absence of agglomeration and congestion forces, the effect of geography on the spatial distribution of economic activity is well behaved and well understood. But as emphasized in Section 3, there is much interesting economics that comes from incorporating agglomeration and congestion forces. We now make a first pass at trying to incorporate these forces. 

# 4.5 Special case #2: Two locations

We now consider the case of two locations, symmetric trade costs, and no innate amenity differences. This is the case that can be most closely connected to the Core-Periphery model that we presented in Section 3. After some tedious derivations, equations (31) and (32) can be written as the following second order equation ratio of wages, ${ w _ { 1 } / w _ { 2 } }$ : 

$$
\frac {A _ {2} ^ {\sigma - 1}}{A _ {1} ^ {\sigma - 1}} \left(\frac {w _ {1} ^ {(\sigma - 1)}}{w _ {2} ^ {(\sigma - 1)}}\right) ^ {2} + \left(1 - \frac {A _ {2} ^ {\sigma - 1}}{A _ {1} ^ {\sigma - 1}}\right) \tau^ {1 - \sigma} \frac {w _ {1} ^ {\sigma - 1}}{w _ {2} ^ {\sigma - 1}} - 1 = 0.
$$

The solution of this equation is the positive solution of the quadratic root of this equation, which is, 

$$
\frac {w _ {1}}{w _ {2}} = \frac {1}{2} \left(1 - \frac {A _ {1} ^ {\sigma - 1}}{A _ {2} ^ {\sigma - 1}}\right) \tau^ {1 - \sigma} + \frac {1}{2} \frac {A _ {1} ^ {\sigma - 1}}{A _ {2} ^ {\sigma - 1}} \sqrt {\left(\frac {A _ {2} ^ {\sigma - 1}}{A _ {1} ^ {\sigma - 1}} - 1\right) ^ {2} (\tau^ {1 - \sigma}) ^ {2} + 4 \frac {A _ {2} ^ {\sigma - 1}}{A _ {1} ^ {\sigma - 1}}}
$$

We assume, without loss of generality that $A _ { 1 } > A _ { 2 }$ . The above expression can be shown to imply $w _ { 2 } / w _ { 1 } < 1 ,$ , i.e. the more productive region offers higher wages. Notice that this ratio is decreasing in $\tau$ so long as $\sigma > 1$ , i.e. the wage advantage of location 1 decreases when trade costs increase. This is intuitive as market access decreases and there is a tendency for workers to endogenously move to the location with the better access (see Allen and Arkolakis, 2018). 

To solve for the population share we need to resort to simulations. Figure 8 in the left panel 


Figure 8: The distribution of spatial economic activity in the Allen and Arkolakis (2014) model


![image](https://cdn-mineru.openxlab.org.cn/result/2026-02-25/4e4c25dd-e494-4502-8cd1-6a5c0f2e691d/74d67f7594e2d816bf61f2e9b19323d61de1641f0d259dc920fbabdce7fb7535.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-02-25/4e4c25dd-e494-4502-8cd1-6a5c0f2e691d/93db9c04c010f2e925d40e2baae4aeb801fabfc6151b0fd0f10631af93edca99.jpg)



Notes: This figure shows the fraction of labor allocated to the more productive region (y-axis) for different values of trade costs (x-axis). The left panel considers an economy without any productivity or amenity externalities; the right panel considers an economy with agglomeration forces arising through productivity externalities.


plots the share of population on location 1 when location 1 is the most productive one, when $\alpha = \beta = 0$ and the productivity in location 1 is $1 0 \%$ higher than the location 2. Naturally, the share of population choosing location 1 is more than half. In addition, the share of population in location 1 increases with trade costs, i.e. trade costs increase concentration of economic activity. Notice that this is in sharp contrast with predictions of the Core-Periphery framework displayed in Figure 7 where the share of the dominant region decreases with higher trade costs. 

It is useful to revisit the prediction of the Core-Periphery model for the share of population as a function of trade costs, the well known Tomahawk graph depicted above in the right panel of Figure 7. To analyze that we plot our model with positive agglomeration externality, $\alpha =$ 0.1, and $\beta = 0$ in the right panel of Figure 8. In the quantitative geography model there are multiple equilibria with large trade costs since this is when more concentration takes place. The stable equilibria are ones of concentration, consistent with the previous discussion, while the symmetric equilibrium is now an unstable one for high trade costs. The new figure resembles a “Pitchfork” rather than a “Tomahawk”. A different prediction, but still quite sharp! 

It is important to point out another important antecedent of our work. In particular, our model under a certain parameterization becomes identical to the Helpman (1998) model, as we will explain below. That model also features entry, as in Krugman (1991), which acts as the positive agglomeration spillover, but there is only a single tradable sector, the manufacturing sector. Consumers consume a second good that is nontradeable, housing, which in the case of the Helpman (1998) framework acts as a negative amenity spillover. The Core-Periphery model 

of Krugman (1991) and the model of Helpman (1998) have different properties regarding the role of transport costs. In Krugman (1991), as transport costs fall, a core-periphery pattern of the concentration of economic activity becomes feasible. In contrast, in Helpman (1998), lower transport costs lead to a dispersion of economic activity, exactly as shown above. Given that our model is formally equivalent to that of Helpman (1998), under certain parameterization, this different role of trade costs is reflected in our framework as well (see, for example, the patterns of concentration when geography is a line in our early work Allen and Arkolakis (2014)). 

# 4.6 Special case #3: Symmetric trade costs

Suppose that trade costs are symmetric, i.e. for all $i \in \mathcal N$ and $j \in \mathcal { N } .$ , we have $\tau _ { i j } = \tau _ { j i }$ . In this special case, it turns out that the inward and outward market access terms are proportional to each other, as noted by Anderson and Van Wincoop (2003) and Donaldson and Hornbeck (2016). (Similar results can also be obtained when trade costs are “quasi-symmetric”, i.e. $\tau _ { i j } =$ $\tau _ { j i } \tau _ { i } ^ { A } \tau _ { j } ^ { B }$ for some exogenous vectors $\big \{ \tau _ { i } ^ { A } , \tau _ { i } ^ { B } \big \} _ { i \in \mathcal { N } }$ ; see Allen, Arkolakis, and Takahashi (2020)). To see this, we can write the inward and outward market access equations as follows: 

$$
P _ {i} ^ {1 - \sigma} = \sum_ {j \in \mathcal {N}} \tau_ {j i} ^ {1 - \sigma} \Pi_ {j} ^ {\sigma - 1} w _ {j} L _ {j} \tag {43}
$$

$$
\Pi_ {i} ^ {1 - \sigma} = \sum_ {j \in \mathcal {N}} \tau_ {i j} ^ {1 - \sigma} P _ {j} ^ {\sigma - 1} w _ {j} L _ {j} \tag {44}
$$

Suppose that $P _ { i } ^ { 1 - \sigma } = \kappa \Pi _ { i } ^ { 1 - \sigma }$ for some scalar $\kappa > 0$ . Then substituting out $P _ { i } ^ { 1 - \sigma }$ in either equation for $\kappa \Pi _ { i } ^ { 1 - \sigma }$ yields: 

$$
\kappa \Pi_ {i} ^ {1 - \sigma} = \sum_ {j \in \mathcal {N}} \tau_ {i j} ^ {1 - \sigma} \Pi_ {j} ^ {\sigma - 1} w _ {j} L _ {j}, \tag {45}
$$

demonstrating that equation (45) and $P _ { i } ^ { 1 - \sigma } = \kappa \Pi _ { i } ^ { 1 - \sigma }$ is indeed a solution to the system of equations (43) and (44). (It also turns out to be the unique solution, which we show in Section 6.5). 

If $P _ { i } ^ { 1 - \sigma } = \kappa \Pi _ { i } ^ { 1 - \sigma }$ for some $\kappa > 0$ , we can derive a simple log-linear relationship between the population $L _ { i }$ and wage $w _ { i }$ in a location. To do so, combine the welfare equalization equilibrium condition (30) with the outward market access equation (33) using the fact that inward and outward market accesses are proportional. This results in the following expression: 

$$
L _ {i} ^ {1 + (\beta - \alpha) (\sigma - 1)} \propto w _ {i} ^ {(1 - 2 \sigma)} \left(\bar {A} _ {i} / \bar {u} _ {i}\right) ^ {\sigma - 1}, \tag {46}
$$

where the proportionality is determined by the aggregate labor market clearing condition (27) and the choice of numeraire. By creating a log linear relationship between the equilibrium population and wages in a location, equation (46) allows us to collapse the 2N system of equi-

librium equations given by equations (31) and (32) into $N$ that depend only on the equilibrium population in each location: 

$$
\left(\frac {W}{\bar {L} ^ {\alpha + \beta}}\right) ^ {\sigma - 1} l _ {i} ^ {\tilde {\sigma} (1 - \alpha (\sigma - 1) - \beta \sigma)} = \sum_ {j \in \mathcal {N}} \tau_ {i j} ^ {1 - \sigma} \left(\bar {A} _ {i} \bar {u} _ {j}\right) ^ {(\sigma - 1) \tilde {\sigma}} \left(\bar {u} _ {i} \bar {A} _ {j}\right) ^ {\sigma \tilde {\sigma}} l _ {j} ^ {\tilde {\sigma} (1 + \sigma \alpha + \beta (\sigma - 1))}, \tag {47}
$$

where $l _ { i } \equiv L _ { i } / \bar { L }$ is the fraction of the total population residing in $i \in \mathcal N$ and $\tilde { \sigma } \equiv \left( \sigma - 1 \right) / \left( 2 \sigma - 1 \right) .$ There are two things to note about equation (47): first, the aggregate labor endowment $\bar { L }$ has no effect on the distribution of labor (it only affects the interpretation of the endogenous scalar); second, the proportionality factor $\kappa$ between inward and outward market access does not affect the equilibrium, as it scales both sides of the equation equally. If we can find the equilibrium population in each location that solves equation (47), we can then use equation (46) to recover the equilibrium wages in each location. It turns out that an equilibrium exists and is unique as long as $\begin{array} { r } { \left| \frac { 1 + \sigma \alpha + \beta ( \bar { \sigma } - 1 ) } { 1 - \alpha ( \sigma - 1 ) - \beta \sigma } \right| \leq 1 . } \end{array}$ , i.e. the agglomeration forces are no stronger than the congestion forces. We provide a sketch of these results in the following subsection. 

# 4.6.1 [Technical] Existence and Uniqueness

To find the equilibrium population in each location, we first note that the mathematical structure of equation (47) is similar to the general system of equations (31) and (32), where the population shares $l _ { i }$ in location $i \in \mathcal N$ (to some power) are proportional to a weighted sum of the populations in all other locations (to some potentially different power), where the weights depend on the geography and the powers are functions of the model parameters. If the exponents on both sides of the equation are equal, i.e. if $\tilde { \sigma } \left( 1 - \alpha \left( \sigma - 1 \right) - \beta \sigma \right) = \tilde { \sigma } \left( 1 + \sigma \alpha + \beta \left( \sigma - 1 \right) \right)$ , or equivalently $\alpha + \beta = 0 .$ , then equation (47) is a linear system of equations and the equilibrium population is (log) proportional to the eigenvector of the geography matrix, as in Section 4.4. 

Suppose instead that $\alpha + \beta \neq 0 , \mathrm { s o }$ that the exponents on the two sides of equation (47) are not equal. Then the system of equations is nonlinear, which requires a new set of tools. We will discuss these tools in depth in Section 6 below, but offer a glimpse here. 

Let us first tackle the question of existence. Here we do not restrict ourselves to focusing only on interior equilibria and we assume $\alpha \left( \sigma - 1 \right) + \beta \sigma < 1$ . (We refer the interested reader to Karlin and Nirenberg (1967) for a more complicated existence proof when $\alpha ( \sigma - 1 ) + \beta \sigma > 1 )$ . Define the function $F : \Delta ^ { N } \to \Delta ^ { N }$ : 

$$
F (l) _ {i} \equiv \frac {\left(\sum_ {j \in \mathcal {N}} T _ {i j} l _ {j} ^ {\tilde {\sigma} (1 + \sigma \alpha + \beta (\sigma - 1))}\right) ^ {\frac {1}{\tilde {\sigma} (1 - \alpha (\sigma - 1) - \beta \sigma)}}}{\sum_ {i \in \mathcal {N}} \left(\sum_ {j \in \mathcal {N}} T _ {i j} l _ {j} ^ {\tilde {\sigma} (1 + \sigma \alpha + \beta (\sigma - 1))}\right) ^ {\frac {1}{\tilde {\sigma} (1 - \alpha (\sigma - 1) - \beta \sigma)}}}, \tag {48}
$$

where we now define $T _ { i j } \equiv \tau _ { i j } ^ { 1 - \sigma } \left( \bar { A } _ { i } \bar { u } _ { j } \right) ^ { ( \sigma - 1 ) \tilde { \sigma } } \left( \bar { u } _ { i } \bar { A } _ { j } \right) ^ { \sigma \tilde { \sigma } }$ . Equation (48) takes as an input a vector 

of population shares $\{ l _ { i } \} _ { i \in \mathcal { N } }$ and returns as an output another vector of population shares. It is continuous (as long as $1 - \alpha \left( \sigma - 1 \right) - \beta \sigma \neq 0 )$ . A fixed point $\boldsymbol { l } ^ { * } = \boldsymbol { F } \left( \boldsymbol { l } ^ { * } \right)$ of equation (48) is a solution to equation (47), where the endogenous scalar is equal to the denominator, i.e. $\begin{array} { r } { \left( \frac { W } { \bar { L } ^ { \alpha + \beta } } \right) ^ { \sigma - 1 } = \sum _ { i \in \mathcal { N } } \left( \sum _ { j \in \mathcal { N } } T _ { i j } l _ { j } ^ { \tilde { \sigma } ( 1 + \sigma \alpha + \beta ( \sigma - 1 ) ) } \right) ^ { \frac { 1 } { \tilde { \sigma } ( 1 - \alpha ( \sigma - 1 ) - \beta \sigma ) } } } \end{array}$ σ−1 . Because is a continuous function $F$ that maps a convex compact space to itself, by the Brouwer’s fixed-point theorem, there exists a solution to equation (47) as long as $1 - \alpha \left( \sigma - 1 \right) - \beta \sigma \neq 0 .$ . 

Now let us tackle the question of uniqueness. We begin with a carefully chosen change of variables. (Re-)define $x _ { i } \equiv \left( W / \bar { L } ^ { \alpha + \beta } \right) ^ { - \frac { ( \sigma ^ { - } 1 ) ( 1 - \alpha ( \sigma - 1 ) - \beta \sigma ) } { ( 2 \sigma - 1 ) ( \alpha + \beta ) } } l _ { i } ^ { \tilde { \sigma } ( 1 - \alpha ( \sigma - 1 ) - \beta \sigma ) } .$ , so that equation (47) can be written as: 

$$
x _ {i} = \sum_ {j \in \mathcal {N}} T _ {i j} x _ {j} ^ {\rho}, \tag {49}
$$

where $\textstyle \rho \equiv { \frac { 1 + \sigma \alpha + \beta ( \sigma - 1 ) } { 1 - \alpha ( \sigma - 1 ) - \beta \sigma } }$ . There are two things to note about the change of variables: first, the endogenous variable $\{ x _ { i } \} _ { i \in \mathcal { N } }$ includes both the population $L _ { i }$ and the global welfare $W$ ; and second, if we can solve for $\{ x _ { i } \} _ { i \in \mathcal { N } }$ using equation (49), then we can separately recover the population and the welfare using the aggregate labor market clearing condition (27) (and then recover the equilibrium wages using equation (46)). 

We now establish sufficient conditions for the uniqueness in a system of equations of the type in equation (47) by providing conditions under which the system is a contraction. The following contradiction argument is loosely based on Karlin and Nirenberg (1967). Define a (new) function $F : \mathbb { R } ^ { N } \to \mathbb { R } ^ { N }$ : 

$$
F (\mathbf {y}) _ {i} \equiv \log \left(\sum_ {j \in \mathcal {N}} T _ {i j} \exp (\rho y _ {j})\right) \tag {50}
$$

and a distance metric $d \left( \mathbf { y } ^ { A } , \mathbf { y } ^ { B } \right) \equiv \mathsf { m a x } _ { i \in \mathcal { N } } \left| y _ { i } ^ { A } - y _ { i } ^ { B } \right|$ . If there exists a $c < 1$ such that for all $\mathbf { y } ^ { A } , \mathbf { y } ^ { B } \in \mathbb { R } ^ { N }$ that $d \left( F \left( \mathbf { y } ^ { A } \right) , F \left( \mathbf { y } ^ { B } \right) \right) \leq c \times d \left( \mathbf { y } ^ { A } , \mathbf { y } ^ { B } \right)$ , then by the Banach fixed-point theorem, there exists a unique vector $\mathbf { y } ^ { * }$ such that $\mathbf { y } ^ { * } = F \left( \mathbf { y } ^ { * } \right)$ , which can be found by starting at any initial guess $\mathbf { y } _ { 0 } \in \mathbb { R } ^ { N }$ and iteratively updating the guess using the function, i.e. if $\mathbf { y } _ { n } = F \left( \mathbf { y } _ { n - 1 } \right)$ then $\operatorname* { l i m } _ { n \to \infty } \mathbf { y } _ { n } = \mathbf { y } ^ { * }$ . Note that once the fixed point $\mathbf { y } ^ { * }$ to equation (50) has been found, then the set $x _ { i } ^ { * } \equiv \exp ( y _ { i } ^ { * } )$ will solve equation (47). 

Choose any $\mathbf { y } ^ { A } , \mathbf { y } ^ { B } \in \mathbb { R } ^ { N }$ . Then we have: 

$$
\left. d \left(F \left(\mathbf {y} ^ {A}\right), F \left(\mathbf {y} ^ {B}\right)\right) = \max  _ {i \in \mathcal {N}} \left| \log \left(\sum_ {j \in \mathcal {N}} C _ {i j} \exp \left(\rho \left(y _ {j} ^ {A} - y _ {j} ^ {B}\right)\right)\right) \right|, \right. \tag {51}
$$

where $C _ { i j } \equiv T _ { i j } \exp \left( \rho y _ { j } ^ { B } \right) / \sum _ { k } T _ { i k } \exp \left( \rho y _ { k } ^ { B } \right)$ . Because $\begin{array} { r } { \sum _ { j \in \mathcal { N } } C _ { i j } = 1 } \end{array}$ , we can bound the distance 

as follows: 

$$
\left. \right. \max  _ {i \in \mathcal {N}} \left| \log \left(\sum_ {j \in \mathcal {N}} C _ {i j} \exp \left(\rho \left(y _ {j} ^ {A} - y _ {j} ^ {B}\right)\right)\right)\right| \leq | \rho | \max  _ {i \in \mathcal {N}} \left| y _ {j} ^ {A} - y _ {j} ^ {B} \right|, \tag {52}
$$

or equivalently: 

$$
d \left(F \left(\mathbf {y} ^ {A}\right), F \left(\mathbf {y} ^ {B}\right)\right) \leq | \rho | d \left(\mathbf {y} ^ {A}, \mathbf {y} ^ {B}\right). \tag {53}
$$

As a result, if $| \rho | < 1 ,$ , then there exists a unique interior equilibrium. 

What does it mean that $| \rho | < 1 ?$ Recall from above that $\textstyle \rho \equiv { \frac { 1 + \sigma \alpha + \beta ( \sigma - 1 ) } { 1 - \alpha ( \sigma - 1 ) - \beta \sigma } }$ . Focusing on the case where $\sigma \geq 1$ and $\alpha , \beta \in [ - 1 , 1 ] ,$ , this condition will hold as long as $\alpha + \beta < 0 ,$ i.e. the externalities are net congestive rather than net agglomerative. Loosely speaking, there is a unique equilibrium as long as individuals would prefer to spread out across locations. 

What happens if $| \rho | > 1 ?$ First, from above, there will still exist interior equilibria (as long as $1 - \alpha \left( \sigma - 1 \right) - \beta \sigma \neq 0 )$ ; however, if the $\rho$ becomes too large (so that $1 - \alpha \left( \sigma - 1 \right) - \beta \sigma < 0 ) .$ , then a “black hole” equilibrium where all agents reside in a single location becomes possible (as noted by Fujita, Krugman, and Venables, 1999) and any interior equilibria become unstable (in the sense that welfare would increase if individuals all moved to the same location). Second, there will exist certain geographies for which there are guaranteed to be multiple interior equilibria; in the words of Allen, Arkolakis, and Li (2024), the condition is sufficient and “globally” necessary. For example, Allen and Arkolakis (2014) show that when $\rho > 1 ,$ , there are an infinite number of equilibria featuring greater economic activity near an arbitrarily chosen location when identical locations are arrayed around a ring. But what remains an open question in the field is how the particular geography of the world determines whether or not there will be multiplicity if $| \rho | > 1$ . Recent analysis by Kucheryavyy, Lyn, and Rodr´ıguez-Clare (2024) of a two location model, for example, has shown that uniqueness can be maintained as long as trade costs are sufficiently low. 

# 5 New models, same old equilibrium

The simple quantitative economic geography model presented in the previous section relies on a number of strong assumptions (e.g. perfect labor mobility, identical agents, perfect competition, etc.). But it turns out to be remarkably general. This is because a number of interesting departures from these strong assumptions are mathematically equivalent to the mathematical structure developed above. This is because these departures act identically to the productivity and/or amenity externalities in the baseline model. We now briefly describe seven alternative models that are isomorphic (and an eighth that is not). The goal is to understand better how different economic assumptions relate to agglomeration or dispersion forces present in the quantitative economic geography model, thereby affecting the labor supply and demand curves determining the spatial equilibrium. 

# 5.1 Ricardian comparative advantage

In the simple model in Section 4, trade arises because consumers have a love of variety and each location produces a differentiated variety. However, one can instead develop a mathematically equivalent model with perfect competition and Ricardian comparative advantage. The key insight—developed by Eaton and Kortum (2002)—is to assume there is a continuum of goods sold competitively by each region. Each good is provided by the least cost supplying country in each destination. If productivities of these products $\omega$ follow a Frechet distribution, i.e. $\varepsilon _ { i } \left( \omega \right) \sim e ^ { - A _ { i } z ^ { - \theta } } .$ , then Eaton and Kortum (2002) show that expenditure share of country $j$ on goods produced in $i$ is given by: 

$$
\lambda_ {i j} = \frac {\left(\frac {w _ {i} \tau_ {i j}}{A _ {i}}\right) ^ {- \theta}}{\sum_ {i} \left(\frac {w _ {i ^ {\prime}} \tau_ {i ^ {\prime} j}}{A _ {i ^ {\prime}}}\right) ^ {- \theta}}. \tag {54}
$$

Equation (54) is isomorphic to the gravity equation (21) with $\theta = \sigma - 1$ . As a result, the Armington justification for trade in the baseline model presented in Section 4 and this Ricardian trade model have the same implications for aggregate trade elasticity. 

# 5.2 Non-tradable sector

In the simple model in Section 4, all goods produced in the economy are tradable. But in reality, some goods cannot be (easily) traded. To allow for this possibility we follow Helpman (1998). The key idea is that the workers spend a constant fraction of their income $\delta$ on the differentiated good and the rest to a locally produced good that is not traded with the rest of the locations, such as services or housing. The non-tradable sector is perfectly competitive and its returns are equally distributed across workers that reside in each location. This change affects only the labor supply of the economy. The non-tradable sector is akin to introducing a congestion amenity externality and is equivalent to our baseline model so long as $\beta = - \left( 1 - \delta \right) / \delta ,$ i.e. the degree of congestion in consumption equals the negative of the ratio of non-tradeable to tradeable good shares. 

# 5.3 Fixed factors of production

In the simple model in Section 4, labor is the only factor of production. To incorporate factors of production such as land and capital, we follow Donaldson and Hornbeck (2016). The key idea is that the production function features an additional fixed factor, e.g capital, $Y _ { i } = A _ { i } K _ { i } ^ { \tilde { \alpha } } L _ { i } ^ { 1 - \tilde { \alpha } }$ , and $\tilde { \alpha }$ is the share of the fixed factor. Under that configuration the new model is isomorphic to the baseline model presented in the previous section if $\tilde { \alpha } = - \alpha ,$ i.e. the fixed factor acts as a negative productivity externality. 

# 5.4 Endogenous labor supply

In the simple model in Section 4, agents inelastically supply the entirety of their labor endowment. To incorporate an endogenous labor supply decision, suppose agents have preferences over both their consumption and their labor supply as follows: 

$$
U \left(C _ {i}, L _ {i}\right) = C _ {i} u _ {i} - \frac {\left(L _ {i}\right) ^ {1 + \eta}}{1 + \eta}, \tag {55}
$$

where the first term is the CES aggregator as specified in equation (17), with budget constraint: 

$$
\sum_ {i \in \mathcal {N}} p _ {i j} c _ {i j} = P _ {j} C _ {j} = w _ {j} L _ {j}. \tag {56}
$$

Agents now decide how many labor units to devote to work. Straightforward application of first order conditions of the constrained maximization problem of equations (55) and (56) yields a labor supply function $L _ { i } = ( w _ { i } u _ { i } / P _ { i } ) ^ { 1 / \eta }$ . Expressing this equation in the form of equation (36), we obtain an exact equivalence if we set $- \beta = \eta ,$ i.e. an endogenous labor supply acts isomorphically to a congestion externality on amenities. 

# 5.5 Idiosyncratic preferences

Rather than assuming that all agents have identical preferences across locations, we can extend the model to allow agents to have heterogeneous preferences across locations, as in Redding (2016). Further discussion of this approach and its use in urban settings can be found in Ahlfeldt, Redding, Sturm, and Wolf (2015) and are discussed in detail in Redding (2025). 

The key idea of incorporating heterogeneity in agent’s preferences is to assume that the welfare of each worker $\omega$ choosing to reside in location $j$ is $\begin{array} { r } { \frac { w _ { j } } { P _ { j } } u _ { j } \times \varepsilon _ { j } \left( \omega \right) . } \end{array}$ wj , where $\varepsilon _ { j } \left( \omega \right)$ is assumed to be an idiosyncratic preference term that is agent-specific. If we assume that distribution is Frechet, then following Eaton and Kortum (2002), it is straightforward to show that the share of workers that choose i as their location is 

$$
L _ {i} / \bar {L} = \frac {\left(w _ {i} u _ {i} / P _ {i}\right) ^ {\theta}}{\sum_ {j} \left(w _ {j} u _ {j} / P _ {j}\right) ^ {\theta}}. \tag {57}
$$

Equation (57) is isomorphic to the labor supply equation (36) in the baseline model with the strength of the amenity externality set to $\beta = - 1 / \theta$ . Intuitively, the more people who reside in a location, the lower on average their idiosyncratic preferences are for that location, which acts in a mathematically equivalent way to a congestion externality in amenities. 

# 5.6 Idiosyncratic productivities

Rather than assuming that all agents are equally productive in all locations, we can also extend the model to allow agents to be heterogeneous in their productivities, as in Bryan and Morten (2019). The approach works in a similar way as incorporating heterogeneous preferences, where we assume that the welfare of each worker $\omega$ choosing to reside in location $j$ is wjP uj × εj (ω), but we now interpret εj (ω) as the idiosyncratic productivity of agent ω in lo- $\begin{array} { r } { \frac { w _ { j } } { P _ { j } } u _ { j } \times \varepsilon _ { j } \left( \omega \right) , } \end{array}$ $\varepsilon _ { j } \left( \omega \right)$ $\omega$ cation $j$ . In this case, the number of workers choosing a location still follows equation (57). However, because the workers who choose a location will tend to be more productive in that location, the effective units of labor in location $i , \tilde { L } _ { i } ,$ are: 

$$
\tilde {L} _ {i} = L _ {i} ^ {\frac {\theta - 1}{\theta}}. \tag {58}
$$

As it is the effective units of labor in a location that determines the amount of goods produced and the payments to labor, we can write the supply and demand curves in terms of effective units of labor. Using equation (58), it is straightforward to show that these labor and supply curves are equivalent to the baseline model when $\beta = - 1 / \left( \theta - 1 \right) .$ , albeit with a modified aggregate effective labor market clearing constraint $\begin{array} { r } { \sum _ { i \in \mathcal { N } } \tilde { L } _ { i } ^ { \frac { \theta } { \theta - 1 } } = \bar { L } ^ { \frac { \theta - 1 } { \theta } } } \end{array}$ . Intuitively, heterogeneity in productivities, like heterogeneity in preferences, acts as a congestion externality in amenities. 

# 5.7 Monopolistic competition with free entry

In the simple model in Section 4, we assumed perfect competition. However, that model is also isomorphic to a model featuring monopolistic competition and free entry, as in the Krugman (1991) model discussed in Section 3.2. Comparing the labor supply and demand equations (14) and (15) from Section 3.2 to the labor supply and demand equations (35) and (36) from Section 4, we see that the former is isomorphic to the latter when $\alpha = 1 / \left( \sigma - 1 \right) .$ , i.e. monopolistic competition with free entry acts equivalently to a productivity agglomeration force. 

# 5.8 Round-about production

The previous models emphasize the flexibility of the simple model in Section 4 in relaxing a variety of assumptions. We conclude with an example of an extension that is not isomorphic. In the simple model, we abstracted from intermediate goods production. Following Krugman and Venables (1995); Eaton and Kortum (2002), suppose we incorporate intermediate goods in the following simple way: we assume that goods are used for final consumption but also as intermediate inputs in the production of other goods. In particular, we assume that firms produce their products using a Cobb-Douglas combination of labor and intermediate goods with share $\gamma$ going to labor and share $1 - \gamma$ going to intermediate goods, which we further 

assume are aggregated in the same way consumers aggregate their final goods. 

Working through the derivations of Section 4.2.2 in this alternative framework yields the following system of equations that define the equilibrium: 

$$
W ^ {\gamma (\sigma - 1)} w _ {i} ^ {\sigma} L _ {i} = \sum_ {j} \tau_ {i j} ^ {1 - \sigma} A _ {i} ^ {\sigma - 1} \left(u _ {i} ^ {1 - \gamma}\right) ^ {1 - \sigma} u _ {j} ^ {\sigma - 1} w _ {j} ^ {\sigma} L _ {j} \tag {59}
$$

$$
W ^ {\gamma (\sigma - 1)} w _ {i} ^ {1 - \sigma} = \sum_ {j} \tau_ {j i} ^ {1 - \sigma} A _ {j} ^ {\sigma - 1} u _ {i} ^ {\sigma - 1} u _ {j} ^ {(1 - \gamma) 1 - \sigma} w _ {j} ^ {1 - \sigma}. \tag {60}
$$

A comparison of equations (59) and (60) to the equilibrium conditions (31) and (32) of the baseline model shows two differences. First, the elasticity of welfare to the endogenous scalar changes due to the presence of roundabout production. This change is to be expected, as the presence of intermediate inputs magnifies the welfare gains from trade (see for example Arkolakis, Costinot, and Rodr´ıguez-Clare (2012); Ossa (2015)). The second difference is a bit more surprising: now the amenities in both locations affect both equilibrium conditions. Intuitively, this is because round-about production, combined with welfare equalization, introduces a new way through which amenities elsewhere affect the price index in a location. Hence, a formal isormorphism to the baseline model only occurs in the special case where amenities are identical in all locations. 

# 6 A workhorse quantitative economic geography framework

We now present a workhorse quantitative economic geography framework, special cases of which include the models presented above. 

# 6.1 Setup

As in Section 4, we consider a world comprising a compact set N of locations and inhabited by workers who consume and produce in the location in which they reside. 

There are two components of the framework. The first component comprises the spatial linkages between locations. We assume that locations are linked together through the flow of goods. Consistent with the evidence presented in Section 2.2 and 2.3, we assume that the value of trade flow from $i \in \mathcal N$ to $j \in \mathcal N$ follow a gravity equation: 

$$
X _ {i j} = T _ {i j} \times \frac {Y _ {i}}{M A _ {i} ^ {\text {o u t}}} \times \frac {E _ {j}}{M A _ {j} ^ {\text {i n}}}, \tag {61}
$$

where $T _ { i j } \le 1$ is an (exogenous) trade friction, $Y _ { i }$ is the aggregate value of economic output in $i \in \mathcal { N } , E _ { j }$ is the aggregate value of economic consumption in $j \in \mathcal { N } , M A _ { i } ^ { o u t }$ is the outward 

market access in $i \in \mathcal N$ and $M A _ { j } ^ { i n }$ is the inward market access in $j \in \mathcal N$ . 

The second component comprises the labor markets in each location. We assume that in each location $i \in \mathcal { N }$ , the (inverse) labor supply curve can be written as the following log-linear function: 

$$
\ln w _ {i} = \varepsilon_ {\text {l o c a l}} ^ {S} \ln L _ {i} - \varepsilon_ {\text {g l o b a l}} ^ {S} \ln M A _ {i} ^ {i n} - \ln C _ {i} ^ {S} - \ln \phi^ {S}, \tag {62}
$$

where $\varepsilon _ { l o c a l } ^ { S }$ is the elasticity of wages to the labor supplied, εSglobal is the elasticity of wages to the $\varepsilon _ { g l o b a l } ^ { S }$ inward market access, $C _ { i } ^ { S }$ is an (exogenous) labor supply factor, and $\phi ^ { S }$ is an endogenous scalar common to all locations to ensure an aggregate labor market condition (or aggregate welfare normalization) is satisfied. Note that our sign convention is such that each parameter would typically be positive, e.g. the supply curve is normally upward sloping and better inward market access or higher $C _ { i } ^ { S }$ would normally imply that labor would be willing to accept a lower wage to work in a location. 

Similarly, we assume that in each location $i \in \mathcal N$ ,the (inverse) labor demand curve can be written as: 

$$
\ln w _ {i} = - \varepsilon_ {\text {l o c a l}} ^ {D} \ln L _ {i} + \varepsilon_ {\text {g l o b a l}} ^ {D} \ln M A _ {i} ^ {\text {o u t}} + \ln C _ {i} ^ {D} + \ln \phi^ {D}, \tag {63}
$$

where $\varepsilon _ { l o c a l } ^ { D }$ is the elasticity of wages to labor demanded, εDglobal is the elasticity of wages to the $\varepsilon _ { g l o b a l } ^ { D }$ outward market access, $C _ { i } ^ { D }$ is an (exogenous) labor demand factor, and $\phi ^ { D }$ is an endogenous scalar common to all locations to ensure the choice of numeraire is satisfied. Again, we maintain a sign convention that each parameter should typically be positive. 

It is immediately evident that the Rosen (1979)-Roback (1982) model of location choice presented in Section (3.1), the Krugman (1991) core-periphery model presented in Section (3.2) and the Allen and Arkolakis (2014) quantitative economic geography model presented in Section 4 are all special cases of this workhorse quantitative economic geography framework. Table 1 shows how the parameters in each of these models (and others) map to the local and global supply and demand elasticities in the workhorse framework. 

A couple of notes are necessary. First, note that when the global elasticities are absent, i.e. εSglobal $\varepsilon _ { g l o b a l } ^ { S } = \varepsilon _ { l o c a l } ^ { D } = 0$ Dloca (as in Rosen (1979)-Roback (1982) model), the equilibrium is trivial to solve, so in what follows we focus on the more interesting case where they are non zero. Second, note that unlike the simple model from Section 4, here we do not need to distinguish between the iceberg trade costs $\tau _ { i j }$ and the trade elasticity $\sigma - 1$ : for the purposes of determining how the equilibrium distribution of economic activity is affected by a change in the underlying geography, only their composite $T _ { i j } = \tau _ { i j } ^ { 1 - \sigma }$ is what matters. This is not to say, however, that the trade elasticity does not matteit affects the model elasticities $( \varepsilon _ { l o c a l } ^ { S } , \varepsilon _ { g l o b a l } ^ { S } , \varepsilon _ { l o c a l } ^ { D } , \varepsilon _ { g l o b a l } ^ { D } ) ;$ ，mS. sis. Indeed, it matters in two ways: first, second, it affects the mapping from the endogenous scalars to agent welfare. 


ngs of models to the workhorse quantitative economic geogra


<table><tr><td>Paper</td><td>εSlocal</td><td>εSglobal</td><td>εDlocal</td><td>εDglobal</td><td>Model parameters</td><td>Notes</td></tr><tr><td>Rosen (1979) and Roback (1982)</td><td>-β</td><td>0</td><td>α</td><td>0</td><td>β: Externality of population on utilityα: Capital share in production</td><td>Based on the version developed in Section 3.1.</td></tr><tr><td>Krugman (1991)</td><td>0</td><td>1/σ-1</td><td>0</td><td>1/σ</td><td>σ: Elasticity of substitution</td><td>Based on the extended version developed in Section 3.2, assuming no agricultural sector.</td></tr><tr><td>Helpman (1998)</td><td>1/μ/μ</td><td>1/σ-1</td><td>0</td><td>1/σ</td><td>σ: Elasticity of substitutionμ: Expenditure share of manufacturing goods</td><td>Based on the multi-regional version developed by Redding and Sturm (2008)</td></tr><tr><td>Eaton and Kortum (2002)</td><td>0</td><td>1/θ</td><td>1/1+θ</td><td>1/1+θ</td><td>θ: Frechet shape parameter</td><td>We differ from Eaton and Kortum (2002) by assuming (1) no round-about production, (2) no non-manufacturing sector, and (3) assuming perfect la-bor mobility across regions.</td></tr><tr><td>Allen and Arkolakis (2014)</td><td>-β</td><td>1/σ-1</td><td>1-α(σ-1)/σ</td><td>1/σ</td><td>σ: Elasticity of substitutionα: Externality of population on productivityβ: Externality of population on utility</td><td></td></tr><tr><td>Donaldson and Hornbeck (2016)</td><td>0</td><td>1/θ</td><td>1+αθ/1+(α+γ)/θ</td><td>1/1+(α+γ)/θ</td><td>θ: Frechet shape parameterα: Expenditure share of landγ: Expenditure share of labor</td><td></td></tr><tr><td>Redding (2016)(with constant return to scale production)</td><td>1+(1-α)ε/αε</td><td>1/θ</td><td>1/1+θ</td><td>1/1+θ</td><td>θ: Frechet shape parameter of productivity drawα: Expenditure share of manufacturing goodsε: Frechet shape parameter of preference draw</td><td rowspan="2">Redding (2016) develops two models, the first with Eaton-Kortum style constant return to scale technol-ogy, and the second with Krugman-style increasing return to scale technology.</td></tr><tr><td>Redding (2016)(with increasing return to scale production)</td><td>1+(1-α)ε/αε</td><td>1/σ-1</td><td>0</td><td>1/σ</td><td>σ: Elasticity of substitutionα: Expenditure share of manufacturing goodsε: Frechet shape parameter of preference draw</td></tr><tr><td>Idiosyncratic productivitya la Bryan and Morten (2019)</td><td>1-βθ/θ-1</td><td>0</td><td>1-α(σ-1)/θ</td><td>0</td><td>θ: θ/1-ρ, β the shape parameter and ρ the correlation parameter of the multi-variate Frechet distributionσ: Elasticity of substitutionα: Externality of population on productivityβ: Externality of population on utility</td><td>We differ from Bryan and Morten (2019) by assum-ing there is no migration cost.</td></tr></table>


y  pp  y ws how the workhorse gravity framework comprises a number of notable economi 


# 6.2 Equilibrium

Analogously to Section 4, we define the geography of the world as the set $\left\{ \left\{ T _ { i j } \right\} _ { i , j \in \mathcal { N } } , \left\{ C _ { i } ^ { S } , C _ { i } ^ { D } \right\} _ { i \in \mathcal { N } } \right\} ,$ the model elasticities as the set $\left\{ \varepsilon _ { l o c a l } ^ { S } , \varepsilon _ { g l o b a l } ^ { S } , \varepsilon _ { l o c a l } ^ { D } , \varepsilon _ { g l o b a l } ^ { D } \right\} ,$ , and the distribution of economic activity as the set $\left\{ \left\{ L _ { i } , w _ { i } , M A _ { i } ^ { i n } , M A _ { i } ^ { o u t } \right\} _ { i \in \mathcal { N } } , \phi ^ { S } , \phi ^ { D } \right\}$ . 

Given any geography and model parameters, an equilibrium is a distribution of economic activity such that the following four conditions hold: 

1. Goods markets clear, i.e. $\textstyle Y _ { i } = \sum _ { j \in { \mathcal { N } } } X _ { i j }$ and $\textstyle E _ { i } = \sum _ { j \in { \mathcal { N } } } X _ { j i }$ for all $i \in \mathcal N$ . 

2. Labor markets clear, i.e. labor supply (defined by equation (62)) and labor demand (defined by equation (63)) are equalized for all $i \in \mathcal N$ . 

3. Output and expenditure is equal to the total payments to labor, i.e. $w _ { i } L _ { i } = Y _ { i } = E _ { i }$ for all $i \in \mathcal N$ . 

4. The endogenous scalars $\phi ^ { S }$ and $\phi ^ { D }$ are determined by a choice of numeraire and a labor market clearing condition, i.e. $\begin{array} { r } { \sum _ { i \in \mathcal { N } } L _ { i } = \bar { L } } \end{array}$ . 

Combining the equilibrium condition 1 with the gravity equation (61) for bilateral trade flows yields the following system of equations relating inward and outward market access to economic activity in all locations: 

$$
M A _ {i} ^ {\text {o u t}} = \sum_ {j \in \mathcal {N}} T _ {i j} \times \frac {E _ {j}}{M A _ {j} ^ {\text {i n}}} \tag {64}
$$

$$
M A _ {i} ^ {i n} = \sum_ {j \in \mathcal {N}} T _ {j i} \times \frac {Y _ {j}}{M A _ {j} ^ {\text {o u t}}}, \tag {65}
$$

which is the workhorse form of equations (43) and (44) in Section 4.6. 

In the following three subsections, we summarize the state of knowledge about the equilibrium properties of this framework. If you would prefer to skip them, the punchline is as follows: (1) under very mild conditions, an equilibrium will exist; (2) as long as agglomeration forces are not “too” strong, the equilibrium will be unique; and (3) these results hold for any number of locations and any geography (although we still do not know if more precise conditions can be given for uniqueness that depend on the particular geography being considered). 

# 6.3 [Technical] Re-defining the equilibrium system

We can then solve for the inward market access in the labor supply equation (62) and the outward market access in the labor demand equation (63), impose equilibrium condition (3), and 

substitute into the market access equations (64) and (65) to yield: 

$$
\begin{array}{l} \tilde {\lambda} w _ {i} ^ {\frac {1}{\varepsilon^ {D} \text {g l o b a l}}} L _ {i} ^ {\frac {\varepsilon^ {D} \text {l o c a l}}{\varepsilon^ {D} \text {g l o b a l}}} = \sum_ {j \in \mathcal {N}} K _ {i j} w _ {j} ^ {\frac {\varepsilon^ {S} \text {g l o b a l} + 1}{\varepsilon^ {S} \text {g l o b a l}}} L _ {j} ^ {\frac {\varepsilon^ {S} \text {g l o b a l} - \varepsilon^ {S} \text {l o c a l}}{\varepsilon^ {S} \text {g l o b a l}}} (66) \\ \tilde {\lambda} w _ {i} ^ {- \frac {1}{\varepsilon_ {g l o b a l} ^ {S}} \frac {\varepsilon_ {l o c a l} ^ {S}}{\varepsilon_ {g l o b a l} ^ {S}}} L _ {i} ^ {- \frac {\varepsilon_ {g l o b a l} ^ {D} - 1}{\varepsilon_ {g l o b a l} ^ {D}} \frac {\varepsilon_ {g l o b a l} ^ {D} - \varepsilon_ {l o c a l} ^ {D}}{\varepsilon_ {g l o b a l} ^ {D}}} = \sum_ {j \in \mathcal {N}} K _ {j i} w _ {j} ^ {- \frac {\varepsilon_ {j i} ^ {D} - 1}{\varepsilon_ {j i} ^ {D}} \frac {\varepsilon_ {j i} ^ {D} - \varepsilon_ {j i} ^ {D}}{\varepsilon_ {j i} ^ {D}}} L _ {j} ^ {- \frac {\varepsilon_ {j i} ^ {D} - 1}{\varepsilon_ {j i} ^ {D}}} (67) \\ \end{array}
$$

where λ˜ ≡  ϕD −1/εDglob $\tilde { \lambda } \equiv \left( \phi ^ { D } \right) ^ { - 1 / \varepsilon _ { g l o b a l } ^ { D } } \left( \phi ^ { S } \right) ^ { - 1 / \varepsilon _ { g l o b a l } ^ { S } }$ −1/εSglo is an endogenous scalar and $K _ { i j } \equiv T _ { i j } \left( C _ { i } ^ { D } \right) ^ { 1 / \varepsilon _ { g l o b a l } ^ { D } } \left( C _ { j } ^ { S } \right) ^ { 1 / \varepsilon _ { g l o b a l } ^ { S } }$ 1/εD 1/εSglobal is a function of the geography and the model parameters. It is not surprising that equations (66) and (67) bear a strong resemblance to the equilibrium system of equations (31) and (32) for the simple quantitative economic geography framework presented in Section 4, as that simple framework is a special case of this more general framework. Both systems have a similar mathematical structure, where a particular log linear combination of wages and population in one location are equal to a weighted sum of a different log linear combination of wages and population in all other location. In both systems, the weights of the sum depend on the geography of the world, whereas the exponents depend on the labor supply and demand elasticities. In the simple quantitative economic geography framework presented in Section 4, these supply and demand elasticities were constructed from particular assumptions regarding consumer preferences and production and amenity externalities. In contrast, the workhorse framework considers the supply and demand elasticities directly. 

As in Section (4.6), it turns out to be helpful to rewrite the equilibrium system in terms of shares rather than levels. Define $y _ { i } \equiv w _ { i } L _ { i } / Y ^ { W }$ as the share of world income in location $i \in \mathcal N$ and recall that $l _ { i } \equiv L _ { i } / \bar { L }$ is the population share. We can then rewrite equations (66) and (67) as follows: 

$$
\begin{array}{l} \lambda y _ {i} ^ {\frac {1}{\varepsilon_ {g l o b a l} ^ {D}} \frac {\varepsilon_ {l o c a l} ^ {D} - 1}{\varepsilon_ {g l o b a l} ^ {D}}} l _ {i} = \sum_ {j \in \mathcal {N}} K _ {i j} y _ {j} ^ {\frac {\varepsilon_ {g l o b a l} ^ {S} + 1}{\varepsilon_ {g l o b a l} ^ {S}}} l _ {j} - \frac {1 + \varepsilon_ {l o c a l} ^ {S}}{\varepsilon_ {g l o b a l} ^ {S}} (68) \\ \lambda y _ {i} ^ {- \frac {1}{\varepsilon_ {\text {g l o b a l}} ^ {S}} \frac {\varepsilon_ {\text {l o c a l}} ^ {S} + 1}{\varepsilon_ {\text {g l o b a l}} ^ {S}}} l _ {i} ^ {} = \sum_ {j \in \mathcal {N}} K _ {j i} y _ {j} ^ {\frac {\varepsilon_ {\text {g l o b a l}} ^ {D} - 1}{\varepsilon_ {\text {g l o b a l}} ^ {D}}} l _ {j} ^ {\frac {1 - \varepsilon_ {\text {l o c a l}} ^ {D}}{\varepsilon_ {\text {g l o b a l}} ^ {D}}}, (69) \\ \end{array}
$$

where λ ≡ λ˜ $\begin{array} { r } { \lambda \equiv \tilde { \lambda } \frac { \frac { \varepsilon _ { l o c a l } ^ { D } - 1 } { \varepsilon _ { g l o b a l } ^ { D } } + \frac { 1 + \varepsilon _ { l o c a l } ^ { S } } { \varepsilon _ { g l o b a l } ^ { S } } } { 1 + \frac { 1 } { \varepsilon _ { g l o b a l } ^ { S } } - \frac { 1 } { \varepsilon _ { g l o b a l } ^ { D } } } } \end{array}$ εDlocal −1 1+εSlocal . This re-definition of the equilibrium in terms of output and population shares emphasizes that the aggregate scale of the economy (in either its nominal terms or in its aggregate labor endowment) has no effect on the equilibrium distribution of economic activity in the economy, as changes in either $\bar { L }$ or $Y ^ { W }$ simply change the endogenous scalar. In 

what follows, we view the endogenous scalar $\lambda$ as a measure of the welfare of the equilibrium, although how exactly it maps into welfare will depend on the details of the micro-foundation being considered. For example, in the simple framework presented in Section 4, $\lambda$ is log proportional to the welfare of the system, with the proportionality equal to the trade elasticity $( \sigma - 1 )$ . We refer the interested reader to Fajgelbaum and Gaubert (2025) for a thorough discussion of the normative properties of spatial models. 

# 6.4 [Technical] Existence

We first consider the question of whether or not there exists a solution to the system of equations defined by equations (68) and (69). Define a (new) function $F ( \pmb { y } , \pmb { l } ) : \Delta ^ { N } + \Delta ^ { N }  \Delta ^ { N } + \Delta ^ { N }$ which maps from any combination of output shares and labor shares into a new set of output 

εSglobal +1 1+εSlocal shares and labor shares as follows. Let F1 (y, l)i ≡ ∑j∈N Kijyj εSglobal l − εSgl o b al and $F _ { 2 } \left( y , l \right) _ { i } \equiv$ D εglobal − 1 1−εDlocal ∑j∈N Kji y εDglobali l εDglobalj D for all $i \in \mathcal N$ . Then define: 

$$
F (\boldsymbol {y}, \boldsymbol {l}) _ {i} = \left\{ \begin{array}{l l} F _ {1} (\boldsymbol {y}, \boldsymbol {l}) _ {i} ^ {b _ {1 1}} F _ {2} (\boldsymbol {y}, \boldsymbol {l}) _ {i} ^ {b _ {1 2}} / \sum_ {i \in \mathcal {N}} F _ {1} (\boldsymbol {y}, \boldsymbol {l}) _ {i} ^ {b _ {1 1}} F _ {2} (\boldsymbol {y}, \boldsymbol {l}) _ {i} ^ {b _ {1 2}} & \text {f o r} i \in \{1, \dots , N \} \\ F _ {1} (\boldsymbol {y}, \boldsymbol {l}) _ {i} ^ {b _ {2 1}} F _ {2} (\boldsymbol {y}, \boldsymbol {l}) _ {i} ^ {b _ {2 2}} / \sum_ {i \in \mathcal {N}} F _ {1} (\boldsymbol {y}, \boldsymbol {l}) _ {i} ^ {b _ {2 1}} F _ {2} (\boldsymbol {y}, \boldsymbol {l}) _ {i} ^ {b _ {2 2}} & \text {f o r} i \in \{N + 1, \dots , 2 N \} \end{array} , \right. \tag {70}
$$

1/εDglobal εDlocal − 1 1 where $b _ { k l }$ is the $( k , l ) ^ { t h }$ element of the $2 \times 2$ matrix B ≡ 1 εSlocal +1 εDglobal . Intuitively, εSglobal εSglobal what the function $F$ does is take any set of output and labor shares as an input and return as an output the implied output and labor shares using equations (68) and (69). As a result, a fixed point of $F \left( y ^ { * } , l ^ { * } \right) = F \left( y ^ { * } , l ^ { * } \right)$ is a set of equilibrium output and labor shares that solves equations (68) and (69). As in Section (4.6), as long as B exists (which allows us to invert the left hand side of equations (68) and (69) to recover the implied output and labor shares), then because $F$ is a continuous operator mapping a convex compact space to itself, Brouwer’s fixed point theorem ensures that there exists an equilibrium. We refer the interested reader to Allen, Arkolakis, and Takahashi (2020) for details, including demonstrating that the equilibrium is an interior one. 

# 6.5 [Technical] Uniqueness

We now provide sufficient conditions for the uniqueness of the equilibrium. We provide a sketch of the proof below and refer the interested reader to Allen, Arkolakis, and Li (2024) for 

details. We begin with a change of variables. Define $x _ { i } \equiv y _ { i } ^ { ~ } \stackrel { \varepsilon _ { l o c a l } ^ { D } - 1 } { \varepsilon _ { g l o b a l } ^ { D } } l _ { i } ^ { ~ \varepsilon _ { g l o b a l } ^ { D } }$ and zi ≡ y− εSglobal l εSglobal εSlocal +1 so that equations (68) and (69) can be written as: 

$$
\lambda x _ {i} = \sum_ {j \in \mathcal {N}} K _ {i j} x _ {j} ^ {a _ {1 1}} z _ {j} ^ {a _ {1 2}} \tag {71}
$$

$$
\lambda z _ {i} = \sum_ {j \in \mathcal {N}} K _ {j i} x _ {j} ^ {a _ {2 1}} z _ {j} ^ {a _ {2 2}}, \tag {72}
$$

where A = a11 a12! ≡  εSglobal +1εSglobalεD −1 εSglobal1−εDlocal 1+εSlocal 1εDglobal εDlocal −1εDglobalεSlocal +1  − 1 a21 a22 1 1 εDglobal εDglobal εSglobal εSglobal 

ticities. In Section 4.6, we showed that in a non-linear system of equations featuring a single equation for each location, a sufficient condition for uniqueness was that the exponent on the right hand side endogenous variable was less than or equal to one in magnitude. It turns out that a similar sufficient condition for uniqueness can be derived in this workhorse framework featuring two equations for each location. Let $| \mathbf { A } | \equiv { \binom { \left| a _ { 1 1 } \right| \left| a _ { 1 2 } \right| } { \left| a _ { 2 1 } \right| \left| a _ { 2 2 } \right| } }$ be the absolute value of the matrix of the model elasticities. If the spectral radius of $| \mathbf { A } |$ (i.e. if the largest eigenvalue in terms of magnitude) is less than or equal to one, i.e. $\rho \left( \left| \mathbf { A } \right| \right) \leq 1 ,$ then the equilibrium is unique. 

To see this, we proceed by contradiction. Suppose that $\rho \left( | \mathbf { A } | \right) \leq 1$ and that there are distinct solutions $\left\{ \lambda ^ { A } , \big \{ \bar { x _ { i } ^ { A } } , z _ { i } ^ { A } \big \} _ { i \in \mathcal { N } } \right\}$ and $\left\{ \lambda ^ { B } , \left\{ x _ { i } ^ { B } , z _ { i } ^ { B } \right\} _ { i \in \mathcal { N } } \right\}$ that both solve the system of equations (68) and (69). Let a hatted variable indicate the ratio of the two solutions, e.g. $\hat { x } _ { i } \equiv x _ { i } ^ { B } / x _ { i } ^ { A }$ , so that we can write equations (68) and (69) as follows: 

$$
\hat {\lambda} \hat {x} _ {i} = \sum_ {j \in \mathcal {N}} C _ {i j} \hat {x} _ {j} ^ {a _ {1 1}} \hat {z} _ {j} ^ {a _ {1 2}} \tag {73}
$$

$$
\hat {\lambda} \hat {z} _ {i} = \sum_ {j \in \mathcal {N}} D _ {i j} \hat {x} _ {j} ^ {a _ {2 1}} \hat {z} _ {j} ^ {a _ {2 2}}, \tag {74}
$$

where Cij $\begin{array} { r } { C _ { i j } \equiv \frac { K _ { i j } \left( x _ { j } ^ { A } \right) ^ { a _ { 1 1 } } \left( z _ { j } ^ { A } \right) ^ { a _ { 1 2 } } } { \sum _ { k \in \mathcal { N } } K _ { i k } \left( x _ { k } ^ { A } \right) ^ { a _ { 1 1 } } \left( z _ { k } ^ { A } \right) ^ { a _ { 1 2 } } } } \end{array}$ ≡ j j∑k∈N Kik (x Ak )a11 (zAk )a12 and Di j ≡ $\begin{array} { r } { D _ { i j } \equiv \frac { K _ { j i } \left( x _ { j } ^ { A } \right) ^ { a _ { 2 1 } } \left( z _ { j } ^ { A } \right) ^ { a _ { 2 2 } } } { \sum _ { k \in \mathcal { N } } K _ { k i } \left( x _ { j } ^ { A } \right) ^ { a _ { 2 1 } } \left( z _ { k } ^ { A } \right) ^ { a _ { 2 2 } } } } \end{array}$ 21 ( z Ak ) a 22 . Define $\begin{array} { r } { \mu _ { 1 } \equiv \frac { \operatorname* { m a x } _ { i \in \mathcal { N } } \hat { x } _ { i } } { \operatorname* { m i n } _ { i \in \mathcal { N } } \hat { x } _ { i } } } \end{array}$ and $\mu _ { 2 } \equiv$ $\begin{array} { r } { \frac { \operatorname* { m a x } _ { i \in \mathcal { N } } \hat { z } _ { i } } { \operatorname* { m i n } _ { i \in \mathcal { N } } \hat { z } _ { i } } } \end{array}$ . Because $\begin{array} { r } { \sum _ { j \in \mathcal { N } } C _ { i j } = \sum _ { j \in \mathcal { N } } D _ { i j } = 1 , } \end{array}$ , we can use equations (73) and (74) to construct the following bounds: 

$$
\mu_ {1} \leq \mu_ {1} ^ {\left| a _ {1 1} \right|} \mu_ {1} ^ {\left| a _ {1 2} \right|} \tag {75}
$$

$$
\mu_ {2} \leq \mu_ {1} ^ {| a _ {2 1} |} \mu_ {2} ^ {| a _ {2 2} |}, \tag {76}
$$

with the inequality strict for at least one of the two equations because the two solutions are 


Figure 9: Uniqueness of the equilibrium



(a) With small spatial linkages


![image](https://cdn-mineru.openxlab.org.cn/result/2026-02-25/4e4c25dd-e494-4502-8cd1-6a5c0f2e691d/adf14e2b18b3753f8dd037a608440d19bba3be60e0a8c7e891f9df9ae7af5451.jpg)



(b) With moderate spatial linkages


![image](https://cdn-mineru.openxlab.org.cn/result/2026-02-25/4e4c25dd-e494-4502-8cd1-6a5c0f2e691d/5c91d61369dc18bc6386cf61d5778f8bd53d518ae5cc1eb200db6b98c51570b8.jpg)



(c) With large spatial linkages


![image](https://cdn-mineru.openxlab.org.cn/result/2026-02-25/4e4c25dd-e494-4502-8cd1-6a5c0f2e691d/9e2460caee0fe23cbacf7a9dfd9f447f0a4d27eb0b6b046bd03bac2456e36cd7.jpg)


Notes: This figure displays the range of values for the local labor supply elasticity $\varepsilon _ { l o c a l } ^ { S }$ and demand elasticity $\varepsilon _ { l o c a l } ^ { D }$ for which a unique equilibrium is assured. Yellow indicates uniqueness. 

assumed to be distinct. In matrix notation, this becomes: 

$$
\binom {\ln \mu_ {1}} {\ln \mu_ {2}} \leq | \mathbf {A} | \binom {\ln \mu_ {1}} {\ln \mu_ {2}}, \tag {77}
$$

with at least one of the two inequalities strict. Given their definitions, $\ln \mu _ { 1 } \geq 0$ and $\ln \mu _ { 2 } \geq 0$ with at least one inequality strict again because the two solutions are distinct. According to the Collatz-Wielandt formula (i.e. that ρ (A) = maxln µ∈R2 ,ln µ̸=0, minh∈{1,2},ln µh̸=0 $\begin{array} { r } { \rho \left( \mathbf { A } \right) \ = \ \operatorname* { m a x } _ { \ln \mu \in \mathbb { R } _ { + } ^ { 2 } , \ln \mu \neq 0 } \operatorname* { m i n } _ { h \in \{ 1 , 2 \} , \ln \mu _ { h } \neq 0 } \frac { ( | \mathbf { A } | \ln \bar { \mu } ) _ { h } } { \ln \mu _ { h } } ) , } \end{array}$ (|A| ln µ)h ), equation (77) then implies that $\rho \left( \left| \mathbf { A } \right| \right) > 1 ,$ a contradiction. Hence there is a unique solution. 

Note that the argument above generalizes readily to non-linear systems of equations with any number of equations in each location, as shown in Allen, Arkolakis, and Li (2024). Allen, Arkolakis, and Li (2024) also show that a similar argument can be applied to non log-linear systems of equations by considering a matrix of the bounds of the elasticities of the system. Note too that this result establishes the claim made in Section 4.6 that when trade costs are symmetric, the unique inward and outward market access terms are those that are equal up to scale. 

What does the condition that $\rho \left( \left| \mathbf { A } \right| \right) \leq 1$ mean in economic terms? Panel (a) of Figure 9 are very small (e.guniqueness is that εSglobal $\varepsilon _ { g l o b a l } ^ { S }$ $\varepsilon _ { g l o b a l } ^ { D }$ are close to zero). In this case, a sufficient condition forhis condition holds if the supply curve is upward sloping $\varepsilon _ { l o c a l } ^ { S } + \varepsilon _ { l o c a l } ^ { D } \geq 0 .$ $( \varepsilon _ { l o c a l } ^ { S } > 0 )$ and the demand curve is downward sloping $( \varepsilon _ { l o c a l } ^ { D } > 0 )$ ) but also holds if the supply curve is upward sloping and the demand curve is not “too” upward sloping (or vice versa). As the spatial linkages get stronger (as illustrated in Panels (b) and (c)), the range of acceptable supply and demand elasticities for which uniqueness is guaranteed shrinks. 

It is important to emphasize that the condition $\rho \left( \left| \mathbf { A } \right| \right) ~ \leq ~ 1$ is a sufficient condition for 

uniqueness. Although Allen, Arkolakis, and Li (2024) show that there will always exist geographies for which there are multiple equilibria if the condition is not satisfied (i.e. the condition is “globally” necessary), given a particular geography, the equilibrium may still be unique even if $\rho \left( | \mathbf { A } | \right) > 1$ . It remains an outstanding question how in these cases the particular geography interacts with the model elasticities to determine the number of possible equilibria; exciting recent research on uniqueness and multiplicity in related spatial economic models includes Bifulco, Gluck, Krebs, and Kukharskyy (2022), Ouazad (2024), and Garg (2025). ¨ 

# 6.6 Counterfactual analysis

We now turn to the question of how to use this quantitative economic geography framework to conduct counterfactual analyses. 

In what follows, we will assume that the geography of the world was initially given by the $N \times N$ matrix $\mathbf { K } ^ { A } = \left[ K _ { i j } ^ { A } \right]$ and we are interested in how the equilibrium distribution of economic activity would change if the geography changed to a new $N \times N$ matrix $\mathbf { K } ^ { B } = \left[ K _ { i j } ^ { B } \right]$ . 

We proceed using the “exact hat” algebra introduced by Dekle, Eaton, and Kortum (2008) for its use in quantitative trade models and subsequently popularized by Costinot and Rodr´ıguez-Clare (2014) (although similar techniques were present in Karlin and Nirenberg (1967)). We assume that the initial economy with geography $\mathbf { K } ^ { A }$ is in equilibrium $A$ with equilibrium distribution of economic activity $\big \{ l _ { i } ^ { A } , y _ { i } ^ { A } \big \} _ { i \in \mathcal { N } }$ and subsequently with geography $\mathbf { K } ^ { B }$ would be in equilibrium B with equilibrium distribution of economic activity $\big \{ l _ { i } ^ { B } , y _ { i } ^ { B } \big \} _ { i \in \mathcal { N } }$ . As in Section 6.5, we let a hatted variable indicate the ratio of the two equilibria, e.g. $\hat { x } _ { i } \equiv x _ { i } ^ { B } / x _ { i } ^ { A }$ . We can then re-write the equilibrium system of equations (68) and (69) in changes as follows: 

$$
\begin{array}{l} \hat {\lambda} \hat {y} _ {i} ^ \frac {1}{\varepsilon_ {g l o b a l} ^ {D}} \hat {l} _ {i} ^ {\frac {\varepsilon_ {l o c a l} ^ {D} - 1}{\varepsilon_ {g l o b a l} ^ {D}}} = \sum_ {j \in \mathcal {N}} C _ {i j} \hat {K} _ {i j} \hat {y} _ {j} ^ {\frac {\varepsilon_ {g l o b a l} ^ {S} + 1}{\varepsilon_ {g l o b a l} ^ {S}}} \hat {l} _ {j} ^ {- \frac {1 + \varepsilon_ {l o c a l} ^ {S}}{\varepsilon_ {g l o b a l} ^ {S}}} (78) \\ \hat {\lambda} \hat {y} _ {i} ^ {- \frac {1}{\varepsilon_ {g l o b a l} ^ {S}} \frac {\varepsilon_ {l o c a l} ^ {S} + 1}{\varepsilon_ {g l o b a l} ^ {S}}} \hat {l} _ {i} ^ {- \frac {\varepsilon_ {g l o b a l} ^ {D} - 1}{\varepsilon_ {g l o b a l} ^ {D}}} = \sum_ {j \in \mathcal {N}} D _ {j i} \hat {K} _ {j i} \hat {y} _ {j} ^ {- \frac {\varepsilon_ {g l o b a l} ^ {D} - 1}{\varepsilon_ {g l o b a l} ^ {D}}} \hat {l} _ {j} ^ {- \frac {1 - \varepsilon_ {l o c a l} ^ {D}}{\varepsilon_ {g l o b a l} ^ {D}}}, (79) \\ \end{array}
$$

$$
\text {w h e r e w (n o w) d e f i n e} C _ {i j} \equiv \frac {K _ {i j} ^ {A} \left(y _ {j} ^ {A}\right) ^ {\frac {\varepsilon^ {S} _ {g l o b a l} + 1}{\varepsilon^ {S} _ {g l o b a l}}} \left(l _ {j} ^ {A}\right) ^ {- \frac {1 + \varepsilon^ {S} _ {l o c a l}}{\varepsilon^ {S} _ {g l o b a l}}}}{\sum_ {k \in \mathcal {N}} K _ {i k} ^ {A} \left(y _ {k} ^ {A}\right) ^ {\frac {\varepsilon^ {S} _ {g l o b a l} + 1}{\varepsilon^ {S} _ {g l o b a l}}} \left(l _ {k} ^ {A}\right) ^ {- \frac {1 + \varepsilon^ {S} _ {l o c a l}}{\varepsilon^ {S} _ {g l o b a l}}}} \text {a n d} D _ {j i} \equiv \frac {K _ {j i} ^ {A} \left(y _ {j} ^ {A}\right) ^ {\frac {\varepsilon^ {D} _ {g l o b a l} - 1}{\varepsilon^ {D} _ {g l o b a l}}} \left(l _ {j} ^ {A}\right) ^ {\frac {1 - \varepsilon^ {D} _ {l o c a l}}{\varepsilon^ {D} _ {g l o b a l}}}}{\sum_ {k \in \mathcal {N}} K _ {k i} ^ {A} \left(y _ {k} ^ {A}\right) ^ {\frac {\varepsilon^ {D} _ {g l o b a l} - 1}{\varepsilon^ {D} _ {g l o b a l}}} \left(l _ {k} ^ {A}\right) ^ {\frac {1 - \varepsilon^ {D} _ {l o c a l}}{\varepsilon^ {D} _ {g l o b a l}}}}.
$$

The key insight of Dekle, Eaton, and Kortum (2008) (and certainly not present in Karlin and Nirenberg (1967)) is that $C _ { i j }$ and $D _ { j i }$ are observable in the data. Indeed, it is straightforward to show that $\begin{array} { r } { C _ { i j } = \frac { X _ { i j } ^ { A } } { Y _ { i } ^ { A } } } \end{array}$ are the initial export shares of $i$ to $j$ and and $\begin{array} { r } { D _ { j i } = \frac { X _ { j i } ^ { A } } { E _ { i } ^ { A } } } \end{array}$ E Ai are the initial import 

shares of i from $j ,$ allowing us to write equations (78) and (79) as follows: 

$$
\hat {\lambda} \hat {y} _ {i} ^ \frac {1}{\varepsilon_ {g l o b a l} ^ {D}} \hat {l} _ {i} ^ {\frac {\varepsilon_ {l o c a l} ^ {D} - 1}{\varepsilon_ {g l o b a l} ^ {D}}} = \sum_ {j \in \mathcal {N}} \left(\frac {X _ {i j} ^ {A}}{Y _ {i} ^ {A}}\right) \hat {K} _ {i j} \hat {y} _ {j} ^ {\frac {\varepsilon_ {g l o b a l} ^ {S} + 1}{\varepsilon_ {g l o b a l} ^ {S}} \hat {l} _ {j} ^ {- \frac {1 + \varepsilon_ {l o c a l} ^ {S}}{\varepsilon_ {g l o b a l} ^ {S}}}} \tag {80}
$$

$$
\hat {\lambda} \hat {y} _ {i} ^ {- \frac {1}{\varepsilon_ {g l o b a l} ^ {S}} \frac {\varepsilon_ {l o c a l} ^ {S} + 1}{\varepsilon_ {g l o b a l} ^ {S}}} \hat {l} _ {i} ^ {\frac {\varepsilon_ {l o c a l} ^ {S} + 1}{\varepsilon_ {l o c a l} ^ {S}}} = \sum_ {j \in \mathcal {N}} \left(\frac {X _ {j i} ^ {A}}{E _ {i} ^ {A}}\right) \hat {K} _ {j i} \hat {y} _ {j} ^ {\frac {\varepsilon_ {g l o b a l} ^ {D} - 1}{\varepsilon_ {g l o b a l} ^ {D}} \frac {1 - \varepsilon_ {l o c a l} ^ {D}}{\varepsilon_ {g l o b a l} ^ {D}}}. \tag {81}
$$

There are three important things to note about equations (80) and (81). First, they can be used to calculate how any potential change in geography $\big \{ \hat { K } _ { i j } \big \} _ { i , j \in \mathcal { N } }$ affects the entire distribution of economic activity while requiring only knowledge of the initial equilibrium trade flows and model elasticities. Second, because the counterfactual equations (80) and (81) share an identical mathematical structure to the equilibrium equations (68) and (69), they inherent all the mathematical properties of that system. In particular, if the sufficient conditions for uniqueness are satisfied for the initial equilibrium, then any counterfactual being considered will also be unique. Third, we can also use equations (80) and (81) to perform normative analysis by examining how the endogenous scalar changes from equilibrium $A$ to equilibrium $B$ , i.e. $\hat { \lambda }$ . To solve for $\hat { \lambda } ,$ we use the fact that output and population shares sum to one to derive the constraints $\begin{array} { r } { \sum _ { i \in \mathcal { N } } y _ { i } ^ { A } \hat { y } _ { i } = \sum _ { i \in \mathcal { N } } l _ { i } ^ { A } \hat { l } _ { i } = 1 . } \end{array}$ 

In the next (technical) subsection, we turn to understanding the comparative statics of the quantitative economic geography framework. In particular, we highlight how counterfactuals can be interpreted as shocks propagating through the trading network with the observed trade data reflecting the strength of network ties. 

# 6.7 [Technical] Comparative statics

Section 6.6 showed how to calculate the change in the equilibrium distribution of economic activity for any change in the geography of the world. By focusing on small changes in geography, we can derive new insights into the economic mechanisms at work. 

To begin, we follow Allen, Arkolakis, and Takahashi (2020), Kleinman, Liu, and Redding (2024a), Adao, Arkolakis, and Esposito, 2019, and Kleinman, Liu, and Redding (2024b) by log differentiating the equilibrium system of equations (68) and (69) for some small change 

$\{ d \ln K _ { i j } \} _ { i , j \in \mathcal { N } }$ to yield: 

$$
\begin{array}{l} d \ln \lambda + \frac {1}{\varepsilon_ {\text {g l o b a l}} ^ {D}} d \ln y _ {i} - \frac {1 - \varepsilon_ {\text {l o c a l}} ^ {D}}{\varepsilon_ {\text {g l o b a l}} ^ {D}} d \ln l _ {i} = \sum_ {j \in \mathcal {N}} \left(\frac {X _ {i j} ^ {A}}{Y _ {i} ^ {A}}\right) \left(d \ln K _ {i j} + \frac {\varepsilon_ {\text {g l o b a l}} ^ {S} + 1}{\varepsilon_ {\text {g l o b a l}} ^ {S}} d \ln y _ {j} - \frac {1 + \varepsilon_ {\text {l o c a l}} ^ {S}}{\varepsilon_ {\text {g l o b a l}} ^ {S}} d \ln l _ {j}\right) \tag {82} \\ d \ln \lambda - \frac {1}{\varepsilon_ {g l o b a l} ^ {S}} d \ln y _ {i} + \frac {1 + \varepsilon_ {l o c a l} ^ {S}}{\varepsilon_ {g l o b a l} ^ {S}} d \ln l _ {i} = \sum_ {j \in \mathcal {N}} \left(\frac {X _ {j i} ^ {A}}{E _ {i} ^ {A}}\right) \left(d \ln K _ {j i} + \frac {\varepsilon_ {g l o b a l} ^ {D} - 1}{\varepsilon_ {g l o b a l} ^ {D}} d \ln y _ {j} + \frac {1 - \varepsilon_ {l o c a l} ^ {D}}{\varepsilon_ {g l o b a l} ^ {D}} d \ln l _ {j}\right), \\ \end{array}
$$

where the fact that output and population shares sum to one additionally implies $\begin{array} { r l } { \sum _ { i \in \mathcal { N } } y _ { i } ^ { A } d \ln y _ { i } = } \end{array}$ 0 and $\Sigma _ { i \in \mathcal { N } } l _ { i } ^ { A } d \ln l _ { i } = 0$ . We can re-write equations (82) and (83) in matrix notation as follows: 

$$
\begin{array}{l} d \ln \lambda + \frac {1}{\varepsilon_ {\text {g l o b a l}} ^ {D}} d \ln \mathbf {y} - \frac {1 - \varepsilon_ {\text {l o c a l}} ^ {D}}{\varepsilon_ {\text {g l o b a l}} ^ {D}} d \ln \boldsymbol {l} = \operatorname {d i a g} \left(\mathbf {T} d \ln \mathbf {K} ^ {T}\right) + \mathbf {T} \left(\frac {\varepsilon_ {\text {g l o b a l}} ^ {S} + 1}{\varepsilon_ {\text {g l o b a l}} ^ {S}} d \ln \mathbf {y} - \frac {1 + \varepsilon_ {\text {l o c a l}} ^ {S}}{\varepsilon_ {\text {g l o b a l}} ^ {S}} d \ln \boldsymbol {l}\right) (84) \\ d \ln \lambda - \frac {1}{\varepsilon_ {\text {g l o b a l}} ^ {S}} d \ln y + \frac {1 + \varepsilon_ {\text {l o c a l}} ^ {S}}{\varepsilon_ {\text {g l o b a l}} ^ {S}} d \ln l = \operatorname {d i a g} \left(\mathbf {S} ^ {T} d \ln \mathbf {K}\right) + \mathbf {S} ^ {T} \left(\frac {\varepsilon_ {\text {g l o b a l}} ^ {D} - 1}{\varepsilon_ {\text {g l o b a l}} ^ {D}} d \ln y - \frac {1 - \varepsilon_ {\text {l o c a l}} ^ {D}}{\varepsilon_ {\text {g l o b a l}} ^ {D}} d \ln l\right), (85) \\ \end{array}
$$

with adding up constraints $\begin{array} { r } { { \pmb y } ^ { T } d \ln { \pmb y } = 0 } \end{array}$ and $l ^ { T } d \ln l = 0$ , where $d \ln { \boldsymbol { y } } = \left[ { \boldsymbol { y } } _ { i } \right]$ and $d \ln l = [ l _ { i } ]$ are $N \times 1$ vectors of the respective changes in output and population shares, $d \ln { \bf K } = \left[ d \ln K _ { i j } \right]$ is an $N \times N$ matrix of the changes in geography, $\mathbf { T } = \left[ X _ { i j } ^ { A } / Y _ { i } ^ { A } \right]$ is an $N \times N$ matrix of export shares in the initial equilibrium, and $\mathbf { S } \equiv \left[ X _ { i j } ^ { A } / E _ { j } ^ { A } \right]$ is an $N \times N$ matrix of import shares in the initial equilibrium. 

To proceed, let us consider a change of variables equivalent to the one used when proving uniqueness in Section 6.5. Define: 

$$
d \ln x \equiv \frac {1}{\varepsilon_ {\text {g l o b a l}} ^ {D}} d \ln y - \frac {1 - \varepsilon_ {\text {l o c a l}} ^ {D}}{\varepsilon_ {\text {g l o b a l}} ^ {D}} d \ln l \tag {86}
$$

$$
d \ln z \equiv - \frac {1}{\varepsilon_ {\text {g l o b a l}} ^ {S}} d \ln y + \frac {1 + \varepsilon_ {\text {l o c a l}} ^ {S}}{\varepsilon_ {\text {g l o b a l}} ^ {S}} d \ln l \tag {87}
$$

We can then re-write equations (84) and (85) as the following $2 N$ system of equations: 

$$
\left( \begin{array}{l} d \ln x \\ d \ln z \end{array} \right) + d \ln \lambda = \left( \begin{array}{l} \operatorname {d i a g} \left(\mathbf {T} d \ln \mathbf {K} ^ {T}\right) \\ \operatorname {d i a g} \left(\mathbf {S} ^ {T} d \ln \mathbf {K}\right) \end{array} \right) + \left( \begin{array}{l l} a _ {1 1} \mathbf {T} & a _ {1 2} \mathbf {T} \\ a _ {2 1} \mathbf {S} ^ {T} & a _ {2 2} \mathbf {S} ^ {T} \end{array} \right) \left( \begin{array}{l} d \ln x \\ d \ln z \end{array} \right), \tag {88}
$$

where the $2 \times 2$ matrix A was defined in Section 6.5. Re-writing the adding up constraints in terms of the new variables yields: 

$$
\boldsymbol {y} ^ {T} \left(b _ {1 1} d \ln x + b _ {1 2} d \ln z\right) = 0 \tag {89}
$$

$$
\boldsymbol {l} ^ {T} \left(b _ {2 1} d \ln x + b _ {2 2} d \ln z\right) = 0, \tag {90}
$$

where $\begin{array} { r } { \mathbf { B } = \left( \begin{array} { l l } { b _ { 1 1 } } & { b _ { 1 2 } } \\ { b _ { 2 1 } } & { b _ { 2 2 } } \end{array} \right) \equiv \left( \begin{array} { l l } { \frac { 1 } { \varepsilon _ { g l o b a l } ^ { D } } } & { \frac { \varepsilon _ { l o c a l } ^ { D } - 1 } { \varepsilon _ { g l o b a l } ^ { D } } } \\ { - \frac { 1 } { \varepsilon _ { g l o b a l } ^ { S } } } & { \frac { \varepsilon _ { l o c a l } ^ { S } + 1 } { \varepsilon _ { g l o b a l } ^ { S } } } \end{array} \right) ^ { - } } \end{array}$ εDglobalεSlocal +1 1 . Applying either constraint to equation (88) b21 εSglobal 

allows us to solve for the change in the endogenous scalar as follows: 

$$
d \ln \lambda = \frac {\boldsymbol {y} ^ {T} \left(b _ {1 1} \operatorname {d i a g} \left(\mathbf {T} d \ln \mathbf {K} ^ {T}\right) + b _ {1 1} \mathbf {T} \left(a _ {1 1} d \ln x + a _ {1 2} d \ln z\right) + b _ {1 2} \operatorname {d i a g} \left(\mathbf {S} ^ {T} d \ln \mathbf {K}\right) + b _ {1 2} \mathbf {S} ^ {T} \left(a _ {2 1} d \ln x + a _ {2 2} d \ln z\right)\right)}{b _ {1 1} + b _ {1 2}} \tag {91}
$$

$$
d \ln \lambda = \frac {\boldsymbol {l} ^ {T} \left(b _ {2 1} \operatorname {d i a g} \left(\mathbf {T} d \ln \mathbf {K} ^ {T}\right) + b _ {2 1} \mathbf {T} \left(a _ {1 1} d \ln x + a _ {1 2} d \ln z\right) + b _ {2 2} \operatorname {d i a g} \left(\mathbf {S} ^ {T} d \ln \mathbf {K}\right) + b _ {2 2} \mathbf {S} ^ {T} \left(a _ {2 1} d \ln x + a _ {2 2} d \ln z\right)\right)}{b _ {2 1} + b _ {2 2}}, \tag {92}
$$

which in turn results (finally) in the following linear system of equations which summarizes how the distribution of economic activity responds to any change in the underlying geography: 

$$
\left( \begin{array}{l} d \ln x \\ d \ln z \end{array} \right) = \left(\mathbf {I} - (\mathbf {I} - \mathbf {Q}) \left( \begin{array}{l l} a _ {1 1} \mathbf {T} & a _ {1 2} \mathbf {T} \\ a _ {2 1} \mathbf {S} ^ {T} & a _ {2 2} \mathbf {S} ^ {T} \end{array} \right)\right) ^ {- 1} (\mathbf {I} - \mathbf {Q}) \left( \begin{array}{l} \operatorname {d i a g} \left(\mathbf {T} d \ln \mathbf {K} ^ {T}\right) \\ \operatorname {d i a g} \left(\mathbf {S} ^ {T} d \ln \mathbf {K}\right) \end{array} \right), \tag {93}
$$

where $\mathbf { Q } \equiv \left( \begin{array} { c c } { \frac { b _ { 1 1 } } { b _ { 1 1 } + b _ { 1 2 } } \mathbf { y } ^ { T } } & { \frac { b _ { 1 2 } } { b _ { 1 1 } + b _ { 1 2 } } \mathbf { y } ^ { T } } \\ { \frac { b _ { 2 1 } } { b _ { 2 1 } + b _ { 2 2 } } \mathbf { l } ^ { T } } & { \frac { b _ { 2 2 } } { b _ { 2 1 } + b _ { 2 2 } } \mathbf { l } ^ { T } } \end{array} \right)$ T and $\mathbf { y } ^ { T }$ and $1 ^ { T }$ are $N \times N$ matrices with the output and b21+b22 population shares $1 \times N$ vectors, respectively, stacked across $N$ rows. The Q matrix ensures the underlying output and population shares sum to one. Given $d \ln x$ and $d \ln z ,$ , the distribution of economic activity can be immediately recovered by inverting equations (86) and (87). 

There are three things to note about equation (93). First, with the system of equations equations (80) and (81) governing how arbitrary changes in geography affect the equilibrium spatial distribution of economic activity, equation (93) shows that for any small change in the geography, the resulting change in the distribution of economic activity can be calculating knowing only the model elasticities, the observed trade flows (which allow one to calculate the income shares T and expenditure shares S), and the observed distribution of economic activity $\{ y _ { i } , l _ { i } \} _ { i \in \mathcal { N } }$ . 

Second, equation (93) provides some intuition about where the sufficient condition for uniqueness $\rho \left( \left| \mathbf { A } \right| \right) \leq 1$ comes from. Note that if $\rho \left( | \mathbf { A } | \right) \leq 1 ,$ , we can follow Allen, Arkolakis, and Takahashi (2020) and Kleinman, Liu, and Redding (2024a) and consider a Neumann series ex-

pansion of equation (93): 

$$
\binom {d \ln x} {d \ln z} = \sum_ {k \geq 0} \left((\mathbf {I} - \mathbf {Q}) \left( \begin{array}{c c} a _ {1 1} \mathbf {T} & a _ {1 2} \mathbf {T} \\ a _ {2 1} \mathbf {S} ^ {T} & a _ {2 2} \mathbf {S} ^ {T} \end{array} \right)\right) ^ {k} (\mathbf {I} - \mathbf {Q}) \binom {\operatorname {d i a g} \left(\mathbf {T} d \ln \mathbf {K} ^ {T}\right)} {\operatorname {d i a g} \left(\mathbf {S} ^ {T} d \ln \mathbf {K}\right)}. \tag {94}
$$

We can consider each element of this infinite sum as the degree of its effect. For example, the zeroth degree (direct) effect (corresponding to $k = 0$ ) is $\begin{array} { r } { \mathbf { \Gamma } \left( \mathbf { I } - \mathbf { Q } \right) \left( \underset { \mathrm { d i a g } } { \mathrm { d i a g } } \left( \mathbf { T } d \ln \mathbf { K } ^ { T } \right) \right) } \end{array}$ ,i.e. how a change in geography affects each location directly through export shares and import shares, holding constant the economic activity in all other locations. The first degree effect (corresponding to $k = 1$ ) then captures how the direct effect of the geography shock indirectly affects a location through its trading partners via the matrix $\left( \mathbf { I } - \mathbf { Q } \right) \left( \begin{array} { c c } { { { } } } & { { { } } } \\ { { { a _ { 1 1 } } { \mathbf { \bar { T } } } } } & { { { a _ { 1 2 } } { \mathbf { T } } } } \\ { { { a _ { 2 1 } } { \mathbf { S } } ^ { T } } } & { { { a _ { 2 2 } } { \mathbf { S } } ^ { T } } } \end{array} \right) ,$ a22ST . And in general, the $k ^ { t h }$ degree effect captures how a location is affected by the $k - 1 ^ { t h }$ degree effect on its trading partners. As long as $\rho \left( \left| \mathbf { A } \right| \right) \leq 1 ,$ the indirect effects will get smaller and smaller as we consider higher and higher degree effects. But if $\rho \left( \left| \mathbf { A } \right| \right) > 1 ,$ , then it is possible that small shocks to the geography can become increasingly amplified as they percolate throughout the economy. As we mentioned in Section (6.5), an unresolved question is how the underlying geography affects whether or not there are multiple equilibria when $\rho \left( \left| \mathbf { A } \right| \right) > 1$ . It is possible that further analysis of the properties of equation (93) will help shed light on this question. 

The third thing to note about equation (93) is that, when combined with equations (91) and/or (92), it allows us to efficiently calculate the first order welfare changes from any small changes in the underlying geography. This can assist in the search for globally optimal spatial policies, which is an active and exciting area of inquiry, see e.g. the recent work of Fajgelbaum and Gaubert (2020) and Donald, Fukui, and Miyauchi (2023). We refer the interested reader to Fajgelbaum and Gaubert (2025) for more details on the efficiency properties of spatial models. 

# 7 Bringing the model to the data

One of the major advantages of the quantitative economic geography framework is its tight link with empirical data, a point we now discuss in detail. 

Suppose you are interested in assessing how a particular policy changes the geography of the world and, in particular, how it affects the spatial distribution of economic activity and the equilibrium welfare of agents. This policy may have already occurred, in which case your analysis would be conducted ex post, or it may be a policy that may occur in the future, in which case your analysis would be conducted ex ante. In either case, the quantitative economic geography framework is a powerful tool for such an analysis. 

Recall from Section 6.6 that the quantitative economic geography framework can be used 

to assess how any change to the underlying geography affects the equilibrium spatial distribution of economic activity and the associated welfare using equations (80) and (81). These equations also highlight the three necessary inputs for such analysis: (1) observed spatial data, including trade flows $\big \{ X _ { i j } \big \} _ { i , j \in \mathcal { N } ^ { \prime } }$ income $\{ Y _ { i } \} _ { i \in \mathcal { N } } ,$ expenditure $\{ E _ { i } \} _ { i \in \mathcal { N } } ,$ and the distribution of population derlying geog $\{ L _ { i } \} _ { i \in { \mathcal { N } } }$ wledge of how the policy or eve.; and (3) the model parameters un-. As $\big \{ \hat { K } _ { i j } \big \} _ { i , j \in \mathcal { N } }$ $\left\{ \varepsilon _ { l o c a l } ^ { S } , \varepsilon _ { g l o b a l } ^ { S } , \varepsilon _ { l o c a l } ^ { D } , \varepsilon _ { g l o b a l } ^ { D } \right\}$ highlighted in Section 2, these spatial data are readily available for many countries, so in what follows, we will focus on how to determine how a policy affects the underlying geography (Fact #2) and the estimation of the model parameters (Fact #3). 

# 7.1 Estimating the effect of spatial policies on the geography

Recall from Section 6 that Kij ≡ Tij  CDi 1/εDglobal  $K _ { i j } \equiv T _ { i j } \left( C _ { i } ^ { D } \right) ^ { 1 / \varepsilon } g _ { g l o b a l } ^ { D } \left( C _ { j } ^ { S } \right) ^ { 1 / \varepsilon } g _ { g l o b a l } ^ { S }$ C Sj , i.e. there are three components of geography that a potential policy may affect: the bilateral trade frictions $\big \{ T _ { i j } \big \} _ { i , j \in \mathcal { N } ^ { \prime } }$ an exogenous factor that affects the demand for labor in a location $\{ C _ { i } ^ { D } \} _ { i \in \mathcal { N } } ,$ and an exogenous factor that affects the supply for a labor in a location $\left\{ C _ { j } ^ { S } \right\} _ { j \in \mathcal { N } }$ . There are many examples of policies that affect bilateral trade frictions, e.g. infrastructure improvements to highways (as in Allen and Arkolakis (2014) and Jaworski and Kitchens (2019)), railroads (as in Donaldson (2018) and Donaldson and Hornbeck (2016)), or subways (as in Severen (2023)); we refer the interested reader to Donaldson (2025) for a much more thorough analysis. Here, we begin by positing that bilateral trade frictions are some function $f$ —parameterized by a vector $\beta$ —of a vector of observables $\mathbf { z } ,$ i.e. $T _ { i j } = f _ { i j } ( \mathbf { z } ; \beta ) ,$ , so that we can write the gravity equation (61) as follows: 

$$
X _ {i j} = f _ {i j} (\mathbf {z}; \boldsymbol {\beta}) \times \frac {Y _ {i}}{M A _ {i} ^ {o u t}} \times \frac {E _ {j}}{M A _ {j} ^ {i n}}. \tag {95}
$$

Following Silva and Tenreyro (2006), equation (95) can be estimated by assuming a particular functional for $f \left( \cdot \right)$ (oftentimes log-linear) and then using a Poisson psuedo-maximumlikelihood (PPML) estimator to estimate: 

$$
X _ {i j} = \exp \left(\log f _ {i j} (\mathbf {z}; \boldsymbol {\beta}) + \gamma_ {i} + \delta_ {j}\right) \epsilon_ {i j}, \tag {96}
$$

where the origin fixed effect $\begin{array} { r } { \gamma _ { i } \equiv \ln { \frac { Y _ { i } } { M A _ { i } ^ { o u t } } } } \end{array}$ and destination fixed effect $\begin{array} { r } { \delta _ { j } \equiv \ln { \frac { E _ { j } } { M A _ { j } ^ { i n } } } } \end{array}$ and $\epsilon _ { i j }$ is an error term. Because $Y _ { i }$ and $E _ { i }$ are observed, an advantage of this approach is that we recover the inward and outward market access terms from the estimated origin and destination fixed effects. Moreover, as Fally (2015) shows, the estimated inward and outward market access terms from the PPML estimation will satisfy the system of equations (64) and (65) for the observed incomes and expenditures and estimated bilateral trade frictions. As a re-

sult, this procedure yields estimated trade frictions $\tilde { T } _ { i j } \equiv f _ { i j } \left( \mathbf { z } ; \tilde { \boldsymbol { \beta } } \right)$ and market access measures $\left\{ \widetilde { M A _ { i } ^ { o u t } } , \widetilde { M A _ { i } ^ { i n } } \right\} _ { i \in \mathcal { N } } ,$ where we use a tilde to indicate a variable has been estimated. 

One limitation to recovering the market access directly from the regression of equation (96) is that the estimated origin and destination fixed effects cannot be separately identified from any origin- or destination-specific components of the trade frictions function. This is a symptom of a broader issue that the supply and demand shifters $\left\{ C _ { i } ^ { S } , C _ { i } ^ { D } \right\}$ cannot be separately identified from the trade frictions $\left\{ T _ { i j } \right\}$ given observed trade data $\left\{ X _ { i j } \right\}$ . Fortunately, as we will discuss below, this does not affect the ability of the framework to perform counterfactuals, so typically, researchers impose a normalization that there are no own trade frictions, i.e. $T _ { i i } = 1$ . 

An alternative way of recovering the market access given estimated trade frictions is to directly rely on the system of equations (64) and (65). Given income and expenditure, the results from Section 6.5 imply that for any set of trade frictions, there exists unique (to-scale) inward and outward market accesses consistent with equations (64) and (65) (as it is straightforward to show that the spectral radius of the absolute value of the elasticities of the system is one). This alternative procedure is helpful in situations where bilateral trade flows are not observed for all $i , j \in \mathcal { N }$ pairs. For example, this may occur when locations are U.S. counties, but trade frictions are estimated based on state-to-state trade flows. 

# 7.1.1 Algorithms for calculating changes in travel times

Once the function $f \left( \mathbf { z } ; { \tilde { \boldsymbol { \beta } } } \right)$ is estimated, it is in principle straightforward to determine how a policy affects trade frictions: simply determine how the observables z change in response to the policy and calculate the resulting change in the trade frictions. In practice, however, this process can be a bit more involved. A common assumption is that trade frictions are a (log linear) function of the travel time (or cost) along the fastest (or least costly) route from an origin to a destination. Under this assumption, z is a high dimensional set of observables including the time (or cost) of travel on every part of the entire infrastructure network. A policy that, say, builds a new link on the infrastructure network could affect the trade frictions between many origins and destinations simultaneously. While this seems like it would be quite difficult to calculate, there fortunately exist efficient algorithms for determining the fastest travel time (or lowest costs). Loosely speaking, the Dijkstra algorithm (Dijkstra, 1959)—and its continuous space generalization, the Fast Marching Method (Tsitsiklis, 1995; Sethian, 1996)—search outward from a given origin, keeping track of the least cost route so far, so that a single iteration of the algorithm calculates the least cost route to all possible destinations. Such procedures can be employed to calculate how large scale infrastructure improvements affect trade frictions between all origins and destinations. For example, Donaldson (2018) applies the Dijkstra algorithm to study the impact of the Indian railroad expansion, and Allen and Arkolakis (2014) apply the Fast Marching Method to study the impact of the U.S. interstate highway system. 

# 7.1.2 An analytical solution for changes in travel times

One limitation with the Dijkstra and Fast Marching Method algorithms is that their black box nature limits the ability to analytically assess how improvements to the transportation infrastructure network affect the spatial distribution of economic activity. An alternative method proposed by Allen and Arkolakis (2022) explicitly derives the the least cost route between any two locations along the transportation network. The basic idea is as follows. Suppose that there is a unit mass of perfectly competitive transportation firms $\nu \in [ 0 , 1 ]$ trying to ship goods from origin $i \in \mathcal N$ to destination $j \in \mathcal N$ over a transportation network summarized by the $N \times N$ matrix $\mathbf { t } \equiv [ t _ { k l } ]$ , where $t _ { k l }$ is the (iceberg) trade cost incurred from moving a good directly along a link between $k$ and l. Each transportation firm $\nu$ will then choose a route—i.e. a sequence of links beginning in $i$ and ending in $j .$ —in order to minimize the total transportation costs incurred from $i$ to $j \mathrm { : }$ : 

$$
r _ {i j} (\nu) = \min  _ {r \in \mathcal {R} _ {i j} ^ {K}, K \geq 1} \left(\prod_ {k = 1} ^ {K} t _ {k, k + 1}\right) \times \varepsilon_ {i j, r} (\nu), \tag {97}
$$

where $\mathcal { R } _ { i j } ^ { K }$ is the set of all possible routes from $i$ to $j$ of length $K$ and $\varepsilon _ { i j , r } \left( \nu \right)$ is an idiosyncratic preference of route $r$ for the transportation firm $\nu$ , which we assume is Frechet distributed with shape parameter $\theta$ . Note that as $\theta  \infty$ , all transportation firms will choose the lowest cost route, so the solution to equation (97) will coincide exactly with the output of the Dijkstra algorithm. 

Let us define $\tau _ { i j } \equiv \mathbb { E } \left[ r _ { i j } \left( \nu \right) \right]$ as the expected total iceberg trade cost from $i$ to $j$ across all transportation firms. The Frechet distributional assumption on the idiosyncratic term allows us to derive an analytical expression for $\tau _ { i j }$ by explicitly enumerating all possible routes from i to $j$ : 

$$
\tau_ {i j} ^ {- \theta} = \left(\sum_ {K = 0} ^ {\infty} \sum_ {k _ {1} = 1} ^ {N} \sum_ {k _ {2} = 1} ^ {N} \dots \sum_ {k _ {K - 1} = 1} ^ {N} \sum_ {k _ {K} = 1} ^ {N} \left(t _ {i, k _ {1}} \times t _ {k _ {1}, k _ {2}} \times \dots \times t _ {k _ {K - 1}, k _ {K}} \times t _ {k _ {K}, j}\right) ^ {- \theta}\right), \tag {98}
$$

where we note that pairs of locations that are not directly connected in the transportation network will have infinitely high direct iceberg costs, thereby not affecting the sum in equation (98). We can write the expression for the expected transportation cost more succinctly in matrix notation as a geometric sum of the $N \times N$ adjacency matrix $\mathbf { A } \equiv \left[ t _ { i j } ^ { - \theta } \right]$ : 

$$
\tau_ {i j} ^ {- \theta} = \sum_ {K = 0} ^ {\infty} \mathbf {A} _ {i j} ^ {K}. \tag {99}
$$

As long as the spectral radius of A is less than one, this geometric sum will be equal to its 

Leontief inverse $\mathbf { B } = \left[ b _ { i j } \right] \equiv \left( \mathbf { I } - \mathbf { A } \right) ^ { - 1 } .$ , so that: 

$$
\tau_ {i j} ^ {- \theta} = b _ {i j}. \tag {100}
$$

Finally, if we assume that the shape parameter governing the idiosyncratic component of the transportation firms is equal to the trade elasticity (justified, for example, by the setup presented in Allen and Arkolakis (2022) where the producers receive productivity draws over productroute pairs), then the trade friction $T _ { i j }$ from equation (61) is equal to the Leontief inverse, i.e. $T _ { i j } = b _ { i j }$ . This, in turn, allows us to write equilibrium equations (64) and (65) for market access directly as functions of the infrastructure network directly as follows: 

$$
M A _ {i} ^ {\text {o u t}} - \sum_ {j \in \mathcal {N}} t _ {i j} ^ {- \theta} M A _ {j} ^ {\text {o u t}} = \frac {E _ {i}}{M A _ {i} ^ {\text {i n}}} \tag {101}
$$

$$
M A _ {i} ^ {i n} - \sum_ {j \in \mathcal {N}} t _ {j i} ^ {- \theta} M A _ {i} ^ {i n} = \frac {Y _ {i}}{M A _ {i} ^ {\text {o u t}}}. \tag {102}
$$

Equations (101) and (102) highlight that market access of a location can be written recursively as only dependent on the market accesses of locations to which it is directly linked. When combined with the labor supply and demand equations (62) and (63), one can then directly calculate how a spatial policy which changes the transport costs along certain segments of the transportation network affects the equilibrium market access in all locations without the need to calculate the change in bilateral trade frictions. Moreover, despite the difference in its mathematical structure, Allen and Arkolakis (2022) shows that one can use similar mathematical arguments to those sketched above to establish conditions for the existence and uniqueness of an equilibrium where the transportation network (rather than the bilateral trade frictions) comprise the underlying geography, even in the presence of traffic congestion externalities (in which case the equilibrium distribution of economic activity will depend on the aggregate labor endowment). It also turns out that incorporating optimal routing into a quantitative economic geography framework is quite flexible; for example, recent work by Fan, Lu, and Luo (2023) and Fuchs and Wong (2022) has extended this optimal routing setup to incorporate multiple modes of transit and flexible transshipment costs. 

# 7.2 Estimation of model elasticities

We now turn to the estimation of the market elasticities. We offer two alternative estimation procedures. 

# 7.2.1 Estimation of the labor supply and demand

The basic idea is straightforward: we simultaneously estimate the labor supply and demand equations (62) and (63), which we write here in terms of the observed output $Y _ { i }$ and population $L _ { i }$ : 

$$
\ln Y _ {i} = \left(1 + \varepsilon_ {\text {l o c a l}} ^ {S}\right) \ln L _ {i} - \varepsilon_ {\text {g l o b a l}} ^ {S} \ln M A _ {i} ^ {\text {i n}} - \ln C _ {i} ^ {S} - \ln \phi^ {S} \tag {103}
$$

$$
\ln Y _ {i} = \left(1 - \varepsilon_ {\text {l o c a l}} ^ {D}\right) \ln L _ {i} + \varepsilon_ {\text {g l o b a l}} ^ {D} \ln M A _ {i} ^ {\text {o u t}} + \ln C _ {i} ^ {D} + \ln \phi^ {D}. \tag {104}
$$

Equations (103) and (104) look like regressions you could run (and indeed, that is the direction we are going), but with several major issues. The first issue is that the inward and outward market access terms $\ln M A _ { i } ^ { i n }$ and $\ln M A _ { i } ^ { o u t }$ are not directly observed in the data. Recall from Section 7.1 that we can estimate inward and outward market access from either a gravity regression or from the inversion of equations (64) and (65), yielding estimates $\ln M A _ { i } ^ { i n }$ and $\ln M A _ { i } ^ { o u t }$ , which can be used as proxies for the true market accesses. 

The second issue is a one of simultaneity: both the inward and outward market accesses and the equilibrium population will be correlated with the supply and demand shifters, so it would be inappropriate to use ordinary least squares to estimate equations (103) and (104). Instead, we need to construct instruments for the market access and population distribution in each regression. Valid instruments would have to be correlated with the population and market access but uncorrelated with unobserved supply shifters in equation (103) and unobserved demand shifters in equation (104). If there are origin-specific and/or destination-specific components of the trade costs, the instruments would have to be uncorrelated with them as well, as such origin-specific and destination-specific components act isomorphically to the demand and supply shifters, respectively. 

What are potentially valid instruments? The standard strategy employed in the estimation of supply and demand curves is to estimate the labor supply equation (103) by using an observable that shifts labor demand but that is uncorrelated with labor supply as an instrument. Similarly, one could estimate the labor demand equation (104) by using an observable that shifts labor supply but that is uncorrelated with labor demand as an instrument. The difficulty, of course, is identifying such observables. For example, Allen and Donaldson (2018), following Glaeser and Gottlieb (2009), propose using climate variation as an instrument for estimating labor demand equation (104), arguing that the advent of air conditioning shifted labor supply in locations with hot climates by improving the amenities there but did not directly affect labor demand. They also propose using crop suitability as an instrument for estimating labor supply equation (103), arguing that changes in international demand for certain crops were unlikely to be correlated with local labor supply shocks. 

One difficulty with this approach is that estimation of either equation (103) or (104) requires 

two instruments: one for the population of a location and one for its inward or outward market access. An option for constructing two instruments is as follows. Suppose you are estimating the labor supply equation (103) and have an instrument $\ln z _ { i }$ that is uncorrelated with the supply shifter $\ln C _ { i } ^ { S }$ . If the instrument is also uncorrelated with the supply shifter in all other locations, then you can construct a second instrument as a trade-friction weighted sum of the instrument in all other locations, i.e. $\begin{array} { r } { \ln M A _ { i } ^ { i n , Z } \equiv \ln \sum _ { j \in \mathcal { N } } \tilde { T } _ { j i } z _ { j } } \end{array}$ . Intuitively, the local variation in the instrument will be used to identify the local labor supply elasticity $\varepsilon _ { l o c a l } ^ { S }$ whereas the tradefriction weighted variation in the instrument elsewhere will be used to identify the global labor supply elasticity $\varepsilon _ { g l o b a l } ^ { S }$ . 

Once the model elasticities have been estimated, one can recover the supply and demand shifters $\left\{ C _ { i } ^ { S } , C _ { i } ^ { D } \right\}$ (to-scale) directly from equations (103) and (104). Intuitively, the supply and demand shifters act as residuals so that the equilibrium distribution of economic activity in the theory matches exactly the observed distribution of economic activity. 

# 7.2.2 The market access approach

A second approach is to use an existing policy change as an instrument. Suppose that you are conducting ex-post evaluation of a policy which affects trade frictions but does not affect the supply or demand shifters. For example, Ahlfeldt, Redding, Sturm, and Wolf (2015) analyze the impact of the construction and destruction of the Berlin Wall on the spatial distribution of economic activity in Berlin and identify model elasticities under the assumption that the distance of a location to the wall is uncorrelated with changes in its amenities or productivities. By equating equations (103) and (104), we can solve for $L _ { i }$ and $Y _ { i }$ as functions of inward and outward market accesses and the supply and demand shifters, as in the “market access” approach of Donaldson and Hornbeck (2016). Writing these expressions in first differences to compare variables before and after the policy was implemented yields: 

$$
\Delta \ln L _ {i} = \frac {\varepsilon_ {g l o b a l} ^ {S}}{\varepsilon_ {l o c a l} ^ {S} + \varepsilon_ {l o c a l} ^ {D}} \Delta \ln M A _ {i} ^ {i n} + \frac {\varepsilon_ {g l o b a l} ^ {D}}{\varepsilon_ {l o c a l} ^ {S} + \varepsilon_ {l o c a l} ^ {D}} \Delta \ln M A _ {i} ^ {o u t} + \Delta \ln \epsilon_ {i} ^ {L} \tag {105}
$$

$$
\Delta \ln Y _ {i} = \frac {\left(1 - \varepsilon_ {l o c a l} ^ {D}\right) \varepsilon_ {g l o b a l} ^ {S}}{\varepsilon_ {l o c a l} ^ {S} + \varepsilon_ {l o c a l} ^ {D}} \Delta \ln M A _ {i} ^ {i n} + \frac {\left(1 + \varepsilon_ {l o c a l} ^ {S}\right) \varepsilon_ {g l o b a l} ^ {D}}{\varepsilon_ {l o c a l} ^ {S} + \varepsilon_ {l o c a l} ^ {D}} \Delta \ln M A _ {i} ^ {o u t} + \Delta \ln \epsilon_ {i} ^ {Y}, \tag {106}
$$

where $\begin{array} { r } { \mathrm { l n } \epsilon _ { i } ^ { L } \equiv \frac { 1 } { \varepsilon _ { l o c a l } ^ { S } + \varepsilon _ { l o c a l } ^ { D } } \Delta \ln C _ { i } ^ { S } + \frac { 1 } { \varepsilon _ { l o c a l } ^ { S } + \varepsilon _ { l o c a l } ^ { D } } \Delta \ln C _ { i } ^ { D } + \frac { 1 } { \varepsilon _ { l o c a l } ^ { S } + \varepsilon _ { l o c a l } ^ { D } } \ln \phi ^ { D } + \frac { 1 } { \varepsilon _ { l o c a l } ^ { S } + \varepsilon _ { l o c a l } ^ { D } } \Delta \ln C _ { i } ^ { S } } \end{array}$ ln ϕS and ∆ ln ϵYi $\begin{array} { r } { \equiv \frac { 1 - \varepsilon _ { l o c a l } ^ { D } } { \varepsilon _ { l o c a l } ^ { S } + \varepsilon _ { l o c a l } ^ { D } } \Delta \ln C _ { i } ^ { S } + \frac { 1 + \varepsilon _ { l o c a l } ^ { S } } { \varepsilon _ { l o c a l } ^ { S } + \varepsilon _ { l o c a l } ^ { D } } \Delta \ln C _ { i } ^ { D } + \frac { 1 - \varepsilon _ { l o c a l } ^ { D } } { \varepsilon _ { l o c a l } ^ { S } + \varepsilon _ { l o c a l } ^ { D } } \Delta \ln \phi ^ { S } + \frac { 1 + \varepsilon _ { l o c a l } ^ { S } } { \varepsilon _ { l o c a l } ^ { S } + \varepsilon _ { l o c a l } ^ { D } } \Delta \ln \phi ^ { D } } \end{array}$ εSloc al +εDloca εSlocal εSlocal 1−εDlocal 1+εSlocal are residuals that depend on the changes in both the supply and demand shifters. One can instrument for the changes in inward and outward market access with exposure to the policy in equations (105) and (106), which are valid instruments as long as the policy exposure is not correlated with the supply and demand shifters. The four resulting coefficient estimates can then be used 

to recover the four model elasticities. (For example, the local demand elasticity $\varepsilon _ { l o c a l } ^ { D }$ is one minus the ratio of the inward market access coefficients in the output and labor regressions). How should one construct the policy exposure instruments? The simplest way would be to use equations (64) and (65) to calculate the market access terms before and after the policy shock, holding fixed income and expenditure before the policy shock to ensure that the changes in market access captured by the instruments is only arising from the effect of the policy shock on trade frictions. 

There is one important (but subtle) caveat to the market access approach: when trade frictions are symmetric, the inward and outward market access terms will be equal up to scale (see Section 4.6). As a result, this approach will not be able to recover the full set of model elasticities. Nevertheless, one can still determine the resulting change in the equilibrium distribution of economic activity, a point which we show in the next (technical) subsection. 

# 7.2.3 [Technical] The market access approach with symmetric trade costs

Recall that when trade costs are symmetric, equations (64) and (65) become: 

$$
\kappa M A _ {i} ^ {i n} = \sum_ {j \in \mathcal {N}} T _ {i j} \times \frac {E _ {j}}{M A _ {j} ^ {i n}}, \tag {107}
$$

where $M A _ { i } ^ { o u t } = \kappa M A _ { i } ^ { i n }$ for some scalar $\kappa _ { \iota }$ , which we can write in changes as follows: 

$$
\widehat {\kappa M A _ {i} ^ {i n}} = \sum_ {j \in \mathcal {N}} \left(\frac {X _ {i j} ^ {A}}{Y _ {i} ^ {A}}\right) \hat {T} _ {i j} \frac {\hat {E} _ {j}}{M \hat {A} _ {j} ^ {i n}}. \tag {108}
$$

Equation the supply and demand equations (103) and (104) and solving for the change in output holding constant the supply and demand shifters yields: 

$$
\hat {Y _ {i}} \propto \widehat {M A _ {i} ^ {i n}} \frac {\left(1 - \varepsilon_ {l o c a l} ^ {D}\right) \varepsilon_ {g l o b a l} ^ {S} + \left(1 + \varepsilon_ {l o c a l} ^ {S}\right) \varepsilon_ {g l o b a l} ^ {D}}{\varepsilon_ {l o c a l} ^ {S} + \varepsilon_ {l o c a l} ^ {D}}, \tag {109}
$$

which, when combined with equation (108) and imposing that $E _ { i } = Y _ { i }$ yields: 

$$
\widehat {\nu M A _ {i} ^ {i n}} = \sum_ {j \in \mathcal {N}} \left(\frac {X _ {i j} ^ {A}}{Y _ {i} ^ {A}}\right) \hat {T} _ {i j} \left(\widehat {M A _ {j} ^ {i n}}\right) ^ {- 1 + \frac {\left(1 - \varepsilon_ {l o c a l} ^ {D}\right) \varepsilon_ {g l o b a l} ^ {S} + \left(1 + \varepsilon_ {l o c a l} ^ {S}\right) \varepsilon_ {g l o b a l} ^ {D}}{\varepsilon_ {l o c a l} ^ {S} + \varepsilon_ {l o c a l} ^ {D}}}, \tag {110}
$$

for an endogenous scalar $\nu$ that depends on $\kappa$ and the endogenous scalars $\phi ^ { D }$ and $\phi ^ { S }$ . Equation (110) shares the same mathematical structure as equation (47) in Section 4.6 and hence its properties follow from the discussion there. In particular, as long as $\begin{array} { r } { \left| 1 - \frac { \left( 1 - \varepsilon _ { l o c a l } ^ { D } \right) \varepsilon _ { g l o b a l } ^ { S } + \left( 1 + \varepsilon _ { l o c a l } ^ { S } \right) \varepsilon _ { g l o b a l } ^ { \bar { D } } } { \varepsilon _ { l o c a l } ^ { S } + \varepsilon _ { l o c a l } ^ { D } } \right| < } \end{array}$ 

1,there exists a unique (to-scale) set of inward market access changes that solve equation (110). Given these $\left\{ \widetilde { M A _ { i } ^ { i n } } \right\} _ { i \in \mathcal { N } } ,$ the equilibrium changes in output $\big \{ \hat { Y } _ { i } \big \} _ { i \in \mathcal { N } }$ can be calculated from equation (109) (and the choice of numeraire requiring that $\Sigma _ { j \in \mathcal { N } } Y _ { i } ^ { A } \hat { Y } _ { i } = 1$ . An equivalent derivation as in equation (109) can be used to solve for the equilibrium change in population, yielding $\hat { L } _ { i } \propto \widetilde { M A _ { i } ^ { i n } } ^ { \frac { \varepsilon _ { g l o b a l } ^ { s } + \varepsilon _ { g l o b a l } ^ { L } } { \varepsilon _ { l o c a l } ^ { S } + \varepsilon _ { l o c a l } ^ { D } } }$ εSglobal +εDglobal where the scale is determined by $\begin{array} { r } { \sum _ { j \in \mathcal { N } } \left( L _ { i } ^ { A } / \bar { L } \right) \hat { L } _ { i } = 1 } \end{array}$ . 

# 7.3 Testing the Models

So far we have discussed how to generate predictions from the quantitative geography framework and how to fit it to the data. Because fitting the model to the data involves exactly matching the observed distribution of economic activity, this leaves little room for testing the model. How can we then test whether the framework yields reliable predictions? 

Testing a multi-region model is a significant challenge: a precursor of quantitative spatial models has been the large literature of Computable General Equlibrium (CGE) models that have been used to predict what happens in the economy when tariffs in certain sectors decrease (Brown and Stern (1989); Brown, Deardorff, and Stern (1991)); within and across countries. Despite nearly five decades of development CGE models have faced criticism for failing to accurately predict what happened to output and sectoral trade after tariff changes (see e.g. Kehoe (2005) and Kehoe, Pujolas, and Rossbach (2017)). We believe that the quantitative geography model hold more promise, for two reasons: first, the frameworks, as we showed in Sections 4-6, they are grounded in modern economic theory. Second, as we showed in Sections 7.1 and 7.2, their model parameters are estimated using modern techniques (unlike CGE models where parameters are typically externally calibrated). 

We now discuss a few ways to test these models. The first way to test quantitative economic geography framework is to look at aspects of its cross-sectional predictions that are not related to the data we aim to fit. This is difficult, as the model is calibrated to exactly match observed incomes and populations, but one way to do so is to see if the calibrated supply and demand recovered “make sense.” For example, Desmet and Rossi-Hansberg (2013) show that the amenities recovered to calibrate the model indeed correlate well with observed characteristics that likely correlate strongly with amenities such as climate and life-quality measures. 

A second way to test the quantitative economic geography framework is to see how well it is able to predict how the spatial distribution of economic activity actually responds to an observed shock. Typically, this is done by holding the calibrated supply and demand shifters fixed at their pre-shock levels, simulating the effect of the spatial shock, and comparing the predicted changes to the observed changes. For example, Donaldson (2018) shows that the observed change in real rental rates in response to railroad construction in Colonial India closely matches the model predicted changes. Similarly, Heblich, Redding, and Sturm (2020) and Nagy 

(2023) show that quantitative models are able to closely match the observed change in the spatial equilibrium resulting from the expansion of the London rail and underground network and the U.S. railroad expansion, respectively. Adao, Arkolakis, and Esposito (2019) and Allen, Fuchs, Ganapati, Graziano, Madera, and Montoriol-Garriga (2020) use linearized versions of spatial models to look at whether more contemporary changes in cross-sectional income are predicted by the ones in the model, when taking into account indirect general equilibrium effects discussed in Section 6.7. The two papers simulate an arguably exogenous shock, the China shock (Autor, Dorn, and Hanson (2013)) and tourism shock in Barcelona, respectively, after controlling for variation orthogonal to the shock. In fact, Allen, Fuchs, Ganapati, Graziano, Madera, and Montoriol-Garriga (2020) show that adding higher order effects, as in equation (94), offers a successively better approximation to the full equilibrium effects of tourism shock in Barcelona. One issue with these tests is that it is not clear what constitutes a “good” model fit. 

Recently, Adao, Costinot, and Donaldson (2025) offer a formal econometric test for quan- ˜ titative models based on a similar approach. In their test, one calibrates the supply and demand shifters to match the observed distribution of economic activity both before and after the shock. If the quantitative framework correctly predicts the causal changes in the distribution of economic activity from some exogenous shock, then the changes in the supply and demand shifters should themselves be uncorrelated with that shock, which offers moment restrictions that can be tested. They illustrate their approach by testing the setup of Fajgelbaum, Goldberg, Kennedy, and Khandelwal (2020) that is used to evaluate the impact of Trump tariffs across regions in the United States. 

# 8 Lessons learned

Over the past decade, there have been many successful examples of papers that combine quantitative regional frameworks with detailed spatial data to perform policy-relevant counterfactual analyses. In this section, we briefly highlight several of the lessons learned from this fruitful approach. 

# 8.1 The importance of market access

It has been long recognized that a location’s connection to other markets can have important impacts on firms, residents and workers, see e.g. Redding and Venables (2004) and Hanson (2005). And, as discussed above, market access plays a crucial role in quantitative regional models as the measure for how the rest of the economy affects economic outcomes in a given location. This insight has enabled the tractable analyses of some of the largest infrastructure investments in history. For example, Donaldson and Hornbeck (2016) examine how the con-

struction of the U.S. railroad network affected the spatial economy through its impact on each location’s market access. They verify empirically that improvements in market access arising from the railroad network expansion causally increase the agricultural land values in those locations and then use this estimated elasticity to quantify the aggregate welfare gains that arise from the complete railroad expansion. Allen and Arkolakis (2014) employ a similar approach to estimate the welfare gains from the construction of the U.S. interstate highway system. In both cases, the authors estimate substantial economic gains from these large infrastructure projects through increased economic integration. Subsequent research has confirmed that improvements in market access lead to substantial impacts on the local economy in empirical contexts as varied as highway construction in India (Alder, 2016; Vogel, Hanson, Khandelwal, Liu, and Park, 2024), high speed rail in China (Egger, Loumeau, and Loumeau, 2023), and land-mine removal in Mozambique (Chiovelli, Michalopoulos, and Papaioannou, 2018). 

# 8.2 The importance of agglomeration forces

A second lesson from the advent of quantitative regional frameworks is the crucial role for agglomeration forces play when estimating the impact of economic policies. For example, the Allen and Arkolakis (2014) estimates of the welfare gains from the construction of the interstate highway system vary by $3 0 \%$ depending on the chosen strength of agglomeration forces. This is roughly the same order of magnitude seen by Ahlfeldt, Redding, Sturm, and Wolf (2015) and Heblich, Redding, and Sturm (2020) when varying the strength of agglomeration forces to assess the impacts of the Berlin Wall and the London Underground, respectively, in related quantitative urban frameworks (see Redding (2025) for details). It is not just the quantitative predictions of the model that depend on the chosen strength of agglomeration forces, however: the qualitative predictions matter too. As emphasized in Section above (see equations (35) and (36)), the sign of the supply and demand labor elasticities depends on the strength of the agglomeration forces (and hence so too can the sign of comparative statics). And, as noted in Section 6.7, counterfactual analyses may no longer be well defined if agglomeration forces are so strong so as to admit multiplicity. Relatedly, Allen and Donaldson (2018) show that agglomeration forces can also result in path dependence, whereby the initial conditions of an economy can lead it toward different steady state spatial economies, holding constant the underlying geography. This phenomenon admits the possibility of that short-run economic policies can have long-run impacts by moving an economy from one steady state path to another.Heterogeneous effects of homogeneous policies 

A third lesson from quantitative regional models is that the same policy can have very different impacts on different locations depending on the underlying geography. This heterogeneity arises naturally from the theory; for example, equation (6.7)—showing that the welfare gains of reducing a bilateral trade cost is equal to its share of trade flows—highlights that the same 

reduction in trade costs will have much larger effects when the partners trade more. But what is especially impressive is that these model predictions about the heterogeneity of policy impacts, when tested empirically, perform quite well. Recall from above that Donaldson (2018) documented that a quantitative spatial model does a good job explaining the observed variation in welfare gains from railroads across colonial India. Similarly, in their pioneering work, Redding and Sturm (2008) analyzed the impact of the fall of the Iron Curtain on population growth in West Germany, documenting that population growth was disproportionately slower for West German cities nearer to East Germany. Less obviously, they showed that the population declines were relatively greater for smaller cities nearer the East German border, which is consistent with the idea that the Iron curtain resulted in a larger decline in market access for smaller cities (who traded less with themselves). To give a third example, Monte, Redding, and Rossi-Hansberg (2018) show that there are greater increases in employment in locations that “win” a million-dollar plant (c.f. Greenstone, Hornbeck, and Moretti (2010)) when their “commuting market access” is greater. Each of these examples highlights not only the fact that geography matters when designing economic policy but also that the quantitative regional framework does a good job explaining exactly how. 

# 9 Exciting extensions

In this section, we briefly describe recent advances that have extended the workhorse quantitative regional economic framework presented in the previous section in fruitful directions. 

# 9.1 Multiple sectors

The framework presented previously considers a single sector. Extending the framework to incorporate multiple sectors is a main motivation for quantitative work: since the pioneering work of Ricardo in 1817 there is a long tradition in international economics to consider multisector models and a literature in quantitative macroeconomics that considers multi-sector models as important for addressing the observed growth patterns (Long Jr and Plosser (1983)). More recently Anderson and Yotov (2010), Costinot, Donaldson, and Komunjer (2012), Ossa (2015), Caliendo and Parro (2015), Levchenko and Zhang (2016), among others, have discussed an extension of a quantitative international trade model with multiple sectors and Caliendo, Parro, Rossi-Hansberg, and Sarte (2018), Faber and Gaubert (2019), Arkolakis and Walsh (2023) have considered such an extension in a quantitative economic geography framework. Incorporating multiple sectors in the framework of Section 6 is straightforward. We now sketch out the procedure. 

Suppose that there are $s \in \{ 1 , . . . , S \} \equiv \mathcal { S }$ sectors and that the utility of a representative 

consumer residing in location $j \in \mathcal N$ is given by: 

$$
W _ {j} = \left(\sum_ {s \in \mathcal {S}} \beta_ {j} ^ {s} \left(\sum_ {i \in \mathcal {N}} c _ {i j} ^ {\frac {\sigma^ {s} - 1}{\sigma^ {s}}}\right) ^ {\frac {\sigma^ {s}}{\sigma^ {s} - 1} \frac {\sigma - 1}{\sigma}}\right) ^ {\frac {\sigma}{\sigma - 1}} u _ {j}, \tag {111}
$$

where $\begin{array} { r } { \beta _ { s } ^ { j } , \sum _ { s } \beta _ { s } ^ { j } = 1 , } \end{array}$ , are consumption weights that a consumer in location $j$ assigns to the consumption of sector $s \in S$ . The elasticity of substitution $\sigma _ { s }$ is now across locations at the sector level and captures the varying degrees of substitutability that each goods in every sector may experience, whereas the elasticity of substitution $\sigma$ is now across sectors. For example, estimates for petroleum indicate that this commodity is very substitutable, while advanced semiconductors are clearly less. 

Maximizing the consumer’s demand from equation (111) yields the following sectoral gravity equation: 

$$
X _ {i j} ^ {s} = \frac {\left(p _ {i j} ^ {s}\right) ^ {1 - \sigma^ {s}}}{\left(P _ {j} ^ {s}\right) ^ {1 - \sigma^ {s}}} e _ {j} ^ {s} w _ {j} L _ {j}, \tag {112}
$$

where $\begin{array} { r } { e _ { j } ^ { s } \equiv \frac { \left( P _ { j } ^ { s } \right) ^ { 1 - \sigma } } { P _ { j } ^ { 1 - \sigma } } \mathbf { j } } \end{array}$  P sj  1 − σ is the consumption share of sector $s$ in location $\begin{array} { r } { j , \left( P _ { j } ^ { s } \right) ^ { 1 - \sigma ^ { s } } \equiv \sum _ { i \in \mathcal { N } } \left( p _ { i j } ^ { s } \right) ^ { 1 - \sigma ^ { s } } } \end{array}$ is a sector-specific price index, and $\begin{array} { r } { P _ { j } ^ { 1 - \sigma } \equiv \sum _ { s \in \mathcal { S } } \left( \beta _ { j } ^ { s } \right) ^ { \sigma } \left( P _ { j } ^ { s } \right) ^ { 1 - \sigma } } \end{array}$ is the composite price index. 

The rest of the multi-sector framework bear a strong resemblance to the single sector baseline framework. Welfare of agents residing in location $j$ can still be written in indirect utility form as $W _ { i } = w _ { i } u _ { i } / P _ { i }$ . Goods market clearing now simply requires summing over sectors as well, i.e. $\begin{array} { r } { Y _ { i } = \sum _ { s \in \mathcal { S } } \sum _ { j \in \mathcal { N } } X _ { i j } ^ { s } } \end{array}$ and $E _ { i } = \textstyle \sum _ { s \in { \mathcal { S } } } \sum _ { j \in { \mathcal { N } } } X _ { j i } ^ { s }$ . Labor market clearing remains unchanged, i.e. $w _ { i } L _ { i } = E _ { i } = Y _ { i }$ . 

For example, consider the simple case where labor is the only factor of production (so that $p _ { i j } ^ { s } = \tau _ { i j } ^ { s } w _ { i } / A _ { i } ^ { s } )$ and the productivity and amenity externalities follow equations (22) and (23). Then following the derivations in Section 4.3, we arrive at the following (implicit) labor demand curve: 

$$
w _ {i} L _ {i} = \sum_ {s \in \mathcal {S}} L _ {i} ^ {\alpha (\sigma^ {s} - 1)} w _ {i} ^ {1 - \sigma_ {s}} \sum_ {j \in \mathcal {N}} \left(\tau_ {i j} ^ {s} / \bar {A} _ {i} ^ {s}\right) ^ {1 - \sigma_ {s}} P _ {j} ^ {\sigma^ {s} - \sigma} w _ {j} L _ {j}, \tag {113}
$$

while the labor supply curve is identical to equation (36) (just the new definition of the consumer price index). 

The addition of multiple sectors can be further extended to include round-about production (as in Krugman and Venables (1995) and Eaton and Kortum (2002) and discussed in Section 5.8) and input-output linkages between sectors, as in Caliendo and Parro (2015) and Caliendo, Dvorkin, and Parro (2019). While the techniques for establishing uniqueness discussed above do not directly apply to settings with multiple sectors with rich input-output structures, Bi-

fulco, Gluck, Krebs, and Kukharskyy (2022) show how an alternative methodology based on ¨ the monotonicity of the partial cross-derivatives of the system can be used to establish the (toscale) uniqueness of the equilibrium. 

# 9.2 Migration costs

The framework presented in Section 6 abstracts from any costs associated with movement across space. Yet there is substantial evidence that workers and the local economy responds slowly to changes to economic shocks, see e.g. Autor, Dorn, and Hanson (2013), Dix-Carneiro and Kovak (2017), and Yagan (2019), suggesting that dynamics play an important role in understanding the evolution of the spatial economy. There has been substantial recent progress extending the quantitative economic geography framework to incorporate such migration costs, including Tombe and Zhu (2019), Allen, de Castro Dobbin, and Morten (2018), Desmet, Nagy, and Rossi-Hansberg (2018), Caliendo, Dvorkin, and Parro (2019), Allen and Donaldson (2018) and Kleinman, Liu, and Redding (2023) Here we sketch out a simple static quantitative economic geography model with migration costs, and we refer the interested reader to Desmet and Parro (2025) for more details of recent dynamic advances. 

Suppose that each location $i \in \mathcal N$ is endowed with an initial allocation of agents $\bar { L } _ { i }$ . Each agent $\omega$ initially allocated to $i \in \mathcal N$ then choose in which location to live in order to solve: 

$$
\max  _ {j \in \mathcal {N}} \frac {W _ {j}}{\mu_ {i j}} \times \varepsilon_ {i j} (\omega), \tag {114}
$$

where $W _ { j }$ is the welfare (real wage and amenity value) of residing in location $j \in \mathcal { N } , \mu _ { i j } \geq 1$ is a bilateral migration friction, and $\varepsilon _ { i j } \left( \omega \right)$ is an idiosyncratic preference term, which we assume is Frechet distributed, i.e. $\varepsilon _ { i j } \left( \omega \right) \sim e ^ { - z ^ { - \theta } }$ . As in Sections 5.1 and 5.5, maximization of (114) leads to a probabilistic decision where the fraction of workers initially allocated to $i \in \mathcal N$ that migrate to $j \in \mathcal N$ is: 

$$
L _ {i j} / \bar {L} _ {i} = \left(W _ {j} / \mu_ {i j}\right) ^ {\theta} / \sum_ {k \in \mathcal {N}} \left(W _ {k} / \mu_ {i k}\right) ^ {\theta}. \tag {115}
$$

Summing across all possible origin locations then allows us to construct the following modified labor supply curve: 

$$
\ln w _ {j} = \frac {1}{\theta} \ln L _ {j} + \frac {1}{1 - \sigma} \ln P _ {j} ^ {1 - \sigma} - \ln u _ {i} - \frac {1}{\theta} \ln \sum_ {i \in \mathcal {N}} \frac {\mu_ {i j} ^ {- \theta} \bar {L} _ {i}}{\sum_ {k \in \mathcal {N}} \left(W _ {k} / \mu_ {i k}\right) ^ {\theta}}. \tag {116}
$$

The first three terms of equation (116) are identical to the standard labor supply curve in equation (62), but the fourth term is new. This fourth term is a migration market access term that summarizes how costly it is for migrants to reach location $j$ relative to all other locations. Intuitively, the labor supply curve will shift outward for locations to which it is easier to migrate. 

# 9.3 Commuting

The seminal work of Ahlfeldt, Redding, Sturm, and Wolf (2015) studies the linkages within cities and is discussed in detail in Redding (2025). Here we make the connection with the workhorse setup and discuss associated extensions. We assume that each agent $\omega$ that chooses to reside in $i \in \mathcal N$ and then chooses to commute in location $j$ solves: 

$$
\max  _ {i, j \in \mathcal {N}} \frac {u _ {i} w _ {j}}{\mu_ {i j}} \times \varepsilon_ {i j} (\omega), \tag {117}
$$

where $\mu _ { i j } \geq 1$ is a bilateral commuting friction, and $\varepsilon _ { i j } \left( \omega \right)$ is an idiosyncratic preference term, which again is Frechet distributed, i.e. $\varepsilon _ { i j } \left( \omega \right) \sim e ^ { - z ^ { - \theta } }$ . The total number of workers residing in $i$ and working in $j$ is given by 

$$
L _ {i j} = \mu_ {i j} ^ {- \theta} \times u _ {i} ^ {\theta} \times w _ {j} ^ {\theta} \times \frac {\bar {L}}{W ^ {\theta}}
$$

where $\begin{array} { r } { W ^ { \theta } = \sum \mu _ { i j } ^ { - \theta } \times u _ { i } ^ { \theta } \times w _ { j } ^ { \theta } , } \end{array}$ , and $\bar { L } = \sum _ { i , k \in \mathcal { N } } L _ { i k }$ is the total number of workers in location $j$ . The assumptions of the workhorse framework continue to hold but we replace the goods market clearing condition with the labor adding up constraints so that $\begin{array} { r } { L _ { i } ^ { R } = \sum _ { k \in \mathcal { N } } L _ { i k } } \end{array}$ , is the total number of residents in location i, $\begin{array} { r } { , L _ { j } ^ { F } = \sum _ { k \in \mathcal { N } } L _ { k j } } \end{array}$ . 

We can specify equation (61) as 

$$
L _ {i j} = T _ {i j} \times \frac {L _ {i} ^ {R}}{M A _ {i} ^ {o u t}} \times \frac {L _ {j} ^ {F}}{M A _ {j} ^ {i n}} \tag {118}
$$

where now $L _ { i j }$ is the number of commuters from location $i$ to location $j$ . Thus, equation 118 apportions the share of residents from i that chose $j$ as their place to commute. Due to gravity, this depends on the size of each commuting destination, thus the presence of $L _ { j } ^ { F }$ , but also on the access of both locations. Equations (64), (65) continue to hold with $Y _ { i } = L _ { i } ^ { R }$ and $E _ { j } = L _ { j } ^ { F }$ and $M A _ { i } ^ { o u t }$ and $M A _ { j } ^ { i n }$ are interpreted as outward and inward commuting access, respectively. 

To derive the equilibrium system notice that the adding up constraints provide a link between the model equilibrium variables, $L _ { i } ^ { R }$ , and $L _ { i } ^ { F }$ and $M A _ { i } ^ { o u t } , M A _ { i } ^ { i n } .$ $M A _ { i } ^ { i n }$ , respectively, analogous to the labor demand and supply curves discussed above. Combining these equations with (118), and 64, 65 provides a system of equation parallel to (66) and (67). However, there two important final things to be determined, in particular, the relationship of $u _ { i }$ and $w _ { j }$ to other endogenous variables. Ahlfeldt, Redding, Sturm, and Wolf (2015) assume that there is a single homogeneous freely traded good, so that $w _ { j } = A _ { j } ,$ where $A _ { j }$ is the productivity of the location of work. They then specify a relation between $u _ { i }$ and residential population in nearby regions, and $A _ { j }$ and worker population nearby. Allen and Arkolakis (2022) specify local spillovers in the form of equations (22) and (23). Monte, Redding, and Rossi-Hansberg (2018) assume each region produces a differentiated good and that there are trade costs across regions specify a 

relationship of these functions with nearby population. This implies that wage is given by equation (35), which adds further richness to the equilibrium system. 

# 9.4 Spatial spillovers

So far we have discussed spatial connections through labor as a factor of production and trade in goods. An intriguing possibility of spatial spillovers is the spread of knowledge, i.e. knowledge spillovers. A way to model knowledge spillovers is through the dependence of productivity, $A _ { i } ,$ on knowledge in nearby regions. Arkolakis, Lee, and Peters (2020); Pellegrina and Sotelo (2021); Cai, Caliendo, Parro, and Xiang (2022); Crews (2023) discuss knowledge spillovers through the movement of skilled labor across countries and regions. An alternative possibility is knowledge spillover through firms. Giroud, Lenzu, Maingi, and Mueller (2024) and Liang (2024) discuss within-firm across-plants knowledge spillovers in productivity and intangible capital, respectively, using multi-location multi-plant frameworks. 

Additional factors of production may be mobile and generate the possibility of spatial spillovers. For example, Conte, Desmet, and Rossi-Hansberg (2022) introduce energy as a factor of production that is freely traded across the world. In addition, Arkolakis and Walsh (2023) extend the workhorse setup to incorporate trade in electricity through transmission lines. 

A particularly important spatial aspect is reserved for environmental economics topics due to the environmental diversity across space and the role of environmental externalities. For example, Kang (2024) develops a model of fisheries where fishing takes place in potentially remote locations and fish is traded across the world. Cruz and Rossi-Hansberg (2024) study negative externalities in productivities and amenities across regions that arise from the rise in temperature due to global carbon emissions. This line of research and its relation to spatial economics is discussed in detail in Balboni and Shapiro (2025). 

# 10 The next generation

In the previous sections, we have discussed the major advances made in turning the theory of earlier generations “quantitative” so that it can be combined with detailed spatial data to perform real world analyses. We have also highlighted current research extending and expanding the scope of quantitative economic geography frameworks. Here, we conclude by offering our best guess of what the next generation of quantitative models may bring. To do so, we highlight four possible paths forward for future research. 

# 10.1 Taking inequality and agent heterogeneity more seriously

In the quantitative economic geography framework developed above, all agents were identical (except, perhaps, up to an idiosyncratic component). But, of course, in reality people differ along a wide variety of dimensions, including their education, skills, race, sex, gender, age, etc. These differences affect agent’s preferences of where they live, where they work, and how they interact with others. And an important concern when making public policy is how such policies differentially impact different groups of people. 

Extending the quantitative regional framework to incorporate such a rich heterogeneity is necessary to tackle issues related to inequality head on, and we think this is a fruitful direction for the next generation of research. Indeed, there are a number of recent papers that have made important strides in this direction. For example, Tsivanidis (2022) incorporates multiple types of workers endowed with varying skills and non-homothetic preferences to evaluate the differential impacts of rapid bus transit in Bogota. Similarly, Zarate (2022) examines the differential ´ impacts of subway expansions in Mexico City using a model where agents are ex-ante identical but ex-post sort into different sectors and locations using a Roy (1951)-like framework. 

One important policy-relevant topic that requires a richer treatment of agent heterogeneity is segregation. Here too there are several recent exciting advances. Hoelzlein (2019) introduces agents who differ in their skill endowments and non-homothetic preferences to study the sorting of firms and residents and the impact of Empowerment Zones on segregation in Los Angeles. Recent papers by Bagagli (2023) and Weiwu (2023) specifically study racial segregation, incorporating homophily into agents’ preferences and used such a framework to understand how highway construction has amplified segregation in cities. Allen, Arkolakis, and Li (2024) characterize the theoretical properties of quantitative spatial frameworks (i.e. including regional, urban, and trade) with many types of agents and arbitrary cross-group preferences. 

Another dimension on which agents can differ is in their family structure. Recent work has emphasized how family ties affect preferences through remittances (see Albert and Monras (2022)), but there are number of interesting questions relating to the spatial implications of “family economics” that are yet to be addressed. How do families with and without children differ in their preferences and choices of location? How should couples jointly choose their location? How do people’s preferences over locations and different types of amenities change over their life cycle? 

There are also interesting and important issues related to inequality across generations. How do parents’ choices of the children’s human capital acquisition affect the intergenerational transfer of wealth, long-run migration patterns, and ultimately spatial inequality? To date, most dynamic quantitative economic geography models have assumed an infinitely lived agent (see Desmet and Parro (2025) ). While there has been some early work on the topic (see e.g. Eckert and Kleineberg (2024) for an application to U.S. school funding and Pellegrina and 

Sotelo (2021) for the implications for the evolution of comparative advantage), extending the dynamic economic geography frameworks to more seriously model inter-generational transfers of wealth and capital can better guide policies meant to ameliorate systemic inequalities and help us better understand the long-term effects of spatial policies. 

# 10.2 Quantum quantitative economic geography

One of the central achievements of quantitative economic geography is its ability to incorporate a rich geography: locations vary in their innate productivities, amenities, and costs in interacting with other locations. Yet while the geography being considered is incredibly rich in these dimensions, it is also surprisingly poor across others. Here are two dimensions in which we think future work could make substantial improvements. 

First, in the workhorse framework presented above, all locations produce differentiated varieties that enter consumer demand in the same way. While it is straightforward to extend this framework to multiple sectors (see Section 9.1), currently absent is the idea that locations may play fundamentally different roles. For example, the “Central Place Theory” of Christaller (1933) proposed a hierarchy of locations, where large central cities are surrounded by smaller towns, which are in turn surrounded by agricultural areas. The central cities are able to sustain more specialized businesses (e.g. sports stadiums) that surrounding residents also enjoy, whereas firms in smaller surrounding towns are only those necessary for more every day needs (e.g. grocery stores). Such concepts of spatial specialization are absent from current quantitative frameworks. However, there have been some recent advances in this direction: for example, Nagy (2023) develops a model featuring cities surrounded by agricultural hinterlands. 

A necessary step toward allowing different locations to host distinctly different sets of firms is to “quantize” quantitative economic geography by considering individual firms rather than a continuum of firms. Determining how a discrete number of firms optimally choose their locations is a difficult problem, but there have been exciting steps forward in recent work by Oberfield, Rossi-Hansberg, Sarte, and Trachter (2024). Technical advances in combinatorial discrete choice problems by Arkolakis, Eckert, and Shi (2017) and applications of clustering techniques (see e.g. Allen (2023)) may offer other potential paths forward. 

# 10.3 Optimal policy design

One of the most exciting components of quantitative economic geography is its ability to offer normative prescriptions to real world spatial policies. Existing research has mostly used the framework to evaluate the welfare effects of economic policies ex post. But we think there is an enormous potential to use the same frameworks to design better spatial policies ex ante. 

As summarized in Fajgelbaum and Gaubert (2025), there has been substantial recent progress 

in characterizing the welfare implications of this class of quantitative economic geography models, and we refer the interested reader to that chapter for a much more thorough treatment of the related research. Yet while the work by Fajgelbaum and Gaubert (2020), Fajgelbaum and Schaal (2020), and Donald, Fukui, and Miyauchi (2023) have clarified the normative properties of this broad class of models and recent work has started to assess political economy considerations in spatial policies (see e.g. Bordeu (2023) and Fajgelbaum, Gaubert, Gorton, Morales, and Schaal (2023)), there remains the opportunity to apply these models to policy design. Pressing questions in this vein include: How should we best design congestion tolls on city streets? How should a government best allocate its infrastructure budget across its road network? What is the best spatial transfer scheme to address both efficiency and equity considerations? 

Of particular interest to policy design is the possibility of multiplicity of equilibria described above, as temporary policies can have long-run or even permanent effects. For example, Bleakley and Lin (2012) find evidence of long-run impacts of the presence of “portage” sites long after river trade had declined in importance; consistent with this empirical result, Allen and Donaldson (2018) show that dynamic extensions of the framework above can generate path dependence where small shocks lead to permanently different steady-state equilibria. Owens III, Rossi-Hansberg, and Sarte (2020) convincingly demonstrate that coordination failures can yield to multiple equilibria that can help help explain why parts of Detroit have not been able to redevelop. And Monte, Porcher, and Rossi-Hansberg (2023) show that the COVID pandemic shifted larger cities like New York and San Francisco to a new equilibrium in which most workers commute to one in which most workers work partially from home. An outstanding question is how one can fully enumerate all possible spatial equilibria given the underlying geography; however, exciting recent research by Garg (2025) and Ouazad (2024) using homotopic methods has made substantial progress on this question. 

# 10.4 Market power

A fourth exciting direction for future research would be to explore the implications of incorporating firm market power. The existing workhorse models described above typically assume perfect or monopolistic competition in both product and labor markets, but in reality firms likely exhibit substantial market power in both. While there has been substantial number of papers examining the role of market power in international trade (see e.g. Arkolakis, Costinot, Donaldson, and Rodr´ıguez-Clare (2019) for the study of imperfect competition in product markets and Felix (2021) for the study of imperfect competition in labor markets), much less has been done studying the role of market power in the spatial distribution of economic activity. Notable exceptions include Gutierrez (2022), who develops a model featuring imperfect com- ´ petition in both product and labor markets to understand how market power affects the welfare gains from trade, and Hong (2024), who develops a model with monopsony power of firms to 

consider the two-sided matching of heterogeneous firms and workers and its implications for the equilibrium distribution of economic activity across space. 

Another potential source of market power is in the transportation sector. Typically, quantitative economic models assume that trade frictions are exogenous or, equivalently, that the transportation sector is perfectly competitive. But in reality, the transportation sector is oftentimes quite concentrated, and the imperfect competition in the sector can affect the equilibrium spatial distribution of economic activity. There have been several recent studies examining how imperfect competition in the transportation sector affects trade costs, including Brancaccio, Kalouptsidi, and Papageorgiou (2020) (bulk shipping), Wong (2022) (container shipping), and Allen, Atkin, Cleves, and Hernandez (2024) (trucking), but it remains an unresolved question about how such endogenous trade costs affects the equilibrium distribution of economic activity. 

# Conclusion

Over the past decade, there has been a remarkable advance in the understanding of how geography shapes the spatial distribution of the economy. In this handbook chapter, we have shown how a new quantitative economic geography framework has been able to unify the economic mechanisms from two seminal models of the previous generation in such a way that they can be combined with detailed spatial economic data to better understand real world economic policy. There are many exciting aspects of the framework. It is tractable: it is mathematically wellbehaved and easy to solve. It is also powerful, allowing one to evaluate the spatial impacts of nearly any change in the underlying geography. And it is flexible: it has many different micro-foundations and it can be extended in many directions. Most exciting, while some of those extensions have been explored, most have not, suggesting that the future of economic geography is bright. 

# References



ADAO, R., C. ARKOLAKIS, AND F. ESPOSITO (2019): “General equilibrium effects in space: Theory and measurement,” Discussion paper, National Bureau of Economic Research. 43, 55 





ADAO˜ , R., A. COSTINOT, AND D. DONALDSON (2025): “Putting Quantitative Models to the Test: An Application to Trump’s Trade War [Practitioner’s Guide],” Quarterly Journal of Economics, Forthcoming. 55 





AHLFELDT, G. M., S. J. REDDING, D. M. STURM, AND N. WOLF (2015): “The economics of density: Evidence from the Berlin Wall,” Econometrica, 83(6), 2127–2189. 32, 52, 56, 60 





ALBERT, C., AND J. MONRAS (2022): “Immigration and spatial equilibrium: the role of expenditures in the country of origin,” American Economic Review, 112(11), 3763–3802. 62 





ALDER, S. (2016): “Chinese roads in India: The effect of transport infrastructure on economic development,” Discussion paper, Work. Pap., Univ. North Carolina, Chapel Hill. 56 





ALLEN, T. (2023): “The topography of nations,” Discussion paper, National Bureau of Economic Research. 63 





ALLEN, T., AND C. ARKOLAKIS (2014): “Trade and the Topography of the Spatial Economy,” The Quarterly Journal of Economics, 129(3), 1085–1140. 16, 26, 27, 30, 35, 36, 47, 48, 56 





(2018): “Modern spatial economics: a Primer,” World Trade Evolution, p. 435. 25 





(2022): “The welfare effects of transportation infrastructure improvements,” The Review of Economic Studies, 89(6), 2911–2957. 49, 50, 60 





ALLEN, T., C. ARKOLAKIS, AND X. LI (2024): “On the Equilibrium Properties of Spatial Models,” American Economic Review: Insights, 6(4), 472–89. 30, 39, 41, 42, 62 





ALLEN, T., C. ARKOLAKIS, AND Y. TAKAHASHI (2020): “Universal gravity,” Journal of Political Economy, 128(2), 393–433. 27, 39, 43, 45 





ALLEN, T., D. ATKIN, S. C. CLEVES, AND C. E. HERNANDEZ (2024): “Trucks,” Discussion paper, Working paper. 65 





ALLEN, T., C. DE CASTRO DOBBIN, AND M. MORTEN (2018): “Border walls,” Discussion paper, National Bureau of Economic Research. 59 





ALLEN, T., AND D. DONALDSON (2018): “The geography of path dependence,” Unpublished manuscript. 51, 56, 59, 64 





ALLEN, T., S. FUCHS, S. GANAPATI, A. GRAZIANO, R. MADERA, AND J. MONTORIOL-GARRIGA (2020): “Is tourism good for locals? Evidence from Barcelona,” Dartmouth College, Mimeograph. 55 





ANDERSON, J. E. (1979): “A theoretical foundation for the gravity equation,” The American economic review, 69(1), 106–116. 5, 16 





ANDERSON, J. E., AND E. VAN WINCOOP (2003): “Gravity with Gravitas: A Solution to the Border Puzzle,” American Economic Review, 93(1), 170–192. 5, 22, 27 





ANDERSON, J. E., AND Y. V. YOTOV (2010): “The changing incidence of geography,” American Economic Review, 100(5), 2157–2186. 57 





ARKOLAKIS, C., A. COSTINOT, D. DONALDSON, AND A. RODR´IGUEZ-CLARE (2019): “The elusive procompetitive effects of trade,” The Review of Economic Studies, 86(1), 46–80. 64 





ARKOLAKIS, C., A. COSTINOT, AND A. RODR´IGUEZ-CLARE (2012): “New trade models, same old gains?,” American Economic Review, 102(1), 94–130. 34 





ARKOLAKIS, C., F. ECKERT, AND R. SHI (2017): “Combinatorial discrete choice: A Quantitative Model of Multinational Location,” Available at SSRN 3455353. 63 





ARKOLAKIS, C., S. K. LEE, AND M. PETERS (2020): “European immigrants and the United States’ rise to the technological frontier,” in 2019 Meeting Papers, vol. 1420. 61 





ARKOLAKIS, C., AND C. WALSH (2023): “Clean Growth,” Discussion paper, National Bureau of Economic Research. 57, 61 





ARMINGTON, P. S. (1969): “A Theory of Demand for Products Distinguished by Place of Production,” International Monetary Fund Staff Papers, 16, 159–178. 16 





AUTOR, D. H., D. DORN, AND G. H. HANSON (2013): “The China syndrome: Local labor market effects of import competition in the United States,” American economic review, 103(6), 2121–2168. 55, 59 





BAGAGLI, S. (2023): “The (Express) Way to Segregation: Evidence from Chicago,” . 62 





BALBONI, C., AND J. S. SHAPIRO (2025): “Spatial Environmental Economics,” in Handbook of Regional and Urban Economics, ed. by D. Donaldson, and S. J. Redding, vol. 6. Elsevier North Holland, Amsterdam. 61 





BIFULCO, P., J. GLUCK ¨ , O. KREBS, AND B. KUKHARSKYY (2022): “Single and attractive: Uniqueness and stability of economic equilibria under monotonicity assumptions,” arXiv preprint arXiv:2209.02635. 42, 58 





BLEAKLEY, H., AND J. LIN (2012): “Portage and path dependence,” Quarterly Journal of Economics, 127(2), 587–644. 64 





BORDEU, O. (2023): “Commuting infrastructure in fragmented cities,” Job Market Paper, University of Chicago Booth School of Business, 6. 64 





BRANCACCIO, G., M. KALOUPTSIDI, AND T. PAPAGEORGIOU (2020): “Geography, transportation, and endogenous trade costs,” Econometrica, 88(2), 657–691. 65 





BROWN, D., A. DEARDORFF, AND R. STERN (1991): “A North American free trade agreement: Analytical issues and a computational assessment,” World Scientific Studies in International Economics, p. 393. 54 





BROWN, D. K., AND R. M. STERN (1989): “US-Canada bilateral tariff elimination: The role of product differentiation and market structure,” in Trade policies for international competitiveness, pp. 217–254. University of Chicago Press. 54 





BRYAN, G., AND M. MORTEN (2019): “The aggregate productivity effects of internal migration: Evidence from Indonesia,” Journal of Political Economy, 127(5), 2229–2268. 33, 36 





CAI, S., L. CALIENDO, F. PARRO, AND W. XIANG (2022): “Mechanics of spatial growth,” Discussion paper, National Bureau of Economic Research. 61 





CALIENDO, L., M. DVORKIN, AND F. PARRO (2019): “Trade and labor market dynamics: General equilibrium analysis of the china trade shock,” Econometrica, 87(3), 741–835. 58, 59 





CALIENDO, L., AND F. PARRO (2015): “Estimates of the Trade and Welfare Effects of NAFTA,” The Review of Economic Studies, 82(1), 1–44. 57, 58 





CALIENDO, L., F. PARRO, E. ROSSI-HANSBERG, AND P.-D. SARTE (2018): “The impact of regional and sectoral productivity changes on the US economy,” The Review of economic studies, 85(4), 2042–2096. 57 





CHIOVELLI, G., S. MICHALOPOULOS, AND E. PAPAIOANNOU (2018): “Landmines and spatial development,” Discussion paper, National Bureau of Economic Research. 56 





CHRISTALLER, W. (1933): Die zentralen Orte in S ¨uddeutschland. 63 





CONTE, B., K. DESMET, AND E. ROSSI-HANSBERG (2022): “On the geographic implications of carbon taxes,” Discussion paper, National Bureau of Economic Research. 61 





COSTINOT, A., D. DONALDSON, AND I. KOMUNJER (2012): “What goods do countries trade? A quantitative exploration of Ricardo’s ideas,” The Review of economic studies, 79(2), 581–608. 57 





COSTINOT, A., AND A. RODR´IGUEZ-CLARE (2014): “Trade theory with numbers: Quantifying the consequences of globalization,” in Handbook of international economics, vol. 4, pp. 197–261. Elsevier. 42 





CREWS, L. G. (2023): “A dynamic spatial knowledge economy,” . 61 





CRUZ, J.-L., AND E. ROSSI-HANSBERG (2024): “The economic geography of global warming,” Review of Economic Studies, 91(2), 899–939. 61 





DEKLE, R., J. EATON, AND S. KORTUM (2008): “Global Rebalancing with Gravity: Measuring the Burden of Adjustment,” IMF Staff Papers, 55(3), 511–540. 20, 42 





DESMET, K., D. K. NAGY, AND E. ROSSI-HANSBERG (2018): “The Geography of Development,” Journal of Political Economy, 126(3), 903–983. 59 





DESMET, K., AND F. PARRO (2025): “Spatial Dynamics,” in Handbook of Regional and Urban Economics, ed. by D. Donaldson, and S. J. Redding, vol. 6. Elsevier North Holland, Amsterdam. 20, 59, 62 





DESMET, K., AND E. ROSSI-HANSBERG (2013): “Urban accounting and welfare,” American Economic Review, 103(6), 2296–2327. 54 





DIJKSTRA, E. W. (1959): “A note on two problems in connexion with graphs,” Numerische Mathematik, 1(1), 269–271. 48 





DIX-CARNEIRO, R., AND B. K. KOVAK (2017): “Trade liberalization and regional dynamics,” American Economic Review, 107(10), 2908–2946. 59 





DONALD, E., M. FUKUI, AND Y. MIYAUCHI (2023): “Unpacking Aggregate Welfare in a Spatial Economy,” Discussion paper, Tech. rep. 46, 64 





DONALDSON, D. (2018): “Railroads of the Raj: Estimating the impact of transportation infrastructure,” American Economic Review, 108(4-5), 899–934. 47, 48, 54, 57 





(2025): “Transportation,” in Handbook of Regional and Urban Economics, ed. by D. Donaldson, and S. J. Redding, vol. 6. Elsevier North Holland, Amsterdam. 47 





DONALDSON, D., AND R. HORNBECK (2016): “Railroads and American Economic Growth: A “Market Access” Approach,” Quarterly Journal of Economics, 131(2), 799–858. 27, 31, 36, 47, 52, 55 





EATON, J., AND S. KORTUM (2002): “Technology, geography, and trade,” Econometrica, 70(5), 1741–1779. 5, 31, 32, 33, 36, 58 





ECKERT, F., AND T. KLEINEBERG (2024): “The geography of opportunity: Education, work, and intergenerational mobility across us counties,” Discussion paper, Working Paper. 62 





EGGER, P. H., G. LOUMEAU, AND N. LOUMEAU (2023): “China’s dazzling transport-infrastructure growth: Measurement and effects,” Journal of International Economics, 142, 103734. 56 





FABER, B., AND C. GAUBERT (2019): “Tourism and economic development: Evidence from Mexico’s coastline,” American Economic Review, 109(6), 2245–2293. 57 





FAJGELBAUM, P. D., AND C. GAUBERT (2020): “Optimal spatial policies, geography, and sorting,” The Quarterly Journal of Economics, 135(2), 959–1036. 46, 64 





(2025): “Optimal Spatial Policies,” in Handbook of Regional and Urban Economics, ed. by D. Donaldson, and S. J. Redding, vol. 6. Elsevier North Holland, Amsterdam. 39, 46, 63 





FAJGELBAUM, P. D., C. GAUBERT, N. GORTON, E. MORALES, AND E. SCHAAL (2023): “Political Preferences and the Spatial Distribution of Infrastructure: Evidence from California’s High-Speed Rail,” Discussion paper, National Bureau of Economic Research. 64 





FAJGELBAUM, P. D., P. K. GOLDBERG, P. J. KENNEDY, AND A. K. KHANDELWAL (2020): “The return to protectionism,” The Quarterly Journal of Economics, 135(1), 1–55. 55 





FAJGELBAUM, P. D., AND E. SCHAAL (2020): “Optimal transport networks in spatial equilibrium,” Econometrica, 88(4), 1411–1452. 64 





FALLY, T. (2015): “Structural gravity and fixed effects,” Journal of international economics, 97(1), 76–85. 47 





FAN, J., Y. LU, AND W. LUO (2023): “Valuing domestic transport infrastructure: A view from the route choice of exporters,” Review of Economics and Statistics, 105(6), 1562–1579. 50 





FELIX, M. (2021): “Trade, labor market concentration, and wages,” Job Market Paper. 64 





FOGEL, R. W. (1964): Railroads and American Economic Growth: Essays in Econometric History. Baltimore: Johns Hopkins University Press. 25 





FROBENIUS, FERDINAND GEORG, G. (1912): “Uber Matrizen aus nicht negativen Elementen,” . 24 ¨ 





FUCHS, S., AND W. F. WONG (2022): “Multimodal transport networks,” Discussion paper, Working Paper. 50 





FUJITA, M., P. KRUGMAN, AND A. J. VENABLES (1999): The Spatial Economy: Cities, Regions, and International Trade. MIT Press. 10, 14, 15, 30 





GARG, T. (2025): “Can Industrial Policy Overcome Coordination Failures? Theory and Evidence,” Discussion paper. 42, 64 





GIROUD, X., S. LENZU, Q. MAINGI, AND H. MUELLER (2024): “Propagation and amplification of local productivity spillovers,” Econometrica, 92(5), 1589–1619. 61 





GLAESER, E. L., AND J. D. GOTTLIEB (2009): “The wealth of cities: Agglomeration economies and spatial equilibrium in the United States,” Journal of economic literature, 47(4), 983–1028. 8, 51 





GREENSTONE, M., R. HORNBECK, AND E. MORETTI (2010): “Identifying agglomeration spillovers: Evidence from winners and losers of large plant openings,” Journal of political economy, 118(3), 536–598. 57 





GROSSMAN, G. M., AND E. ROSSI-HANSBERG (2010): “External economies and international trade redux,” The Quarterly Journal of Economics, 125(2), 829–858. 16 





GUTIERREZ ´ , A. (2022): “Labor market power and the pro-competitive gains from trade,” Unpublished Working Paper. 64 





HANSON, G. H. (2005): “Market potential, increasing returns and geographic concentration,” Journal of international economics, 67(1), 1–24. 55 





HEAD, K., AND T. MAYER (2014): “Gravity equations: Workhorse, toolkit, and cookbook,” in Handbook of international economics, vol. 4, pp. 131–195. Elsevier. 7 





HEBLICH, S., S. J. REDDING, AND D. M. STURM (2020): “The making of the modern metropolis: evidence from London,” The Quarterly Journal of Economics, 135(4), 2059–2133. 54, 56 





HELPMAN, E. (1998): “The Size of Regions,” in Topics in Public Economics, ed. by E. S. David Pines, and I. Zilcha, pp. 33–54. Cambridge University Press. 26, 27, 31, 36 





HENDERSON, J. V. (1974): “The Sizes and Types of Cities,” The American Economic Review, 64(4), 640–656. 8 





HOELZLEIN, M. (2019): “Two-sided sorting and spatial inequality in cities,” Unpublished manuscript. 62 





HONG, G. (2024): “Two-sided sorting of workers and firms: Implications for spatial inequality and welfare,” Discussion paper, Canadian Labour Economics Forum (CLEF), University of Waterloo. 64 





HULTEN, C. R. (1978): “Growth accounting with intermediate inputs,” The Review of Economic Studies, 45(3), 511–518. 25 





JACKSON, M. O., ET AL. (2008): Social and economic networks, vol. 3. Princeton university press Princeton. 24 





JAWORSKI, T., AND C. T. KITCHENS (2019): “National policy for regional development: Historical evidence from Appalachian highways,” Review of Economics and Statistics, 101(5), 777–790. 47 





KANG, E. (2024): “Global Fisheries: Quantifying the Externalities from Open Access,” Job Market Paper, University of Michigan. 61 





KARLIN, S., AND L. NIRENBERG (1967): “On a theorem of P. Nowosad,” Journal of Mathematical Analysis and Applications, 17(1), 61–67. 28, 29, 42 





KEHOE, T. J. (2005): “An Evaluation of the Performance of Applied General Equilibrium Models on the Impact of NAFTA,” in Frontiers in Applied General Equilibrium Modeling: Essays in Honor of Herbert Scarf, ed. by T. J. Kehoe, T. Srinivasan, and J. Whalley, chap. 13, pp. 344–77. Cambridge University Press. 54 





KEHOE, T. J., P. S. PUJOLAS, AND J. ROSSBACH (2017): “Quantitative trade models: Developments and challenges,” Annual Review of Economics, 9(1), 295–325. 54 





KLEINMAN, B., E. LIU, AND S. J. REDDING (2023): “Dynamic Spatial General Equilibrium,” Econometrica, 91(2), 385–424. 59 





(2024a): “International Friends and Enemies,” American Economic Journal: Macroeconomics, 16(4), 350–85. 43, 45 





KLEINMAN, B., E. LIU, AND S. J. REDDING (2024b): “The linear algebra of economic geography models,” in AEA Papers and Proceedings, vol. 114, pp. 328–333. American Economic Association 2014 Broadway, Suite 305, Nashville, TN 37203. 43 





KRUGMAN, P. (1991): “Increasing returns and economic geography,” Journal of political economy, 99(3), 483–499. 2, 7, 10, 12, 16, 26, 27, 33, 35, 36 





KRUGMAN, P., AND A. J. VENABLES (1995): “Globalization and the Inequality of Nations,” Quarterly Journal of Economics, 110(4), 857–880. 33, 58 





KUCHERYAVYY, K., G. LYN, AND A. RODR´IGUEZ-CLARE (2023): “Grounded by gravity: A well-behaved trade model with industry-level economies of scale,” American Economic Journal: Macroeconomics, 15(2), 372–412. 16 





(2024): “Spatial equilibria: The case of two regions,” Journal of International Economics, 152, 104008. 30 





LEVCHENKO, A. A., AND J. ZHANG (2016): “The evolution of comparative advantage: Measurement and welfare implications,” Journal of Monetary Economics, 78, 96–111. 57 





LIANG, J. (2024): “Knowledge and Firm Growth in Space,” mimeo, Yale University. 61 





LONG JR, J. B., AND C. I. PLOSSER (1983): “Real business cycles,” Journal of political Economy, 91(1), 39–69. 57 





MARSHALL, A. (1890): Principles of Economics, vol. 1. Macmillan, London. 1 





MISES, R., AND H. POLLACZEK-GEIRINGER (1929): “Praktische Verfahren der Gleichungsauflosung.,” ¨ ZAMM-Journal of Applied Mathematics and Mechanics/Zeitschrift f ¨ur Angewandte Mathematik und Mechanik, 9(1), 58–77. 24 





MONTE, F., C. PORCHER, AND E. ROSSI-HANSBERG (2023): “Remote work and city structure,” Discussion paper, National Bureau of Economic Research. 64 





MONTE, F., S. J. REDDING, AND E. ROSSI-HANSBERG (2018): “Commuting, migration and local employment elasticities,” American Economic Review, 108(12), 3855–90. 57, 60 





MORETTI, E. (2011): “Chapter 14 - Local Labor Markets,” vol. 4 of Handbook of Labor Economics, pp. 1237–1313. Elsevier. 8, 16 





NAGY, D. K. (2023): “Hinterlands, city formation and growth: Evidence from the US westward expansion,” Review of Economic Studies, 90(6), 3238–3281. 54, 63 





OBERFIELD, E., E. ROSSI-HANSBERG, P.-D. SARTE, AND N. TRACHTER (2024): “Plants in space,” Journal of Political Economy, 132(3), 867–909. 63 





OSSA, R. (2015): “Why trade matters after all,” Journal of International Economics, 97(2), 266–277. 34, 57 





OUAZAD, A. C. (2024): “Equilibrium Multiplicity: A Systematic Approach using Homotopies, with an Application to Chicago,” arXiv preprint arXiv:2401.10181. 42, 64 





OWENS III, R., E. ROSSI-HANSBERG, AND P.-D. SARTE (2020): “Rethinking detroit,” American Economic Journal: Economic Policy, 12(2), 258–305. 64 





PELLEGRINA, H. S., AND S. SOTELO (2021): “Migration, Specialization, and Trade: Evidence from Brazil’s March to the West,” Discussion paper, National Bureau of Economic Research. 61, 62 





PERRON, O. (1907): “Zur theorie der matrices,” Mathematische Annalen, 64(2), 248–263. 24 





REDDING, S., AND A. J. VENABLES (2004): “Economic geography and international inequality,” Journal of international Economics, 62(1), 53–82. 11, 55 





REDDING, S. J. (2016): “Goods trade, factor mobility and welfare,” Journal of International Economics, 101, 148–167. 32, 36 





(2025): “Quantitative Urban Economics,” in Handbook of Regional and Urban Economics, ed. by D. Donaldson, and S. J. Redding, vol. 6. Elsevier North Holland, Amsterdam. 3, 32, 56, 60 





REDDING, S. J., AND D. M. STURM (2008): “The Costs of Remoteness: Evidence from German Division and Reunification,” American Economic Review, 98(5), 1766–1797. 36, 57 





ROBACK, J. (1982): “Wages, rents, and the quality of life,” The Journal of Political Economy, pp. 1257–1278. 2, 7, 8, 9, 10, 16, 35, 36 





ROSEN, S. (1979): “Wage-Based Indexes of Urban Quality of Life,” in Current Issues in Urban Economics, ed. by P. Mieszkowski, and M. Straszheim, pp. 74–104. Johns Hopkins University Press. 7, 8, 9, 10, 16, 35, 36 





ROY, A. D. (1951): “Some thoughts on the distribution of earnings,” Oxford economic papers, 3(2), 135–146. 62 





SAMUELSON, P. A. (1954): “The transfer problem and transport costs, II: Analysis of effects of trade impediments,” The Economic Journal, 64(254), 264–289. 17 





SETHIAN, J. (1996): “A Fast Marching Level Set Method for Monotonically Advancing Fronts,” Proceedings of the National Academy of Sciences, 93(4), 1591–1595. 48 





SEVEREN, C. (2023): “Commuting, labor, and housing market effects of mass transportation: Welfare and identification,” Review of Economics and Statistics, 105(5), 1073–1091. 47 





SILVA, J. S., AND S. TENREYRO (2006): “The log of gravity,” The Review of Economics and statistics, 88(4), 641–658. 47 





SMITH, A. (1776): “An inquiry into the nature and causes of the wealth of nations: Volume One,” London: printed for W. Strahan; and T. Cadell, 1776. 1 





TOMBE, T., AND X. ZHU (2019): “Trade, migration, and productivity: A quantitative analysis of china,” American Economic Review, 109(5), 1843–1872. 59 





TSITSIKLIS, J. (1995): “Efficient Algorithms for Globally Optimal Trajectories,” Automatic Control, IEEE Transactions on, 40(9), 1528–1538. 48 





TSIVANIDIS, N. (2022): “Evaluating the impact of urban transit infrastructure: Evidence from Bogota’s transmilenio,” Unpublished manuscript, 18. 62 





VAHRENKAMP, R. (1976): “Derivatives of the dominant root,” Applied Mathematics and Computation, 2(1), 29–39. 24 





VOGEL, K. B., G. H. HANSON, A. KHANDELWAL, C. LIU, AND H. PARK (2024): “Using Satellite Imagery to Detect the Impacts of New Highways: An Application to India,” Discussion paper, National Bureau of Economic Research. 56 





VON THUNEN ¨ , J. (1826): Der Isolierte Staat in Beziehung auf Landwirtschaft und National¨okonomie. Puthes, Hamburg, The Isolated State, Oxford: Pergamon, 1966. 1 





WEIWU, L. (2023): “Unequal Access: Racial Segregation and the Distributional Impacts of Interstate Highways in Cities,” Discussion paper, Technical report, MIT. 62 





WONG, W. F. (2022): “The round trip effect: Endogenous transport costs and international trade,” American Economic Journal: Applied Economics, 14(4), 127–166. 65 





YAGAN, D. (2019): “Employment hysteresis from the great recession,” Journal of Political Economy, 127(5), 2505–2558. 59 





ZARATE ´ , R. D. (2022): “Spatial misallocation, informality, and transit improvements: Evidence from mexico city,” . 62 

