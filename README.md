# coral_fish

This project models the functional consequences of range-shifting tropical fish

### core scripts

**fish_trait_master.R** handles formatting and updates from the global tropical fish triat database (M Beger)<br />
**fish_sites_master.R** handles survey data from Japan and Australia: summarising and correcting biomass (then clustering sites based on species biomass); running species accumulation curves and estimating species richness<br />
**optimal_clusters.R** explores optimal clustering soution of trait data calling clVal.R function
**clVal.R** performs cluster validation based on n iteration bootstrapping (set at 95% subsample) or trait dendrogram, selection of optimal clusters based on average silhouette width, Jaccard stability index and Rand matching index<br />
**tropicalization_paper.R** performs core paper analyses
