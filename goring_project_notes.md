# Scott's notes

+ rarefaction of richness in pollen grains counted and pollen taxa found
+ grasses are a problem: can have lots of grass species in one site, e.g.
+ To Do
	+ Build BC phylogenetic tree
	+ comparing trees in space, mantel tests, or clustering trees?
	+ check clustering trees method
	+ 
+ hypothesis: immediating following glaciation, first colonists are a select group, can migrate rapidly, very stress-tolerant, so get things like Pinus, then get infilling (e.g., understorey communities) but also during this time gradually increasing precipitation
	+ early, would have long branch distances among species in the community
	+ then filling in the species in the tree
	+ the big unknown is entomophilous species, without having them in the pollen dataset...
	+ things stabilized at about 6000 years before present
	+ climatic stability should effect species composition and perhaps phylogeneitc structure
		+ 
	+ "test relationships between climatic change and phylogenetic structure"
	+ two different regions in BC, similar dominant canopy species, 
		+ "coastal": should be about 6000 years old
		+ interior wet belt: should be about 1000 years old probably
	+ the two regions probably wouldn't end up the same way
		+ coastal much warmer, interior much colder
		+ longer degree growing days on coast
	+ question: is species richness associated with metrics of phylogenetic structure?
	+ climatic equitability: you get real stratification of ecozones, so more difficult for species to move across them.

+ Questions:
	+ Is species richness related to phylogenetic tree shape metrics?
		+ be careful that metrics aren't scaling with tree size by default
	+ Is continentality (diff. between coldest and warmest month of year) related to phylogenetic measures?
		+ GAM regression fits to data...

+ Phylogeny (DEADLINE: before Oct. 1st meeting)
	+ one phylogeny for whole species pool
	+ then a phylogeny for each community
	+ missing species - change Family names? 
	+ wants branch lengths - use bladj
	+ polytomies - resolve them using Arne code -> ran into problems, Beast keeps erroring saying XML is not formatted correctly, emailed Arne about the problem
	+ phylogeny data: 
		+ just read in whole file instead of line by line

+ Oct. 1st meeting
	+ me: think of how we can use the data
		+ how does balance of trees change through time
		+ does precipitation change tree shape, etc.?
	+ open science, run the paper as totally open science exercise

+ NDVI - simon will get. 

+ Simon will ask Jana about how to define the regional species pool

+ Randomization - attempt to weight by abundance, etc. How do we do the randomizations?

+ Tree shape metrics
	+ Just use the multi2di tree for the ones that need dichotomous trees, fix soon!