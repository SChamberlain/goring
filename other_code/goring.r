########################################################################
########################################################################
# Read in tree with branch lengths
require(ape); require(gdata); library(picante); library(apTreeshape); library(vegan)
setwd("/Mac/Projects/Goring")
tree <- read.tree(file="/Mac/Projects/Goring/goring_tree_withbranchlengths.txt")

# Get data.frame to do name matching between acronyms and full species names
setwd("/Mac/Projects/Goring")
matchspnames <- read.xls("taxon.table.pollenequiv_acronyms.xlsx", sheet="matchto")
matchspnames <- matchspnames[,c("genus_sp","spec.code")]

numrows <- nrow(read.csv("CleanedPresAbs.csv")) # number of rows to go through
header <- names(read.csv("CleanedPresAbs.csv", nrows=1)) # header needed for loop

# Function to calculate phylogeny shape metrics on each community
calcmets <- function(phy, commdat) {

	# Numer of polytomies in tree
	numpolys <- function(x) {
		degree <- tabulate(x$edge[, 1])
		length(which(degree > 2))
	}
	numpols <- numpolys(phy)
	
	# Faith's diversity - using picante fxn pd() has to use a matrix as input for comm. data
	# CAN'T CALCULATE THIS - NO ABUNDANCE DATA
# 	faithdiv <- pd(t(data.frame(commdat)), phy)[,1]
	
	# Taxonomic distinctiveness
	mod <- taxondive(t(data.frame(commdat)), cophenetic.phylo(phy))
	ttd <- mod$Lambda[[1]]
	
	# Tree size
	numtips <- Ntip(phy)
	
	# Number of nodes
	numnodes <- Nnode(phy)
	
	# Colless tree shape metric
	if(is.binary.tree(phy)){phy <- phy} else {phy <- multi2di(phy)}
	phy_ <- as.treeshape(phy) # convert to apTreeshape format
	collessmet <- colless(phy_, norm="yule")
	sackinmet <- sackin(phy_, norm="yule")
	
	# Gamma statistic of Pybus and Harvey
	gamma <- gammaStat(phy)

	# mpd and mntd
  mntd <- mntd(t(commdat), cophenetic(phy))
  mpd <- mpd(t(commdat), cophenetic(phy))
  
  data.frame(numpols=numpols, ttd=ttd, numtips=numtips, numnodes=numnodes, collessmet=collessmet, 
             sackinmet=sackinmet, gamma=gamma, mntd=mntd, mpd=mpd)
}

# Create an empty data.frame to put data in to
outdf <- data.frame(plotname=NA, numpols=NA, ttd=NA, numtips=NA, numnodes=NA, 
										collessmet=NA, sackinmet=NA, gamma=NA, mntd=NA, mpd=NA)

# The for loop
for(i in 1:numrows) {
# 	i <- 1
	dd <- read.csv("/Mac/Projects/Goring/CleanedPresAbs.csv", skip=i, nrows=1, header=F)
	plotname <- as.character(dd[1,1])
	names(dd) <- header
	ee <- dd[,-1]
	if(sum(as.numeric(ee)) == 0){
		message(paste0("no species in ", plotname, " writing NA's to data.frame"))
		dddd <- data.frame(plotname=as.character(plotname), rep(NA, ncol(outdf)-1))
		outdf[i,] <- dddd
	} else
	{
		ff <- t(ee)
		dat_present <- ff[ff > 0 , ]
		if(length(dat_present) < 3) {
			message(paste0("less than 3 species in ", plotname, " inserting NA's to data.frame"))
			dddd <- data.frame(plotname=as.character(plotname), rep(NA, ncol(outdf)-1))
			outdf[i,] <- dddd
		} else{
			sp_present <- names(dat_present)
			getnames <- matchspnames[matchspnames$spec.code %in% sp_present , "genus_sp"]
			names(dat_present) <- getnames
			removetips <- tree$tip.label[!tree$tip.label %in% getnames]
			tree2 <- drop.tip(tree, removetips)
			dat_present <- dat_present[names(dat_present) %in% tree2$tip.label]
			message(paste0("writing ", plotname, " data in to data.frame"))
			ddd <- calcmets(tree2, dat_present)
			dddd <- data.frame(plotname=plotname, ddd)
			dddd$plotname <- as.character(dddd$plotname)
			outdf[i,] <- dddd
		}
	}
}

# The output data.frame
head(outdf); str(outdf)
write.csv(outdf, "goring_tree_metrics.csv")