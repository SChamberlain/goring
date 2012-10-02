#### Calculate phylogenetic tree metrics on trees for each plot
# Read in tree with branch lengths
require(ape); require(gdata); library(picante); library(apTreeshape); library(vegan); library(plyr)
library(foreach); library(multicore); library(doMC)

# Load in tree
tree <- read.tree(file="goring_tree_withbranchlengths.txt")

# Get data.frame to do name matching between acronyms and full species names
matchspnames <- read.xls("taxon.table.pollenequiv_acronyms.xlsx", sheet="matchto")
matchspnames <- matchspnames[,c("genus_sp","spec.code")]

# Load plot data
dat <- read.csv("CleanedPresAbs.csv") # read in data file
numrows <- nrow(dat) # number of rows to go through
header <- names(dat) # header needed for loop

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
  #   faithdiv <- pd(t(data.frame(commdat)), phy)[,1]
  
  # Taxonomic distinctiveness
  # THIS IS FUNKY, NEED TO FIGURE OUT HOW TO USE
  #   taxdis <- taxa2dist(dune.taxon, varstep=TRUE)
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
  
  data.frame(numpols, ttd, numtips, numnodes, collessmet, sackinmet, gamma, mntd, mpd)
}


# function to run metrics calculations on each tree for each plot
runmets <- function(x) {
  dd <- dat[x, ]
  plotname <- as.character(dd[1,1])
  names(dd) <- header
  ee <- dd[,-1]
  if(sum(as.numeric(ee)) == 0){
    message(paste0("no species in ", plotname, " writing NA's to data.frame"))
    dddd <- data.frame(as.character(plotname), t(data.frame(rep(NA, ncol(outdf)-1))))
    names(dddd) <- NULL; row.names(dddd) <- NULL
    write.table(ddd, file = "goringdfout.txt", append=T, row.names=F, col.names=F)
  } else
  {
    ff <- t(ee)
    dat_present <- ff[ff > 0 , ]
    if(length(dat_present) < 3) {
      message(paste0("less than 3 species in ", plotname, " inserting NA's to data.frame"))
      dddd <- data.frame(as.character(plotname), t(data.frame(rep(NA, ncol(outdf)-1))))
      names(dddd) <- NULL; row.names(dddd) <- NULL
      write.table(ddd, file = "goringdfout.txt", append=T, row.names=F, col.names=F)
    } else{
      sp_present <- names(dat_present)
      getnames <- matchspnames[matchspnames$spec.code %in% sp_present , "genus_sp"]
      names(dat_present) <- getnames
      if(!length(unique(getnames)) == length(getnames)){
        bbb <- data.frame(names_ = getnames, dat = dat_present)
        bbb_ <- ddply(bbb, .(names_), summarise, length(names_))
        dat_present <- rep(1, nrow(bbb_))
        names(dat_present) <- bbb_$names_
      } else
      {
        removetips <- tree$tip.label[!tree$tip.label %in% unique(getnames)]
        tree2 <- drop.tip(tree, removetips)
        dat_present <- dat_present[names(dat_present) %in% tree2$tip.label]
        message(paste0("writing ", plotname, " data in to data.frame"))
        ddd <- calcmets(tree2, dat_present)
        dddd <- data.frame(plotname=plotname, ddd)
        dddd$plotname <- as.character(dddd$plotname)
        names(dddd) <- NULL
        write.table(dddd, file = "goringdfout.txt", append=T, row.names=F, col.names=F)
      }
    }
  }
}

# Create an empty data.frame to put data in to
outdf <- data.frame(plotname=NA, numpols=NA, ttd=NA, numtips=NA, numnodes=NA, 
                    collessmet=NA, sackinmet=NA, gamma=NA, mntd=NA, mpd=NA)

# Create empty data.frame to hold results
write.table(outdf, "goringdfout.txt", row.names=F) 

registerDoMC(cores=8)
outout <- llply(1:numrows, runmets, .parallel=T)

dat <- read.table("goringdfout.txt", header=T)[-1,]

write.csv(dat, "goring_metrics.csv")