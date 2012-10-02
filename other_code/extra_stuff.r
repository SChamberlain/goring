############################ Goring paper project
# ## plot with ggphylo
# install_github("ggphylo", "gjuggler")
# require(ggphylo)
# ggphylo(randtree)

### Format names
# library(gdata); library(stringr)
# dat <- read.xls("taxon.table.pollenequiv_acronyms.xlsx", sheet=2)
# head(dat); str(dat)
# dat$genus_sp <- tolower(dat$genus_sp) # lowercase
# dat$genus <- tolower(dat$genus) # lowercase
# dat$family <- tolower(dat$family) # lowercase
# dat$genus_sp <- str_replace(dat$genus_sp, "_$", "") # replace trailing "_"
# write.csv(dat, "dat_fixed.csv")

# ######### Bladj, to get branch lengths
# setwd("/Users/ScottMac/Downloads/phylocom-4.2/mac") # to directory where phylocom is
# write.tree(tree, "phylo") # write tree to directory where Phylocom executable lives
# system("./phylocom bladj > phyloout.txt") # run Bladj
# tree_bladjed <- read.tree("phyloout.txt") # read in tree with branch lengths
# plot(extract.clade(tree_bladjed, "alismatales"), cex=0.6, no.margin=T)
