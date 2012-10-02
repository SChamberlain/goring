### Writing plot and tree data to file for each plot (if needed)


######## Data 
# ## get each communities data
# setwd("/Mac/Projects/Goring")
# # tt <- read.csv("CleanedPresAbs.csv", nrows=10)[,1:2]
# matchspnames <- read.xls("taxon.table.pollenequiv_acronyms.xlsx", sheet="matchto")
# matchspnames <- matchspnames[,c("genus_sp","spec.code")]
# head(matchspnames)
# 
# numrows <- nrow(read.csv("CleanedPresAbs.csv"))
# header <- names(read.csv("CleanedPresAbs.csv", nrows=1))
# setwd("/Mac/Projects/Goring/eachplot")
# for(i in 1:numrows) {
# # 	i <- 1
# 	dd <- read.csv("/Mac/Projects/Goring/CleanedPresAbs.csv", skip=i, nrows=1, header=F)
# 	plotname <- dd[1,1]
# 	names(dd) <- header
# 	ee <- dd[,-1]
# 	if(sum(as.numeric(ee)) == 0){
# 		message(paste0("no species in ", plotname))
# 		NULL
# 	} else
# 	{
# 		ff <- t(ee)
# 		dat_present <- ff[ff > 0 , ]
# 		
# 		if(length(dat_present) < 3) {
# 			message(paste0("less than 3 species in ", plotname, " just wrote data file"))
# 			write.csv(data.frame(X = plotname, t(dat_present)), paste0(plotname, "_data.csv"), row.names=F)
# 		} else{
# 			sp_present <- names(dat_present)
# 			getnames <- matchspnames[matchspnames$spec.code %in% sp_present , "genus_sp"]
# 			removetips <- tree$tip.label[!tree$tip.label %in% getnames]
# 			tree2 <- drop.tip(tree, removetips)	
# 			message(paste0("writing ", plotname))
# 			write.tree(tree2, paste0(plotname, "_tree.txt"))
# 			write.csv(data.frame(X = plotname, t(dat_present)), paste0(plotname, "_data.csv"), row.names=F)
# 		}
# 	}
# }
# 
# # Use the data and/or tree from a plot
# setwd("/Mac/Projects/Goring/eachplot")
# plot(read.tree(file="2883_tree.txt"))
# read.csv("2883_data.csv")