#LOADING LIBRARIES
#################
library(partykit) #conditional inference trees
library(ggplot2)
library(tidyr)


#TIDYING UP DATA
#################
#importing data
basa <- read.table("data/basa_char_par.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
marbore <- read.table("data/marbore_char_par.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

#target columns
target.basa <- c("cal.BP", "CHAR", "Pinus", "Abies", "Betula", "Corylus", "herb", "Dec_Querc")
target.marbore <- c("cal.BP", "CHAR", "Pinus", "Abies", "Betula", "Corylus", "herb", "Dec_Querc")

#checking match between target and actual columns
if(sum(target.basa %in% colnames(basa)) != length(target.basa)){message("some target columns are missing in basa")}
if(sum(target.marbore %in% colnames(marbore)) != length(target.marbore)){message("some target columns are missing in marbore")}

#subsetting data
basa <- basa[, target.basa]
marbore <- marbore[, target.marbore]

#new names
colnames(basa) <- c("age", "char", "pinus", "abies", "betula", "corylus", "herbs", "quercus")
colnames(marbore) <- c("age", "char", "pinus", "abies", "betula", "corylus", "herbs", "quercus")

#plot df
basa.plot <- basa
basa.plot$site <- "Basa de la Mora"
marbore.plot <- marbore
marbore.plot$site <- "Marboré"
plot.df <- rbind(basa.plot, marbore.plot)
rm(basa.plot, marbore.plot)
plot.df <- tidyr::gather(plot.df, variable, value, 2:8)
  ggplot(data = plot.df, aes(x = age, y = value, group=variable, color=variable)) +
  geom_line() +
  facet_wrap("site", ncol = 1)


#ADDING LAG (each sample is aligned with the previous one)
##########################################################

#remove first case of data
basa.lagged <- basa[-1, ]
marbore.lagged <- marbore[-1, ]

#remove last case of data
basa <- basa[-nrow(basa), ]
marbore <- marbore[-nrow(marbore), ]

#new names the lagged data using "antecedent" as suffix
colnames(basa.lagged) <- c("age.antecedent", "char.antecedent", "pinus.antecedent", "abies.antecedent", "betula.antecedent", "corylus.antecedent", "herbs.antecedent", "quercus.antecedent")
colnames(marbore.lagged) <- c("age.antecedent", "char.antecedent", "pinus.antecedent", "abies.antecedent", "betula.antecedent", "corylus.antecedent", "herbs.antecedent")

#merging the data
basa <- cbind(basa, basa.lagged)
marbore <- cbind(marbore, marbore.lagged)

#computing age difference
basa$age.difference <- abs(basa$age - basa$age.antecedent)
marbore$age.difference <- abs(marbore$age - marbore$age.antecedent)

#into list
sites <- list()
sites$basa <- basa
sites$marbore <- marbore

#removing temporary datasets
rm(basa.lagged, marbore.lagged, basa, marbore)


#RECURSIVE PARTITION MODELING
###############################

#modeling formulas
###############################
formulas <- list()
formulas$pinus <- as.formula(paste("pinus ~ pinus.antecedent + char + char.antecedent"))
formulas$abies <- as.formula(paste("abies ~ abies.antecedent + char + char.antecedent"))
formulas$betula <- as.formula(paste("betula ~ betula.antecedent + char + char.antecedent"))
formulas$corylus <- as.formula(paste("corylus ~ corylus.antecedent + char + char.antecedent"))
formulas$herbs <- as.formula(paste("herbs ~ herbs.antecedent + char + char.antecedent"))
formulas$quercus <- as.formula(paste("quercus ~ quercus.antecedent + char + char.antecedent"))

#fitting models for site "basa"
###############################
models <- list()
models$basa.pinus <- ctree(formulas$pinus, data=sites$basa)
models$basa.abies <- ctree(formulas$abies, data=sites$basa)
models$basa.betula <- ctree(formulas$betula, data=sites$basa)
models$basa.corylus <- ctree(formulas$corylus, data=sites$basa)
models$basa.herbs <- ctree(formulas$herbs, data=sites$basa)
models$basa.quercus <- ctree(formulas$quercus, data=sites$basa)

#fitting models for site "marbore"
###############################
models$marbore.pinus <- ctree(formulas$pinus, data=sites$marbore)
models$marbore.abies <- ctree(formulas$abies, data=sites$marbore)
models$marbore.betula <- ctree(formulas$betula, data=sites$marbore)
models$marbore.corylus <- ctree(formulas$corylus, data=sites$marbore)
models$marbore.herbs <- ctree(formulas$herbs, data=sites$marbore)


#PLOTTING MODELS
################################
#name substitutions
site.names <- c("Basa de la Mora", "Marboré")
names(site.names) <- c("basa", "marbore")

taxon.names <- c("Pinus", "Abies", "Corylus", "Betula", "Herbs", "Deciduous Quercus")
names(taxon.names) <- c("pinus", "abies", "corylus", "betula", "herbs", "quercus")

#opening pdf and iterating through models
pdf(file = "output/ctree_marbore_basa.pdf", width=12, height=9, pointsize=15)
for(i in 1:length(models)){

  #model and site name name
  site.taxon <- unlist(strsplit(x=names(models)[i], split=".", fixed = TRUE))
  site <- site.taxon[1]
  taxon <- site.taxon[2]

  #substitution of plain names
  site.name <- site.names[names(site.names)==site]
  taxon.name <- taxon.names[names(taxon.names)==taxon]

  #plot
  plot(models[[i]], main=paste(site.name, taxon.name, sep=" - "))

}
dev.off()


#SAVING OUTPUT
####################################
save(formulas, models, sites, file = "output_marbore_basa.RData")
