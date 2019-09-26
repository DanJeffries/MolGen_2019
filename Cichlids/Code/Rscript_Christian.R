library(hierfstat)
library(vcfR)
library(adegenet)
library(outliers)

setwd("/Users/cdeguttr/Desktop/MolGen_2019-master/Cichlids/Data/")

dta_pop <- read.delim("pop.txt")

### Open the hierfstqt object ###
hierf_stat_cichlid <- read.fstat.data("cichlid_hierfstat.dat")
### Add population info ###
popcodes <- read.csv("pop.txt", sep = "\t", header = F)
hierf_stat_cichlid$Pop <- popcodes$V2
### Compute observed and expected heterozygosity and basic population statistics ###
stat <- basic.stats(hierf_stat_cichlid)

### Calculate the mean Ho and Hs for the 4 populations and plot the results ###
mean1 <- colMeans(stat$Ho,na.rm = T)
mean2 <- colMeans(stat$Hs,na.rm = T)
plot(mean1, mean2)
abline(0,1)

### Fst between pop ###
dta4 <- pairwise.WCfst(hierf_stat_cichlid, diploid = T)

### PCA per individual using Hierfstat ###
pca1 <- indpca(hierf_stat_cichlid, ind.labels=hierf_stat_cichlid[,1])
plot(pca1, col=hierf_stat_cichlid[,1], cex=0.7)

### PCA per individual with ADEGENET ###
### Load data from PLINK .raw format ###
data <- read.PLINK("plink.raw")
data
### Add population names to the data ###
pops <- dta_pop$STRATA
pop(data) <- pops
### Run PCA analysis selecting the first 2 axis given that they explain more than 80% of variability ###
pca1 <- glPca(data, parallel=TRUE, n.cores=5)
### Plot the PCA results ###
myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=1)
abline(h=0,v=0, col="grey")

### Per locus Pairwise Fst ### Optional exercise, we provided the Fst file during the course ###
attach(jerome_stat$hierfstat)
var_dat <- varcomp.glob(data.frame(POP_ID), hierf_stat_cichlid$hierfstat[,-c(1)])
detach(jerome_stat$hierfstat)

### Create a Manhattan plot ###
### Load the dataset with Fst values ###
Fst_cichlids_2_pop <- read.delim("Di_1-Di_2_global_mod.tsv")

plot(Fst_cichlids_2_pop$Fst, pch = 20, col = c("black", "grey70")[(as.numeric(as.factor(Fst_cichlids_2_pop$Chr))%%2)+1],
     xlab= "Genomic location (Scaffold)", ylab = "Fst", cex = 0.8, xaxt = "n")


### Detecting outliers in the Fst output ###
### First we can plot the data and detect them visually. In these case the outliers are expectet on hte top part of the graph ###
plot(Fst_cichlids_2_pop$Fst, type="b")
### This function returns the most extreme value: largest diff. between this point and the mean of the data ###
out <- outlier(Fst_cichlids_2_pop$Fst) 
grubbs.test(Fst_cichlids_2_pop$Fst)### if it is an outlier the test is significant!
### I can remove this individual and keep doing the same analysis until the P-value of the test is not significant



### In this dataset only Fst == 0 are significant ###

### Plot with the outliers test significances line ###
plot(Fst_cichlids_2_pop$Fst, pch = 20, col = c("black", "grey70")[(as.numeric(as.factor(Fst_cichlids_2_pop$Chr))%%2)+1],
     xlab= "Genomic location (Scaffold)", ylab = "Fst", cex = 0.8, xaxt = "n")
abline(h=0.375, lty = 4)












