## Load packages

library(adegenet)
library(pegas)
library(hierfstat)

## set working directory
setwd("~/Data/MolGen_2019/Frogs/Data/")

## import from fstat

### TO DO ### Change the locus names

hyla_genind <- read.fstat.data("hyla_FSTAT.dat")

# add the population info 

popcodes <- read.csv("populations_numeric_codes.txt", sep = "\t", header = F)
hyla_genind$Pop <- popcodes$V2

## summary hierfstat

hyla_summ_hier <- basic.stats(hyla_genind)

## Explore hierfsat summary stats

names(hyla_summ_hier)

hyla_summ_hier$n.ind.samp  ## per population and per locus genotype numbers
hyla_summ_hier$pop.freq  ## per populations and per locus allele frequencies
hyla_summ_hier$Ho  ## per population and per locus observed heterozygosity
hyla_summ_hier$Hs  ## per population and per locus expected heterozygosity
hyla_summ_hier$Fis  ## per pop and per locus Fis
hyla_summ_hier$perloc  ## per locus summary stats across all pops
hyla_summ_hier$overall  ## Overall stats for the whole dataset

## Plot Ho vs Hs

plot(hyla_summ_hier$perloc$Ho, 
     hyla_summ_hier$perloc$Hs,
     pch = 16,
     cex = 0.5,
     xlab="Ho",
     ylab="He",
     main="Ho vs. He - Per locus")

abline(0,1)

## Plot by population

plot(colMeans(hyla_summ_hier$Ho, na.rm = T), 
     colMeans(hyla_summ_hier$Hs, na.rm = T),
     pch = 16,
     cex = 0.75,
     xlab="Ho",
     ylab="He",
     main="Ho vs. He - populations")

text(colMeans(hyla_summ_hier$Ho, na.rm = T), 
     colMeans(hyla_summ_hier$Hs, na.rm = T), 
     unique(hyla_genind$Pop), 
     pos = 1)

abline(0,1)

?plot

## Test isolation by distance

## adegenet?

hyla_genind_adegenet <- adegenet::read.fstat("hyla_FSTAT.dat")
pop(hyla_genind_adegenet) <- popcodes$V2

Hyla_lat_long <- read.csv("Hyla_coordinates.tsv", sep = "\t")

hyla_genpop <- genind2genpop(hyla_genind_adegenet)

Dgen <- dist.genpop(hyla_genpop, method = 2)
Dgeo <- dist(Hyla_lat_long)  ## can I use this with hierfstat?

ibd <- mantel.randtest(Dgen, Dgeo)

ibd

plot(ibd)

plot(Dgeo, Dgen)
abline(lm(Dgen~Dgeo))



## First calclate the pairwise FSTs between each combination of populations.

ppfst.hyla_genind<-pairwise.neifst(hyla_genind)


## make dataframe with populations as rownames

row.names(hyla_genind$tab)

## PCA

popcols = rep(rainbow(12))[hyla_genind$Pop]

pcaSNP<-indpca(hyla_genind,ind.labels=hyla_genind$Pop)
plot(pcaSNP, col=popcols,eigen=F, pch = 16)




