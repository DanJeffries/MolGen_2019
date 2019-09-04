## Load packages

library(adegenet)
library(pegas)
library(hierfstat)

## set working directory
setwd("~/Data/MolGen_2019/Frogs/Data/")

## import vcf
hyla_vcf <- read.vcf("batch_1.vcf", )

hyla_vcf$`1288_41`

## convert to genind
hyla_genind <- loci2genind(hyla_vcf, c("NA", "."), ploidy = 2) #, na.alleles = "NA")

hyla_genind$tab[1:30,1:10]

# add the population info 

popcodes <- read.csv("popmap_kept_largepops_codes.txt", sep = "\t", header = F)
pop(hyla_genind) <- popcodes$V2

## explore genind object
pop(hyla_genind)
indNames(hyla_genind)

## summary adegenet
hyla_summ_adeg<-summary(hyla_genind)

hyla_genind$tab[1:20,1:10]


## explore adegenet summary

hyla_summ$n   ## total number of samples
hyla_summ$n.by.pop  ## number of samples per population
hyla_summ$loc.n.all  ## number of alleles per locus
hyla_summ$pop.n.all  ## number of alleles per population
hyla_summ$NA.perc   ## percentage of missing data
hyla_summ$Hobs  ## observed heterozygosity for all loci
hyla_summ$Hexp  ## expected heterozygosity for all loci

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

plot(hyla_summ_hier$perloc$Ho, hyla_summ_hier$perloc$Hs)
abline(0,1)

## Plot by population

plot(colMeans(hyla_summ_hier$Ho, na.rm = T), colMeans(hyla_summ_hier$Hs, na.rm = T))
text(colMeans(hyla_summ_hier$Ho, na.rm = T), colMeans(hyla_summ_hier$Hs, na.rm = T), levels(pop(hyla)), pos = 1)
abline(0,1)


## Test isolation by distance

## make dataframe with populations as rownames

row.names(hyla_genind$tab)

## PCA

pcaSNP<-indpca(hyla_genind,ind.labels=pop(hyla_genind))
plot(pcaSNP) #, col=popcols,eigen=F)



