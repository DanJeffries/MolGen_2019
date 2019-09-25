library(adegenet)
#for hw test
library(pegas)
#for pairwise fsts
library(hierfstat)
setwd("~/Desktop/mol_course/")

eel_snps<-read.genepop("eel_snps_subset.gen")
eel_lat_lon<-read.table("eel_pops_lat_lon.tsv", header=TRUE)
levels(eel_snps@pop)<-eel_lat_lon$Pop
eels_summ_hier<-basic.stats(eel_snps)

#explore hierfstat summstats
names(eels_summ_hier)

eels_summ_hier$n.ind.samp #per population and per locus genotype numbers
eels_summ_hier$pop.freq #per population and per locus allele freqs
eels_summ_hier$Ho #per population and per locus observed heterozygosity
eels_summ_hier$Hs #per population and per locus expected heterozygosity
eels_summ_hier$Ho #per population and per locus observed heterozygosity
eels_summ_hier$Fis #per population and per locus Fis
eels_summ_hier$perloc #per locus summary stats across all pops
eels_summ_hier$overall #per population and per locus observed heterozygosity

plot(eels_summ_hier$perloc$Ho, eels_summ_hier$perloc$Hs, pch=16, cex=0.5, xlab="Ho", ylab="He", main="Ho vs. He - per locus")
abline(0, 1)

plot(colMeans(eels_summ_hier$Ho, na.rm=T), colMeans(eels_summ_hier$Hs, na.rm=T), pch=16, cex=0.5, xlab="Ho", ylab="He", main="Ho vs. He - per locus")

text(colMeans(eels_summ_hier$Ho, na.rm=T), colMeans(eels_summ_hier$Hs, na.rm=T), unique(eel_snps$pop), pos=2)
#change this?
abline(0,1)

#test IBD
eel_snps_pop<-genind2genpop(eel_snps)
Dgen<-dist.genpop(eel_snps_pop, method=2)
Dgeo<-as.dist(geod(lon=eel_lat_lon$Lon, lat=eel_lat_lon$Lat, R=6371))
ibd<-mantel.randtest(Dgen, Dgeo)
plot(ibd)

plot(Dgeo, Dgen)
ibd_lm<-lm(Dgen~Dgeo)
abline(ibd_lm)

#do pairwise fsts

eels_ppfst<-pairwise.fst(eel_snps)



#do PCA
popcols=rep(rainbow(12))[eel_snps$pop]
eel_pca<-indpca(eel_snps, ind.labels = eel_snps$pop)
plot(eel_pca, col=popcols, eigen=F, pch=16)

#alternative
eel_pca_scale<-scaleGen(eel_snps, NA.method="mean")
eel_pca_2<-dudi.pca(eel_pca_scale, cent=FALSE, scale=FALSE,scannf=FALSE, nf=3)
s.class(eel_pca_2$li, eel_snps@pop)
add.scatter.eig(eel_pca_2$eig[1:20],nf=2,xax=1,yax=2)
colorplot(eel_pca_2$li, eel_pca_2$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
