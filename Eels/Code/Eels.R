## Install necassary packages (make this a separate script?)

install.packages("devtools")
library("devtools")

install_github("thierrygosselin/radiator")  ## took a good 10-15 mins to install

install_github("zhengxwen/gdsfmt")
install_github("zhengxwen/SNPRelate")


## Filters using Radiator?

library(radiator)

setwd("~/Data/MolGen_2019/Frogs/")

Hyla <- read_vcf("batch_1.vcf", strata = "popmap_kept_largepops_latitude.txt")

Hyla_conversion <- genomic_converter("batch_1.vcf", 
                  strata = "popmap_kept_largepops_latitude.txt", 
                  output = "snprelate", 
                  filename = "Hyla_SNPrelate_format")

## Use SNPRelate for the popgen

library(SNPRelate)

GDS_file = "01_radiator_genomic_converter_20190829@1553/Hyla_SNPrelate_format_snprelate_20190829@1554.gds.rad"

snpgdsSummary("01_radiator_genomic_converter_20190829@1553/Hyla_SNPrelate_format_snprelate_20190829@1554.gds.rad")

Hyla_data <- snpgdsOpen("01_radiator_genomic_converter_20190829@1553/Hyla_SNPrelate_format_snprelate_20190829@1554.gds.rad")

pops <- read.csv("popmap_kept_largepops_latitude.txt", sep = "\t")
pops$STRATA <- as.character(pops$STRATA)

add.gdsn(GDS_file, "sample.annot", pops)

