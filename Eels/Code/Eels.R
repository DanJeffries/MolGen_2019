## Install necassary packages (make this a separate script?)

install.packages("devtools")
library("devtools")

devtools::install_github("thierrygosselin/radiator")  ## took a good 10-15 mins to install

install_github("zhengxwen/gdsfmt")
install_github("zhengxwen/SNPRelate")


## Filters using Radiator?

library(radiator)

setwd("~/Data/MolGenMethods_2019/Frogs/")

Hyla <- read_vcf("batch_1.vcf", strata = "popmap_kept_largepops_latitude.txt")

Hyla_conversion <- genomic_converter("batch_1.vcf", 
                  strata = "popmap_kept_largepops_latitude.txt", 
                  output = "snprelate", 
                  filename = "Hyla_SNPrelate_format")

## Use SNPRelate for the popgen

library(SNPRelate)

snpgdsSummary("01_radiator_genomic_converter_20190829@1553/Hyla_SNPrelate_format_snprelate_20190829@1554.gds.rad")

Hyla_data <- snpgdsOpen()

  
  
  "01_radiator_genomic_converter_20190829@1553/Hyla_SNPrelate_format_snprelate_20190829@1554.gds.rad")
