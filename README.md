
## Recommendations of how to improve this procedure are welcome.

``` r
#install.packages("tidyverse")
#install.packages("rsnps")
#install.packages("rentrez")
#install.packages("pbapply")
#install.packages("protr")

#BiocManager::install("Biostrings")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("biomaRt")
#BiocManager::install("DECIPHER")
#BiocManager::install("ggmsa")

#install.packages("devtools")
#devtools::install_github("Close-your-eyes/igsc")

library(magrittr)

## FCMR (FC receptor for IgM, human); EntrezID: 9214
id <- 9214
```

Single nucleotide polymorphisms (snp) can be acquired from NCBI. They
are mapped to genomic nucleotide positions. But how to map a snp to the
corresponding position on protein level in order to find out its effect
on the protein (e.g. silent mutation yes/no; premature stop codon
induced yes/no, change to a similar aa or to one with very different
chemical properties, etc.).  
Below snps within exons of FCMR (FAIM3, Toso), the human Fc receptor for
IgM, are mapped to protein level and their effects are plotted in a
multiple alignment.  
Data for FCMR are acquired manually and with biomaRt. The latter just to
demonstrate that doing this job purely from R is feasible.

``` r
## manual data acquisition from NCBI and UCSC
# this is always the plus-strand (gene for FCMR is on minus strand)
genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
# gene boundaries from UCSC genome browser
fcmr_end_gen <- 206903317
fcmr_start_gen <- 206921941
# genomic
gen_man_plus <- Biostrings::DNAStringSet(stats::setNames(as.character(genome$chr1[fcmr_end_gen:fcmr_start_gen]), "gen_man_plus"))
gen_man_minus <- Biostrings::reverseComplement(gen_man_plus)
names(gen_man_minus) <- "gen_man_minus"
# exons
exons_man_minus <- Biostrings::DNAStringSet(unlist(protr::readFASTA("NM_005449.5.exons.fa")))
exons_man_plus <- Biostrings::reverseComplement(exons_man_minus)
# cds
cds_man_minus <- Biostrings::DNAStringSet(stats::setNames(protr::readFASTA("huFcmr.fasta")[[1]], "cds_man_minus"))
```

This kind of plotting is a sanity check to make sure that we know which
sequences we have in our variables.

``` r
igsc::MultiplePairwiseAlignmentsToOneSubject(subject = gen_man_plus,
                                             patterns = exons_man_plus,
                                             order.patterns = T,
                                             type = "local")[["match.plot"]]
igsc::MultiplePairwiseAlignmentsToOneSubject(subject = gen_man_minus,
                                             patterns = exons_man_minus,
                                             order.patterns = T,
                                             type = "local")[["match.plot"]]
igsc::MultiplePairwiseAlignmentsToOneSubject(subject = cds_man_minus,
                                             patterns = exons_man_minus,
                                             type = "local")[["match.plot"]]
```

![alt
text](20220806_snp_to_protein_fcmr_files/figure-html/plot%20exons%20against%20genomic%20seq%20and%20cds-1.png)
![alt
text](20220806_snp_to_protein_fcmr_files/figure-html/plot%20exons%20against%20genomic%20seq%20and%20cds-2.png)
![alt
text](20220806_snp_to_protein_fcmr_files/figure-html/plot%20exons%20against%20genomic%20seq%20and%20cds-3.png)
