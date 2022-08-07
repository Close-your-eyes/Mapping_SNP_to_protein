## ---- setup --------
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

## ---- get genomic and exon data for FCMR manually --------
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

## ---- plot exons against genomic seq and cds --------
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

## ---- get genomic and exon data for FCMR with bioMart --------
## data acquisition with bioMart
# same information in the end; just to demonstrate how to purely work from R
# queries for this FCMR (different isoforms are documented) only work without trailing .8 even though this is correct identifier at ensemble
# also choosing the 'correct' identifier (which is the desired isoform) was not trivial and had to be checked manually
# bioMart retrieves the gene information (cds; genomic; exons+introns) on the genes strand; in case of FCMR it is the minus strand
mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# genomic
gen_mart_minus <-  biomaRt::getSequence(id = "ENST00000367091", #.8
                                        type = "ensembl_transcript_id",
                                        seqType = "gene_exon_intron",
                                        mart = mart)
gen_mart_minus <- Biostrings::DNAStringSet(stats::setNames(gen_mart_minus$gene_exon_intron, "gen_mart_minus"))
gen_mart_plus <- Biostrings::reverseComplement(gen_mart_minus)
names(gen_mart_plus) <- "gen_mart_plus"

# exons
exons_mart_minus_df <- biomaRt::getBM(attributes = c("gene_exon", "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end"),
                                      filters = "ensembl_transcript_id",
                                      values = "ENST00000367091", #.8
                                      mart = mart)
exons_mart_minus_df <- exons_mart_minus_df[order(exons_mart_minus_df$exon_chrom_start, decreasing = TRUE),]
exons_mart_minus <- Biostrings::DNAStringSet(stats::setNames(exons_mart_minus_df$gene_exon, exons_mart_minus_df$ensembl_exon_id))
exons_mart_plus <- Biostrings::reverseComplement(exons_mart_minus)

# cds
cds_mart <- biomaRt::getSequence(id = "ENST00000367091", #.8
                                 type = "ensembl_transcript_id",
                                 seqType = "coding",
                                 mart = mart)
cds_mart <- Biostrings::DNAStringSet(stats::setNames(cds_mart$coding, "cds_mart"))

## ---- plot exons against genomic seq and cds 2 --------
# plot genomic seq from manual gathering and bioMart query against each other
igsc::MultiplePairwiseAlignmentsToOneSubject(subject = gen_man_minus,
                                             patterns = gen_mart_minus,
                                             type = "local")[["match.plot"]]
# align exons against genomic seq and cds
igsc::MultiplePairwiseAlignmentsToOneSubject(subject = gen_mart_minus,
                                             patterns = exons_mart_minus,
                                             type = "local",
                                             order.patterns = T)[["match.plot"]]
igsc::MultiplePairwiseAlignmentsToOneSubject(subject = gen_mart_plus,
                                             patterns = exons_mart_plus,
                                             type = "local",
                                             order.patterns = T)[["match.plot"]]
igsc::MultiplePairwiseAlignmentsToOneSubject(subject = cds_mart,
                                             patterns = exons_mart_minus,
                                             type = "local",
                                             order.patterns = T)[["match.plot"]]

## ---- get genomic start and end positions for each exon -------------
# get genomic start and end positions for each exon by aligning exons to DNA minus strand
exon_in_genome <- purrr::map_df(.x = as.character(exons_man_minus), .id = "exon", .f = function(x) {
  pa <- Biostrings::pairwiseAlignment(subject = gen_man_minus,
                                      pattern = x,
                                      type = "local")
  data.frame(start_al = pa@subject@range@start,
             end_al = pa@subject@range@start + pa@subject@range@width - 1,
             start_genome = fcmr_start_gen - pa@subject@range@start + 1,
             end_genome = fcmr_start_gen - (pa@subject@range@start + pa@subject@range@width - 1) + 1)
})
exon_in_genome

## ---- get snp data --------
# get snp identifiers from NCBI
snp_links <- rentrez::entrez_link(dbfrom = "gene", id = id, db = "snp")
fcmr_snps <- snp_links[["links"]][["gene_snp"]]

# get details for each snp
ncbi_snps <- parallel::mclapply(fcmr_snps, function(x) {
  tryCatch({
    # pulls from assembly 38
    rsnps::ncbi_snp_query(paste0("rs", x))
  }, error = function(e) {
    NULL
  })
}, mc.cores = 32)

#attr <- biomaRt::listAttributes(mart)
#marts <- biomaRt::listMarts()
#snps <- biomaRt::getBM( attributes=c("snp_chromosome_strand","chr_name","chrom_start"),filters=c("chr_name","chrom_start","chrom_end"), values=list(c(1),c(fcmr_start_gen), c(fcmr_end_gen)), mart)

# make a data frame; filter for snv* only; filter for snps within exons only
# *snv = single nucleotide variant
# https://stackoverflow.com/questions/15917233/elegant-way-to-vectorize-seq
seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
ncbi_snps_df <-
  dplyr::bind_rows(ncbi_snps) %>%
  dplyr::filter(bp %in% unlist(seq2(exon_in_genome[,"end_genome"], exon_in_genome[,"start_genome"]))) %>%
  dplyr::filter(class == "snv") %>%
  dplyr::arrange(-bp)
head(ncbi_snps_df,10)

## ---- align cds to exons -------------
# align cds to exons
# exons have more nts as they also comprise UTRs
# used alignment below to identify where cds starts
cds_exon_al <- Biostrings::pairwiseAlignment(pattern = cds_man_minus, subject = do.call(Biostrings::xscat, exons_man_minus), type = "local")
igsc::MultiplePairwiseAlignmentsToOneSubject(patterns = cds_man_minus, subject = do.call(Biostrings::xscat, exons_man_minus), type = "local")[["match.plot"]]

## ---- create a data frame with one row for each exonic nucleotide --------
# create a data frame with one row for each exonic nucleotide
df_exon_gen_pos <- purrr::map_df(.x = names(exons_man_minus), .f = function(x) {
  data.frame(exon = x,
             nt = strsplit(as.character(exons_man_minus[[x]]), "")[[1]],
             gen_pos = exon_in_genome[which(exon_in_genome$exon == x), "start_genome"]:exon_in_genome[which(exon_in_genome$exon == x), "end_genome"])
}) %>%
  dplyr::mutate(row_num = dplyr::row_number()) %>%
  dplyr::mutate(cds = dplyr::between(row_num, cds_exon_al@subject@range@start, cds_exon_al@subject@range@start + cds_exon_al@subject@range@width - 1))
head(df_exon_gen_pos,10)

## ---- translate snps to their effect on protein level -------
# translate snps to their effect on protein level
range <- c(1:30) # arbitrary range of aa to consider below
fcmr_prot <- lapply(1:nrow(ncbi_snps_df), function(y) {
  bp <- ncbi_snps_df[y,"bp",drop=T]

  if (as.character(Biostrings::complement(Biostrings::DNAString(ncbi_snps_df[y,"ancestral_allele",drop=T]))) != df_exon_gen_pos[which(df_exon_gen_pos$gen_pos == bp),"nt",drop=T]) {
    # test if ancestral_allele from ncbi snp database matches the information in our data frame
    # turn ncbi info into complement as they come from plus strand, but FCMR is on minus strand
    stop("mismatch")
  }

  protein_with_snp <- lapply(strsplit(ncbi_snps_df[y,"variation_allele",drop=T], ",")[[1]], function(x) {
    df_exon_gen_pos[which(df_exon_gen_pos$gen_pos == bp),"nt"] <- as.character(Biostrings::complement(Biostrings::DNAString(x)))
    Biostrings::translate(Biostrings::DNAString(paste0(df_exon_gen_pos[which(df_exon_gen_pos$cds),"nt"], collapse = "")))[range]
  })
  names(protein_with_snp) <- rep(ncbi_snps_df[y,"rsid"], length(protein_with_snp))
  return(protein_with_snp)
})
names(fcmr_prot) <- ncbi_snps_df$rsid
fcmr_prot <- purrr::flatten(fcmr_prot)
# fcmr_prot is a list of protein sequences with modification according to the snp
# so one list entry for each snv-snp with exons
fcmr_prot[1:5]

## ---- exclude snps which cause a silent mutation ---------
# exclude snps which cause a silent mutation
cds_aa_range <- Biostrings::translate(DNAStringSet(cds_man_minus))[[1]][range]
fcmr_prot2 <- lapply(fcmr_prot, function(x) {
  if (as.character(x) == as.character(cds_aa_range)) {
    return(NULL)
  } else {
    return(x)
  }
})
fcmr_prot2 <- fcmr_prot2[which(!sapply(fcmr_prot2, is.null))]

# add cds at the first index
fcmr_prot2 <- c(list(cds_aa_range), fcmr_prot2)
names(fcmr_prot2)[1] <- "cds"

## ---- mutiple alignment of proteins with snps and plot -------------
# mutiple alignment of proteins with snps and plot
fcmr_prot2_mult_aln <- DECIPHER::AlignSeqs(Biostrings::AAStringSet(fcmr_prot2), verbose = F)
names(fcmr_prot2_mult_aln) <- make.unique(names(fcmr_prot2_mult_aln))

# print in easy-to-read-format like with printPairwiseAlignment
# http://emboss.sourceforge.net/docs/themes/AlignFormats.html

# or use DECIPHER (which opens the browser)
# DECIPHER::BrowseSeqs(fcmr_prot2_mult_aln, highlight = 1)

# or use ggmsa to get a ggplot2 object; note: stop-codons are plotted as empty positions (not so nice currently; DECIPHER plots asterisk (*))
ggmsa::ggmsa(msa = fcmr_prot2_mult_aln,
             seq_name = T,
             consensus_views = T,
             ref = "cds",
             use_dot = T)


