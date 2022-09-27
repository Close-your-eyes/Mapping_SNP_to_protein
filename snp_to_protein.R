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
<<<<<<< HEAD
}, mc.cores = 32)
=======
}, mc.cores = 16)
#saveRDS(ncbi_snps, file.path(getwd(), "ncbi_snps.RDS"))
#ncbi_snps <- readRDS(file.path(getwd(), "ncbi_snps.RDS"))
>>>>>>> c21779332bbe2b41c16da63def038d361a828c3c

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
cds_exon_al <- Biostrings::pairwiseAlignment(pattern = cds_man_minus, subject = do.call(Biostrings::xscat, as.list(exons_man_minus)), type = "local")
## ---- align cds to exons plot -------------
igsc::MultiplePairwiseAlignmentsToOneSubject(patterns = cds_man_minus, subject = do.call(Biostrings::xscat, as.list(exons_man_minus)), type = "local")[["match.plot"]]

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
        Biostrings::translate(Biostrings::DNAString(paste0(df_exon_gen_pos[which(df_exon_gen_pos$cds),"nt"], collapse = "")))#[range]
    })
    names(protein_with_snp) <- paste0(ncbi_snps_df[y,"rsid"], "_", strsplit(ncbi_snps_df[y,"variation_allele",drop=T], ",")[[1]])
    return(protein_with_snp)
})
#names(fcmr_prot) <- ncbi_snps_df$rsid
fcmr_prot <- purrr::flatten(fcmr_prot)
# fcmr_prot is a list of protein sequences with modification according to the snp
# so one list entry for each snv-snp with exons
fcmr_prot[1:5]

## ---- exclude snps which cause a silent mutation ---------
# exclude snps which cause a silent mutation
cds_aa_range <- Biostrings::translate(Biostrings::DNAStringSet(cds_man_minus))[[1]]#[range]
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
# create a mutiple alignment of proteins which have non-silent snps incorporated (snv only, see above)
fcmr_prot2_mult_aln <- DECIPHER::AlignSeqs(Biostrings::AAStringSet(fcmr_prot2), verbose = F)
names(fcmr_prot2_mult_aln) <- make.unique(names(fcmr_prot2_mult_aln))

# print in easy-to-read-format like with printPairwiseAlignment
# http://emboss.sourceforge.net/docs/themes/AlignFormats.html

# or use DECIPHER (which opens the browser); this looks nicer and allows scrolling but cannot produce a figure
# DECIPHER::BrowseSeqs(fcmr_prot2_mult_aln, highlight = 1)

# or use ggmsa to get a ggplot2 object; note: stop-codons are plotted as empty positions (not so nice currently; DECIPHER plots asterisk (*) for stop-codons)
gg_aln <- ggmsa::ggmsa(msa = fcmr_prot2_mult_aln,
                       seq_name = T,
                       consensus_views = T,
                       ref = "cds",
                       use_dot = T)

## manual ggplotting
library(ggplot2)
out3 <- purrr::map(as.list(fcmr_prot2_mult_aln), as.character)
out4 <- purrr::flatten(purrr::map(out3, strsplit, split = ""))
out5 <- purrr::map_dfr(out4, function(x) {
    stack(setNames(x, seq(1, length(x))))
}, .id = "rsid")
out6 <-
    out5 %>%
    tidyr::pivot_wider(names_from = rsid, values_from = values) %>%
    dplyr::mutate(dplyr::across(.cols = unique(out5$rsid)[which(unique(out5$rsid) != "cds")],
                                .fns = ~ ifelse(.x == cds, ".", .x))) %>%
    tidyr::pivot_longer(cols = unique(out5$rsid), names_to = "rsid", values_to = "aa") %>%
    dplyr::mutate(rsid = factor(rsid, levels = unique(out5$rsid)))


AA_cols <- setNames(ggmsa:::scheme_AA$Chemistry_AA, rownames(ggmsa:::scheme_AA))
AA_cols <- c(AA_cols, "*" = "grey80")
p <- ggplot(out6, aes(x = ind, y = rsid, label = aa)) +
    geom_tile(aes(fill = aa), color = "black") +
    theme_bw() +
    theme(axis.title = element_blank(), legend.position = "none", panel.border = element_blank(),
          axis.ticks = element_blank()) +
    geom_text() +
    scale_x_discrete(breaks = c(10,20,30)) +
    scale_y_discrete(limits = rev) +
    scale_fill_manual(values = AA_cols, na.value = "white") +
    coord_fixed(ratio = nlevels(out6$ind)/nlevels(out6$rsid)*0.8)
#ggsave(p, filename = paste0("alignment_30aa.png"), device = "png", path = getwd(), dpi = "retina", width = 5*nlevels(out6$ind)/nlevels(out6$rsid)*0.8, height = 5)
ggsave(p, filename = paste0("alignment_all_aa.png"), device = "png", path = getwd(), dpi = "retina", width = 30*nlevels(out6$ind)/nlevels(out6$rsid)*0.8, height = 30)



## ---- for github (ggmsa only) -------
aln <- DECIPHER::AlignSeqs(Biostrings::AAStringSet(named_protein_chr_vector), verbose = F)
out <- purrr::map(as.list(aln), as.character)
out <- purrr::flatten(purrr::map(out, strsplit, split = ""))
out <- purrr::map_dfr(out, function(x) {
    stack(setNames(x, seq(1, length(x))))
}, .id = "protein_name")

ref_protein_name <- "ref"
out <-
    out %>%
    tidyr::pivot_wider(names_from = protein_name, values_from = values) %>%
    dplyr::mutate(dplyr::across(.cols = unique(out$protein_name)[which(unique(out$protein_name) != ref_protein_name)],
                                .fns = ~ ifelse(.x == !!sym(ref_protein_name), ".", .x))) %>%
    tidyr::pivot_longer(cols = unique(out$protein_name), names_to = "protein_name", values_to = "aa") %>%
    dplyr::mutate(protein_name = factor(protein_name, levels = unique(out$protein_name)))

AA_cols <- setNames(ggmsa:::scheme_AA$Chemistry_AA, rownames(ggmsa:::scheme_AA))
AA_cols <- c(AA_cols, "*" = "grey80")
p <- ggplot(out, aes(x = ind, y = protein_name, label = aa)) +
    geom_tile(aes(fill = aa), color = "black") +
    theme_bw() +
    theme(axis.title = element_blank(), legend.position = "none", panel.border = element_blank(),
          axis.ticks = element_blank()) +
    geom_text() +
    scale_x_discrete(breaks = c(10,20,30)) +
    scale_y_discrete(limits = rev) +
    scale_fill_manual(values = AA_cols, na.value = "white") +
    coord_fixed(ratio = nlevels(out$ind)/nlevels(out$protein_name)*0.8)
