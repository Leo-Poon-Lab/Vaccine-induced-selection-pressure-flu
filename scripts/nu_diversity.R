library(tidyverse)
library(Biostrings)
library(writexl)

# read reference seq
ref_seq <- readDNAStringSet("../data/reference_HK68_seg.fasta")
stop_t <- cumsum(width(ref_seq))
df_pos <- tibble(segment = names(ref_seq), start = c(1,stop_t[-length(stop_t)]+1), stop = stop_t)
df_pos$width <- df_pos$stop - df_pos$start + 1
stopifnot(df_pos$width == width(ref_seq))

# use snpgenie to estimat gene diversity
dir.create("../results/snpgenie")
## prepare GTF
df_gtf <- read_tsv("../data/reference_seg.GTF", col_names = F)
df_gtf <- df_gtf %>% filter(X3 == "CDS")
df_gtf <- df_gtf %>% filter(X1>"KY33")
unique(df_gtf$X1)
df_gtf$X1[df_gtf$X1=="KY348535.1"] <- "NS"
df_gtf$X1[grepl("gene-M", df_gtf$X9)] <- "M"
df_gtf$X1[grepl("gene-NA", df_gtf$X9)] <- "NA"
df_gtf$X1[grepl("gene-HA", df_gtf$X9)] <- "HA"
df_gtf$X1[grepl("gene-NP", df_gtf$X9)] <- "NP"
df_gtf$X1[grepl("gene-PA", df_gtf$X9)] <- "PA"
df_gtf$X1[grepl("gene-PB1", df_gtf$X9)] <- "PB1"
df_gtf$X1[grepl("gene-PB2", df_gtf$X9)] <- "PB2"
check <- sapply(df_gtf$X1, function(x){which(df_pos$segment==x)})
df_gtf$X4 <- df_gtf$X4+df_pos$start[check]-1
df_gtf$X5 <- df_gtf$X5+df_pos$start[check]-1
df_gtf$X1 <- "HK68"
df_gtf$X9
write_tsv(df_gtf, "../data/reference_mod.GTF", col_names = F, na="", quote_escape="none")

## prepare vcf
files_vcf <- list.files("../results/ivar_snvs/", "tsv", full.names = T)
lapply(files_vcf, function(x){
	# x <- files_vcf[1]
	sample_name_t <- strsplit(x, "ivar_")[[1]][3]
	sample_name_t <- strsplit(sample_name_t, ".tsv", fixed=T)[[1]][1]
	df_tmp <- read_tsv(x)
	df_vcf <- tibble('#CHROM' = df_tmp$REGION)
	df_vcf$POS <- df_tmp$POS
	df_vcf$ID <- "."
	df_vcf$REF <- df_tmp$REF
	df_vcf$ALT <- df_tmp$ALT
	df_vcf$QUAL <- df_tmp$ALT_QUAL
	df_vcf$FILTER <- "PASS"
	df_vcf$INFO <- paste0("DP=", df_tmp$REF_DP+df_tmp$ALT_DP, ";", "AF=", df_tmp$ALT_FREQ)
	# df_vcf$FORMAT <- 
	df_vcf$`<SAMPLE>` <- sample_name_t
	dir.create(paste0("../results/snpgenie/", sample_name_t))
	write_tsv(df_vcf, paste0("../results/snpgenie/", sample_name_t, "/", sample_name_t, ".vcf"))
	file.copy("../data/reference.fasta", paste0("../results/snpgenie/", sample_name_t, "/reference_HK68_seg.fasta"), overwrite=T)
	file.copy("../data/reference_mod.GTF", paste0("../results/snpgenie/", sample_name_t, "/reference_mod.gtf"), overwrite=T)
	cur_wd <- getwd()
	setwd(paste0("../results/snpgenie/", sample_name_t))
	system("rm -r SNPGenie_Results")
	system("~/softwares/SNPGenie/snpgenie.pl --vcfformat=2")
	setwd(cur_wd)
})

# summerise diversity
file_pi <- list.files("../results/", "product_results", recursive = T, full.names = T)
df_pi <- lapply(file_pi, function(x){
	df_tmp <- read_tsv(x, col_types = cols(.default = "c"))
})
df_pi <- bind_rows(df_pi)

df_pi_N <- df_pi %>% pivot_wider("file", names_from = product, values_from = piN)
write_xlsx(df_pi_N, "../results/data_pi_N.xlsx")
df_pi_S <- df_pi %>% pivot_wider("file", names_from = product, values_from = piS)
write_xlsx(df_pi_S, "../results/data_pi_S.xlsx")

df_pi$`piN/piS` <- as.numeric(df_pi$piN)/as.numeric(df_pi$piS)
df_piS_over_piN <- df_pi %>% pivot_wider("file", names_from = product, values_from = `piN/piS`)
write_xlsx(df_piS_over_piN, "../results/data_piS_over_piN.xlsx")

file_pi <- list.files("../results/", "population_summary", recursive = T, full.names = T)
df_pi <- lapply(file_pi, function(x){
	df_tmp <- read_tsv(x, col_types = cols(.default = "c"))
})
df_pi <- bind_rows(df_pi)
write_xlsx(df_piS_over_piN, "../results/data_pi_sample_summary.xlsx")

