library(tidyverse)
library(Biostrings)
library(ggsci)
library(ggpval)
library(writexl)
source("./helper/save_pptx.r")

# read reference seq
ref_seq <- readDNAStringSet("../data/reference_HK68_seg.fasta")
stop_t <- cumsum(width(ref_seq))
df_pos <- tibble(segment = names(ref_seq), start = c(1,stop_t[-length(stop_t)]+1), stop = stop_t)
df_pos$width <- df_pos$stop - df_pos$start + 1
stopifnot(df_pos$width == width(ref_seq))

# use snpgenie to estimat gene diversity
# https://github.com/chasewnelson/SNPGenie#output
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
## each gene
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

## full genome
file_pi <- list.files("../results/", "population_summary", recursive = T, full.names = T)
df_pi_pop <- lapply(file_pi, function(x){
	df_tmp <- read_tsv(x, col_types = cols(.default = "c"))
})
df_pi_pop <- bind_rows(df_pi_pop)
write_xlsx(df_pi_pop, "../results/data_pi_sample_summary.xlsx")

## ploting and testing
parse_data <- function(df){
	df$sample <- toupper(gsub(".vcf", "", df$file, fixed=T))
	df$passage <- ifelse(grepl("19H3P18", df$sample), "P18", "P5")
	df$replicate <- gsub("19H3P5", "", df$sample, fixed=T)
	df$replicate <- gsub("19H3P18", "", df$replicate, fixed=T)
	df$type <- sapply(strsplit(df$replicate, ""), function(x){x[1]})
	df$type2 <- ifelse(df$type=="N", "Naive", "Vaccinated")
	return(df)
}

df_pi <- parse_data(df_pi)
df_pi_pop <- parse_data(df_pi_pop)
df_plot_pop <- df_pi_pop %>% select(piN, piS, sample:type2) %>% pivot_longer(piN:piS)
df_plot_pop$product <- "Overall"

df_plot_pi <- df_pi %>% select(product, piN, piS, sample:type2) %>% pivot_longer(piN:piS)
df_plot_pi$product <- gsub("gene-", "", df_plot_pi$product)

df_plot_all_raw <- bind_rows(df_plot_pi, df_plot_pop)
df_plot_all <- df_plot_all_raw %>% group_by(product, name, passage, type2) %>% summarise(mean=mean(as.numeric(value)), sd=sd(as.numeric(value)))
df_plot_all$group <- paste(df_plot_all$name, df_plot_all$type2)
df_plot_all$product <- factor(df_plot_all$product, levels=c("PB2", "PB1", "PB1-F2", "PA", "PA-X", "HA", "NP", "NA", "M1", "M2", "NS1", "NEP", "Overall"))

df_plot_all$group <- factor(df_plot_all$group, levels = c("piN Naive","piN Vaccinated","piS Naive","piS Vaccinated"))
color_t <- c("#ca0020","#0571b0","#f4a582","#92c5de")
names(color_t) <- levels(df_plot_all$group)

label_t <- c("piN Naive"=expression(pi[N]~(Naive)),
	"piN Vaccinated"=expression(pi[N]~(Vaccinated)),
	"piS Naive"=expression(pi[S]~(Naive)),
	"piS Vaccinated"=expression(pi[S]~(Vaccinated)))

df_stat <- df_plot_all_raw %>% group_by(product, name, passage) %>% mutate(value=as.numeric(value)) %>% 
	summarise(p_value_VACCINATEDvsNAIVE = t.test(value[type2=="Vaccinated"], value[type2=="Naive"])$p.value) %>% 
	ungroup()

df_plot_all <- left_join(df_plot_all, df_stat)
df_plot_all$sig_label <- sapply(df_plot_all$p_value_VACCINATEDvsNAIVE, function(x){
	if(is.na(x)){return(NA)}
	if(x>0.05){return(NA)}
	if(x>0.01){return("*")}
	if(>0.001){return("**")}
	if(x<=0.001){return("***")}
})
df_t <- df_plot_all %>% group_by(product, name, passage) %>% summarise(y=max(mean+sd)) 
df_plot_all <- left_join(df_plot_all, df_t)

annotation_df <- df_plot_all %>% filter(!is.na(sig_label)) %>% mutate(start=group[1], end=group[2]) %>% ungroup() %>% select(sig_label:end, product, passage, group) %>% unique()

nudge_y <- 0.0003
nudge_y_text <- nudge_y+0.0001
ggplot(df_plot_all, aes(x=group, y=mean, fill = group))+
	geom_bar(stat="identity", position=position_dodge())+
  	geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
	geom_segment(
		data = annotation_df,
		mapping = aes(x=start, xend=end, y=y+nudge_y, yend=y+nudge_y))+
	geom_text(
		data = annotation_df,
		mapping = aes(x=start, y=y+nudge_y_text, label=sig_label),
		nudge_x=0.5)+
	facet_grid(cols=vars(product), rows = vars(passage), scales="free_x")+ 
	# scale_fill_jco(name="Diversity")+
	scale_fill_manual(name="Diversity", values = color_t, labels=label_t)+
	theme_bw()+
	xlab("Gene")+
	ylab("Nucleotide diversity")+
	theme(legend.title = element_blank(),
		axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
	NULL

ggsave("../results/diversity.pdf", width=8*sqrt(2), height=8)
save_pptx("../results/diversity.pptx", width=8*sqrt(2), height=8)
