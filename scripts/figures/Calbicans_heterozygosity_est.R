## Purpose: Calculate genome heterozygosity % and plot distribution
## Author: Nancy Scott
## Email: scot0854@umn.edu

## Get phenotyping data----
# *C. albicans* specific data from candidemia phenotyping manuscript
source("scripts/figures/Calbicans_redcap_summary.R")

sample_naming <- sample_info %>% 
  mutate(join_id = case_when(!is.na(secondary_id) ~ secondary_id, .default = primary_id)) %>%
  select(primary_id, join_id)

## Load packages----
library(correlation)
library(patchwork)

## Variables----
spreadsheet_list <- "data/metadata/2024_Calbicans_snp_depth_paths.txt"
in_patient_data <- "data/metadata/2024_Calbicans_sorted_patient_info.xlsx"
ordered_patient_data <-  "data/metadata/2025_Calbicans_MEC_ropars_hirakawa_raxml_midpoint_tips.txt"
chr_size <- "scripts/figures/plotting_data/Calbicans_chr_lengths.csv"
save_dir <- "images/Calbicans/"
genome_size <- 14324315
color_file <- read.table("scripts/figures/plotting_data/clade_colors.txt")

# Metadata
pop_data <- read_excel(in_patient_data)
pop_data[pop_data=="NA"] <- NA

pt_order <- read.csv(ordered_patient_data, header = FALSE, col.names = "sample")
pt_order <- pt_order %>% 
  filter(str_detect(sample, "^AMS|MEC")) %>% 
  mutate(sample = replace(sample, sample=="MEC103", "MEC103-2"))
pt_order <- pt_order %>%
  left_join(pop_data %>% filter(study=="MEC") %>% select(sample, mec_pt_code, mec_isolate_code))

# Plotting info
chrs <- read.csv(chr_size, header = TRUE)

chr_labels <- as_labeller(c("1"="Chr1","2"="Chr2", "3"="Chr3", "4"="Chr4",
                            "5"="Chr5", "6"="Chr6", "7"="Chr7", "8"="ChrR"))

clade_colors <- color_file[,2]
names(clade_colors) <- color_file[,1]

## Read in SNP counts for all samples----
snp_files <- scan(spreadsheet_list, what=character())

genome_snp <- read_xlsx(snp_files[1]) %>%
  select(index, pos, snp_count)
names(genome_snp)[names(genome_snp)=="snp_count"] <-str_extract(snp_files[1], "AMS[:digit:]+|MEC[:digit:]+")

for(i in 2:length(snp_files)){
  new_snp <- read_xlsx(snp_files[i]) %>%
    select(index, pos, snp_count)
  names(new_snp)[names(new_snp)=="snp_count"] <-str_extract(snp_files[i], "AMS[:digit:]+|MEC[:digit:]+")

  genome_snp <- genome_snp %>%
    left_join(new_snp, by = join_by(index,pos))
}

snp_again <- genome_snp %>%
  pivot_longer(names_to = "sample", values_to = "snp_count", cols=-c(index,pos)) %>% 
  mutate(sample = replace(sample, sample=="MEC103", "MEC103-2")) %>% 
  left_join(sample_naming, by = c("sample" = "join_id"))

## Estimate het % for genome and chrs----
snp_total <- snp_again %>% 
  group_by(sample) %>%
  summarize(all_snps = sum(snp_count, na.rm = TRUE)) %>%
  mutate(het_percentage = all_snps/genome_size*100)

snp_total <- snp_total %>% 
  left_join(pt_order)%>% 
  left_join(sample_naming, by = c("sample" = "join_id")) %>% 
  left_join(pop_data %>% select(sample, ploidy))

snp_total <- snp_total %>%
  mutate(outlier = ifelse(sample %in% c("AMS5232", "MEC352"), mec_pt_code, NA))

chr_total <- snp_again %>% 
  group_by(sample, index) %>% 
  summarise(chr_snps = sum(snp_count, na.rm = TRUE)) %>% 
  mutate(chr_het_percent = chr_snps/chrs$length * 100) %>% 
  left_join(pt_order) %>% 
  left_join(sample_naming, by = c("sample" = "join_id")) %>% 
  left_join(pop_data %>% select(sample, ploidy))

chr_summary <- chr_total %>% group_by(index) %>% 
  summarise(min_het = min(chr_het_percent), max_het = max(chr_het_percent),
            mean_het = mean(chr_het_percent), median_het = median(chr_het_percent))

## Correlation tests----
# overall heterozygosity, growth rate, SMG 
genome_het_gc <- snp_total %>% 
  filter(is.na(ploidy)) %>% 
  left_join(gc %>% filter(drug=="fluconazole"), by = "primary_id") %>% 
  select(primary_id, het_percentage, r) %>% 
  left_join(mic_info %>% filter(drug=="fluconazole"), by = "primary_id") 
  
genome_het_gc_corr <- genome_het_gc %>%
  filter(is.na(ploidy)) %>% 
  correlation(select = c("het_percentage", "r", "mean_smg"))

# Per chrom heterozygosity, growth rate, SMG correlation
chr_het_gc_corr <- chr_total %>% 
  filter(is.na(ploidy)) %>% 
  left_join(gc %>% filter(drug=="fluconazole"), by = "primary_id") %>% 
  select(primary_id, chr_het_percent, r, index) %>% 
  left_join(mic_info %>% filter(drug=="fluconazole"), by = "primary_id") %>% 
  group_by(index) %>% 
  select(sample, index, chr_het_percent, r, mean_smg) %>% 
  correlation()

## Plots----
# Overall heterozygosity
het_all <- ggplot(snp_total, aes(het_percentage)) +
  geom_histogram(color = "grey50", fill = "grey65", bins = 40) +
  scale_x_continuous(limits = c(0.0,0.8)) +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black")) +
  xlab("Genome heterozygosity, %") +
  ylab("Isolates, N")

ggsave(paste0(save_dir,format(Sys.Date(), "%Y"),"_Calbicans_MEC_genome_het.pdf"),
       het_all,
       width = 6,
       height = 4.5,
       units = "in")

# Faceted chr heterozygosity
het_chr_facet <- chr_total %>% 
  ggplot(aes(chr_het_percent)) +
  geom_histogram(color = "grey50", fill = "grey65", bins = 40) +
  facet_wrap(~index, nrow = 1, labeller = chr_labels) +
  scale_x_continuous(limits = c(0, 0.8), breaks = c(0.2, 0.4, 0.6)) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 7),
        axis.title = element_text(color = "black", size = 8),
        ) +
  #xlab(NULL) +
  ylab(NULL) +
  xlab("Chromosome heterozygosity")
  #ylab("Isolates, N") 

ggsave(paste0(save_dir,format(Sys.Date(),"%Y"),"_Calbicans_MEC_chr_het.pdf"),
       het_chr_facet,
       width = 6,
       height = 4.5,
       units = "in")

# Stacked plots
gen_chr_het <- (het_all+ theme(plot.margin = margin(b=12)))/
  het_chr_facet + plot_annotation(tag_levels = "A")

ggsave(paste0(save_dir,format(Sys.Date(), "%Y"),"_Calbicans_MEC_combined_hets.pdf"),
       gen_chr_het,
       width = 6,
       height = 6.2,
       units = "in")

# Scatter of phenotypes and heterozygosity
r_genome_het_scatter <- snp_total %>% 
  left_join(gc %>% filter(drug=="fluconazole"), by = "primary_id") %>% 
  #filter(is.na(ploidy)) %>% 
  ggplot(aes(y = r, x = het_percentage, color = ploidy)) + 
  geom_point() +
  scale_color_manual(values = c("#0072B2", "#0072B2", "#F0E442", "#D55E00", "#57606C")) +
  #geom_text(aes(label = outlier), size = 4) +
  theme_bw() +
  scale_x_continuous(limits = c(0,0.8)) +
  theme(axis.text = element_text(color = "black", size = 7),
        axis.title = element_text(color = "black"),
  ) +
  xlab("Genome heterozygosity, %") +
  ylab("Growth rate (hour-1)")

k_chr_het_scatter <- chr_total %>% 
  left_join(gc %>% filter(drug=="fluconazole"), by = "primary_id") %>% 
  filter(is.na(ploidy)) %>% 
  ggplot(aes(x = k, y = chr_het_percent)) + 
  geom_point() +
  facet_wrap(~index, nrow = 2, labeller = chr_labels) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 7),
        axis.title = element_text(color = "black"),
  ) +
  xlab("Chromosome heterozygosity, %") +
  ylab("Carrying capacity (k)")

smg_chr_het_scatter <- chr_total %>% 
  left_join(mic_info %>% filter(drug=="fluconazole"), by = "primary_id") %>% 
  #filter(is.na(ploidy)) %>% 
  filter(index == 4) %>% 
  ggplot(aes(y = mean_smg, x = chr_het_percent)) + 
  geom_point(alpha = 0.5, size = 2) +
  theme_bw() +
  theme(axis.text = element_text(color = "black", size = 9),
        axis.title = element_text(color = "black", size = 10),
        #strip.background = element_blank(),
        #strip.text = element_blank()
  ) +
  scale_x_continuous(limits = c(0.0,0.8))+
  scale_y_continuous(limits = c(0.0, 0.8))+
  #xlab("Chr4 heterozygosity, %") +
  #ylab("Mean supra-MIC growth") +
  xlab(NULL)+
  ylab(NULL)+
  annotate("point", x = 0.71, y = 0.556, color="#931A1D", size = 2.5) +
  annotate("text", x = 0.73, y = 0.61, label = "120")

ggsave(paste0(save_dir,format(Sys.Date(),"%Y"),"_Calbicans_MEC_smg_vs_Chr4_het.pdf"),
       smg_chr_het_scatter,
       width = 4,
       height = 3.5,
       units = "in")

# composite with LOH heatmap from separate script
# read in heatmap plot if needed or else run script to get p 
p/smg_chr_het_scatter/het_chr_facet + plot_layout(heights = c(4,1,1))
ggsave("2025_LOH_smg_scatter_chr_het.pdf", height = 9.5, width = 8, units = "in")

# composite of LOH heatmap, genome het and growth rate scatter
het_all + r_genome_het_scatter
ggsave(paste0(save_dir, format(Sys.Date(), "%Y"), "_Calbicans_MEC_genome_het_histo_scatter.pdf"),
       width = 8,
       height = 3,
       units = "in")
