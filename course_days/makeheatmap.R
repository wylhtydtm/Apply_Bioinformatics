
library(tidyverse)
library(DBI)
library(pheatmap)
library(RColorBrewer)

mydb <- dbConnect(RSQLite::SQLite(), "/shared/projects/form_2022_19/pascal/central_database/sample_collab.db")
abun_long <- as_tibble(dbGetQuery(mydb, 'select* from bracken_abundances_long abu left join sample_annot spl using(run_accession);'))
#get sample annotation
sample_annot <- dbGetQuery(mydb, "select * from sample_annot;")

#load input file from snakemake 

#samples <- as_tibble(dbGetQuery(mydb, 'select * from sample_annot order by patient_code;'))

#In the R script, an S4 object named snakemake is available and allows access to input and output files and other parameters.
input_files <- lapply(snakemake@input, function(x) {
    accession <- str_extract(x, "ERR\\d+")
    infile <- read.delim(x, sep = "\t")
    outfile <- infile %>% mutate(run_accession = accession) %>% relocate(run_accession)
    return(outfile)
})

abun_long %>% select(run_accession, taxon_name, fraction_total_reads) %>% filter(fraction_total_reads>0.5) %>% arrange(desc(fraction_total_reads))

#from supplementary figure4, 
pathogens <- c('Burkholderia multivorans','Chlamydia pneumoniae','Dolosigranulum pigrum','Haemophilus influenzae','Haemophilus parainfluenzae','HHV-6','Human coronavirus HKU1','Human coronavirus NL63','Influenza A virus','Klebsiella pneumoniae','Moraxella catharralis','Mycoplasma pneumoniae','Rhinovirus A','Streptococcus pneumoniae','Severe acute respiratory syndrome-related coronavirus')

abun_long %>% filter((str_detect(taxon_name, 'virus') & str_detect(taxon_name, 'Human'))|taxon_name%in% pathogens) %>%
  ggplot(aes(host_subject_id, taxon_name,fill=log2(fraction_total_reads)))+
  geom_raster(hjust=0, vjust=0) +
  scale_fill_gradient(low='white', high='blue')+
  ggtitle('Viral diversity in patient samples')+
  theme(plot.title =element_text(size=20, face='bold'),
        axis.text.x =element_text(angle=45, vjust=1, size=6, hjust=1),
        axis.text.y= element_text(size=6))
        + ylab('Virus')+xlab('samples')


library('heatmaply')
library(RColorBrewer)

# column name as the run_accession
abun_wide <- abun_long%>% pivot_wider(id_cols =taxon_name, names_from =run_accession, values_from= fraction_total_reads)
abun_wide %>% filter(str_detect(taxon_name, 'virus'))
abun_wide_samples <- colnames(abun_wide)[-1]
sample_annot_rtPCR <- abun_long%>%filter(run_accession == abun_wide_samples)%>%select(host_disease_status)
sample_annot_nuc <- abun_long%>%filter(run_accession == abun_wide_samples)%>%select(nuc)

min_abun  <- abun_long%>%filter(fraction_total_reads >0)%>%summarize(min=min(fraction_total_reads, na.rm=TRUE))%>%as.numeric
View(min_abun)

abun_wide %>% head(100) %>% column_to_rownames(var='taxon_name') %>% replace(is.na(.),0) %>% + (min_abun/2) %>%
  mutate_if(is_double,log2) %>%
  heatmaply(k_col=2, k_row=2, xlab='Samples', ylab='taxa', grip_gap=1,
  colors= colorRampPalette(c('white','orange','red'))(50),
  col_side_colors=c(sample_annot_rtPCR,sample_annot_nuc),
  main = 'Oral microbiome diversity in COVID+/- patients', key.title='Abundances',
  plot_method= 'plotly') %>% layout(height=1000, width=800)
  
##alternatively

library('pheatmap')

