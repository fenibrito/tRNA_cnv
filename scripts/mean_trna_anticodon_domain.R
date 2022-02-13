#calculate each domain mean tRNA gene loci per anticodon per genome 
# this analysis exclude the undetermined tRNAs 
# plot the mean of tRNA gene per anticodon as a box plot
# median test = Wilcoxon rank sum test with continuity correction

# import library
library("tidyverse")
library("RColorBrewer")

# save image directory 
figdir <- '/home/feni/repository/trna_project/figures'

# upload GtRNAdb data
gtrnadb.data <- read_csv('/home/feni/repository/trna_project/data/gtrnadb_data.csv')

# define color palette
domainColor <- c("#049F76","#118AB2","#BF805F")
domainColor2 <- c("#A2E1F6","#AFFDE8","#E8D1C5")
qual_col_pals <-  brewer.pal.info[brewer.pal.info$category == 'qual',]
anticodonpal <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#define domain order 
domain.order <- c('Bacteria', 'Archaea', 'Eukaryota')

# mean tRNA gene loci per anticodon per genome
mean.trna.anticodon.per.genome <- gtrnadb.data %>%
  filter(!str_detect(Anticodon, "N|M|Y")) %>% 
  select(Domain,GenomeID,Anticodon) %>%
  group_by(Domain,GenomeID,Anticodon) %>%
  summarise(n.trna.per.anticodon = n()) %>%
  ungroup() %>% 
  group_by(Domain,Anticodon) %>% 
  mutate(mean.trna.per.anticodon = mean(n.trna.per.anticodon)) %>% 
  select(Domain,mean.trna.per.anticodon) %>% 
  distinct() %>% 
  ungroup()

mean.trna.anticodon.per.genome %>%
  group_by(Domain) %>% 
  mutate(mean = mean(mean.trna.per.anticodon),
         median = median(mean.trna.per.anticodon),
         stdev = sd(mean.trna.per.anticodon),
         error = sd(mean.trna.per.anticodon)/sqrt(length(mean.trna.per.anticodon)),
         interquart = IQR(mean.trna.per.anticodon)) %>% 
  select(Domain,mean,median,stdev,error,interquart) %>% 
  distinct() 

# plot
trna.anticodon.per.genome <- mean.trna.anticodon.per.genome %>% 
  ggplot() +
  geom_boxplot(aes(x = factor(Domain, domain.order), y = mean.trna.per.anticodon, fill = Domain),
               outlier.shape = 8,notch = T,lwd = 0.8,width=0.4,color = domainColor,fill = domainColor)+
  geom_jitter(aes(x = factor(Domain, domain.order), y = mean.trna.per.anticodon, color = Anticodon),width=0.2)+
  scale_color_manual(values = anticodonpal)+
  theme_classic()+
  theme(axis.title = element_text(size = 8, hjust = 0.5),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title.align = 1,
        legend.key.size = unit(0.3, "cm"))+
  ylab("Mean of tRNA gene per isoacceptor per genome") + xlab("")

# save plot
ggplot2::ggsave(paste0(figdir,"/trna_anticodon_per_genome.pdf"),
                plot =  trna.anticodon.per.genome,
                device = "pdf",units = "cm",width = 8,height = 8)

#-------------------------------------------------------------------------------------

# median test
archaea.mean.trna.anticodon<- mean.trna.anticodon.per.genome %>% filter(str_detect(Domain, "Archaea"))
bacteria.mean.trna.anticodon <- mean.trna.anticodon.per.genome %>% filter(str_detect(Domain, "Bacteria"))
eukarya.mean.trna.anticodon <- mean.trna.anticodon.per.genome %>% filter(str_detect(Domain, "Eukaryota"))

# median test = wilcoxon test 
arcbac_wilcox <- wilcox.test(archaea.mean.trna.anticodon$mean.trna.per.anticodon, bacteria.mean.trna.anticodon$mean.trna.per.anticodon, paired = F)
arceuk_wilcox <- wilcox.test(archaea.mean.trna.anticodon$mean.trna.per.anticodon, eukarya.mean.trna.anticodon$mean.trna.per.anticodon, paired = F)
eukbac_wilcox <- wilcox.test(bacteria.mean.trna.anticodon$mean.trna.per.anticodon, eukarya.mean.trna.anticodon$mean.trna.per.anticodon, paired = F)

print("Wilcoxon rank sum test results:")
print(paste0('arc vs bac = ',arcbac_wilcox$p.value))
print(paste0('euk vs bac = ',eukbac_wilcox$p.value))
print(paste0('euk vs arc = ',arceuk_wilcox$p.value))
