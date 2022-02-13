# calculate each domain mean tRNA gene loci per isoacceptor per genome 
# this analysis exclude the undetermined tRNAs 
# plot the mean of tRNA gene per isoacceptor as a box plot
# median test = Wilcoxon rank sum test with continuity correction

#import library
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

# define domain order 
domain.order <- c('Bacteria', 'Archaea', 'Eukaryota')

# mean tRNA gene loci per isotype per genome
mean.trna.isotype.per.genome <- gtrnadb.data %>%
  filter(!str_detect(Anticodon, "N|M|Y"),
         !str_detect(Isotype, "Sup")) %>% 
  select(Domain,GenomeID,Isotype) %>%
  group_by(Domain,GenomeID,Isotype) %>%
  summarise(n.trna.per.isotype = n()) %>%
  ungroup() %>% group_by(Domain,Isotype) %>% 
  mutate(mean.trna.per.isotype = mean(n.trna.per.isotype)) %>% 
  select(Domain,mean.trna.per.isotype) %>% distinct %>% ungroup()

mean.trna.isotype.per.genome %>%
  group_by(Domain) %>% 
  mutate(mean = mean(mean.trna.per.isotype),
         median = median(mean.trna.per.isotype),
         stdev = sd(mean.trna.per.isotype),
         error = sd(mean.trna.per.isotype)/sqrt(length(mean.trna.per.isotype)),
         interquart = IQR(mean.trna.per.isotype)) %>% 
  select(Domain,mean,median,stdev,error,interquart) %>% 
  distinct() 

# plot
trna.isotype.per.genome <- mean.trna.isotype.per.genome %>% 
  ggplot() +
  geom_boxplot(aes(x = factor(Domain, domain.order), y = mean.trna.per.isotype, fill = Domain),
               outlier.shape = 6, 
               notch = T,
               lwd = 0.8,
               width=0.4,
               color = domainColor,
               fill = domainColor)+
  geom_jitter(aes(x = factor(Domain, domain.order), y = mean.trna.per.isotype, color = Isotype),width=.2)+
  scale_color_manual(values = anticodonpal)+
  theme_classic()+
  theme(axis.title = element_text(size = 10, hjust = 0.5),
        legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.title.align = .5,
        legend.key.size = unit(0.3, "cm"))+
  ylab("Mean of tRNA gene per isoacceptor per genome") + xlab("")

# save plot
ggplot2::ggsave(paste0(figdir,"/trna_isotype_per_genome.pdf"),
                plot =  trna.isotype.per.genome,
                device = "pdf",units = "cm",width = 8,height = 12)

#-------------------------------------------------------------------------------------

# median test
archaea.mean.trna.isotype <- mean.trna.isotype.per.genome %>%   filter(str_detect(Domain, "Archaea"))
bacteria.mean.trna.isotype <- mean.trna.isotype.per.genome %>%   filter(str_detect(Domain, "Bacteria"))
eukarya.mean.trna.isotype <- mean.trna.isotype.per.genome %>%   filter(str_detect(Domain, "Eukaryota"))

#median test = wilcoxon test 
arcbac_wilcox <- wilcox.test(archaea.mean.trna.isotype$mean.trna.per.isotype, bacteria.mean.trna.isotype$mean.trna.per.isotype, paired = F)
arceuk_wilcox <- wilcox.test(archaea.mean.trna.isotype$mean.trna.per.isotype, eukarya.mean.trna.isotype$mean.trna.per.isotype, paired = F)
eukbac_wilcox <- wilcox.test(bacteria.mean.trna.isotype$mean.trna.per.isotype, eukarya.mean.trna.isotype$mean.trna.per.isotype, paired = F)

print("Wilcoxon rank sum test results:")
print(paste0('arc vs bac = ',arcbac_wilcox$p.value))
print(paste0('euk vs bac = ',eukbac_wilcox$p.value))
print(paste0('euk vs arc = ',arceuk_wilcox$p.value))
