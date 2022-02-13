#compare tRNA gene length: tRNA genomic, mature tRNA, tRNA pseudogene

#import library
library("tidyverse")
library("reshape2")

# save image directory 
figdir <- '/home/feni/repository/trna_project/figures'

# upload GtRNAdb data
gtrnadb.data <- read_csv('/home/feni/repository/trna_project/data/gtrnadb_data.csv')

# define domain order 
domain.order <- c('Bacteria', 'Archaea', 'Eukaryota')

# define color palette
domainColor <- c("#049F76","#118AB2","#BF805F")
domainColor2 <- c("#A2E1F6","#AFFDE8","#E8D1C5")

# pseudo tRNA mean length
ptrna.len.mean <- gtrnadb.data %>% 
  filter(str_detect(Anticodon, "N|M|Y")) %>% 
  select(Domain,GenomeID,length)%>%
  group_by(Domain, GenomeID) %>%
  summarise(ptrna.len.mean = mean(length, na.rm = T)) %>%
  select(Domain, GenomeID, ptrna.len.mean) %>%
  distinct() %>% 
  ungroup()%>% 
  melt(id = c("Domain","GenomeID"),value.name = "mean.gene.length")

# mature tRNA mean length
mat.tRNA.len <- gtrnadb.data %>% 
  filter(!str_detect(Anticodon, "N|M|Y")) %>% 
  select(Domain,GenomeID,len_mature) %>%
  group_by(Domain,GenomeID) %>%
  summarise(mat.trna.len.mean = mean(len_mature, na.rm = T)) %>%
  select(Domain, GenomeID, mat.trna.len.mean) %>%
  distinct() %>% 
  ungroup()%>% 
  melt(id = c("Domain","GenomeID"),value.name = "mean.gene.length")

# canonical tRNA mean length
trna_len_mean <- gtrnadb.data %>% 
  filter(!str_detect(Anticodon, "N|M|Y")) %>% 
  select(Domain,GenomeID,length) %>%
  group_by(Domain,GenomeID) %>%
  summarise(trna.len.mean = mean(length, na.rm = T)) %>%
  select(Domain,GenomeID, trna.len.mean) %>%
  distinct() %>% 
  ungroup() %>% 
  melt(id = c("Domain","GenomeID"),value.name = "mean.gene.length")

trna.gene.length <- rbind(ptrna.len.mean,mat.tRNA.len,trna_len_mean)

genetypelabell <- as_labeller(c("trna.len.mean"="tRNA (genomic locus)","mat.trna.len.mean"="Mature tRNA","ptrna.len.mean"="tRNA pseudogene"))

compare.trna.gene.len <- trna.gene.length %>%
  # filter(!str_detect(variable, 'ptrna.len.mean')) %>% 
  ggplot(aes(x = factor(Domain, domain.order), y = mean.gene.length, color = variable, fill = variable))+
  geom_violin(scale = "width")+
  geom_boxplot(outlier.size = 0.4,outlier.shape = 23,outlier.colour = "black",outlier.fill = "black",outlier.stroke = 1.5,outlier.alpha = .2,
               notch = T,fill="white",color="grey25",width=0.2,lwd = 0.25)+
  facet_wrap(vars(factor(variable,c("trna.len.mean","mat.trna.len.mean","ptrna.len.mean"))), labeller = genetypelabell)+
  scale_fill_manual(values = c("#264653","#242423","#2A9D8F"))+
  scale_color_manual(values = c("#264653","#242423","#2A9D8F"))+
  scale_y_continuous(limits = c(70,120),breaks = c(70,90,110,120))+
  theme_linedraw(base_line_size = 0.2)+ 
  xlab("none")+ylab('tRNA gene length (nt)')+
  theme(axis.title.x=element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.position = "none")

# save plot
ggplot2::ggsave(paste0(figdir,"/compare_trna_gene_len.pdf"),
                plot =  compare.trna.gene.len,
                device = "pdf",units = "cm",width = 12,height = 10)

compare.trna.gene.vs.pseudo.len <- trna.gene.length %>%
  # filter(str_detect(variable, 'ptrna.len.mean')) %>% 
  ggplot(aes(x = factor(Domain, domain.order), y = mean.gene.length, color = variable, fill = variable))+
  geom_violin(scale = "width")+
  geom_boxplot(outlier.size = 0.4,outlier.shape = 23,outlier.colour = "black",outlier.fill = "black",outlier.stroke = 1.5,outlier.alpha = .2,
               notch = T,fill="white",color="grey25",width=0.2,lwd = 0.25)+
  facet_wrap(vars(factor(variable,c("trna.len.mean","mat.trna.len.mean","ptrna.len.mean"))),scales = 'free', labeller = genetypelabell)+
  scale_fill_manual(values = c("#264653","#242423","#2A9D8F"))+
  scale_color_manual(values = c("#264653","#242423","#2A9D8F"))+
  theme_linedraw(base_line_size = 0.2)+ 
  xlab("none")+ylab('tRNA gene length (nt)')+
  theme(axis.title.x=element_blank(),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.position = "none")

#salve plot
ggplot2::ggsave(paste0(figdir,"/compare_trna_gene_vs_pseudo_len.pdf"),
                plot =  compare.trna.gene.vs.pseudo.len,
                device = "pdf",units = "cm",width = 12,height = 10)


#-------------------------------------------------------------------------------------
# pseudo tRNA median test

#average data for each domain separated for the wilcoxon test
archaea.ptrna.len.mean <- ptrna.len.mean %>% filter(str_detect(Domain, "Archaea"))
bacteria.ptrna.len.mean <- ptrna.len.mean %>% filter(str_detect(Domain, "Bacteria"))
eukarya.ptrna.len.mean <- ptrna.len.mean %>% filter(str_detect(Domain, "Eukaryota"))

#median test = wilcoxon test 
arcbac_wilcox <- wilcox.test(archaea.ptrna.len.mean$mean.gene.length, bacteria.ptrna.len.mean$mean.gene.length, paired = F)
arceuk_wilcox <- wilcox.test(archaea.ptrna.len.mean$mean.gene.length, eukarya.ptrna.len.mean$mean.gene.length, paired = F)
eukbac_wilcox <- wilcox.test(bacteria.ptrna.len.mean$mean.gene.length, eukarya.ptrna.len.mean$mean.gene.length, paired = F)

print("Wilcoxon rank sum test results for tRNA-pseudogenes mean length:")
print(paste0('arc vs bac = ',arcbac_wilcox$p.value))
print(paste0('euk vs bac = ',eukbac_wilcox$p.value))
print(paste0('euk vs arc = ',arceuk_wilcox$p.value))


# mature tRNA median test

#average data for each domain separated for the wilcoxon test
archaea.mat.tRNA.len <- mat.tRNA.len %>% filter(str_detect(Domain, "Archaea"))
bacteria.mat.tRNA.len <- mat.tRNA.len %>% filter(str_detect(Domain, "Bacteria"))
eukarya.mat.tRNA.len <- mat.tRNA.len %>% filter(str_detect(Domain, "Eukaryota"))

#median test = wilcoxon test 
arcbac_wilcox <- wilcox.test(archaea.mat.tRNA.len$mean.gene.length, bacteria.mat.tRNA.len$mean.gene.length, paired = F)
arceuk_wilcox <- wilcox.test(archaea.mat.tRNA.len$mean.gene.length, eukarya.mat.tRNA.len$mean.gene.length, paired = F)
eukbac_wilcox <- wilcox.test(bacteria.mat.tRNA.len$mean.gene.length, eukarya.mat.tRNA.len$mean.gene.length, paired = F)

print("Wilcoxon rank sum test results for tRNA maturure mean length:")
print(paste0('arc vs bac = ',arcbac_wilcox$p.value))
print(paste0('euk vs bac = ',eukbac_wilcox$p.value))
print(paste0('euk vs arc = ',arceuk_wilcox$p.value))

# genomic tRNA median test

#average data for each domain separated for the wilcoxon test
archaea.trna_len_mean <- trna_len_mean %>% filter(str_detect(Domain, "Archaea"))
bacteria.trna_len_mean <- trna_len_mean %>% filter(str_detect(Domain, "Bacteria"))
eukarya.trna_len_mean <- trna_len_mean %>% filter(str_detect(Domain, "Eukaryota"))

#median test = wilcoxon test 
arcbac_wilcox <- wilcox.test(archaea.trna_len_mean$mean.gene.length, bacteria.trna_len_mean$mean.gene.length, paired = F)
arceuk_wilcox <- wilcox.test(archaea.trna_len_mean$mean.gene.length, eukarya.trna_len_mean$mean.gene.length, paired = F)
eukbac_wilcox <- wilcox.test(bacteria.trna_len_mean$mean.gene.length, eukarya.trna_len_mean$mean.gene.length, paired = F)

print("Wilcoxon rank sum test results for tRNA gene mean length:")
print(paste0('arc vs bac = ',arcbac_wilcox$p.value))
print(paste0('euk vs bac = ',eukbac_wilcox$p.value))
print(paste0('euk vs arc = ',arceuk_wilcox$p.value))
