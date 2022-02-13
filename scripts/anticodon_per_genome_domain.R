# calculate each domain mean tRNA gene loci per anticodon per genome 
# this analysis exclude the undetermined tRNAs 
# plot the number of tRNA gene per anticodon as a box plot
# median test = Wilcoxon rank sum test with continuity correction

#import library
library("tidyverse")
library('ggplot2')

#save image directory 
figdir <- '/home/feni/repository/trna_project/figures'
#import trna data
gtrnadb.data <- read_csv('/home/feni/repository/trna_project/data/gtrnadb_data.csv')

# count the tRNA gene loci per anticodon per genome
anticodon.set.per.genome <- gtrnadb.data %>% 
  filter(!str_detect(Anticodon, "N|M|Y")) %>% 
  select(Domain,GenomeID,Anticodon) %>% 
  group_by(Domain, GenomeID) %>%
  distinct() %>%
  mutate(n.anticodon = n()) %>%
  select(Domain,GenomeID,n.anticodon) %>%
  distinct() %>% 
  ungroup() %>% group_by(Domain) %>%
  mutate(mean.trna.anticodon = mean(n.anticodon),
         median.trna.anticodon = median(n.anticodon),
         stdev = sd(n.anticodon),
         error = sd(n.anticodon)/sqrt(length(n.anticodon))) %>%
  ungroup()

# mean and median of tRNA anticodon per domain
anticodon.set.per.genome %>%
  group_by(Domain) %>% 
  mutate(interquart = IQR(n.anticodon)) %>% 
  select(Domain,mean.trna.anticodon,median.trna.anticodon,stdev,error,interquart) %>% 
  distinct() 

#define palette color
domainColor <- c("#049F76","#118AB2","#BF805F")
domainColor2 <- c("#A2E1F6","#AFFDE8","#E8D1C5")

#define domain order 
domain.order <- c('Bacteria', 'Archaea', 'Eukaryota')

anticodon.per.genome.domain <- anticodon.set.per.genome %>% 
  ggplot(aes(x = factor(Domain, domain.order), y = n.anticodon)) +
  geom_boxplot(lwd = 1,width=0.6,notch = F,color = domainColor,fill = domainColor,outlier.shape = NA)+
  geom_pointrange(aes(y = mean.trna.anticodon, ymin = mean.trna.anticodon - stdev, ymax = mean.trna.anticodon + stdev), color = "black", lwd = 1.5, fatten = 2)+
  geom_text(aes(y = 20, label = paste0('*',round(mean.trna.anticodon,digits = 1))), size=3, color = "grey10", stat = "unique")+
  scale_y_continuous(breaks=c(0,10,20,30,40,50,60,70),limits = c(20,70))+
  theme_classic()+
  # theme_linedraw(base_line_size = 0.2)+
  theme(axis.title = element_text(size = 10, hjust = 0.5), 
        axis.text = element_text(size = 8, colour = "grey15"))+
  ylab("Mean of tRNA anticodon per genome\n(*mean per genome)") + xlab("")

#save box plot
ggplot2::ggsave(paste0(figdir,"/anticodon_per_genome_domain.pdf"), 
                plot =  anticodon.per.genome.domain, 
                device = "pdf",units = "cm",width = 5.5,height = 8) 


#-------------------------------------------------------------------------------------

# median test
archaea.mean.anticodon <- anticodon.set.per.genome %>% filter(str_detect(Domain, "Archaea"))
bacteria.mean.anticodon <- anticodon.set.per.genome %>% filter(str_detect(Domain, "Bacteria"))
eukarya.mean.anticodon<- anticodon.set.per.genome %>% filter(str_detect(Domain, "Eukaryota"))

#median test = wilcoxon test 
arcbac_wilcox <- wilcox.test(archaea.mean.anticodon$n.anticodon, bacteria.mean.anticodon$n.anticodon, paired = F)
arceuk_wilcox <- wilcox.test(archaea.mean.anticodon$n.anticodon, eukarya.mean.anticodon$n.anticodon, paired = F)
eukbac_wilcox <- wilcox.test(bacteria.mean.anticodon$n.anticodon, eukarya.mean.anticodon$n.anticodon, paired = F)

print("Wilcoxon rank sum test results:")
print(paste0('arc vs bac = ',arcbac_wilcox$p.value))
print(paste0('euk vs bac = ',eukbac_wilcox$p.value))
print(paste0('euk vs arc = ',arceuk_wilcox$p.value))
