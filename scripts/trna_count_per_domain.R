#calculate the number of tRNA loci per genome per domain
# this analysis exclude the undetermined tRNAs 
# plot the mean of tRNA gene count as a bar plot
# plot the distribution of tRNA gene count as a box plot
# mean test = student t test
# median test = Wilcoxon rank sum test with continuity correction

# import library
library("tidyverse")
library("scales")

# set image directory 
figdir <- '/home/feni/repository/trna_project/figures'

# upload GtRNAdb data
gtrnadb.data <- read_csv('/home/feni/repository/trna_project/data/gtrnadb_data.csv')

# count tRNA gene loci per genome 
trna.per.domain <- gtrnadb.data %>% 
  filter(!str_detect(Anticodon, "N|M|Y"))%>% 
  select(Domain,GenomeID) %>% 
  group_by(GenomeID) %>% 
  mutate(n.trna.per.genome = n())  %>% 
  distinct() %>% 
  filter(n.trna.per.genome >= 10) %>% 
  ungroup() %>% 
  group_by(Domain) %>% 
  mutate(mean.trna.per.genome = mean(n.trna.per.genome),
         sdtrna = sd(n.trna.per.genome),
         erro = sd(n.trna.per.genome)/sqrt(length(n.trna.per.genome))) %>%
  ungroup()

trna.per.domain %>% 
  group_by(Domain) %>% 
  mutate(medi = median(n.trna.per.genome),
         interquart = IQR(n.trna.per.genome)) %>% 
  select(Domain,mean.trna.per.genome,sdtrna,erro,medi,interquart) %>% 
  distinct()

# define color palette
domainColor <- c("#049F76","#118AB2","#BF805F")
domainColor2 <- c("#A2E1F6","#AFFDE8","#E8D1C5")

# define domain order 
domain.order <- c('Bacteria', 'Archaea', 'Eukaryota')

# plot mean tRNA gene per genome per domain
tibble.trna.per.domain <- tibble(trna.per.domain)
tibble.trna.per.domain$Domain <- factor(trna.per.domain$Domain, levels = domain.order)

mean.trna.domain.barplot <- tibble.trna.per.domain %>% 
  select(Domain, mean.trna.per.genome, sdtrna, erro) %>%
  distinct() %>% 
  ggplot() +
  geom_bar(ggplot2::aes(x = factor(Domain, domain.order), y = mean.trna.per.genome, fill = Domain),
           stat="identity",width=0.7)+
  geom_errorbar(ggplot2::aes(x = Domain, y = mean.trna.per.genome, ymin = mean.trna.per.genome - erro, ymax = mean.trna.per.genome + erro),
                width = .2, size = .6)+
  scale_fill_manual(values = domainColor)+
  theme_classic()+
  theme(axis.title = element_text(size = 10, hjust = 0.5),
        legend.position = "bottom",legend.text = element_text(size = 8),legend.title.align = .5,
        legend.spacing.x = unit(.2, "char"),legend.key.size = unit(0.3, "cm"))+
  ylab("Mean tRNA gene count per genome\n") + xlab("")

# boxplot trna count distribution with boxplot
trna.count.per.domain.boxplot <- trna.per.domain %>% 
  ggplot() +
  geom_violin(aes(x = factor(Domain, domain.order), y = n.trna.per.genome, fill = Domain), scale = "width", color = NA)+
  scale_fill_manual(values = domainColor2)+
  geom_boxplot(aes(x = factor(Domain, domain.order), y = n.trna.per.genome),
               notch = T,lwd = 0.3,width=0.4,color = domainColor,fill = domainColor)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  theme_classic()+
  theme(axis.title = element_text(size = 10,hjust = 0.5), 
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.title.align = .5,
        legend.spacing.x = unit(.2, "char"),
        legend.key.size = unit(0.3, "cm"))+
  ylab("Log10 (tRNA gene count per genome)\n") + xlab("")

#save barplot mean.trna.domain.barplot
ggplot2::ggsave(paste0(figdir,"/mean_trna_domain_barplot.pdf"), 
                plot =  mean.trna.domain.barplot, 
                device = "pdf",units = "cm",width = 8.5,height = 8) 

#save box plot trna.count.per.domain.boxplot
ggplot2::ggsave(paste0(figdir,"/trna_count_per_domain_boxplot.pdf"), 
                plot =  trna.count.per.domain.boxplot, 
                device = "pdf",units = "cm",width = 8.5,height = 8) 

#-------------------------------------------------------------------------------------

# student t test
hist(trna.per.domain$n.trna.per.genome,main = "average tRNA count per genome",xlab = "average tRNA count per genome")


archaea.mean.trna.count <- trna.per.domain %>% filter(str_detect(Domain, "Archaea"))
bacteria.mean.trna.count <- trna.per.domain %>% filter(str_detect(Domain, "Bacteria"))
eukarya.mean.trna.count <- trna.per.domain %>% filter(str_detect(Domain, "Eukaryota"))

#mean test = student t test
arcbac_ttest <- t.test(archaea.mean.trna.count$n.trna.per.genome, bacteria.mean.trna.count$n.trna.per.genome, paired = F)
arceuk_ttest <- t.test(archaea.mean.trna.count$n.trna.per.genome, eukarya.mean.trna.count$n.trna.per.genome, paired = F)
eukbac_ttest <- t.test(bacteria.mean.trna.count$n.trna.per.genome, eukarya.mean.trna.count$n.trna.per.genome, paired = F)

print("student t test results:")
print(paste0('arc vc bac = ',arcbac_ttest$p.value))
print(paste0('arc vc euk = ',arceuk_ttest$p.value))
print(paste0('euk vc bac = ',eukbac_ttest$p.value))

#median test = wilcoxon test 
arcbac_wilcox <- wilcox.test(archaea.mean.trna.count$n.trna.per.genome, bacteria.mean.trna.count$n.trna.per.genome, paired = F)
arceuk_wilcox <- wilcox.test(archaea.mean.trna.count$n.trna.per.genome, eukarya.mean.trna.count$n.trna.per.genome, paired = F)
eukbac_wilcox <- wilcox.test(bacteria.mean.trna.count$n.trna.per.genome, eukarya.mean.trna.count$n.trna.per.genome, paired = F)

print("Wilcoxon rank sum test results:")
print(paste0('arc vc bac = ',arcbac_wilcox$p.value))
print(paste0('arc vc euk = ',arceuk_wilcox$p.value))
print(paste0('euk vc bac = ',eukbac_wilcox$p.value))

