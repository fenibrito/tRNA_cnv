
#import library
library("tidyverse")
library("RColorBrewer")

# save image directory 
figdir <- '/home/feni/repository/trna_project/figures'

# upload GtRNAdb data
gtrnadb.data <- read_csv('/home/feni/repository/trna_project/data/gtrnadb_data.csv')

#define color palette 
domainColor <- c("#049F76","#118AB2","#BF805F")
domainColor2 <- c("#A2E1F6","#AFFDE8","#E8D1C5")
qual_col_pals <-  brewer.pal.info[brewer.pal.info$category == 'qual',]
anticodonpal <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# define domain order 
domain.order <- c('Bacteria', 'Archaea', 'Eukaryota')

# count n of genomes per species
genomes.by.species <- gtrnadb.data %>% 
  select(Domain, Genome, GenomeID) %>% 
  # filter(str_detect(Genome, '\\w \\w')) %>%
  mutate(myspecie = str_extract(Genome, '\\w* \\w*')) %>% 
  group_by(myspecie) %>%
  distinct() %>% 
  mutate(nGenomes = n()) %>%
  select(Domain,myspecie,GenomeID,nGenomes)

# select all tRNA isotypes from the Archaean species with more than 7 genomes on the database
speciesDF_arc <- genomes.by.species %>% 
  filter(str_detect(Domain, 'Archaea'),
         nGenomes >= 7)

# select all tRNA isotypes from the Bacteria species with more than 48 genomes on the database
speciesDF_bac <- genomes.by.species %>% 
  filter(str_detect(Domain, 'Bacteria'),
         nGenomes >= 48) 

# select all tRNA isotypes from the Eukaryota species with more than 12 genomes on the database
speciesDF_euk <- genomes.by.species %>% 
  filter(str_detect(Domain, 'Eukaryota'),
         nGenomes >= 12)

# data_frame with the species trna data
trna.species <- rbind(speciesDF_arc,speciesDF_bac,speciesDF_euk) %>% 
  left_join(select(gtrnadb.data, GenomeID, Isotype, Anticodon,length),by = 'GenomeID') 

##############################################################################################################
# count tRNA gene per genome per specie
CNV_Species <- trna.species %>% 
  filter(!str_detect(Anticodon, "N|M|Y"))%>% 
  select(Domain,myspecie,GenomeID,nGenomes) %>% 
  group_by(GenomeID) %>% 
  mutate(ntrna_gen = n())  %>% 
  distinct()%>% 
  ungroup() %>% 
  group_by(myspecie) %>% 
  mutate(count_perSpecie = sum(ntrna_gen),
         mean_trna = mean(ntrna_gen), 
         sdtrna = sd(ntrna_gen), 
         erro = sd(ntrna_gen)/sqrt(length(ntrna_gen))) %>%
  ungroup()


# barplot mean trna count per specie
trna.count.per.sp <- CNV_Species %>% 
  select(Domain,myspecie,nGenomes,mean_trna, sdtrna, erro) %>%
  group_by(Domain,myspecie,nGenomes,mean_trna, sdtrna, erro) %>% 
  distinct() %>% 
  ungroup() %>% 
  ggplot() +
  geom_bar(aes(y = myspecie, x = mean_trna, fill = Domain),stat="identity", width=0.7)+
  scale_fill_manual(values = domainColor)+
  geom_errorbar(aes(y = myspecie, x = mean_trna, xmin = mean_trna - sdtrna, xmax = mean_trna + sdtrna),width = .2, size = .6)+
  geom_text(aes(x = 0,y = myspecie, label = paste0("n = ",nGenomes)),hjust = 0,size=2,color = "grey10",stat = "unique")+
  facet_grid(rows = vars(factor(Domain, domain.order)),scales = 'free_y',space = 'free_y',switch = 'y')+
  theme_classic()+
  theme_linedraw(base_line_size = 0.2)+
  theme(axis.title = element_text(size = 10,hjust = 0.5),
        axis.text.y = element_text(angle = 0,hjust = 1),
        legend.position = "none",legend.text = element_text(size = 8),legend.title.align = .5,legend.spacing.y = unit(.2, "char"),legend.key.size = unit(0.3, "cm"))+
  xlab("Mean tRNA gene count per genome\n(n = number of genomes)") + ylab("")

#box and violin plot with trna count distribution per specie
trna.count.per.sp.box <- CNV_Species %>% 
  ggplot() +
  geom_violin(aes(y = myspecie, x = ntrna_gen, fill = Domain, color = Domain), scale = "width")+
  scale_fill_manual(values = domainColor)+
  scale_color_manual(values = domainColor)+
  geom_boxplot(aes(y = myspecie, x = ntrna_gen, fill = Domain),
               outlier.size = .6, outlier.alpha = .8,notch = F,lwd = 0.3,width=0.4)+
  geom_text(aes(y =myspecie, x = 29, label = nGenomes), size = 2.5, stat = 'unique')+
  facet_grid(rows = vars(factor(Domain, domain.order)),scales = 'free_y',space = 'free_y',switch = 'y')+
  # facet_wrap(vars(factor(Domain, domain.order)),scales = 'free',nrow = 3,switch = 'y')+
  theme_classic()+
  theme_linedraw(base_line_size = 0.2)+
  theme(axis.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 0,hjust = 1),
        legend.position = "none",legend.text = element_text(size = 8),legend.title.align = .5,legend.spacing.x = unit(.2, "char"),legend.key.size = unit(0.3, "cm"))+
  xlab("tRNA gene count per genome") + ylab("")


#####################################################################################################
# isoacceptor divesity per specie
# count tRNA isoacceptor set per genome per specie
iscpt.per.genomeSP <- trna.species %>% 
  filter(!str_detect(Isotype, "Sup")) %>%
  select(Domain, myspecie, GenomeID,Anticodon) %>% 
  filter(!str_detect(Anticodon, "N|M|Y")) %>% 
  group_by(Domain,myspecie,GenomeID) %>%
  distinct() %>%  
  mutate(isoacceptors_count = n()) %>%
  select(Domain, myspecie, GenomeID, isoacceptors_count) %>% 
  distinct() %>% ungroup() %>% 
  group_by(myspecie) %>% 
  mutate(isoacceptorsPERsp = sum(isoacceptors_count),
         iqr.isoacceptorsPERsp = IQR(isoacceptors_count),
         mean.isoacceptorsPERsp = mean(isoacceptors_count),
         median.isoacceptorsPERsp = median(isoacceptors_count),
         stdev.isoacceptorsPERsp = sd(isoacceptors_count),
         sterror.isoacceptorsPERsp = sd(isoacceptors_count)/sqrt(length(isoacceptors_count))) %>% 
  ungroup()

#plot
#isoacceptor diversity per specie
isoacc.per.sp <- iscpt.per.genomeSP %>%
  ggplot(aes(y = myspecie, x = isoacceptors_count)) +
  geom_boxplot(outlier.shape = 8,outlier.alpha = .3, outlier.size = 1,notch = F,lwd = .2,width = .8,color = 'black',fill = 'black')+
  geom_pointrange(aes(x = mean.isoacceptorsPERsp,xmin = mean.isoacceptorsPERsp - stdev.isoacceptorsPERsp,xmax = mean.isoacceptorsPERsp + stdev.isoacceptorsPERsp), 
                  color = "#e28413", lwd = .2, fatten = 2)+
  geom_text(aes(x = 15, label = paste0("¨",round(mean.isoacceptorsPERsp, digits = 1))),
            hjust = 0,size=3,color = "grey10",stat = "unique")+
  facet_grid(rows = vars(factor(Domain, domain.order)),scales = 'free_y',space = 'free_y',switch = 'y')+
  scale_x_continuous(breaks=c(0,10,20,30,40,50),limits = c(22,51))+
  theme_classic()+
  theme_linedraw(base_line_size = 0.2)+
  theme(axis.title = element_text(size = 12,hjust = 0.5), 
        axis.text = element_text(size = 10, colour = "grey15"),
        axis.text.y = element_text(angle = 0,hjust = 1))+
  xlab("Number of tRNA isoacceptor per genome") + ylab("")
# ggtitle("tRNA isoacceptor diversity per specie")


##############################################################################################################
# tRNA gene count per isotype per genome per specie
isotype_PERsp <- trna.species %>%
  select(Domain,myspecie,GenomeID,nGenomes,Isotype,Anticodon) %>%
  filter(!str_detect(Anticodon, "N|M|Y")) %>% 
  filter(!str_detect(Isotype, "Sup")) %>% 
  group_by(Domain,myspecie,Isotype) %>%
  mutate(isotype_count = n(),
         mean.iso.sp = isotype_count/nGenomes) %>% 
  select(Domain,myspecie,Isotype,isotype_count,mean.iso.sp) %>% 
  distinct() 

#plot
isotype.sp <- isotype_PERsp %>% 
  ggplot() +
  geom_boxplot(aes(y = myspecie, x = mean.iso.sp), 
               color = 'grey60',fill = 'grey20',outlier.shape = NA,notch = F,lwd = 0.4,width=0.4)+
  geom_jitter(aes(y = myspecie, x = mean.iso.sp, color = Isotype),size = .7,width=.1, alpha = .8)+
  scale_color_manual(values = anticodonpal)+
  facet_grid(rows = vars(factor(Domain, domain.order)),scales = 'free_y',space = 'free_y',switch = 'y')+
  theme_classic()+
  theme_linedraw(base_line_size = 0.2)+
  theme(axis.title = element_text(size = 10,hjust = 0.5),
        axis.text.y = element_text(angle = 0, hjust = 1),
        legend.position = "right",legend.text = element_text(size = 8),legend.title.align = .5,legend.spacing.y = unit(1, "char"),legend.key.size = unit(0.3, "cm"))+
  xlab("Mean tRNA gene count per isotype per genome") + ylab("")

##############################################################################################################
#mean anticodon count per genome block

anticodon.meanPERsp <- trna.species %>%
  select(Domain,myspecie,GenomeID,nGenomes,Isotype,Anticodon) %>%
  filter(!str_detect(Anticodon, "N|M|Y")) %>% 
  filter(!str_detect(Isotype, "Sup")) %>% 
  group_by(Domain,myspecie,Anticodon) %>%
  mutate(isoacc_count = n(),
         mean.isoacc.sp = isoacc_count/nGenomes)%>% 
  select(Domain,myspecie,Anticodon,isoacc_count,mean.isoacc.sp) %>% 
  distinct() %>% ungroup()


anticodon.sp <- anticodon.meanPERsp %>% 
  ggplot() +
  geom_boxplot(aes(y = myspecie, x = mean.isoacc.sp), 
               color = 'grey60',
               fill = 'grey20',
               outlier.shape = NA,  
               notch = F,
               lwd = 0.4,
               width=0.4)+
  geom_jitter(aes(y = myspecie, x = mean.isoacc.sp, color = Anticodon),width=.2)+
  scale_color_manual(values = anticodonpal)+
  facet_grid(rows = vars(factor(Domain, domain.order)),scales = 'free_y',space = 'free_y',switch = 'y')+
  theme_classic()+
  theme_linedraw(base_line_size = 0.2)+
  theme(axis.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.y = element_text(angle = 0, hjust = 1),
        legend.position = "right",
        legend.text = element_text(size = 8),
        legend.title.align = .5,
        legend.spacing.y = unit(1, "char"),
        legend.key.size = unit(0.3, "cm"))+
  xlab("Mean tRNA gene count per isoacceptor per genome") + ylab("")


specie.color <- c("#7e637a","#9d4edd","#e06100","#38574d","#e07a5f","#b1c5fc","#9A031E","#adc178","#F5A3C1","#0d42d3","#643047","#8ac926","#0d3b66","#f4e285")
specie.list <- c("Staphylococcus aureus","Salmonella enterica","Mycobacterium tuberculosis","Listeria monocytogenes","Helicobacter pylori","Escherichia coli","Chlamydia trachomatis","Burkholderia pseudomallei","Sulfolobus islandicus","Methanosarcina mazei","Saccharomyces cerevisiae","Fusarium oxysporum","Cryptococcus gattii","Candida albicans")
abrv.list <- c("Sa","Se","Mt","Lm","Hp","Ec","Ct","Bp","Si","Mm","Sc","Fo","Cg","Ca")
cnv_intra.sp <- trna.species %>% 
  select(Domain,myspecie,GenomeID,nGenomes,Isotype,Anticodon) %>%
  filter(!str_detect(Anticodon, "N|M|Y"),
         !str_detect(Isotype, "Sup")) %>%
  group_by(Domain,myspecie,GenomeID,Isotype) %>%
  mutate(isotype_count = n()) %>% 
  select(Domain,myspecie,GenomeID,Isotype,isotype_count) %>% 
  distinct() %>%
  mutate(spname = myspecie) %>% 
  separate(spname, into = c('name','last'),sep=' ') %>%
  mutate(abrv = paste0(str_extract(name,'^\\w'),str_extract(last,'^\\w'))) %>%
  select(Domain,myspecie,GenomeID,Isotype,isotype_count,abrv) 
  # filter(str_detect(myspecie, 'Escherichia coli|Saccharomyces cerevisiae')) %>%

cnv_intra_sp_plot <- cnv_intra.sp %>% 
  ggplot()+
  geom_violin(aes(x = factor(abrv,abrv.list), y = isotype_count, fill = factor(myspecie,specie.list), color = factor(myspecie,specie.list)), scale = "width")+
  geom_jitter(aes(x = factor(abrv,abrv.list), y = isotype_count, fill = factor(myspecie,specie.list)), shape = 21, color = "white", position = position_dodge(1.0))+
  scale_color_manual(values = specie.color, name = '')+
  scale_fill_manual(values = specie.color, name = '')+
  facet_wrap(vars(Isotype),scales = 'free',ncol = 3)+
  theme_minimal()+
  ylab("tRNA gene count per isotype per genome") + xlab("")+
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(angle = 90),
        legend.position = 'bottom')


# save plots
ggsave(paste0(figdir,"/trna_count_per_sp.pdf"),plot = trna.count.per.sp, device = "pdf",units = "cm",width = 14,height = 10)
ggsave(paste0(figdir,"/trna_count_per_sp_box.pdf"),plot = trna.count.per.sp.box, device = "pdf",units = "cm",width = 14,height = 10)
ggsave(paste0(figdir,"/isoacc_per_sp.pdf"),plot = isoacc.per.sp, device = "pdf",units = "cm",width = 14,height = 10)
ggsave(paste0(figdir,"/isotype_sp.pdf"),plot = isotype.sp, device = "pdf",units = "cm",width = 18,height = 12)
ggsave(paste0(figdir,"/anticodon_sp.pdf"),plot = anticodon.sp, device = "pdf",units = "cm",width = 18,height = 12)
ggsave(paste0(figdir,"/cnv_intra_sp_plot.pdf"),plot =  cnv_intra_sp_plot,device = "pdf",units = "cm",width = 25,height = 35)
