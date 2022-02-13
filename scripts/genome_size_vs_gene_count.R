#compare the number of tRNA gene and protein-coding gene with the genome size 


library("tidyverse")
library("reshape2")
library("scales")
library("ggpmisc")

# save image directory 
figdir <- '/home/feni/repository/trna_project/figures'

# upload GtRNAdb data
gtrnadb.taxa.data <- read_csv('/home/feni/repository/trna_project/data/gtrnadb_taxonomy_data.csv')

# number of canonical tRNA
canon.trna <- gtrnadb.taxa.data %>% 
  filter(!str_detect(Anticodon, "N|M|Y"))%>% 
  select(Domain,GenomeID) %>% 
  group_by(GenomeID) %>% 
  mutate(n.canon.trna.per.genome = n())  %>% 
  distinct() %>% 
  # filter(n.trna.per.genome >= 10) %>% 
  ungroup() 

# number if pseudo tRNAs
ntrna.per.genome <- gtrnadb.taxa.data %>% 
  # filter(str_detect(Anticodon, "N|M|Y"))%>% 
  select(Domain,GenomeID) %>% 
  group_by(GenomeID) %>% 
  mutate(n.trna.per.genome = n())  %>% 
  distinct() %>% 
  # filter(n.trna.per.genome >= 10) %>% 
  ungroup() %>% 
  left_join(canon.trna,by=c("Domain","GenomeID")) %>% 
  mutate(pseudo.trna.per.genome = n.trna.per.genome - n.canon.trna.per.genome) %>% 
  left_join(select(gtrnadb.taxa.data,Domain,GenomeID,kingdom,phylum),by=c('Domain','GenomeID')) %>% 
  distinct()

# import genome information
data_genomes <- read_tsv('/home/feni/repository/trna_project/data/genomas_info_edit.txt', na = "NA",
                         col_select = c('Domain','SubGroup','GenomeID','size_gen','Genes','ntrna_gen','GC','filt'))

data_genomes <- data_genomes %>% filter(str_detect(filt, 'ok'))

trna.data <- left_join(ntrna.per.genome,data_genomes,by = c("GenomeID", "Domain")) %>%
  filter(!is.na(ntrna_gen))


# assigning taxonomic groups on the trna dataframe
trna.data$levels <- trna.data$phylum
trna.data$levels[trna.data$kingdom == "Fungi"] <- "Fungi"
trna.data$levels[trna.data$phylum == "Apicomplexa"] <- "'Protozoa'"
trna.data$levels[trna.data$phylum == "Euglenozoa"] <- "'Protozoa'"
trna.data$levels[trna.data$phylum == "Echinodermata"] <- "'Invertebrates'"
trna.data$levels[trna.data$phylum == "Arthropoda"] <- "'Invertebrates'"
trna.data$levels[trna.data$phylum == "Mollusca"] <- "'Invertebrates'"
trna.data$levels[trna.data$phylum == "Nematoda"] <- "'Invertebrates'"
trna.data$levels[trna.data$phylum == "Chordata"] <- "Vertebrates"
trna.data$levels[trna.data$Domain == "Bacteria"] <- "Bacteria"
trna.data$levels[trna.data$Domain == "Archaea"] <- "Archaea"
trna.data$levels[trna.data$phylum == "Streptophyta"] <- "Land Plants"


leg_level <- c("Bacteria", "Archaea", "Fungi","'Protozoa'", "'Invertebrates'", "Vertebrates", "Land Plants")
colorpolot <- c("#9699B0","#27302e","#985277","#63CAE3","#c18c5d","#bb4430","#208b3a")

# scatter plot
# all domains included (Bateria, Archaea and Eukaryota)
p1 <- trna.data %>%  
  select(levels,size_gen,n.canon.trna.per.genome,Genes)%>%
  melt(id = c('size_gen','levels')) %>% 
  ggplot(aes(x = size_gen, y = value, mapping = variable)) +
  geom_point(aes(colour = levels, shape = variable),size = 2, stat = 'identity', alpha = 1)+
  geom_smooth(method = "lm", se=T, size = 0.5, formula = y ~ x,color = "gray30")+
  scale_x_log10(name = "Log10(Genome size[Mb])", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(name = "Log10(Gene count per genome)", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(name = "", breaks = leg_level, labels = leg_level, values = colorpolot)+
  scale_shape_manual(name = "", values=c(18,20), breaks = c("Genes","n.canon.trna.per.genome"), labels = c("Protein-coding genes","tRNA genes"))+
  theme_classic(base_line_size = 0.3) +
  theme(axis.text = element_text(size = 8),axis.title = element_text(size = 10),
        legend.text=element_text(size=10),legend.key.height= unit(0.4, 'cm'),legend.key.width= unit(0.4, 'cm'),legend.background=element_blank())+
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = T,rr.digits = 2,size = 3,label.y = c(0.9,0.95))

p1


# only eukaryotes genomes included
p2 <- trna.data %>%  
  filter(!str_detect(levels,c("Bacteria|Archaea"))) %>% 
  select(levels,size_gen,n.canon.trna.per.genome,Genes)%>%
  melt(id = c('size_gen','levels')) %>% 
  ggplot(aes(x = size_gen, y = value, mapping = variable)) +
  geom_point(aes(colour = levels, shape = variable),size = 3, stat = 'identity', alpha = 1)+
  geom_smooth(method = "lm", se=T, size = 0.5, formula = y ~ x,color = "gray30")+
  scale_x_log10(name = "Log10(Genome size[Mb])", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(name = "Log10(Gene count per genome)", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(name = "", breaks = leg_level, labels = leg_level, values = colorpolot)+
  scale_shape_manual(name = "", values=c(18,20), breaks = c("Genes","n.canon.trna.per.genome"), labels = c("Protein-coding genes","tRNA genes"))+
  theme_classic(base_line_size = 0.3) +
  theme(axis.text = element_text(size = 8),axis.title = element_text(size = 10),
        legend.text=element_text(size=10),legend.key.height= unit(0.4, 'cm'),legend.key.width= unit(0.4, 'cm'),legend.background=element_blank())+
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = T,rr.digits = 2,size = 3,label.y = c(0.9,0.95))
p2


# only Bacteria genomes included
bac_level <- c("Kin_FCB group","Kin_of_Proteobacteria","Kin_Terrabacteria group","Kin_of_Spirochaetes","Kin_unclassified Bacteria","Kin_PVC group","Kin_of_Fusobacteria","Kin_of_Thermotogae","Kin_of_Synergistetes","Kin_of_Aquificae","Kin_of_Acidobacteria")
bac_leg <- c("Bacteroidetes","Proteobacteria","Terrabacteria group","Spirochaetes","Unclassified Bacteria","PVC group","Fusobacteria","Thermotogae","Synergistetes","Aquificae","Acidobacteria")
bac_color <- c("#355070","#f9c74f","#8367c7","#007F5F","#f15bb5","#f94144","#f3722c","#370617","#f8961e","#90be6d","#43aa8b")

p3 <- trna.data %>%  
  filter(str_detect(Domain,c("Bacteria"))) %>% 
  select(kingdom,size_gen,n.canon.trna.per.genome,Genes)%>%
  melt(id = c('size_gen','kingdom')) %>% 
  ggplot(aes(x = size_gen, y = value, mapping = variable)) +
  geom_point(aes(colour = kingdom, shape = variable),size = 3, stat = 'identity', alpha = 1)+
  geom_smooth(method = "lm", se=T, size = 0.5, formula = y ~ x,color = "gray30")+
  scale_x_log10(name = "Log10(Genome size[Mb])", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(name = "Log10(Gene count per genome)", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(name = "", breaks = bac_level, labels = bac_leg, values = bac_color)+
  scale_shape_manual(name = "", values=c(18,20), breaks = c("Genes","n.canon.trna.per.genome"), labels = c("Protein-coding genes","tRNA genes"))+
  theme_classic(base_line_size = 0.3) +
  theme(axis.text = element_text(size = 8),axis.title = element_text(size = 10),
        legend.text=element_text(size=10),legend.key.height= unit(0.4, 'cm'),legend.key.width= unit(0.4, 'cm'),legend.background=element_blank())+
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = T,rr.digits = 2,size = 3,label.y = c(0.9,0.95))
p3

# only archaea genomes included
arc_level <- c("Methanomada group","DHVE2 group","Methanopyri","Thermoplasmata","Thermococci","Stenosarchaea group","Crenarchaeota","Candidatus Korarchaeota","Thaumarchaeota")
arc_color <- c("#f9c74f","#90be6d","#43aa8b","#007f5f","#577590","#355070","#f15bb5","#f94144","#f3722c")

p4 <- trna.data %>%  
  filter(str_detect(Domain,c("Archaea"))) %>% 
  select(SubGroup,size_gen,n.canon.trna.per.genome,Genes) %>% 
  melt(id = c("size_gen","SubGroup")) %>% 
  ggplot(aes(x = size_gen, y = value, mapping = variable)) +
  geom_smooth(color = "black", method = "lm", se=T, size = 0.5, formula = y ~ x)  +
  geom_point(aes(colour = SubGroup, shape = variable),size = 3, stat = 'identity', alpha = 1)+
  scale_x_log10(name = "Log10(Genome size[Mb])", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(name = "Log10(Gene count)", limits=c(15,11000), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(name = "", breaks = arc_level, labels = arc_level, values = bac_color)+
  scale_shape_manual(name = "", values=c(18,20), breaks = c("Genes","n.canon.trna.per.genome"), labels = c("Protein-coding genes","tRNA genes"))+
  theme_classic(base_line_size = 0.3) +
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text=element_text(size=8),
        legend.key.height= unit(0.4, 'cm'),
        legend.key.width= unit(0.4, 'cm'),
        legend.background=element_blank())+
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..rr.label.., sep = "~~~")),
               parse = T,
               rr.digits = 2,
               size = 3,
               label.y = c(0.9,0.95))
p4

ntrna.vs.protcodgenes <- trna.data %>%  
  select(levels,n.trna.per.genome,Genes) %>% 
  ggplot(aes(x = Genes, y = n.trna.per.genome)) +
  geom_smooth(color = "black", method = "lm", se=T, size = 0.5, formula = y ~ exp(x))  +
  geom_point(aes(colour = levels),size = 2, stat = 'identity', alpha = 1)+
  scale_x_log10(name = "Log10(Number of protein-coding genes per genome)", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(name = "Log10(Number of tRNA genes per genome)", limits=c(15,11000), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(name = "", breaks = leg_level, labels = leg_level, values = colorpolot)+
  theme_classic(base_line_size = 0.3) +
  theme(axis.text = element_text(size = 8),axis.title = element_text(size = 10),
        legend.text=element_text(size=8),legend.key.height= unit(0.4, 'cm'),legend.key.width= unit(0.4, 'cm'),legend.background=element_blank())+
  stat_poly_eq(formula = y ~ exp(x), 
               aes(label = paste(..rr.label.., sep = "~~~")),
               parse = T,rr.digits = 2,size = 4,label.y = c(0.9,0.95))
ntrna.vs.protcodgenes

domainColor <- c("#049F76","#1FB7EA","#BF805F")
domainColor2 <- c("#A2E1F6","#AFFDE8","#E8D1C5")
domain.order <- c('Bacteria', 'Archaea', 'Eukaryota')

# gene density
# Protein coding gene / Mb
pcgene.mb <- trna.data %>% 
  select(Domain,GenomeID,n.canon.trna.per.genome,size_gen,Genes,levels) %>% 
  mutate(trna_dens = n.canon.trna.per.genome/size_gen,
         gene_dens = Genes/size_gen) %>% 
  ggplot(aes(x=gene_dens, group=Domain, color = Domain, fill=Domain)) +
  geom_density(adjust=1.5, alpha=.5) +
  scale_color_manual(name = "", values = domainColor)+
  scale_fill_manual(name = "", values = domainColor)+
  theme_minimal()+ylab('')+xlab('Protein coding gene/Mb')

# tRNA / Mb
trna.mb <- trna.data %>% 
  select(Domain,GenomeID,n.canon.trna.per.genome,size_gen,Genes,levels) %>% 
  mutate(trna_dens = n.canon.trna.per.genome/size_gen,
         gene_dens = Genes/size_gen) %>%
  ggplot(aes(x=trna_dens, group=Domain, fill=Domain, color = Domain)) +
  geom_density(adjust=1.5, alpha=.5) +
  scale_color_manual(name = "", values = domainColor)+
  scale_fill_manual(name = "", values = domainColor)+
  theme_minimal()+ylab('')+xlab('tRNA/Mb')

euk.trna.mb <- trna.data %>% 
  select(Domain,GenomeID,n.canon.trna.per.genome,size_gen,Genes,levels) %>% 
  filter(!str_detect(levels,"Bacteria|Archaea")) %>%
  mutate(trna_dens = n.canon.trna.per.genome/size_gen,
         gene_dens = Genes/size_gen) %>% 
  group_by(levels) %>% 
  mutate(mediana = median(trna_dens)) %>% ungroup() %>% 
  ggplot(aes(x=trna_dens, y = reorder(levels,-mediana), group=levels, fill=levels, color = levels)) +
  geom_jitter(height = .25)+
  geom_boxplot(aes(x = trna_dens),fill = 'NA',
               outlier.size = 0.4,outlier.shape = 23,outlier.colour = "black",outlier.fill = "black",outlier.stroke = 1.5,outlier.alpha = .2,
               notch = F,color='grey10',width=0.8,lwd = .5)+
  scale_color_manual(name = "", breaks = leg_level, labels = leg_level, values = colorpolot)+
  theme(legend.title = element_blank())+
  theme_minimal()+
  theme(legend.position = 'none')+
  xlab('Eukariotes tRNA/Mb') + ylab('')

# save plots
ggsave(paste0(figdir,"/ntrna_vs_genomeSize_domains_.pdf"),plot =  p1,device = "pdf",units = "cm",width = 15,height = 10)
ggsave(paste0(figdir,"/ntrna_vs_genomeSize_euk.pdf"),plot =  p2,device = "pdf",units = "cm",width = 15,height = 10)
ggsave(paste0(figdir,"/ntrna_vs_genomeSize_bac.pdf"),plot =  p3,device = "pdf",units = "cm",width = 15,height = 10)
ggsave(paste0(figdir,"/ntrna_vs_genomeSize_arc.pdf"),plot =  p4,device = "pdf",units = "cm",width = 15,height = 10)
ggsave(paste0(figdir,"/ntrna_vs_protcodgenes.pdf"),plot =  ntrna.vs.protcodgenes,device = "pdf",units = "cm",width = 15,height = 10)
ggsave(paste0(figdir,"/pcgene_mb.pdf"),plot =  pcgene.mb,device = "pdf",units = "cm",width = 10,height = 3)
ggsave(paste0(figdir,"/trna_mb.pdf"),plot =  trna.mb,device = "pdf",units = "cm",width = 10,height = 3)
ggsave(paste0(figdir,"/euk_trna_mb.pdf"),plot =  euk.trna.mb,device = "pdf",units = "cm",width = 10,height = 5)
