# compare the number of tRNA anticodon with other enomic features such as number of protein-coding genes, genome size and CG content


# import library
library("tidyverse")
library("reshape2")
library("scales")
library("ggpmisc")

# save image directory 
figdir <- '/home/feni/repository/trna_project/figures'

# upload GtRNAdb data
gtrnadb.taxa.data <- read_csv('/home/feni/repository/trna_project/data/gtrnadb_taxonomy_data.csv')

# count the number tRNA anticodon per genome
trna.anticodon <- gtrnadb.taxa.data %>%
  filter(!str_detect(Anticodon, "N|M|Y")) %>% 
  select(Domain, Genome, GenomeID, Isotype, Anticodon)%>% 
  group_by(Domain, GenomeID) %>%
  distinct() %>%  
  mutate(n.anticodon = n()) %>%
  select(Domain,Genome,GenomeID,n.anticodon) %>%
  distinct()


ntrna.per.genome <- gtrnadb.taxa.data %>% 
  filter(!str_detect(Anticodon, "N|M|Y"))%>%
  select(Domain,GenomeID) %>% 
  group_by(GenomeID) %>% 
  mutate(n.trna.per.genome = n())  %>% 
  distinct() %>% 
  ungroup() %>% 
  left_join(trna.anticodon,by=c("Domain","GenomeID")) %>% 
  left_join(select(gtrnadb.taxa.data,Domain,GenomeID,kingdom,phylum),by=c('Domain','GenomeID')) %>% 
  distinct()


# import genome information
data_genomes <- read_tsv('/home/feni/repository/trna_project/data/genomas_info_edit.txt', na = "NA",
                         col_select = c('Domain','SubGroup','GenomeID','size_gen','Genes','ntrna_gen','GC','filt'))

data_genomes <- data_genomes %>% filter(str_detect(filt, 'ok'))

trna.data <- left_join(ntrna.per.genome,data_genomes,by = c("GenomeID", "Domain")) %>%
  filter(!is.na(ntrna_gen))


#assigning the taxonomic groups on the trna dataframe
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

# scatter plot of tRNA anticodon count vs tRNA gene count per genome
anticodon.vs.ntrna <- trna.data %>%  
  select(levels,size_gen,n.trna.per.genome,n.anticodon,Genes) %>% 
  ggplot(aes(x = n.trna.per.genome, y = n.anticodon))+
  geom_point(aes(y = n.anticodon, color = levels), size = 1, stat = 'identity', alpha = 0.8)+
  geom_smooth(aes(y = n.anticodon), color = 'black', method="lm", se=T, fullrange=F, formula = y ~ x, size = 0.5)+
  scale_x_log10(name = "Log10(tRNA gene count per genome)", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(name = "", breaks = leg_level, labels = leg_level, values = colorpolot)+
  theme_classic(base_line_size = 0.3) +ylab("Number of tRNA anticodon per genome")+
  theme(axis.text = element_text(size = 8),axis.title = element_text(size = 10),
        legend.text=element_text(size=10),legend.key.height= unit(0.4, 'cm'),legend.key.width= unit(0.4, 'cm'),legend.background=element_blank())+
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = T,rr.digits = 2,size = 3,label.y = c(0.9,0.95))


# scatter plot of tRNA anticodon count vs genome size
anticodon.vs.genomesize <- trna.data %>%  
  select(levels,size_gen,n.trna.per.genome,n.anticodon,Genes) %>% 
  ggplot(aes(x = size_gen, y = n.anticodon))+
  geom_point(aes(y = n.anticodon, color = levels), size = 1, stat = 'identity', alpha = 0.8)+
  geom_smooth(aes(y = n.anticodon), color = 'black', method="lm", se=T, fullrange=F, formula = y ~ x, size = 0.5)+
  scale_x_log10(name = "Log10(Genome size[Mb]", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(name = "", breaks = leg_level, labels = leg_level, values = colorpolot)+
  theme_classic(base_line_size = 0.3) +ylab("Number of tRNA anticodon per genome")+
  theme(axis.text = element_text(size = 8),axis.title = element_text(size = 10),
        legend.text=element_text(size=10),legend.key.height= unit(0.4, 'cm'),legend.key.width= unit(0.4, 'cm'),legend.background=element_blank())+
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = T,rr.digits = 2,size = 3,label.y = c(0.9,0.95))

# scatter plot of tRNA anticodon count vs protein-coding genes
anticodon.vs.protcodgenes <- trna.data %>%  
  select(levels,size_gen,n.trna.per.genome,n.anticodon,Genes) %>% 
  ggplot(aes(x = Genes, y = n.anticodon))+
  geom_point(aes(y = n.anticodon, color = levels), size = 1, stat = 'identity', alpha = 0.8)+
  geom_smooth(aes(y = n.anticodon), color = 'black', method="lm", se=T, fullrange=F, formula = y ~ exp(x), size = 0.5)+
  scale_x_log10(name = "Log10(Number of protein-coding genes per genome)", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(name = "", breaks = leg_level, labels = leg_level, values = colorpolot)+
  theme_classic(base_line_size = 0.3) + ylab("Number of tRNA isoacceptor per genome")+
  theme(axis.text = element_text(size = 8),axis.title = element_text(size = 10),
        legend.text=element_text(size=10),legend.key.height= unit(0.4, 'cm'),legend.key.width= unit(0.4, 'cm'),legend.background=element_blank())+
  stat_poly_eq(formula = y ~ exp(x), 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = T,rr.digits = 2,size = 3,label.y = c(0.9,0.95))



ggsave(paste0(figdir,"/anticodon_vs_ntrna.pdf"),plot =  anticodon.vs.ntrna,device = "pdf",units = "cm",width = 13,height = 10)
ggsave(paste0(figdir,"/anticodon_vs_genomesize.pdf"),plot =  anticodon.vs.genomesize,device = "pdf",units = "cm",width = 13,height = 10)
ggsave(paste0(figdir,"/anticodon_vs_protcodgenes.pdf"),plot =  anticodon.vs.protcodgenes,device = "pdf",units = "cm",width = 13,height = 10)
