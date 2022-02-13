#compare the canonical tRNA gene count vs the tRNA-peseudo gene per genome size

#import library
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
  ungroup() 

# number if pseudo tRNAs
ntrna.per.genome <- gtrnadb.taxa.data %>% 
  # filter(str_detect(Anticodon, "N|M|Y"))%>% 
  select(Domain,GenomeID) %>% 
  group_by(GenomeID) %>% 
  mutate(n.trna.per.genome = n())  %>% 
  distinct() %>% 
  ungroup() %>% 
  left_join(canon.trna,by=c("Domain","GenomeID")) %>% 
  mutate(pseudo.trna.per.genome = n.trna.per.genome - n.canon.trna.per.genome) %>% 
  left_join(select(gtrnadb.taxa.data,Domain,GenomeID,kingdom,phylum),by=c('Domain','GenomeID')) %>% 
  distinct()

# import genome information
data_genomes <- read_tsv('/home/feni/repository/trna_project/data/genomas_info_edit.txt', na = "NA",
                         col_select = c('Domain','GenomeID','size_gen','Genes','ntrna_gen','GC','filt'))

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


# tRNA gene count per gene type (canonical tRNA and pseudo-tRNA)
n.trna.type.per.genome <- trna.data %>%
  select(size_gen,n.trna.per.genome,n.canon.trna.per.genome,pseudo.trna.per.genome) %>% 
  reshape2::melt(id = "size_gen") %>% 
  filter(str_detect(variable, "n.canon.trna.per.genome|pseudo.trna.per.genome")) %>% 
  filter(!value == 0) %>% 
  ggplot(aes(x = size_gen, y = value, colour = variable, shape = variable)) +
  geom_smooth(method = "lm", se=T, size = 0.5, formula = y ~ x)+
  geom_point(aes(colour = variable, shape = variable),size = 1.5, stat = 'identity', alpha = 1)+
  scale_x_log10(name = "Log10(Genome size[Mb])", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(name = "Log10(tRNA genes per genome)", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(name = "", breaks = c("n.canon.trna.per.genome","pseudo.trna.per.genome"), labels = c("tRNA", "pseudo_tRNAs") , values = c("#293241","#2A9D8F","#a4161a"))+
  scale_shape_manual(name = "", values=c(18,20), breaks = c("n.canon.trna.per.genome","pseudo.trna.per.genome"), labels = c("tRNA", "pseudo_tRNAs"))+
  theme_classic(base_line_size = 0.3) +
  theme(axis.text = element_text(size = 8),axis.title = element_text(size = 10),
        legend.text=element_text(size=10),legend.key.height= unit(0.4, 'cm'),legend.key.width= unit(0.4, 'cm'),legend.background=element_blank(),legend.position = c(0.2, 0.78))+
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = T,
               rr.digits = 3,
               size = 4,
               label.y = c(0.9,0.85))

# save scatter plot
ggsave(paste0(figdir,"/n_trna_type_per_genome.pdf"), 
       plot =  n.trna.type.per.genome, 
       device = "pdf",units = "cm",width = 10,height = 10)