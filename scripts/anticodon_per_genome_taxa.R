# calculate each taxonomic group mean tRNA gene loci per anticodon per genome 
# this analysis exclude the undetermined tRNAs 
# plot the number of tRNA gene per anticodon as a box plot
# median test = Wilcoxon rank sum test with continuity correction

# import library
library("tidyverse")

# save image directory 
figdir <- '/home/feni/repository/trna_project/figures'

# upload GtRNAdb data
gtrnadb.taxa.data <- read_csv('/home/feni/repository/trna_project/data/gtrnadb_taxonomy_data.csv')

# assigning taxonomic groups on the trna dataframe
gtrnadb.taxa.data$levels <- "level"
gtrnadb.taxa.data$levels[gtrnadb.taxa.data$phylum == "Streptophyta"] <- "Land_plants"
gtrnadb.taxa.data$levels[gtrnadb.taxa.data$kingdom == "Fungi"] <- "Fungi"
gtrnadb.taxa.data$levels[gtrnadb.taxa.data$phylum == "Apicomplexa"] <- "'Protozoa'"
gtrnadb.taxa.data$levels[gtrnadb.taxa.data$phylum == "Euglenozoa"] <- "'Protozoa'"
gtrnadb.taxa.data$levels[gtrnadb.taxa.data$phylum == "Echinodermata"] <- "'Invertebrates'"
gtrnadb.taxa.data$levels[gtrnadb.taxa.data$phylum == "Arthropoda"] <- "'Invertebrates'"
gtrnadb.taxa.data$levels[gtrnadb.taxa.data$phylum == "Mollusca"] <- "'Invertebrates'"
gtrnadb.taxa.data$levels[gtrnadb.taxa.data$phylum == "Nematoda"] <- "'Invertebrates'"
gtrnadb.taxa.data$levels[gtrnadb.taxa.data$phylum == "Chordata"] <- "Vertebrates"
gtrnadb.taxa.data$levels[gtrnadb.taxa.data$phylum == "Fusobacteria"] <- "Fusobacteria"
gtrnadb.taxa.data$levels[gtrnadb.taxa.data$phylum == "Tenericutes"] <- "Tenericutes"
gtrnadb.taxa.data$levels[gtrnadb.taxa.data$phylum == "Proteobacteria"] <- "Proteobacteria"
gtrnadb.taxa.data$levels[gtrnadb.taxa.data$phylum == "Acinobacteria"] <- "Acinobacteria"
gtrnadb.taxa.data$levels[gtrnadb.taxa.data$phylum == "Chloroflexi"] <- "Chloroflexi"
gtrnadb.taxa.data$levels[gtrnadb.taxa.data$phylum == "Euryarchaeota"] <- "Euryarchaeota"
gtrnadb.taxa.data$levels[gtrnadb.taxa.data$phylum == "Thaumarchaeota"] <- "Thaumarchaeota"


# mean tRNA gene loci per anticodon per genome
anticodon.set.per.taxa <- gtrnadb.taxa.data%>%
  filter(!str_detect(levels, "level")) %>% 
  select(Domain,levels,GenomeID,Isotype,Anticodon) %>% 
  filter(!str_detect(Anticodon, "N|M|Y")) %>%
  group_by(Domain, levels, GenomeID) %>%
  distinct()%>%  
  mutate(n.anticodon = n()) %>% 
  select(Domain,levels,GenomeID,n.anticodon)%>% 
  distinct() %>% 
  ungroup() %>% group_by(levels)%>% 
  mutate(mean.trna.anticodon = mean(n.anticodon),
         median.trna.anticodon = median(n.anticodon),
         stdev = sd(n.anticodon),
         error = sd(n.anticodon)/sqrt(length(n.anticodon))) %>% 
  ungroup()

anticodon.set.per.taxa %>%
  group_by(levels) %>% 
  mutate(interquart = IQR(n.anticodon)) %>% 
  select(Domain,mean.trna.anticodon,median.trna.anticodon,stdev,error,interquart) %>% 
  distinct() 

# genome count per level
anticodon.set.per.taxa <- gtrnadb.taxa.data %>% 
  select(Domain,levels,GenomeID) %>%
  distinct() %>% 
  group_by(levels) %>% 
  summarise(n.genomes = n()) %>% 
  ungroup() %>% 
  left_join(anticodon.set.per.taxa, by = "levels")

# box plot 
anticodon.per.genome.taxa <- anticodon.set.per.taxa %>%
  filter(!str_detect(levels, "level")) %>%
  ggplot(aes(x = reorder(levels, mean.trna.anticodon), y = n.anticodon)) +
  geom_boxplot(lwd = .25,width=.5,notch = F,color = "grey80",fill = 'grey20',outlier.shape = 8, outlier.alpha = .5)+
  geom_pointrange(aes(y = mean.trna.anticodon, ymin = mean.trna.anticodon - stdev, ymax = mean.trna.anticodon + stdev), color = "#e28413", lwd = 0.5, fatten = 1)+
  geom_text(aes(y = 21, label = paste0("n = ",n.genomes)), size=3, color = "grey10",stat = "unique")+
  scale_y_continuous(breaks=c(0,10,20,30,40,50,60),limits = c(20,65))+
  theme_classic()+
  # theme_linedraw(base_line_size = 0.2)+
  theme(plot.title = element_text(size = 10,  hjust = 0),
        plot.title.position = "plot", 
        axis.title = element_text(size = 10, hjust = 0.5), 
        axis.text.x = element_text(size = 10, colour = "grey15", angle = 35, hjust = 1),
        axis.text.y = element_text(size = 10, colour = "grey15"))+
  ylab("Mean of tRNA anticodon per genome\n") + xlab("")

# save box plot
ggsave(paste0(figdir,"/isoacceptor_per_genome_taxa.pdf"), 
                plot =  isoacceptor.per.genome.taxa, 
                device = "pdf",units = "cm",width = 15,height = 8)
