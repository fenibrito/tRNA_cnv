#tRNA codon usage per domain

#import library
library("tidyverse")
library("scales")
library("ggpmisc")

# save image directory 
figdir <- '/home/feni/repository/trna_project/figures'

# upload GtRNAdb data
gtrnadb.data <- read_csv('/home/feni/repository/trna_project/data/gtrnadb_data.csv')
domain.order <- c('Eukaryota','Archaea','Bacteria')

# number of genomes per domain
genomes.per.domain <- gtrnadb.data %>% 
  select (GenomeID, Domain) %>%
  distinct() %>%
  group_by(Domain) %>%
  summarise(ngenomes = n())

# anticodon frequency per domain
anticodon.genome.freq <- gtrnadb.data %>% 
  filter(!str_detect(Anticodon, "N|M|Y")) %>% 
  mutate(isoacceptor = paste0(Anticodon,'_',Isotype)) %>% 
  select(Domain, GenomeID, isoacceptor) %>% 
  group_by(Domain, GenomeID, isoacceptor) %>% 
  distinct() %>% 
  ungroup() %>% 
  group_by(Domain, isoacceptor) %>% 
  summarise(genomes.per.isoaccp = n()) %>% 
  left_join(genomes.per.domain, by = 'Domain') %>% 
  mutate(aticodon.per.genomes = (genomes.per.isoaccp/ngenomes)*100) %>% 
  ungroup()


# bar plot
anticodons.per.genome <- anticodon.genome.freq %>% 
  group_by(Domain) %>% 
  summarise(n=n()) %>% 
  ggplot(aes(x = factor(Domain, domain.order), y = n)) +
  geom_bar(stat="identity")+
  geom_text(aes(label = round(n, digits=2)))+
  xlab('')+ ylab('')+
  theme_minimal()

# balloon plot
roworder <- c("TGC_Ala","TCT_Arg","CCT_Arg","CCG_Arg","GTT_Asn","GTC_Asp","GCA_Cys","TTG_Gln","TTC_Glu","TCC_Gly","GCC_Gly","GTG_His","TAA_Leu","TAG_Leu","CAA_Leu","CAG_Leu","TTT_Lys","CAT_Met","GAA_Phe","TGG_Pro","TGA_Ser","GCT_Ser","TGT_Thr","CCA_Trp","GTA_Tyr","TAC_Val","CGC_Ala","TCG_Arg","CTG_Gln","CTC_Glu","CCC_Gly","CTT_Lys","CGG_Pro","CGA_Ser","CGT_Thr","CAC_Val","GGC_Ala","GAT_Ile","CAT_Ile","GAG_Leu","GGG_Pro","GGA_Ser","GGT_Thr","GAC_Val","AGC_Ala","AAT_Ile","TAT_Ile","AAG_Leu","AGG_Pro","AGA_Ser","AGT_Thr","AAC_Val","ACG_Arg","GCG_Arg","ACC_Gly","TCA_SeC","TCA_Sup","ATT_Asn","ATC_Asp","ACA_Cys","ATG_His","AAA_Phe","GCA_SeC","TCC_SeC","CAA_SeC","ACA_SeC","CTA_SeC","TAA_SeC","CCG_SeC","AAA_SeC","ACT_Ser","TTA_Sup","CTA_Sup","ATA_Tyr")
df <- anticodon.genome.freq %>% select(Domain, isoacceptor, aticodon.per.genomes)
domain.order <- c('Eukaryota','Archaea','Bacteria')
balloon <- ggplot(df, aes(y = factor(isoacceptor, rev(roworder)), x = factor(Domain,domain.order), color = aticodon.per.genomes)) +
  geom_point(aes(size=aticodon.per.genomes))+
  # geom_text(aes(label=round(aticodon.per.genomes, digits=2)), alpha=1.0, size=2, color = 'grey') +
  scale_color_gradientn(name = '',
                        colours = c("#264653","#2a9d8f","#e9c46a","#FAD3B3","#EE9781","#902C14"))+
  theme_minimal(base_line_size = 0.2) + xlab("")+ ylab("") +
  theme(axis.text.x = element_text(angle = 90, size = 10, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(angle = 0, size = 8, vjust = 0.5, hjust = 1),
        legend.position = "right", legend.title = element_blank())
balloon

# save plots
ggsave(paste0(figdir,"/anticodons_per_genome_domain_bar.pdf"),plot =  anticodons.per.genome,device = "pdf",units = "cm",width = 3,height = 4)
ggsave(paste0(figdir,"/anticodons_per_genome_domain_balloon.pdf"),plot =  balloon,device = "pdf",units = "cm",width = 6.5,height = 40)

