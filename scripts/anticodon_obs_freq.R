#tRNA codon usage per domain

#import library
# if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn", force = T)

# install.packages("easypackages")
library(easypackages)
source("~/repository/trna_project/scripts/functions.R")
packs <- c("tidyverse","scales","ggpmisc","pheatmap","RColorBrewer","ggvenn")
libraries(packs)


# save image directory 
figdir <- '/home/feni/repository/trna_project/figures'

# upload GtRNAdb data
gtrnadb.data <- read_csv('/home/feni/repository/trna_project/data/gtrnadb_data.csv')

# define color palette
domain.order <- c('Eukaryota','Archaea','Bacteria')
domainColor <- c("#049F76","#118AB2","#BF805F")

# number of genomes per domain
genomes.per.domain <- gtrnadb.data %>% 
  select (GenomeID, Domain) %>%
  distinct() %>%
  group_by(Domain) %>%
  summarise(ngenomes = n())

# anticodon frequency per genome
anticodon.genome.freq <- gtrnadb.data %>% 
  filter(!str_detect(Anticodon, "N|M|Y")) %>% 
  mutate(comb = paste0(Anticodon,'_',Isotype)) %>% 
  select(Domain, GenomeID, comb) %>% 
  group_by(Domain, GenomeID, comb) %>% 
  distinct() %>% 
  ungroup() %>% 
  group_by(Domain, comb) %>% 
  summarise(genomes.per.isoaccp = n()) %>% 
  left_join(genomes.per.domain, by = 'Domain') %>% 
  mutate(anticodon.per.genomes = (genomes.per.isoaccp/ngenomes)*100) %>% 
  ungroup()

# count the total anticodon set per domain
anticodons.per.genome <- anticodon.genome.freq %>% 
  group_by(Domain) %>% 
  summarise(total.antic.per.domain=n())

# count anticodon freq per genome per domain
anticodon_freq <- gtrnadb.data %>% 
  select(Domain, Isotype, Anticodon) %>% 
  filter(!str_detect(Anticodon, "N|M|Y")) %>% 
  mutate(comb = paste0(Anticodon,'_',Isotype)) %>% 
  group_by(Domain, comb) %>% 
  summarise(anticodon.per.domain = n()) %>% 
  ungroup() %>% 
  group_by(Domain) %>% 
  distinct() %>% 
  mutate(total_trna_perDomain = sum(anticodon.per.domain)) %>% 
  ungroup() %>% 
  select(Domain, comb, anticodon.per.domain,total_trna_perDomain) %>% 
  left_join(genomes.per.domain, by = c('Domain')) %>%
  left_join(anticodons.per.genome, by = c('Domain'))

freq <- anticodon_freq %>% 
  mutate(media.trna.genoma = total_trna_perDomain/ngenomes) %>%
  mutate(fe = 1/total.antic.per.domain) %>%
  mutate(fo = anticodon.per.domain/total_trna_perDomain) %>%
  mutate(cont_obs = media.trna.genoma*fo) %>%
  mutate(cont_esp = media.trna.genoma*fe) %>%
  mutate(fc = cont_obs/cont_esp) %>%
  mutate(lg2fc = log2(fc))


mat_fro <- matrix(data=0, nrow = length(unique(freq$comb)), ncol = length(unique(freq$Domain)))
rownames(mat_fro)<-unique(freq$comb)
colnames(mat_fro)<-unique(freq$Domain)
mat_fro

# fill the mat_fro
for(i in 1:nrow(freq)){
  Domain <- freq$Domain[i]
  comb <- freq$comb[i]
  combs_count_index_i <- which(rownames(mat_fro)==comb)
  combs_count_index_j <- which(colnames(mat_fro)==Domain)
  mat_fro[combs_count_index_i, combs_count_index_j]<-freq$cont_obs[i]
}
mat_fro


# plot frequencies per domain as a heatmap
paletteLength <- 22
quantile_breaks <- function(xs, n) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

myColor <- colorRampPalette(c("black","#40916c","#52b788","#74c69d","#95d5b2","#b7e4c7","#d8f3dc"))(paletteLength)
myBreaks <- seq(0.00000001, max(mat_fro), length.out = 20)
myBreaks <- quantile_breaks(myBreaks, n = 20)
myBreaks <- append(myBreaks, c(0),0)

# set row order
roworder <- c("TGC_Ala","TCT_Arg","CCT_Arg","CCG_Arg","GTT_Asn","GTC_Asp","GCA_Cys","TTG_Gln","TTC_Glu","TCC_Gly","GCC_Gly","GTG_His","TAA_Leu","TAG_Leu","CAA_Leu","CAG_Leu","TTT_Lys","CAT_Met","GAA_Phe","TGG_Pro","TGA_Ser","GCT_Ser","TGT_Thr","CCA_Trp","GTA_Tyr","TAC_Val","CGC_Ala","TCG_Arg","CTG_Gln","CTC_Glu","CCC_Gly","CTT_Lys","CGG_Pro","CGA_Ser","CGT_Thr","CAC_Val","GGC_Ala","GAT_Ile","CAT_Ile","GAG_Leu","GGG_Pro","GGA_Ser","GGT_Thr","GAC_Val","AGC_Ala","AAT_Ile","TAT_Ile","AAG_Leu","AGG_Pro","AGA_Ser","AGT_Thr","AAC_Val","ACG_Arg","GCG_Arg","ACC_Gly","TCA_SeC","TCA_Sup","ATT_Asn","ATC_Asp","ACA_Cys","ATG_His","AAA_Phe","GCA_SeC","TCC_SeC","CAA_SeC","ACA_SeC","CTA_SeC","TAA_SeC","CCG_SeC","AAA_SeC","ACT_Ser","TTA_Sup","CTA_Sup","ATA_Tyr")
mat_fro <- mat_fro[roworder,]

heatmap <- pheatmap(mat_fro,
         color = myColor,
         breaks = myBreaks,
         legend_breaks = c(16,12,8,4,0),
         display_numbers = T,
         cellheight = 7,
         cellwidth = 20,
         number_format = "%.4f",
         number_color = "grey100",
         border_color = "grey60",
         fontsize_row = 6,
         fontsize_col = 8,
         fontsize_number = 4,
         treeheight_row = 40,
         treeheight_col = 20,
         cluster_rows = F,
         #clustering_distance_rows = "euclidean",
         #clustering_distance_cols = "euclidean",
         angle_col = 90,
         filename = paste0(figdir,'/heatmap_obs_freq.pdf'))

qual_col_pals <-  brewer.pal.info[brewer.pal.info$category == 'qual',]
anticodonpal <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
mof <- left_join(freq[,c("Domain","comb","fo")],anticodon.genome.freq,by=c('Domain','comb'))

# pervasiveness
pervasiveness <- mof %>%  
  select(Domain,comb,fo,anticodon.per.genomes)%>%
  ggplot(aes(x = anticodon.per.genomes, y = fo)) +
  geom_smooth(aes(y = fo), color = 'black', method=loess, se=T, fullrange=T, formula = y ~ x, size = 0.3)+
  geom_point(aes(colour = comb, shape = Domain),size = 2, stat = 'identity', alpha = 1)+
  scale_color_manual(values = anticodonpal)+
  theme_classic(base_line_size = 0.3) + 
  ylab("Mean tRNA gene frequency per genome per domain" ) + xlab("Pervasiveness per domain")+
  theme(axis.text = element_text(size = 8),
        axis.title.y = element_text(size = 10),
        legend.position = "right",legend.text=element_text(size=8),legend.key.height= unit(0.4, 'cm'),legend.key.width= unit(0.4, 'cm'))+
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = T,rr.digits = 2,size = 3,label.y = c(0.9,0.95))

# salve plot
ggplot2::ggsave(paste0(figdir,"/trna_pervasiveness_per_domain.pdf"),
                plot =  pervasiveness,device = "pdf",units = "cm",width = 20,height = 12)

vd.70 <- venn70(mof)
plot.vd.70 <- ggvenn(vd.70, 
                     fill_color = domainColor,
                     stroke_size = 0, set_name_size = 5, text_size = 4)

#salve plot
ggplot2::ggsave(paste0(figdir,"/trna_pervasiveness_per_domain.pdf"),plot =  pervasiveness,device = "pdf",units = "cm",width = 20,height = 12)
ggsave(paste0(figdir,"/gvenn70.pdf"),plot =  plot.vd.70,device = "pdf",units = "cm",width = 15,height = 10)

#save table with the isoacceptors shared by all domains
n <- 70
bac <- mof %>% filter(anticodon.per.genomes >= n) %>% filter(str_detect(Domain, 'Bacteria')) %>% mutate(anticodon.per.genomes = round(anticodon.per.genomes,digits = 2))
arc <- mof %>% filter(anticodon.per.genomes >= n) %>% filter(str_detect(Domain, 'Archaea')) %>% mutate(anticodon.per.genomes = round(anticodon.per.genomes,digits = 2))
euk <- mof %>% filter(anticodon.per.genomes >= n) %>% filter(str_detect(Domain, 'Eukaryota')) %>% mutate(anticodon.per.genomes = round(anticodon.per.genomes,digits = 2))

alpha <- inner_join(bac[,c('Domain','comb','anticodon.per.genomes')],arc[,c('Domain','comb','anticodon.per.genomes')],by='comb',suffix = c('.bac','.arc')) %>%
  inner_join(euk[,c('Domain','comb','anticodon.per.genomes')],by='comb') %>%
  dplyr::rename(isoacceptor =comb,
                Bacteria = anticodon.per.genomes.bac,
                Archaea = anticodon.per.genomes.arc,
                Eukaryota = anticodon.per.genomes) %>%
  select(isoacceptor,Bacteria,Archaea,Eukaryota)

alpha <- alpha %>% 
  mutate(isoacc = isoacceptor) %>% 
  separate(isoacc, into = c('anticodon','isotype'), sep = '_') %>% 
  group_by(isotype) %>% 
  mutate(n = n()) %>% ungroup()

alpha <- alpha[order(-alpha$n),]
alpha <- select(alpha,isoacceptor,Bacteria,Archaea,Eukaryota)

write_delim(alpha,"/home/feni/repository/trna_project/tables/most_pervasive_isoacceptors.tab",delim = '\t')
