
# import library
source("~/repository/trna_project/scripts/functions.R")
library("tidyverse")
library("scales")
library("ggpmisc")

# save images directory 
figdir <- '/home/feni/repository/trna_project/figures'

# define color palette
colorpolot <- c("#e27396","#d1ac00","#985277","#63CAE3","#c18c5d","#bb4430","#208b3a","#0f4c5c")

# upload GtRNAdb data
gtrnadb.data <- read_csv('/home/feni/repository/trna_project/data/gtrnadb_data.csv')

#import tRF data and reformat
trf1 <- readr::read_csv('/home/feni/repository/trna_project/data/tRFs/trfdb_trf-1.csv')
trf3 <- readr::read_csv('/home/feni/repository/trna_project/data/tRFs/trfdb_trf-3.csv')
trf5 <- readr::read_csv('/home/feni/repository/trna_project/data/tRFs/trfdb_trf-5.csv')

trfs <- rbind(trf1,trf3,trf5)
rm(trf1,trf3,trf5)


trf.db <- trfs %>% 
  tidyr::separate(`tRF Map Positions`, into = c("s", "e"), sep = " - ") %>% 
  tidyr::separate(s, into = c("stext", "start"), sep = ":") %>%  
  tidyr::separate(e, into = c("etext", "end"), sep = ":") %>% 
  tidyr::separate(`tRNA Name`, into = c("chr", "trna_type"), sep = "-") %>% 
  mutate(Isotype = substr(trna_type,1,3),Anticodon = substr(trna_type,4,6)) %>% 
  mutate(trf.length = str_length(`tRF Sequence`)) %>% 
  select(`tRF ID`,Type,`tRNA Gene Co-ordinates`,`tRNA Sequence`,`Organism`,chr,trna_type,`tRF Sequence`,Isotype,Anticodon,start,end,trf.length) %>% 
  dplyr::rename(trf_id = `tRF ID`,
         type = Type,
         trna_gene_coordinates = `tRNA Gene Co-ordinates`,
         trna_seq = `tRNA Sequence`,
         organism = `Organism`,
         trf_seq = `tRF Sequence`)

trf.db$myspecie <- ""
trf.db$myspecie[trf.db$organism == "human"] <- "Homo sapiens"
trf.db$myspecie[trf.db$organism == "mouse"] <- "Mus musculus"
trf.db$myspecie[trf.db$organism == "drosophila"] <- "Drosophila melanogaster"
trf.db$myspecie[trf.db$organism == "c.elegans"] <- "Caenorhabditis elegans"
trf.db$myspecie[trf.db$organism == "s.pombe"] <- "Schizosaccharomyces pombe"
trf.db$myspecie[trf.db$organism == "r.sphaeroides"] <- "Rhodobacter sphaeroides"
trf.db$myspecie[trf.db$organism == "Xenopus-tropicalis"] <- "Xenopus tropicalis"
trf.db$myspecie[trf.db$organism == "Zebra_fish_Zv9"] <- "Danio rerio"


# count number of tRFs per specie
trf_per_sp <- trf.db %>%
  group_by(myspecie,trf_seq) %>% 
  summarise(trfs.per.seq = n()) %>% 
  distinct() %>% 
  select(myspecie,trfs.per.seq) %>% 
  ungroup() %>% 
  group_by(myspecie) %>% 
  mutate(trf.seq.per.sp = n(),
         trf.per.sp = sum(trfs.per.seq)) %>% 
  select(myspecie,trf.seq.per.sp,trf.per.sp) %>% 
  distinct() %>% 
  ungroup()

# count number of tRNA per specie
sp.trna.data <- gtrnadb.data %>% 
  select(Genome,GenomeID) %>%
  filter(str_detect(Genome,'^Homo sapiens|^Mus musculus|^Drosophila melanogaster|Caenorhabditis elegans|Schizosaccharomyces pombe|Rhodobacter sphaeroides|Xenopus tropicalis|Danio rerio')) %>% 
  mutate(myspecie = str_extract(Genome, '\\w* \\w*')) %>%
  group_by(myspecie) %>%
  distinct() %>% 
  mutate(nGenomes = n())%>%
  ungroup() %>% 
  left_join(select(gtrnadb.data, Domain,GenomeID, Isotype, Anticodon,len_mature),by = 'GenomeID') 

trna.count.sp <- sp.trna.data %>%
  filter(!str_detect(Anticodon, "N|M|Y"))%>% 
  select(Domain,myspecie,GenomeID,nGenomes) %>% 
  group_by(GenomeID) %>% 
  mutate(ntrna_gen = n())  %>% 
  distinct()%>% 
  ungroup() %>% 
  group_by(myspecie) %>% 
  mutate(mean_trna = mean(ntrna_gen), 
         sdtrna = sd(ntrna_gen), 
         erro = sd(ntrna_gen)/sqrt(length(ntrna_gen))) %>%
  ungroup() %>% 
  select(Domain,myspecie,GenomeID,mean_trna,sdtrna,erro)


# tRF count per specie vs tRNA count
trna.trf.per.sp <- trna.count.sp %>% 
  select(Domain,myspecie,mean_trna) %>% 
  distinct() %>% 
  left_join(trf_per_sp,by = c("myspecie"))


# plot mean tRNA vs number of tRF per specie
trf_per_trna <- trna.trf.per.sp %>% 
  ggplot(aes(x = mean_trna, y = trf.per.sp))+
  geom_point(aes(y = trf.per.sp, color = myspecie), size = 2, stat = 'identity')+
  geom_smooth(aes(y = trf.per.sp), color = 'black', method=lm, se=T,fullrange=F, formula = y ~ x, size = 0.5)+
  scale_x_log10(name = "Log10(Mean tRNA gene count per species)", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(name = "Log10(tRF count)", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = colorpolot)+
  theme_classic(base_line_size = 0.3) +
  theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),
        legend.title = element_blank(),legend.text=element_text(size=10,face = 'italic'),legend.background=element_blank())+ #,legend.position = c(0.75,0.25)
  # ylab('Number of tRNA isoacceptor per genome')+
  stat_poly_eq(formula = y ~ x,
               eq.with.lhs = "italic(y)~`=`~",
               aes(label = paste(..rr.label.., sep = "~~")),
               parse = T,rr.digits = 3,size = 4,label.x = 0.05,label.y = 0.95)

# plot tRF count vs number of tRF sequence
seq_per_trf <- trna.trf.per.sp %>% 
  ggplot(aes(x = trf.per.sp, y = trf.seq.per.sp))+
  geom_point(aes(y = trf.seq.per.sp, color = myspecie), size = 2, stat = 'identity')+
  geom_smooth(aes(y = trf.seq.per.sp), color = 'black', method=lm, se=T,fullrange=F, formula = y ~ x, size = 0.5)+
  scale_x_log10(name = "Log10(tRF count)", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(name = "Log10(tRF non-redundant sequence count)", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_color_manual(values = colorpolot)+
  theme_classic(base_line_size = 0.3) +
  theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),
        legend.title = element_blank(),legend.text=element_text(size=10,face = 'italic'),legend.background=element_blank())+ #,legend.position = c(0.75,0.25)
  # ylab('Number of tRNA isoacceptor per genome')+
  stat_poly_eq(formula = y ~ x,
               eq.with.lhs = "italic(y)~`=`~",
               aes(label = paste(..rr.label.., sep = "~~")),
               parse = T,rr.digits = 3,size = 4,label.x = 0.05,label.y = 0.95)

# tRF length and tRNA length
trf_len <- trf.db %>% 
  select(myspecie,trf.length) %>% 
  left_join(sp.trna.data[,c('myspecie','len_mature')],by = 'myspecie') %>% 
  group_by(myspecie) %>% 
  mutate(trf.mean.len = mean(trf.length, na.rm = T),
         trna.mean.len = mean(len_mature, na.rm = T)) %>% 
  select(myspecie,trf.mean.len,trna.mean.len) %>% 
  distinct() %>% 
  ungroup()


# plot tRF length vs tRF CG content
frag.cg.cont <- trf.db %>%
  group_by(myspecie,trf_seq) %>% 
  summarise(trfs.per.seq = n()) %>% 
  ungroup() %>% 
  select(myspecie,trf_seq) %>% 
  mutate(cg_count = cg.content(trf_seq)) %>% 
  group_by(myspecie) %>% 
  summarise(trf.mean.cg = mean(cg_count)) %>% 
  left_join(trf_len[,c('myspecie','trf.mean.len')],by = 'myspecie') %>% 
  ggplot(aes(x = trf.mean.len, y = trf.mean.cg))+
  geom_point(aes(y = trf.mean.cg, color = myspecie), size = 2, stat = 'identity')+
  geom_smooth(aes(y = trf.mean.cg), color = 'black', method=lm, se=T,fullrange=F, formula = y ~ x, size = 0.5)+
  scale_color_manual(values = colorpolot)+
  scale_x_continuous(limits = c(18,24),breaks = c(18,20,22,24))+
  theme_classic(base_line_size = 0.3) +
  theme(axis.text = element_text(size = 10),axis.title = element_text(size = 8),
        legend.title = element_blank(),legend.text=element_text(size=10, face = 'italic'),legend.background=element_blank())+ #,legend.position = c(0.75,0.25)
  ylab('Mean tRF CG content')+xlab('Mean tRF length')+
  stat_poly_eq(formula = y ~ x,
               eq.with.lhs = "italic(y)~`=`~",
               aes(label = paste(..rr.label.., sep = "~~")),
               parse = T,rr.digits = 3,size = 4,label.x = 0.95,label.y = 0.95)

# non-redundant tRF sequence
non_redundant_trfs_sp <- trna.trf.per.sp %>% 
  select(myspecie,trf.seq.per.sp) %>% 
  ggplot(aes(x = trf.seq.per.sp, y = reorder(myspecie, -trf.seq.per.sp), fill = myspecie)) +
  geom_bar(stat = 'identity')+
  geom_text(aes(x = (trf.seq.per.sp + 6), label = trf.seq.per.sp), size = 4, stat = 'unique')+
  scale_fill_manual(values = colorpolot)+
  xlab('Number of non-redundant tRF')+ylab('')+
  theme_minimal()+
  theme(legend.title = element_blank(), legend.position = 'none',
        axis.text.y = element_text(face = 'italic'),
        axis.title.x = element_text(size = 8))

# tRF uniqueness
uniqueness <- trna.trf.per.sp %>% 
  mutate(uniq = trf.seq.per.sp/trf.per.sp) %>% 
  ggplot(aes(x = uniq, y = reorder(myspecie, -uniq), fill = myspecie)) +
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = colorpolot)+
  xlab('Non-redundant tRF frequency')+ylab('')+
  theme_minimal()+
  theme(legend.title = element_blank(), legend.position = 'none',
        axis.text.y = element_text(face = 'italic'),
        axis.title.x = element_text(size = 8))

# boxplot tRF sequence length per specie
sp_lab <- str_sort(unique(trf.db$myspecie))

trf.seq.len.box <- trf.db %>%
  group_by(myspecie,trf_seq) %>% 
  summarise(trfs.per.seq = n()) %>% 
  ungroup() %>% 
  select(myspecie,trf_seq) %>% 
  mutate(frag_len = str_length(trf_seq)) %>%
  group_by(myspecie) %>% 
  mutate(medi.trna.len = median(frag_len)) %>% 
  ungroup() %>% 
  ggplot(aes(y = myspecie, x = frag_len))+
  geom_boxplot(aes(y = reorder(myspecie,-medi.trna.len), fill = myspecie),
                   outlier.size = 0.4,outlier.shape = 23,outlier.colour = "black",outlier.fill = "black",outlier.stroke = 1.5,outlier.alpha = .2,
               notch = F,color='grey10',width=0.6,lwd = 0.25)+
  scale_fill_manual(breaks = sp_lab, labels = sp_lab, values = colorpolot)+
  scale_x_continuous(limits = c(13,38),breaks = c(14,18,22,26,30,34))+
  theme_minimal()+
  theme(axis.title.x = element_text(size = 8),legend.text=element_text(face = 'italic'))+
  xlab('tRF length') + ylab('')


#tRF length vs tRNA length
trna.trf.len <- trf_len %>%
  ggplot(aes(x = trf.mean.len, y = trna.mean.len))+
  geom_point(aes(y = trna.mean.len, color = myspecie), size = 2, stat = 'identity')+
  geom_smooth(aes(y = trna.mean.len), color = 'black', method=lm, se=T,fullrange=F, formula = y ~ x, size = 0.5)+
  scale_color_manual(values = colorpolot)+
  scale_x_continuous(limits = c(18,24),breaks = c(18,20,22,24))+
  theme_classic(base_line_size = 0.3) +
  theme(axis.text = element_text(size = 10),axis.title = element_text(size = 8),
        legend.title = element_blank(),legend.text=element_text(size=10,face = 'italic'),legend.background=element_blank())+ #,legend.position = c(0.75,0.25)
  ylab('Mean tRNA length')+xlab('Mean tRF length')+
  stat_poly_eq(formula = y ~ x,
               eq.with.lhs = "italic(y)~`=`~",
               aes(label = paste(..rr.label.., sep = "~~")),
               parse = T,rr.digits = 3,size = 4,label.x = 0.95,label.y = 0.95)

# tRF pe type
trf.per.type <- trf.db %>%
  select(myspecie,type,trf_seq) %>%
  group_by(myspecie,type,trf_seq) %>%
  mutate(trfs.per.seq = n()) %>%
  distinct() %>%
  ungroup() %>%
  group_by(myspecie,type) %>%
  mutate(trf.seq = n())%>%
  ungroup() %>%
  select(myspecie,type,trf.seq) %>%
  distinct() %>%
  ungroup()

# trF length per type
trf_len_type <- trf.db %>%
  group_by(myspecie,type,trf_seq) %>% 
  summarise(trfs.per.seq = n()) %>% 
  ungroup() %>% 
  select(myspecie,type,trf_seq) %>% 
  mutate(frag_len = str_length(trf_seq)) %>%
  group_by(myspecie,type) %>% 
  mutate(medi.trf.len = median(frag_len),
         meanLen = mean(frag_len)) %>% 
  ungroup()

# tRF length per type
trf_len_per_type <- trf_len_type %>% 
  ggplot(aes(x = type, y = frag_len))+
  geom_jitter(aes(color = type),width = .15)+
  geom_boxplot(aes(x = type, color = type),fill = 'NA',
                 outlier.size = 0.4,outlier.shape = 23,outlier.colour = "black",outlier.fill = "black",outlier.stroke = 1.5,outlier.alpha = .2,
               notch = F,color='grey10',width=0.4,lwd = 0.25)+
  scale_color_manual(values = c('#2ec4b6','#e71d36','#ff9f1c'))+
  scale_y_continuous(limits = c(13,38),breaks = c(14,18,22,26,30,34))+
  theme(legend.title = element_blank())+
  theme_classic()+
  xlab('') + ylab('tRF length')

# tRF length per species
tfr_len_per_sp <- trf.db %>%
  group_by(myspecie,type,trf_seq) %>% 
  summarise(trfs.per.seq = n()) %>% 
  ungroup() %>% 
  select(myspecie,type,trf_seq) %>% 
  mutate(frag_len = str_length(trf_seq)) %>%
  group_by(myspecie,type) %>% 
  mutate(medi.trna.len = median(frag_len)) %>% 
  ungroup()

tfr_len_per_sp_box <- tfr_len_per_sp %>%  
  ggplot(aes(y = myspecie, x = frag_len))+
  geom_jitter(aes(y = reorder(myspecie,-medi.trna.len),color = type),height = .15,size = 1)+
  scale_color_manual(values = c('#2ec4b6','#e71d36','#ff9f1c'))+
  geom_boxplot(aes(y = reorder(myspecie,-medi.trna.len), color = myspecie),fill = 'NA',
               outlier.size = 0.4,outlier.shape = 23,outlier.colour = "black",outlier.fill = "black",outlier.stroke = 1.5,outlier.alpha = .2,
               notch = F,color='grey10',width=0.6,lwd = 0.25)+
  scale_x_continuous(limits = c(13,38),breaks = c(14,18,22,26,30,34))+
  theme_classic()+
  theme(legend.title = element_blank(),
        axis.text.y = element_text(face = 'italic'),
        axis.title.x = element_text(size = 10))+
  xlab('tRF length') + ylab('')

unique_trf_type_seq <- trf.per.type  %>%
  ggplot(aes(y = factor(type, levels = c('trf-5','trf-3','trf-1'))))+
  geom_bar(aes(x = trf.seq,fill = myspecie), width = 0.9,size=0.05,stat="identity", position="dodge")+
  scale_fill_manual(values = colorpolot)+
  facet_grid(rows = vars(factor(myspecie,levels = c('Homo sapiens','Mus musculus','Xenopus tropicalis','Danio rerio','Drosophila melanogaster','Caenorhabditis elegans','Schizosaccharomyces pombe','Rhodobacter sphaeroides'))))+
  scale_x_continuous(limits = c(0,max(trf.per.type$trf.seq)),breaks = seq(0, 80, by=20))+
  theme_classic() + theme(legend.position = 'none')+
  xlab('Non-redundant tRFs per type') + ylab('')

# save plot
ggsave(paste0(figdir,"/trf_per_trna.pdf"),plot =  trf_per_trna, device = "pdf",units = "cm",width = 14,height = 8)
ggsave(paste0(figdir,"/seq_per_trf.pdf"),plot =  seq_per_trf, device = "pdf",units = "cm",width = 14,height = 8)
ggsave(paste0(figdir,"/trf_seq_div.pdf"),plot =  trf_seq_div, device = "pdf",units = "cm",width = 12,height = 6)
ggsave(paste0(figdir,"/uniqueness.pdf"),plot =  uniqueness, device = "pdf",units = "cm",width = 10,height = 6)
ggsave(paste0(figdir,"/trf_seq_len_box.pdf"),plot =  trf.seq.len.box, device = "pdf",units = "cm",width = 10,height = 6)
ggsave(paste0(figdir,"/trna_trf_len.pdf"),plot =  trna.trf.len, device = "pdf",units = "cm",width = 14,height = 8)
ggsave(paste0(figdir,"/frag_cg_cont.pdf"),plot =  frag.cg.cont, device = "pdf",units = "cm",width = 14,height = 8)
ggsave(paste0(figdir,"/non_redundant_trfs_sp.pdf"),plot =  non_redundant_trfs_sp, device = "pdf",units = "cm",width = 12,height = 6)
ggsave(paste0(figdir,"/trf_len_type.pdf"),plot =  trf_len_per_type, device = "pdf",units = "cm",width = 8,height = 8)
ggsave(paste0(figdir,"/tfr_len_per_sp_box.pdf"),plot =  tfr_len_per_sp_box, device = "pdf",units = "cm",width = 14,height = 6)
ggsave(paste0(figdir,"/unique_trf_type_seq.pdf"),plot =  unique_trf_type_seq, device = "pdf",units = "cm",width = 8,height = 8)
