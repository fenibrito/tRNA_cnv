# scatter plot comparing the tRNA mean gene count per domain
#import library
library("tidyverse")
library("RColorBrewer")
library("reshape2")
install.packages("scatterplot3d")
library("scatterplot3d")

#save images directory 
figdir <- '/home/feni/repository/trna_project/figures'

#import trna data
gtrnadb.data <- read_csv('/home/feni/repository/trna_project/data/gtrnadb_data.csv')

#define palette color
domainColor <- c("#049F76","#118AB2","#BF805F")
domainColor2 <- c("#A2E1F6","#AFFDE8","#E8D1C5")
qual_col_pals <-  brewer.pal.info[brewer.pal.info$category == 'qual',]
anticodonpal <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#define domain order 
domain.order <- c('Bacteria', 'Archaea', 'Eukaryota')

# mean tRNA gene loci per isotype per genome
mean.trna.isotype.per.genome <- gtrnadb.data %>%
  filter(!str_detect(Anticodon, "N|M|Y"),
         !str_detect(Isotype, "Sup")) %>% 
  select(Domain,GenomeID,Isotype) %>%
  group_by(Domain,GenomeID,Isotype) %>%
  summarise(n.trna.per.isotype = n()) %>%
  ungroup() %>% group_by(Domain,Isotype) %>% 
  mutate(mean.trna.per.isotype = mean(n.trna.per.isotype)) %>% 
  select(Domain,mean.trna.per.isotype) %>% distinct %>% ungroup()

arc_mean_isot <- mean.trna.isotype.per.genome %>%
  filter(str_detect(Domain,'Archaea')) %>% 
  melt(id = 'Isotype','mean.trna.per.isotype', value.name = 'arc_mean_isotype_per_genome') %>%
  select(Isotype,arc_mean_isotype_per_genome)

euk_mean_isot <- mean.trna.isotype.per.genome %>%
  filter(str_detect(Domain,'Eukaryota')) %>% 
  melt(id = 'Isotype','mean.trna.per.isotype', value.name = 'euk_mean_isotype_per_genome') %>%
  select(Isotype,euk_mean_isotype_per_genome)

bac_mean_isot <- mean.trna.isotype.per.genome %>%
  filter(str_detect(Domain,'Bacteria')) %>% 
  melt(id = 'Isotype','mean.trna.per.isotype', value.name = 'bac_mean_isotype_per_genome') %>%
  select(Isotype,bac_mean_isotype_per_genome)

arc_euk_isot <- left_join(arc_mean_isot,euk_mean_isot,by='Isotype')
mean_isot_dom <- left_join(arc_euk_isot,bac_mean_isot,by='Isotype')

x <-  mean_isot_dom$euk_mean_isotype_per_genome
y <-  mean_isot_dom$arc_mean_isotype_per_genome
z <-  mean_isot_dom$bac_mean_isotype_per_genome

color <- anticodonpal[as.numeric(as.factor(mean_isot_dom$Isotype))]
s3d <- scatterplot3d(x = x, y = y, z = z,
              angle = 100,
              pch = 16,
              scale.y = 0.8,
              color = color,
              main="Mean tRNA gene per\nisotype per genome",
              zlab = "Bacteria",
              xlab = "Eukarya",
              ylab = "Archaea")
myLM <- lm(y ~ x + z,data=mean_isot_dom)
s3d$plane3d(myLM, lty.box = "solid",draw_polygon = T, draw_lines = T)

legend(s3d$xyz.convert(7.5, 3, 4.5), pch = 19, yjust=1,
       legend = levels(as.factor(mean_isot_dom$Isotype)), col = color)


merge(arc_mean_isot,euk_mean_isot,by = 'Isotype') %>%
  ggplot(aes(x = arc_mean_isotype_per_genome, y = euk_mean_isotype_per_genome, mapping = Isotype)) +
  geom_point(aes(colour = Isotype), shape = 20,size = 2, stat = 'identity', alpha = 1)+
  geom_text(aes(label = Isotype,color = Isotype),size = 5,hjust=0,vjust=0)+
  theme_classic()+
  theme(legend.position = 'none') + ylab('Eukarya mean tRNA gene per genome')+ xlab('Achaea mean tRNA gene per genome')

merge(arc_mean_isot,bac_mean_isot,by = 'Isotype') %>%
  ggplot(aes(x = arc_mean_isotype_per_genome, y = bac_mean_isotype_per_genome, mapping = Isotype)) +
  geom_point(aes(colour = Isotype), shape = 20,size = 2, stat = 'identity', alpha = 1)+
  geom_text(aes(label = Isotype,color = Isotype),size = 5)+
  scale_color_manual(values = anticodonpal)+
  theme_classic()+
  theme(legend.position = 'none') + ylab('Bacteria mean tRNA gene per genome')+ xlab('Achaea mean tRNA gene per genome')

merge(euk_mean_isot,bac_mean_isot,by = 'Isotype') %>%
  ggplot(aes(x = euk_mean_isotype_per_genome, y = bac_mean_isotype_per_genome, mapping = Isotype)) +
  geom_point(aes(colour = Isotype), shape = 20,size = 2, stat = 'identity', alpha = 1)+
  geom_text(aes(label = Isotype,color = Isotype),size = 5)+
  scale_color_manual(values = anticodonpal)+
  theme_classic()+
  theme(legend.position = 'none') + ylab('Bacteria mean tRNA gene per genome')+ xlab('Eukarya mean tRNA gene per genome')


# mean tRNA gene loci per anticodon per genome
mean.trna.anticodon.per.genome <- gtrnadb.data %>%
  filter(!str_detect(Anticodon, "N|M|Y")) %>% 
  select(Domain,GenomeID,Anticodon) %>%
  group_by(Domain,GenomeID,Anticodon) %>%
  summarise(n.trna.per.anticodon = n()) %>%
  ungroup() %>% 
  group_by(Domain,Anticodon) %>% 
  mutate(mean.trna.per.anticodon = mean(n.trna.per.anticodon)) %>% 
  select(Domain,mean.trna.per.anticodon) %>% 
  distinct() %>% 
  ungroup()

arc_mean_antc <- mean.trna.anticodon.per.genome %>%
  filter(str_detect(Domain,'Archaea')) %>% 
  melt(id = 'Anticodon','mean.trna.per.anticodon', value.name = 'arc_mean_anticodon_per_genome') %>%
  select(Anticodon,arc_mean_anticodon_per_genome)


euk_mean_antc <- mean.trna.anticodon.per.genome %>%
  filter(str_detect(Domain,'Eukaryota')) %>% 
  melt(id = 'Anticodon','mean.trna.per.anticodon', value.name = 'euk_mean_anticodon_per_genome') %>%
  select(Anticodon,euk_mean_anticodon_per_genome)

bac_mean_antc <- mean.trna.anticodon.per.genome %>%
  filter(str_detect(Domain,'Bacteria')) %>% 
  melt(id = 'Anticodon','mean.trna.per.anticodon', value.name = 'bac_mean_anticodon_per_genome') %>%
  select(Anticodon,bac_mean_anticodon_per_genome)

euk_arc_antc <- left_join(euk_mean_antc,arc_mean_antc,by = 'Anticodon') 
mean_antc_dom <- left_join(euk_arc_antc,bac_mean_antc,by = 'Anticodon') 

# mean_antc_dom %>% 
#   ggplot(aes(y = arc_mean_anticodon_per_genome, x = euk_mean_anticodon_per_genome, mapping = Anticodon)) +
#   geom_point(aes(colour = Anticodon), shape = 20,size = 2, stat = 'identity', alpha = 1)+
#   geom_text(aes(label = Anticodon,color = Anticodon),size = 5)+
#   scale_color_manual(values = anticodonpal)+
#   theme_classic()+
#   theme(legend.position = 'none') + xlab('Eukarya mean tRNA gene per genome')+ ylab('Achaea mean tRNA gene per genome')
# 
# mean_antc_dom %>% 
#   ggplot(aes(y = arc_mean_anticodon_per_genome, x = bac_mean_anticodon_per_genome, mapping = Anticodon)) +
#   geom_point(aes(colour = Anticodon), shape = 20,size = 2, stat = 'identity', alpha = 1)+
#   geom_text(aes(label = Anticodon,color = Anticodon),size = 5)+
#   scale_color_manual(values = anticodonpal)+
#   theme_classic()+
#   theme(legend.position = 'none') + xlab('Bacteria mean tRNA gene per genome')+ ylab('Achaea mean tRNA gene per genome')
# 
# mean_antc_dom %>% 
#   ggplot(aes(x = euk_mean_anticodon_per_genome, y = bac_mean_anticodon_per_genome, mapping = Anticodon)) +
#   geom_point(aes(colour = Anticodon), shape = 20,size = 2, stat = 'identity', alpha = 1)+
#   geom_text(aes(label = Anticodon,color = Anticodon),size = 5)+
#   scale_color_manual(values = anticodonpal)+
#   theme_classic()+
#   theme(legend.position = 'none') + ylab('Bacteria mean tRNA gene per genome')+ xlab('Eukarya mean tRNA gene per genome')

x2 = mean_antc_dom$euk_mean_anticodon_per_genome
y2 = mean_antc_dom$arc_mean_anticodon_per_genome
z2 = mean_antc_dom$bac_mean_anticodon_per_genome

antc_col <- anticodonpal[as.numeric(as.factor(mean_antc_dom$Anticodon))]
antic_s3d <- scatterplot3d(x = x2, y = y2, z = z2,
              angle = 100,
              pch = 16,
              color = antc_col,
              main="Mean tRNA gene per\nanticodon per genome",
              ylab = "Archaea",
              xlab = "Eukarya",
              zlab = "Bacteria")

myLM2 <- lm(y ~ x + z,data=mean_antc_dom)
antic_s3d$plane3d(myLM, lty.box = "solid",draw_polygon = T, draw_lines = T)

legend(antic_s3d$xyz.convert(25, 3, 4.5), pch = 19, yjust=1,
       legend = levels(as.factor(mean_antc_dom$Anticodon)), col = seq_along(levels(as.factor(mean_antc_dom$Anticodon))))
