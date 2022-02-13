# scatter plot comparing the tRNA mean gene count per isoacceptor and anticodon per domain

#import library
# install.packages("scatterplot3d")
library("tidyverse")
library("RColorBrewer")
library("reshape2")
library("scatterplot3d")

# save image directory 
figdir <- '/home/feni/repository/trna_project/figures'

# upload GtRNAdb data
gtrnadb.data <- read_csv('/home/feni/repository/trna_project/data/gtrnadb_data.csv')

#define palette color
domainColor <- c("#049F76","#118AB2","#BF805F")
domainColor2 <- c("#A2E1F6","#AFFDE8","#E8D1C5")
qual_col_pals <-  brewer.pal.info[brewer.pal.info$category == 'qual',]
anticodonpal <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# plot tRNA count per isoacceptor per genome ----------------------------------------------------

# mean tRNA gene loci per isoacceptor per genome
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


# color vector
color_ord <- factor(mean_isot_dom$Isotype) 
color <- anticodonpal[as.numeric(factor(mean_isot_dom$Isotype))]

# define axis
x <-  mean_isot_dom$euk_mean_isotype_per_genome
y <-  mean_isot_dom$arc_mean_isotype_per_genome
z <-  mean_isot_dom$bac_mean_isotype_per_genome

# plot
pdf(file="/home/feni/repository/trna_project/figures/3Dscatter_isotype_per_genome.pdf",width = 5, height = 4)
s3d <- scatterplot3d(x = x, y = y, z = z,
                     angle = 55,
                     pch = 16,
                     scale.y = 1,
                     color = color,
                     type = 'h',
                     main="Mean tRNA gene per isoacceptor per genome",
                     zlab = "Bacteria",
                     xlab = "Eukarya",
                     ylab = "Archaea")

# set legent to plot
legend(horiz = T, inset = 1, xpd = TRUE, 
       s3d$xyz.convert(-1, -3, 0), pch = 19, yjust=0, cex = 0.35,
       legend = levels(factor(mean_isot_dom$Isotype,color_ord)), 
       col = anticodonpal[as.numeric(factor(mean_isot_dom$Isotype))])
dev.off()
# myLM <- lm(y ~ x + z,data=mean_isot_dom)
# s3d$plane3d(myLM, lty.box = "solid",draw_polygon = F, draw_lines = T)



# plot tRNA count per anticodon per genome ----------------------------------------------------
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

#color vector
color_ant_ord <- factor(mean_antc_dom$Anticodon) 
antc_col <- anticodonpal[as.numeric(as.factor(mean_antc_dom$Anticodon))]

# define axis
x2 = mean_antc_dom$euk_mean_anticodon_per_genome
y2 = mean_antc_dom$arc_mean_anticodon_per_genome
z2 = mean_antc_dom$bac_mean_anticodon_per_genome

# plot
pdf(file="/home/feni/repository/trna_project/figures/3Dscatter_anticodon_per_genome.pdf",width = 5, height = 4)
antic_s3d <- scatterplot3d(x = x2, y = y2, z = z2,
                           angle = 55,
                           pch = 16,
                           color = antc_col,
                           main="Mean tRNA gene per anticodon per genome",
                           ylab = "Archaea",
                           xlab = "Eukarya",
                           zlab = "Bacteria")

# myLM2 <- lm(y ~ x + z,data=mean_antc_dom)
# antic_s3d$plane3d(myLM, lty.box = "solid",draw_polygon = F, draw_lines = T)

# set legend to plot
legend(horiz = T, inset = 1, xpd = TRUE, 
       s3d$xyz.convert(-1, -3, 0), pch = 19, yjust=0, cex = 0.35,
       legend = levels(factor(mean_antc_dom$Anticodon,color_ant_ord)), 
       col = anticodonpal[as.numeric(factor(mean_antc_dom$Anticodon))])
dev.off()
