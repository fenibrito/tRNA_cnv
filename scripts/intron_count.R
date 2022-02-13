#plot fraction of intron count per tRNA gene per domain

#import library
library("tidyverse")

figdir <- '/home/feni/repository/trna_project/figures'

# upload GtRNAdb data
gtrnadb.data <- read_csv('/home/feni/repository/trna_project/data/gtrnadb_data.csv')#define palette color

# define domain colors
domainColor <- c("#049F76","#118AB2","#BF805F")

# define domain order 
domain.order <- c('Bacteria', 'Archaea', 'Eukaryota')

# intron count per domain
intron.per.domain <- gtrnadb.data %>%
  select(Domain, IntronCount) %>% 
  group_by(Domain, IntronCount) %>% 
  count() %>%
  ungroup() %>% 
  group_by(Domain) %>% 
  mutate(fraction = prop.table(n))
print(intron.per.domain)

# plot
intron.count <- intron.per.domain %>% 
  ggplot(aes(x = factor(Domain, rev(domain.order)), y = IntronCount)) +
  geom_point(mapping = aes(size = fraction, color = Domain))+
  scale_y_continuous(breaks = c(0,1,2,3), limits=c(0, 3))+
  scale_x_discrete(breaks = c('Eukaryota','Archaea','Bacteria'), labels = c('EUK', 'ARC', 'BAC'))+
  scale_size(name = '',breaks = c(0.01, 0.25, 0.5, 0.75, 0.9), labels = c('1%','25%','50%','75%','100%'))+
  scale_color_manual(values = domainColor)+
  theme_minimal(base_line_size = 0.2) + ylab("Intron Count")+ xlab("") +
  theme(axis.text.x = element_text(angle = 0, size = 10, vjust = 0.5, hjust = .5),
        axis.text.y = element_text(angle = 0, size = 7, vjust = 0.5, hjust = .5),
        legend.position = "right")

# save plot
ggsave(paste0(figdir,"/intron_count.pdf"),plot = intron.count,device = "pdf",units = "cm",width = 7,height = 8)  
  
