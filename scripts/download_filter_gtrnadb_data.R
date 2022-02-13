# upload and modifying GtRNAdb data
library("tidyverse")
library("Biostrings")

# upload GtRNAdb data
gtrnadb.raw.data <- readr::read_delim("/home/feni/GtRNAdb_data/gtrnadb-search7469-17dez19.out",
                  "\t", escape_double = FALSE, trim_ws = TRUE)

# table edit block
# substitute "scaffold:" for "scaffold_"
gtrnadb.raw.data$Locus <- str_replace(gtrnadb.raw.data$Locus, pattern = "scaffold:", replacement = "scaffold_")

# changing columns to obtain the 'start' and 'end' of tRNA locus and the genomic length
gtrnadb <- gtrnadb.raw.data %>%
  separate(Locus, into = c("chr", "position"), sep = ":") %>%
  separate(position, into = c("start", "end"), sep = "-") %>%
  separate(tRNAscanID, into = c("tRNAscanID","tRNAscanIDcode"), sep = "-",extra = "merge") %>%
  mutate(IDs = paste0(GenomeID, "_", GtRNAdbID, "#", tRNAscanID)) %>% 
  mutate(length = (as.numeric(end) - as.numeric(start)) + 1)

# reassign isotype name
gtrnadb.data <- gtrnadb
gtrnadb.data$Isotype <- as.character(gtrnadb.data$Isotype)
gtrnadb.data$Isotype[gtrnadb.data$Isotype == "iMet"] <- "Met"
gtrnadb.data$Isotype[gtrnadb.data$Isotype == "fMet"] <- "Met"
gtrnadb.data$Isotype[gtrnadb.data$Isotype == "Ile2"] <- "Ile"
gtrnadb.data$Isotype <- as.factor(gtrnadb.data$Isotype)

rm(gtrnadb,gtrnadb.raw.data)

# remove genomes with less than 20 different isoacceptors
select.genomes <- gtrnadb.data %>%
  filter(!str_detect(Anticodon, "N|M|Y")) %>% 
  select(Domain, Genome, GenomeID, Isotype, Anticodon)%>% 
  group_by(Domain, GenomeID) %>%
  distinct() %>%  
  mutate(n.isoacceptor = n()) %>% 
  filter(n.isoacceptor >= 20) %>%
  select(Domain,Genome,GenomeID,n.isoacceptor) %>%
  distinct()


gtrnadb.data <- left_join(select.genomes[,c("Domain","Genome","GenomeID")],gtrnadb.data,by = c("Domain","Genome","GenomeID"))

# upload tRNA mature sequence file
maduros <- readRNAStringSet("/home/feni/tRNA/gtrnadb.ucsc.edu/genomes/trnaMature.fa")
IDseq = names(maduros)
seq = (maduros)
dfasta <- data.frame(IDseq,seq, row.names = NULL)

# dataframe with the mature tRNA length
trna.mature.length <- dfasta %>%
  separate(IDseq, into = c('id','waste'), sep = '\\) ', extra = 'drop') %>% 
  mutate(len_mature = str_length(seq)) %>% 
  mutate(IDs = id)%>% 
  select(IDs,len_mature) 

gtrnadb.data <- left_join(gtrnadb.data,trna.mature.length,by = "IDs")
write_csv(gtrnadb.data,'/home/feni/repository/trna_project/data/gtrnadb_data.csv')
