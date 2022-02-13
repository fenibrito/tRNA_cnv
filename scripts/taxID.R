
# add taxonomic information to the GtRNAdb data
library("tidyverse")
library("Biostrings")

# import GtRNAdb data
gtrnadb.data <- read_csv('/home/feni/repository/trna_project/data/gtrnadb_data.csv')

# select list of Genomes
genome.list <- as.data.frame(gtrnadb.data %>% 
                                       select(Genome) %>% 
                                       distinct())

# save genome list
write_csv(genome.list,'/home/feni/repository/trna_project/data/list_of_genomes.csv', col_names = F)


# upload the genome list with the taxID retrieved from NCBI Taxonomy https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi
genomes.taxID <- read_delim("/home/feni/repository/trna_project/data/list_of_genomesTAXID.txt",
                       "\t", escape_double = T, trim_ws = TRUE, col_select = c(name,taxid))

# select the taxID of the species not founded by the genome name
# For those cases we searched for the species taxID and, in some cases, the information was manually curated
cadidatus.sp <- genomes.taxID %>%
  filter(is.na(taxid)) %>%
  select(name) %>% 
  distinct() %>% 
  filter(str_detect(name,"^Candidatus")) %>%
  mutate(myspecie = str_extract(name, '\\w* \\w* \\w*')) %>% 
  select(myspecie) %>% 
  distinct() 

list.no.taxid <- genomes.taxID %>%
  filter(is.na(taxid)) %>%
  filter(!str_detect(name,"^Candidatus")) %>%
  select(name) %>% 
  distinct() %>% 
  mutate(myspecie = str_extract(name, '\\w* \\w*')) %>% 
  select(myspecie) %>% 
  distinct() %>% 
  rbind(cadidatus.sp) %>% 
  mutate(species = ifelse(test = str_detect(myspecie, 'sp$'), yes = paste0(myspecie,'.'), no = myspecie)) %>% 
  select(species)


# save list of species without taxid 
write_csv(list.no.taxid,'/home/feni/repository/trna_project/data/list_of_species_no_taxid.csv', col_names = F)

# upload list of species with its taxid
sp.taxID <- read_delim("/home/feni/repository/trna_project/data/list_of_species_taxid.txt",
           "\t", escape_double = T, trim_ws = TRUE, col_select = c(name,taxid))

genomes.taxID <- sp.taxID %>% right_join(genomes.taxID,by = 'name') %>%
  mutate(taxid = ifelse(test = is.na(taxid.x) & !is.na(taxid.y), yes = taxid <-  taxid.y, no = taxid <-  taxid.x)) %>%
  select(name, taxid)

catidatus.taxid <- genomes.taxID %>%
  filter(str_detect(name,"^Candidatus")) %>% 
  dplyr::rename(species = name) %>% 
  mutate(myspecie = str_extract(species, '\\w* \\w* \\w*')) %>% 
  mutate(name = ifelse(test = str_detect(myspecie, 'sp$'), yes = paste0(myspecie,'.'), no = myspecie)) %>% 
  left_join(sp.taxID, by = 'name') %>% 
  mutate(taxid = ifelse(test = is.na(taxid.x), yes = taxid <-  taxid.y, no = taxid <-  taxid.x)) %>% 
  dplyr::rename(Genome=species)

taxid.list <- genomes.taxID %>%
  filter(!str_detect(name,"^Candidatus")) %>% 
  dplyr::rename(species = name) %>% 
  mutate(myspecie = str_extract(species, '\\w* \\w*')) %>% 
  mutate(name = ifelse(test = str_detect(myspecie, 'sp$'), yes = paste0(myspecie,'.'), no = myspecie)) %>%  
  left_join(sp.taxID, by = 'name') %>% 
  mutate(taxid = ifelse(test = is.na(taxid.x), yes = taxid <-  taxid.y, no = taxid <-  taxid.x)) %>% 
  dplyr::rename(Genome=species) %>% 
  rbind(catidatus.taxid) %>%
  select(Genome,taxid)

unique.taxid.list <- taxid.list %>% 
  select(taxid) %>% 
  distinct()

# export data to recover the taxonomy lineage
write_csv(unique.taxid.list,'/home/feni/repository/trna_project/data/genomes_taxid_list.csv')

# upload taxonomic data retrieved from http://bioinfo.icb.ufmg.br/taxallnomy/
taxallnomy <- read_delim("/home/feni/repository/trna_project/data/taxallnomy_lineage.txt","\t", escape_double = T, trim_ws = TRUE)

gtrnadb.data <- left_join(gtrnadb.data,taxid.list,by='Genome')
gtrnadb.taxonomy <- left_join(gtrnadb.data,taxallnomy,by='taxid')
write_csv(gtrnadb.taxonomy,'/home/feni/repository/trna_project/data/gtrnadb_taxonomy_data.csv')
