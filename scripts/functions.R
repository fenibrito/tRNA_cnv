
cg.content <- function(seq){
  library(stringr)
  sequence <- seq
  num_g <- str_count(sequence, "G")
  num_c <- str_count(sequence, "C")
  gc_content <- (num_g + num_c) / str_length(sequence) * 100 
}

venn50 <- function(x){
  n <- 50
  bac <- x %>% filter(anticodon.per.genomes >= n) %>% filter(str_detect(Domain, 'Bacteria')) 
  arc <- x %>% filter(anticodon.per.genomes >= n) %>% filter(str_detect(Domain, 'Archaea'))
  euk <- x %>% filter(anticodon.per.genomes >= n) %>% filter(str_detect(Domain, 'Eukaryota'))
  
  list.venn.50 <- list(
    Bacteria = bac$comb,
    Archaea = arc$comb,
    Eukaryota = euk$comb)
  
}

venn70 <- function(x){
  n <- 70
  bac <- x %>% filter(anticodon.per.genomes >= n) %>% filter(str_detect(Domain, 'Bacteria')) 
  arc <- x %>% filter(anticodon.per.genomes >= n) %>% filter(str_detect(Domain, 'Archaea'))
  euk <- x %>% filter(anticodon.per.genomes >= n) %>% filter(str_detect(Domain, 'Eukaryota'))
  
  list.venn.70 <- list(
    Bacteria = bac$comb,
    Archaea = arc$comb,
    Eukaryota = euk$comb)
  
}


venn90 <- function(x){
  n <- 90
  bac <- x %>% filter(anticodon.per.genomes >= n) %>% filter(str_detect(Domain, 'Bacteria')) 
  arc <- x %>% filter(anticodon.per.genomes >= n) %>% filter(str_detect(Domain, 'Archaea'))
  euk <- x %>% filter(anticodon.per.genomes >= n) %>% filter(str_detect(Domain, 'Eukaryota'))
  
  list.venn.90 <- list(
    Bacteria = bac$comb,
    Archaea = arc$comb,
    Eukaryota = euk$comb)
  
}

