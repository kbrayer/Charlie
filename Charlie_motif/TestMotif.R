library(MotifDb)
library(seqLogo)
library(motifStack)
library(Biostrings)
library(GenomicFeatures)
library(org.Sc.sgd.db)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)

query(MotifDb, "DAL80")   
pfm.dal80.jaspar <- query(MotifDb,"DAL80")[[1]]
seqLogo(pfm.dal80.jaspar)
dal1 <- "YIR027C"
chromosomal.loc <- 
  transcriptsBy(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, by="gene") [dal1]
promoter.dal1 <- 
  getPromoterSeq(chromosomal.loc, Scerevisiae, upstream=1000, downstream=0)
pcm.dal80.jaspar <- round(100 * pfm.dal80.jaspar)
matchPWM(pcm.dal80.jaspar, unlist(promoter.dal1)[[1]], "90%")

#2nd
query(MotifDb,"DAL80")
dal80.jaspar <- query(MotifDb,"DAL80")[[1]]
dal80.scertf <-query(MotifDb,"DAL80")[[2]]
seqLogo(dal80.jaspar)
seqLogo(dal80.scertf)

pfm.dal80.jaspar <- new("pfm", mat=query(MotifDb, "dal80")[[1]], 
                        name="DAL80-JASPAR")
pfm.dal80.scertf <- new("pfm", mat=query(MotifDb, "dal80")[[2]], 
                        name="DAL80-ScerTF")

#next line doesn't work right
plotMotifLogoStack(DNAmotifAlignment(c(pfm.dal80.scertf, pfm.dal80.jaspar)))

#3rd
query(MotifDb, "gat1")
pfm.gat1.jaspar = new("pfm", mat=query(MotifDb, "gat1")[[1]], 
                      name="GAT1-JASPAR")
pfm.gat1.scertf = new("pfm", mat=query(MotifDb, "gat1")[[2]], 
                      name="GAT1-ScerTF")
pfm.gat1.uniprobe = new("pfm", mat=query(MotifDb, "gat1")[[3]], 
                        name="GAT1-UniPROBE")

#next line doesn't work right
plotMotifLogoStack(c(pfm.gat1.uniprobe, pfm.gat1.scertf, pfm.gat1.jaspar))

#4th

pfm.dal80.scertf <- query(MotifDb, "dal80")[[2]]
pcm.dal80.scertf <- round(100 * pfm.dal80.scertf)

pfm.gat1.jaspar <- query(MotifDb, "gat1")[[1]]
pcm.gat1.jaspar <- round(100 * pfm.gat1.jaspar)

pfm.gat1.scertf <- query(MotifDb, "gat1")[[2]]
pcm.gat1.scertf <- round(100 * pfm.gat1.scertf)

genes <- c("DAL1", "DAL2", "DAL4", "DAL5", "DAL7", "DAL80", "GAP1")
orfs <- as.character(mget(genes, org.Sc.sgdCOMMON2ORF))

grlY <- transcriptsBy(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, by="gene") [orfs]
Ypromoter.seqs <- getPromoterSeq(grlY, Scerevisiae, upstream=1000,
                                downstream=0)

print (class(Ypromoter.seqs))
promoter.seqs <- unlist(Ypromoter.seqs)
print (class(Ypromoter.seqs))
#next line doesn't work
matchPWM(pcm.dal80.scertf, Ypromoter.seqs[[1]], "90%")

Ypwm.hits <- sapply(Ypromoter.seqs[[1]], 
                   function(pseq) 
                     matchPWM(pcm.dal80.scertf, pseq, min.score="90%"))

dal80.scertf.hits <- sapply(promoter.seqs, function(pseq) 
  matchPWM(pcm.dal80.scertf, pseq, min.score="90%"))
gat1.scertf.hits  <- sapply(promoter.seqs, function(pseq) 
  matchPWM(pcm.gat1.scertf, pseq, min.score="90%"))
gat1.jaspar.hits  <- sapply(promoter.seqs, function(pseq) 
  matchPWM(pcm.gat1.jaspar, pseq, min.score="90%"))

dal80.scertf <- sapply(dal80.scertf.hits, length)
gat1.jaspar  <- sapply(gat1.jaspar.hits,  length)
gat1.scertf  <- sapply(gat1.scertf.hits,  length)

tbl.gata     <- data.frame(gene=genes, dal80.scertf, gat1.jaspar, gat1.scertf)

