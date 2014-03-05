##Finding new motifs from DNA sequences

library(MotIV)
library(MotifDb)
library(seqLogo)
library(motifStack)
library(Biostrings)
library(GenomicFeatures)
library(hgu133plus2.db)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(rGADEM)
library(ChIPpeakAnno)



##load file of sequences to be analyzed
FastaFile <- paste("SigGenes_Proms_Fasta.fasta", ".fasta",)
#path <- ~/Volumes/LaCie/ness-code/Charlie_motif/"SigGenes_Proms_Fasta.fasta"
Sequences <- readDNAStringSet(FastaFile, "fasta")
gadam<-GADEM(FastaFile, verbose=1,genome=Hsapiens, pValue=0.02)

##First trying with "seeded" analysis

query(MotifDb, "MYB")
pfm.MYB.jaspar <- query (query (MotifDb, "mmusculus"),"MYB") [[3]]
pcm.MYB.jaspar <- round(100 * pfm.MYB.jaspar)
