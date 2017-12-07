library("GenomicRanges")
library("Biostrings")
library("BSgenome.Hsapiens.UCSC.hg19")

#' chromInfo
#'
#' two-columns dataframe with chromosomes length
#'
#' @export chromInfo
chromInfo <- function(version = "hg19") {
  if (version == "hg19") {
    chrom <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", 
               "10", "11", "12", "13", "14", "15", "16", "17", "18", 
               "19", "20", "21", "22", "X", "Y")
    chrom <- paste("chr", chrom, sep = "")
    size <- c(249250621, 243199373, 198022430, 191154276, 
              180915260, 171115067, 159138663, 146364022, 141213431, 
              135534747, 135006516, 133851895, 115169878, 107349540, 
              102531392, 90354753, 81195210, 78077248, 59128983, 
              63025520, 48129895, 51304566, 155270560, 59373566)
    data.frame(chrom = chrom, size = size, stringsAsFactors = FALSE)
  }
}

#' convertToGRanges
#'
#' high-level interface for GRanges object creation
#'
#' @export convertToGRanges
convertToGRanges <- function (chr, startpos, endpos, strands, chromLevels = NULL, chromSizes = NULL, ...) {
  chrChar <- as.character(chr)
  chrChar[chrChar == "23"] <- "X"
  chrChar[chrChar == "24"] <- "Y"
  if (!is(strands, "Rle")) {
    strands <- Rle(strands)
  }
  gr <- GRanges(seqnames = Rle(chrChar), ranges = IRanges(start = as.numeric(startpos),  end = as.numeric(endpos)), strand = strands, ...)
  if (!is.null(chromLevels) & !is.null(chromSizes)) {
    seqlevels(gr, pruning.mode = "coarse") <- as.character(chromLevels)
    seqlengths(gr) <- chromSizes
    gr <- sort(gr)
  }
  gr
}

#' PFEAT.exons
#'
#' @param H.data GRanges object, H3K4me3 data, unstranded
#' @param T.data GRanges object, TFBS data, unstranded
#' @param F.data GRanges object, CAGE peaks data, stranded
#' @param S.data GRanges object, Splicing Site data, stranded
#'
#' @export PFEAT.exons
PFEAT.exons <- function(H.data, T.data, F.data, S.data, ncores=1) {
  # remove meta information so to speed up the calculations
  mcols(H.data) <- NULL
  mcols(T.data) <- NULL
  mcols(F.data) <- NULL
  mcols(S.data) <- NULL
  
  # HTFS regions
  predictor <- H.data
  predictor <- subsetByOverlaps(predictor, T.data, ignore.strand=TRUE)  # because strands of H3K4me3 data is *, 'ignore.strand' value is not relevant
  predictor <- subsetByOverlaps(predictor, F.data, ignore.strand=TRUE)  # because strands of H3K4me3 data is *, 'ignore.strand' value is not relevant
  predictor <- subsetByOverlaps(predictor, S.data, ignore.strand=TRUE)  # because strands of H3K4me3 data is *, 'ignore.strand' value is not relevant
  
  # HTFS merged regions
  predictor <- reduce(predictor)
  
  
  predictor.hits.ss20 <- findOverlaps(predictor, S.data, ignore.strand=TRUE)
  predictor.hits.cage <- findOverlaps(predictor, F.data, ignore.strand=TRUE)
  predictor.hits.tfbs <- findOverlaps(predictor, T.data, ignore.strand=TRUE)
  
  predictor.df <- do.call("rbind", mclapply(1:length(predictor), function(h.idx) {
    h <- predictor[h.idx]
    
    s.idx <- subjectHits(predictor.hits.ss20[queryHits(predictor.hits.ss20) == h.idx])
    f.idx <- subjectHits(predictor.hits.cage[queryHits(predictor.hits.cage) == h.idx])
    t.idx <- subjectHits(predictor.hits.tfbs[queryHits(predictor.hits.tfbs) == h.idx])
    
    s <- S.data[s.idx]
    f.pos <- follow(s, F.data[f.idx])
    t.pos <- follow(s, T.data[t.idx])
    
    bool <- !is.na(f.pos) & !is.na(t.pos)
    
    data.frame(hidx=rep(h.idx, length(which(bool))), sidx=s.idx[bool], fidx=f.idx[f.pos[bool]], tidx=t.idx[t.pos[bool]])
  }, mc.cores=ncores))
  
  predictor.df$chr     <- as.character(seqnames(predictor[predictor.df$hidx]))
  predictor.df$strand  <- as.character(strand(S.data[predictor.df$sidx]))
  
  # For strand +, I define the exon from the beginning of the CAGE peak to the beginning of the Splicing Site
  predictor.df$start   <- start(F.data[predictor.df$fidx])
  predictor.df$end     <- start(S.data[predictor.df$sidx])
  
  # For strand -, I define the exon from the end of the Splicing Site to the end of the CAGE peak
  neg <- predictor.df$strand == "-"
  predictor.df$start[neg] <- end(S.data[predictor.df$sidx])[neg] 
  predictor.df$end[neg]   <- end(F.data[predictor.df$fidx])[neg] 
  
  gr.novel <- convertToGRanges(chr=predictor.df$chr,
                               startpos=predictor.df$start,
                               endpos=predictor.df$end,
                               strands=predictor.df$strand,
                               chromLevels=chromInfo()$chrom,
                               chromSizes=chromInfo()$size,
                               hidx=predictor.df$hidx,
                               sidx=predictor.df$sidx,
                               fidx=predictor.df$fidx,
                               tidx=predictor.df$tidx)
  gr.novel
}

# load the preprocessed data. They are GRanges objects, metadata info not relevant here
load("H3K4me3.brain.RData")
load("cage.brain.RData")
load("tfbs.brain.RData")
load("ss.brain.d20.f05.RData")

# find exons
gr.htfs.exons.all <- HTFS.exons(H.data=H3K4me3.brain, T.data=tfbs.brain, F.data=cage.brain, S.data=ss.brain.d20.f05, ncores=4)

# merge ovelapping exons
gr.htfs.exons <- reduce(gr.htfs.exons.all)

# get sequences
gr.htfs.exons.seqs <- getSeq(x=BSgenome.Hsapiens.UCSC.hg19, gr.htfs.exons)
names(gr.htfs.exons.seqs) <- paste(seqnames(gr.htfs.exons), ":", start(gr.htfs.exons), "-", end(gr.htfs.exons), ":", strand(gr.htfs.exons), sep="")

# write out
write.table(as.data.frame(gr.htfs.exons)[,c(1,2,3,5)], file="PFEAT_exons.bed", append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
writeXStringSet(x=gr.htfs.exons.seqs, filepath="PFEAT_exons.fasta")
