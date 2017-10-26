library(data.table)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(Hmisc)
library(vegan)
library(Biostrings)
# Global variables -----
mapq.min=20
lstart.min=Inf

# Functions ------
# Calculate Quantile Coefficient of Dispersion QCD
qcd.fun <- function(x){
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr=(q3-q1)/2
  midhinge=(q1+q3)/2
  qcd <- iqr/midhinge
  return(qcd)
}

# downamples to values in sample size, as well as to the min sample size. 
# ed.count.matr is a simple matrix of counts. samples shows the sample sizes desired. seed sets the random seed
downsample.reads <- function(ed.count.matr, sample.sizes=NA, seed){
  min.sample.size <- min(colSums(ed.count.matr))
  set.seed(seed)
  samples <- c(min.sample.size, sample.sizes)
  samples <- samples[!is.na(samples)]
  draref.list <- sapply(seq_along(samples), simplify=FALSE, FUN=function(x){
    size <- samples[x]
    message(paste0("starting on draref ", size))
    rf <- t(drarefy(t(ed.count.matr), sample = size))
    message(paste0("done draref ", size))
    return(rf)  
  })
  names(draref.list) <- paste0("draref.", samples)
  rraref.list <- sapply(seq_along(samples), simplify=FALSE, FUN=function(x){
    size <- samples[x]
    message(paste0("starting on rraref ", size))
    rf <- t(rrarefy(t(ed.count.matr), sample = size))
    message(paste0("done rraref ", size))
    return(rf)  
  })
  names(rraref.list) <- paste0("rraref.", samples)
  raref.list <- c(draref.list, rraref.list)
  df.list <- names(raref.list)
  read.count.combined <- Reduce(function(...) merge(..., all=TRUE), sapply(seq_along(df.list), simplify=FALSE, FUN=function(x){
    nm <- df.list[x]
    d.dt <- data.table(melt(raref.list[[x]]))
    setnames(d.dt, c("miRNA.ID2", "fileName", nm))
    d.dt[, c("miRNA.ID", "editing"):=tstrsplit(miRNA.ID2, split=".", fixed=TRUE)]
    setkeyv(d.dt, c("miRNA.ID2", "miRNA.ID", "editing", "fileName"))
  }))
  
  return(read.count.combined)
  
  
}

# Set custom plot theme ------------
my_theme_bw <- theme_bw(base_size=8)
theme_set(my_theme_bw)
my_theme_bw <- theme_update(text=element_text(color="black"), axis.text=element_text(color="black", size=rel(1.0)))
# Editing pool annotation files -----
edited.mirs.bedfile <- "data/editing/edited_miRs.bed"
editing.pool.info.file <- "data/editing/editingPoolIDInfo.txt"
editing.count.Trim.file.gz <- "data/editing/ALL_MT_DG_KW_q10_from_trim4NBowtie.txt.gz"
editing.count.Trim.file <- "data/editing/ALL_MT_DG_KW_q10_from_trim4NBowtie.txt"
pool.sequences.file <- "data/editing/editing_pool_sequences.fa"

if(!file.exists(editing.count.Trim.file)){
  warning(paste0("Decompressing: ", editing.count.Trim.file.gz))
  R.utils::gunzip(editing.count.Trim.file.gz, remove=false)
}

# OUTPUT DIRECTORIES --------
# Same as from the main rebuttal code, but assumes that code was run first, becasue directories are not created here
top.output.directories <- c("main_figures", "supplemental_figures", "tables", "geo_submission")
main.figures=paste0("FIG", 2:7)
geo.dirs=c("geo_processed_data_files", "geo_metadata_files")
to.create <- unlist(sapply(top.output.directories, USE.NAMES = FALSE, simplify=FALSE, function(x){
  if(x=="main_figures"){
    these.dirs <- paste0("./output/", x, "/", main.figures)
  } else if(x=="geo_submission"){
    these.dirs <- paste0("./output/", x, "/", geo.dirs)
  } else{
    these.dirs <- paste0("./output/", x)
  }
}), recursive = TRUE)
# Provides path to output directories for the 7 main figures and any supplemental tables and figures
outdirs <- to.create
names(outdirs) <- basename(outdirs)

# Import editing data ----
edited.mirs <- fread(edited.mirs.bedfile)
setnames(edited.mirs, c("miRNA.ID", "lStart.miR", "lEnd.miR"))
editing.pool.info <- fread(editing.pool.info.file)

# Editing sequences
pool.sequences <- readDNAStringSet(pool.sequences.file, use.names = TRUE)
pool.sequences.rna <- RNAStringSet(pool.sequences)
edited.sequences.names <- c(edited.mirs$miRNA.ID, paste0(edited.mirs$miRNA.ID, "-edited"))
edited.sequences <- pool.sequences.rna[names(pool.sequences.rna)%in%edited.sequences.names]
edited.sequences.orig <- edited.sequences[!grepl("edited", names(edited.sequences))]
edited.sequences.edit <- edited.sequences[grepl("edited", names(edited.sequences))]
n.seqs <- length(edited.sequences.orig)
editing.sites <- sapply(seq(n.seqs), function(x) compareStrings(edited.sequences.orig[x], edited.sequences.edit[x]))
editing.sites.nt <- sub("?", "[A/I]", x = editing.sites, fixed=TRUE)
edited.miR.seq.table <- data.frame(miRNA.ID=names(edited.sequences.orig), sequence=editing.sites, seq.ai=editing.sites.nt)
# Table created with positions of editing events. Left here for convenience, but not included in final publication data.
# write.table(edited.miR.seq.table, "Edited_miRs_with_sequences.txt", row.names=FALSE, col.names=TRUE,sep="\t", quote=FALSE) 

# Import editing count data
import.editing.dat <- function(count.file, pool.info.dt){
  editing.counts.local.dt <- fread(count.file)
  editing.counts.local.dt[, fileName:=sub("TruSeq-", "TS-", fileName)]
  setnames(editing.counts.local.dt, "miR.ID", "miRNA.ID")
  editing.counts.local.dt[, c("pool.ID", "lib.method", "lab", "replicate"):=tstrsplit(fileName, split="_|-")[1:4]]
  editing.counts.local.dt[, replicate:=as.character(ifelse(replicate=="A", "1", ifelse(replicate=="B", "2", ifelse(replicate=="C", "3", replicate))))]
  editing.counts.local.dt[ , editing:=ifelse(grepl("edit", miRNA.ID)==TRUE, "edited", "unedited")]
  editing.counts.local.dt[, miRNA.ID:=sub("-edited", "", miRNA.ID)]
  lab.id <- c(DG="LAB2", MT="LAB5", KW="LAB6")
  pool.info.dt[, lab2:=lab.id[lab]]
  editing.counts.local.dt[, lib.method2:=ifelse(lib.method=="4N", "4N_B", ifelse(lib.method=="TS", "TruSeq", ifelse(lib.method=="NEB", "NEBNext", NA)))]
  setkeyv(pool.info.dt, c("lab", "pool.ID"))
  setkeyv(editing.counts.local.dt, c("lab", "pool.ID"))
  merge(editing.counts.local.dt, pool.info.dt, all.x = TRUE)
  
}

filter.and.summarize <- function(count.dt, edited.mirs.dt, mapq.cutoff=0, lStart.min=Inf){
  count.dt[, tot.reads:=sum(readCount), by=.(lab, pool.ID, percent.edited.expected, lib.method, replicate)]
  count.dt.filt <- subset(count.dt, MAPQ>=mapq.cutoff & lStart.miR<=lStart.min)
  count.dt.filt[, tot.reads:=sum(readCount), by=.(lab, pool.ID, percent.edited.expected, lib.method, replicate)]
  
  count.dt.summary <- count.dt[, .(meanQual=wtd.mean(MAPQ, readCount),
                                   readCount=sum(readCount)), 
                               by=.(fileName,lab, pool.ID, percent.edited.expected, lib.method, replicate, tot.reads, miRNA.ID, editing)]
  count.dt.filt.summary <- count.dt.filt[, .(meanQual=wtd.mean(MAPQ, readCount),
                                             readCount=sum(readCount)), 
                                         by=.(fileName, lab, pool.ID, percent.edited.expected, lib.method, replicate, tot.reads, miRNA.ID, editing)]
  count.dt.filt.summary[, `:=`(filtered.set="FILTERED",
                               pool.ID2=paste0("POOL_Edit_", percent.edited.expected*100))]
  count.dt.summary[, `:=`(filtered.set="UNFILTERED",
                          pool.ID2=paste0("POOL_Edit_", percent.edited.expected*100))]
  
  summary.merge <- rbind(count.dt.summary, count.dt.filt.summary)
  summary.merge.edited <- subset(summary.merge, miRNA.ID%in%edited.mirs.dt$miRNA.ID)
  return(list(all=summary.merge, edited=summary.merge.edited))
}
filter.and.summarize2 <- function(count.dt, edited.mirs.dt, mapq.cutoff=0, lStart.min=Inf){
  count.dt[, tot.reads:=sum(readCount), by=.(lab, pool.ID, percent.edited.expected, lib.method, replicate)]
  sample.info <- count.dt[, .N, by=.(lab, pool.ID, fileName, lib.method, replicate, percent.edited.expected, tot.reads)][, N:=NULL]
  setkey(sample.info, fileName)
  merge.cols <- colnames(sample.info)
  count.dt.filt <- subset(count.dt, MAPQ>=mapq.cutoff & lStart.miR<=lStart.min)
  count.dt.filt[, tot.reads:=sum(readCount), by=.(lab, pool.ID, percent.edited.expected, lib.method, replicate)]
  sample.info.filt <- count.dt.filt[, .N, by=.(lab, pool.ID, fileName, lib.method,  replicate, percent.edited.expected, tot.reads)][, N:=NULL]
  setkey(sample.info.filt, fileName)
  count.dt.summary <- count.dt[, .(readCount=sum(readCount)), 
                               by=.(fileName,lab, pool.ID, percent.edited.expected, lib.method, replicate, tot.reads, miRNA.ID, editing)]
  count.dt.summary.cast.remelt <- melt(dcast.data.table(count.dt.summary, miRNA.ID+editing~fileName, value.var = "readCount", fill=0, fun.aggregate = sum), id.vars = c("miRNA.ID", "editing"), variable.name = "fileName", value.name = "readCount")
  setkey(count.dt.summary.cast.remelt, fileName)
  count.dt.summary2 <- merge(sample.info, count.dt.summary.cast.remelt)
  
  count.dt.filt.summary <- count.dt.filt[, .(readCount=sum(readCount)), 
                                         by=.(fileName, lab, pool.ID, percent.edited.expected, lib.method, replicate, tot.reads, miRNA.ID, editing)]
  count.dt.summary.cast.filt.remelt <- melt(dcast.data.table(count.dt.filt.summary, miRNA.ID+editing~fileName, value.var = "readCount", fill=0, fun.aggregate = sum), id.vars = c("miRNA.ID", "editing"), variable.name = "fileName", value.name = "readCount")
  setkey(count.dt.summary.cast.filt.remelt, fileName)
  count.dt.filt.summary2 <- merge(sample.info.filt, count.dt.summary.cast.filt.remelt)
  count.dt.filt.summary2[, `:=`(filtered.set="FILTERED",
                               pool.ID2=paste0("POOL_Edit_", percent.edited.expected*100))]
  count.dt.summary2[, `:=`(filtered.set="UNFILTERED",
                          pool.ID2=paste0("POOL_Edit_", percent.edited.expected*100))]
  
  summary.merge <- rbind(count.dt.summary2, count.dt.filt.summary2)
  summary.merge.edited <- subset(summary.merge, miRNA.ID%in%edited.mirs.dt$miRNA.ID)
  return(list(all=summary.merge, edited=summary.merge.edited))
}

editing.count.4N.trim <- import.editing.dat(count.file = editing.count.Trim.file, pool.info.dt = editing.pool.info)
editing.count.4N.trim[, `:=`(lab=lab2, lib.method=lib.method2)][, `:=`(lab2=NULL, lib.method2=NULL)]
editing.count.4N.trim.summary.list <- filter.and.summarize2(count.dt = editing.count.4N.trim, edited.mirs.dt = edited.mirs, mapq.cutoff = mapq.min, lStart.min = lstart.min)
editing.count.4N.trim.summary.all <- editing.count.4N.trim.summary.list[[1]]
editing.count.4N.trim.summary.edited.miRs <- editing.count.4N.trim.summary.list[[2]]
editing.count.4N.trim.summary.cast.edited.miRs <- dcast.data.table(editing.count.4N.trim.summary.edited.miRs, ... ~ editing, value.var = c("readCount"), fun.aggregate = mean)
setnames(editing.count.4N.trim.summary.cast.edited.miRs, c("edited", "unedited"), c("readCount_edited", "readCount_unedited"))
editing.count.4N.trim.summary.cast.edited.miRs[, percent.edited.observed:=readCount_edited/sum(readCount_edited, readCount_unedited), by=1:nrow(editing.count.4N.trim.summary.cast.edited.miRs)]
editing.count.4N.trim.summary.cast.filtered.set <- dcast.data.table(editing.count.4N.trim.summary.edited.miRs, ... ~ filtered.set, value.var = c("readCount", "tot.reads"), fun.aggregate = max)

# make QC Summary table ----
editing.count.4N.trim.summary.cast.filtered.set.summary <- editing.count.4N.trim.summary.cast.filtered.set[, .(tot.reads_UNFILTERED=max(tot.reads_UNFILTERED), tot.reads_FILTERED=max(tot.reads_FILTERED)), by=.(fileName, pool.ID2, lab, lib.method, replicate, percent.edited.expected)]
setnames(editing.count.4N.trim.summary.cast.filtered.set.summary, c("fileName", "pool.ID2", "lib.method", "tot.reads_FILTERED", "tot.reads_UNFILTERED"), c("fileBaseName", "Editing Pool ID", "lib.method.detail", "Total Reads (Filtered)", "Total Reads(Unfiltered)"))
this.outfile <- paste0(outdirs["tables"], "/FOR_TableS8_Editing_Sample_QC.xlsx")
# Table S8 Fodder ---- 
# More manipulation was done in Excel to get final table S8
write.xlsx(editing.count.4N.trim.summary.cast.filtered.set.summary, this.outfile,  col.names=TRUE, row.names=FALSE)

# GEO Submission metadata ----
lab.run.read.length <- c(75, 69, 85)
names(lab.run.read.length) <- c("Lab2", "Lab5", "Lab6")
fun.paste <- function(x){
  paste(unique(x), collapse=";")
}
editing.count.4N.trim.summary.cast.filtered.set.summary.geo <- copy(editing.count.4N.trim.summary.cast.filtered.set.summary)
setnames(editing.count.4N.trim.summary.cast.filtered.set.summary.geo, "Editing Pool ID", "Editing.Pool.ID")
editing.count.4N.trim.summary.cast.filtered.set.summary.geo[, `:=`(pool.id.geo=paste0("SynthEdit", 100*percent.edited.expected)), by=1:nrow(editing.count.4N.trim.summary.cast.filtered.set.summary.geo)]
editing.count.4N.trim.summary.cast.filtered.set.summary.geo[, `:=`(fastq.file=sub("_EDITINGPOOL.bam", ".fastq.gz", as.character(fileBaseName)),
                                                                   source=paste0(sub("POOL_Edit_", "Synthetic smallRNA; Editing Pool ", Editing.Pool.ID), "%"),
                                                                   organism="Synthetic",
                                                                   molecule="total RNA",
                                                                   description="",
                                                                   lab=sub("LAB", "Lab", lab),
                                                                   lab.libMethod.Pool=paste(lib.method.detail, lab, pool.id.geo, sep=".")
                                                                   ), by=1:nrow(editing.count.4N.trim.summary.cast.filtered.set.summary.geo)]
editing.count.4N.trim.summary.cast.filtered.set.summary.geo.cast <- dcast.data.table(editing.count.4N.trim.summary.cast.filtered.set.summary.geo, lab.libMethod.Pool + lab.libMethod.Pool + source + organism + lab + lib.method.detail + pool.id.geo + percent.edited.expected + molecule + description ~ replicate, value.var=c("fileBaseName", "fastq.file"))
editing.count.4N.trim.summary.cast.filtered.set.summary.geo.processed.files <- editing.count.4N.trim.summary.cast.filtered.set.summary.geo[, .(fileBaseName=fileBaseName, file.type="bam", md5="x")]
editing.count.4N.trim.summary.cast.filtered.set.summary.geo.raw.files <- editing.count.4N.trim.summary.cast.filtered.set.summary.geo[, .(fastq.file=fastq.file, file.type="fastq", file.checksum="x", instrument.model="Illumina NextSeq 500", read.length=lab.run.read.length[lab], end="single")]
this.outfile <- paste0(outdirs["geo_metadata_files"], "/Geo_metadata_EditingPools.xlsx")
write.xlsx(editing.count.4N.trim.summary.cast.filtered.set.summary.geo.cast, this.outfile, sheetName = "EDITING.METADATA", row.names=FALSE)
write.xlsx(editing.count.4N.trim.summary.cast.filtered.set.summary.geo.processed.files, this.outfile, sheetName = "EDITING.PROCESSED.FILES", row.names=FALSE, append = TRUE)
write.xlsx(editing.count.4N.trim.summary.cast.filtered.set.summary.geo.raw.files, this.outfile, sheetName = "EDITING.RAW.FILES", row.names=FALSE, append = TRUE)

# Subset edited miRs ----
editing.count.4N.trim.edited.miRs <- subset(editing.count.4N.trim, miRNA.ID%in%edited.mirs$miRNA.ID)
editing.count.4N.trim.edited.miRs.FPs <- subset(editing.count.4N.trim.edited.miRs, (editing=="edited" & percent.edited.expected==0.0) | (editing=="unedited" & percent.edited.expected==100))

order.pool.id <- unique(editing.count.4N.trim.summary.cast.edited.miRs, by=c("pool.ID2", "percent.edited.expected"))
order.pool.id.s <- order.pool.id$percent.edited.expected
names(order.pool.id.s) <- order.pool.id$pool.ID2
order.pool.id.s <- order.pool.id.s[order(order.pool.id.s)]

# FIG6B Editing %Observed ----
editing.count.4N.trim.summary.cast.edited.miRs[, lib.method:=factor(lib.method, levels=c("TruSeq", "NEBNext", "4N_B"))]
observed.edit.ratio.g <- ggplot(editing.count.4N.trim.summary.cast.edited.miRs, aes(x=miRNA.ID, color=lab, y=percent.edited.observed)) +
  stat_summary(fun.data=median_hilow, size=0.25, pos=position_dodge(width = 0.5)) + 
  facet_wrap(~lib.method, nrow=3) +
  scale_y_continuous(labels = scales::percent) +
  theme(strip.background = element_rect(size=0.5), 
        strip.text = element_text(size = 8), 
        axis.text.x = element_text(angle=40, color="black", size=8, hjust=1),
        legend.position = "top",
        axis.text.y=element_text(color="black", size=8)) + labs(y="% Observed", x=NULL) +
  geom_hline(aes(yintercept = percent.edited.expected), lty=2)
legend <- get_legend(observed.edit.ratio.g)
observed.edit.ratio.g2 <- observed.edit.ratio.g + theme(legend.position = "none")

observed.edit.ratio.g.grid.filt <- plot_grid(plotlist = sapply(names(order.pool.id.s), simplify=FALSE, FUN=function(x){observed.edit.ratio.g2 %+% subset(editing.count.4N.trim.summary.cast.edited.miRs, pool.ID2 == x & filtered.set=="FILTERED")}), ncol=2, hjust = -0.2,labels = paste0(order.pool.id.s*100, "% edited")); observed.edit.ratio.g.grid.filt
observed.edit.ratio.g.grid.l <- plot_grid(legend, observed.edit.ratio.g.grid.filt , rel_heights = c(0.1, 3), ncol=1); observed.edit.ratio.g.grid.l
this.outfile <- paste0(outdirs["FIG6"], "/FIG6B_Editing_analysis_plot.pdf")
save_plot(this.outfile, plot = observed.edit.ratio.g.grid.l, base_height = 8, base_width = 7)


editing.count.4N.trim.summary.cast.edited.miRs[, percent.obs.minus.exp:=percent.edited.observed-percent.edited.expected]

# Try a mean abs deviation from expected using edited and unedited values
editing.count.4N.trim.summary.cast.edited.miRs[, read.count.both:=readCount_edited+readCount_unedited]
editing.count.4N.trim.summary.cast.edited.miRs[, `:=`(expected.count.edited=read.count.both*percent.edited.expected, expected.count.unedited=read.count.both*(1-percent.edited.expected))]
editing.count.4N.trim.summary.cast.edited.miRs[, abs.dev.expect:=abs(percent.edited.observed-percent.edited.expected)/((percent.edited.expected+(1-percent.edited.expected))/2)]

editing.count.4N.trim.summary.cast.edited.miRs[, `:=`(CPM_edited=(readCount_edited*10^6)/tot.reads, CPM_unedited=(readCount_unedited*10^6)/tot.reads)]
editing.count.4N.trim.summary.cast.edited.miRs.intralab.summary <- editing.count.4N.trim.summary.cast.edited.miRs[filtered.set=="FILTERED", .(
  intralab.mean.percent.edited=mean( percent.edited.observed, na.rm = TRUE), 
  intralab.SD.percent.edited=sd(percent.edited.observed), 
  intralab.qcd.percent.edited=qcd.fun(percent.edited.observed), 
  intralab.SD.CPM.edited=sd(CPM_edited,na.rm = TRUE), 
  intralab.SD.CPM.unedited=sd(CPM_unedited, na.rm=TRUE),
  intralab.mean.CPM.edited=mean(CPM_edited),
  intralab.mean.CPM.unedited=mean(CPM_unedited),
  intralab.qcd.CPM.edited=qcd.fun(CPM_edited),
  intralab.qcd.CPM.unedited=qcd.fun(CPM_unedited),
  intralab.min.abs.dif.expected=min(abs(percent.obs.minus.exp), na.rm = TRUE),
  intralab.max.abs.dif.expected=max(abs(percent.obs.minus.exp), na.rm=TRUE)), by=.(lab, pool.ID2, filtered.set, percent.edited.expected, lib.method, miRNA.ID)]

editing.count.4N.trim.summary.cast.edited.miRs.interlab.summary <- editing.count.4N.trim.summary.cast.edited.miRs.intralab.summary[filtered.set=="FILTERED", .(
  all.mean.percent.edited=mean(intralab.mean.percent.edited, na.rm = TRUE), 
  all.SD.percent.edited=mean(intralab.SD.percent.edited),
  all.qcd.percent.edited=qcd.fun(intralab.mean.percent.edited), 
  all.SD.CPM.edited=sd(intralab.mean.CPM.edited,na.rm = TRUE), 
  all.SD.CPM.unedited=sd(intralab.mean.CPM.unedited, na.rm=TRUE),
  all.mean.CPM.edited=mean(intralab.mean.CPM.edited),
  all.mean.CPM.unedited=mean(intralab.mean.CPM.unedited),
  all.qcd.CPM.edited=qcd.fun(intralab.mean.CPM.edited),
  all.qcd.CPM.unedited=qcd.fun(intralab.mean.CPM.unedited),
  all.min.abs.dif.expected=min(intralab.min.abs.dif.expected, na.rm = TRUE),
  all.max.abs.dif.expected=max(intralab.max.abs.dif.expected, na.rm=TRUE)), by=.(pool.ID2, percent.edited.expected, filtered.set, lib.method, miRNA.ID)]
editing.count.4N.trim.summary.cast.edited.miRs.interlab.summary[, all.CV.percent.edited:=all.SD.percent.edited/all.mean.percent.edited]


# Downsampling ----
sample.sizes=c(10^6)
seed.val=1000

editing.count.4N.trim.summary.all[, miRNA.ID2:=paste0(miRNA.ID, ".", editing), by=1:nrow(editing.count.4N.trim.summary.all)]
editing.count.4N.trim.df <- as.matrix(data.frame(dcast.data.table(editing.count.4N.trim.summary.all[filtered.set=="FILTERED"], miRNA.ID2~fileName, fill=0, value.var = "readCount", fun.aggregate = sum), row.names=1))

editing.count.4N.trim.summary.filt <- editing.count.4N.trim.summary.all[filtered.set=="FILTERED"]
editing.count.4N.trim.summary.filt[, fileName:=make.names(fileName)]
read.count.combined.Trim <- downsample.reads(ed.count.matr = editing.count.4N.trim.df, sample.sizes = sample.sizes, seed = seed.val)

mergecols <- intersect(colnames(editing.count.4N.trim.summary.filt), colnames(read.count.combined.Trim))
setkeyv(editing.count.4N.trim.summary.filt, mergecols)
setkeyv(read.count.combined.Trim, mergecols)
read.count.combined.Trim.si <- merge(editing.count.4N.trim.summary.filt, read.count.combined.Trim)
read.count.combined.Trim.si.edited.miRs <- subset(read.count.combined.Trim.si, miRNA.ID%in%edited.mirs$miRNA.ID)
read.count.combined.Trim.si.edited.miRs[, tot.counts.miR:=sum(readCount), by=.(fileName, miRNA.ID)]

this.outfile <- paste0(outdirs["tables"], "/FOR_TableS8_Editing_Sample_QC.xlsx")
# Table S8 Fodder DOWNSAMPLING---- 
# More manipulation was done in Excel to get final table S8. This adds the downsampling information
write.xlsx(read.count.combined.Trim.si.edited.miRs, this.outfile,  sheetName = "DOWNSAMPLING", append = TRUE, col.names=TRUE, row.names=FALSE)
write.xlsx(editing.count.4N.trim.summary.cast.edited.miRs.interlab.summary, this.outfile,  sheetName = "INTERLAB_Variation", append = TRUE, col.names=TRUE, row.names=FALSE)


