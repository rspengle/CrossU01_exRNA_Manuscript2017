# Load packages -------------------------------------------------------------------
library(data.table)
library(ggplot2)
library(Hmisc)
library(gplots)
library(edgeR)
library(stringr)
library(metaseqR)
library(DESeq2)
library(cqn)
library(RColorBrewer)
library(scales)
library(MASS)
library(GGally)
library(pheatmap)
library(variancePartition)
library(cowplot)
library(KernSmooth)
library(dendsort)
library(vegan)
library(ggbeeswarm)
library(xlsx)

# Load previously saved data -------------------------------------------------------------------
rdata.file="./data/20170901_EQ_Ratio_Plasma_MANUSCRIPT_DATA.RData"
load(rdata.file)
# Set custom plot theme -------------------------------------------------------------------
my_theme_bw <- theme_bw(base_size=8)
theme_set(my_theme_bw)
my_theme_bw <- theme_update(text=element_text(color="black"), axis.text=element_text(color="black", size=rel(1.0)))

# Subset of NEB labs that were done using identical protocols
NEB.subset <- c("NEBNext.Lab1", "NEBNext.Lab2", "NEBNext.Lab4", "NEBNext.Lab5")

# Set variable ordering -----------------------------------
pools <- c("SynthEQ", "Synth[AB]", "PlasmaPool", "*")

new.order <- sapply(pools, simplify=FALSE, function(x){
  lab.libs <- CROSS.U01.METADATA.ALL[grepl(x, pool), unique(as.character(lab.libMethod))]
  lab.libs <- lab.libs[order(lab.libs)]
  n.order <- c(
    grep("TruSeq.?", lab.libs, value=TRUE, ignore.case=TRUE), 
    grep("CleanTag.?", lab.libs, value=TRUE, ignore.case=TRUE),
    grep("NebNext.?", lab.libs, value=TRUE, ignore.case=TRUE),
    grep("4N_A", lab.libs, value=TRUE, ignore.case=TRUE),
    grep("4N_B", lab.libs, value=TRUE, ignore.case=TRUE),
    grep("4N_C", lab.libs, value=TRUE, ignore.case=TRUE),
    grep("4N_D", lab.libs, value=TRUE, ignore.case=TRUE),
    grep("4N_Xu", lab.libs, value=TRUE, ignore.case=TRUE),
    grep("4N_NEXTflex", lab.libs, value=TRUE, ignore.case=TRUE)
  )
})
names(new.order) <- c("equimolar", "ratiometric", "plasma", "all")
lib.method.simple.names <- unique(as.character(CROSS.U01.METADATA.ALL$lib.method.simple))
lib.method.simple.ordered <- c(
  grep("TruSeq", lib.method.simple.names, value=TRUE, ignore.case=TRUE), 
  grep("CleanTag", lib.method.simple.names, value=TRUE, ignore.case=TRUE),
  grep("NebNext", lib.method.simple.names, value=TRUE, ignore.case=TRUE),
  grep("4N", lib.method.simple.names, value=TRUE, ignore.case=TRUE)
)
lib.method.detail.names <- unique(as.character(CROSS.U01.METADATA.ALL$lib.method.detail))
lib.method.detail.ordered <- c(
  grep("TruSeq.?", lib.method.detail.names, value=TRUE, ignore.case=TRUE), 
  grep("CleanTag.?", lib.method.detail.names, value=TRUE, ignore.case=TRUE),
  grep("NebNext.?", lib.method.detail.names, value=TRUE, ignore.case=TRUE),
  grep("4N_A", lib.method.detail.names, value=TRUE, ignore.case=TRUE),
  grep("4N_B", lib.method.detail.names, value=TRUE, ignore.case=TRUE),
  grep("4N_C", lib.method.detail.names, value=TRUE, ignore.case=TRUE),
  grep("4N_D", lib.method.detail.names, value=TRUE, ignore.case=TRUE),
  grep("4N_Xu", lib.method.detail.names, value=TRUE, ignore.case=TRUE),
  grep("4N_NEXTflex", lib.method.detail.names, value=TRUE, ignore.case=TRUE)
)



# Custom color pallate for lib methods --------------------
# Custom palatte.. Does Truseq, CleanTag, NEB and 4N as normal. Then adds a greenish color for DOWNBIAS, and then takes the hue_pal_4 value and shades them for the different 4N methods
custom.lib.method.pallate <- c(hue_pal()(4)[1:3], hcl(285, c = 100, l = 65, alpha = 0.9*seq(4,1)/4), hue_pal(h=c(300,260))(2))
n.labs <- max(as.numeric(sub("Lab", "", unique(CROSS.U01.METADATA.ALL$Lab))))
ann_colors = list(
  lib.method.simple = custom.lib.method.pallate[1:4],
  lib.method.detail = custom.lib.method.pallate,
  Lab = brewer.pal(n.labs, "Blues")
)
lib.simple.levels <- lib.method.simple.ordered
lib.detail.levels <- lib.method.detail.ordered  

names(ann_colors$lib.method.simple) <- lib.simple.levels
names(ann_colors$lib.method.detail) <- lib.detail.levels
names(ann_colors$Lab) <- paste0("Lab", 1:n.labs)




# Set global variables -----------------------------------------------------------
# Set min and max lenghts to use (based on expected size of sequences in pool)
min.seqLen=16
max.seqLen=25
min.adjusted.count=1 # minimum adjusted count for plasma pool filtering

# MAKE OUTPUT DIRECTORIES -----------------------------------------------------------
top.output.directories <- c("main_figures", "supplemental_figures", "tables")
main.figures=paste0("FIG", 2:7)
to.create <- unlist(sapply(top.output.directories, USE.NAMES = FALSE, simplify=FALSE, function(x){
  if(x=="main_figures"){
    these.dirs <- paste0("./output/", x, "/", main.figures)
  } else{
    these.dirs <- paste0("./output/", x)
  }
}), recursive = TRUE)
sapply(to.create, function(x){
  if(!dir.exists(x)){
    dir.create(x, recursive = TRUE)
  }
})
# Provides path to output directories for the 7 main figures and any supplemental tables and figures

outdirs <- to.create
names(outdirs) <- basename(outdirs)
# The main figures have subdirectories for each figure
#outdirs
#FIG2                            FIG3                            FIG4                            FIG5                            FIG6 
#"./output/main_figures/FIG2"    "./output/main_figures/FIG3"    "./output/main_figures/FIG4"    "./output/main_figures/FIG5"    "./output/main_figures/FIG6" 
#FIG7            supplemental_figures                          tables 
#"./output/main_figures/FIG7" "./output/supplemental_figures"               "./output/tables" 


# FUNCTIONS -------------------------------------------------------------------
# Some functions for later
# make a function for averaging technical replicates
FUNCTION.aveTechReps <- function (x, ID = colnames(x), ...){
  if (is.null(x)) 
    return(NULL)
  x <- as.matrix(x)
  if (is.null(ID)) 
    stop("No sample IDs")
  t.x <- t(apply(x, 1, FUN=function(y){
    t(tapply(y, list(ID), mean))
  }))
  colnames(t.x) <- if(class(ID)=="factor"){
    levels(ID)
  } else{
    levels(as.factor(ID))
  }
  
  return(t.x)
  
}
FUNCTION.medianTechReps <- function (x, ID = colnames(x), ...){
  if (is.null(x)) 
    return(NULL)
  x <- as.matrix(x)
  if (is.null(ID)) 
    stop("No sample IDs")
  t.x <- t(apply(x, 1, FUN=function(y){
    t(tapply(y, list(ID), median))
  }))
  colnames(t.x) <- if(class(ID)=="factor"){
    levels(ID)
  } else{
    levels(as.factor(ID))
  }
  
  return(t.x)
  
}

# Function to collapse lane replicates from nextseq
# For character columns, collapses unique values
# char.override indicates column names that should be treated as a character and pasted
# Numeric values are summed
# group.by.var indicates the variables that other variables are collapsed on
# paste.collapse.character decides whether or not to paste the replicate characters. By default, it's off and only numeric values are returned. 
## This is helpful for when the tables are quite large and the pasted duplicates contain superfluous info that can be pulled more efficiently later
FUNCTION.collapse.lanes <- function(dt, char.override=NA, group.by.var, paste.collapse.character=FALSE){
  is.num.cols <- unlist(dt[,lapply(.SD, is.numeric)])
  num.cols <- names(is.num.cols[is.num.cols==TRUE])
  group.by.var <- group.by.var[group.by.var%in%colnames(dt)]
  if(is.na(char.override)){
    sum.cols <- num.cols
  } else{
    sum.cols <- num.cols[!num.cols%in%char.override]
  }
  paste.cols <- names(is.num.cols[!names(is.num.cols)%in%sum.cols]) # These will get pasted
  paste.cols <- unique(paste.cols[!paste.cols%in%group.by.var])
  
  # Collapse the metadata -- only important for the NextSeq files that came with 4 lanes per sample
  if(paste.collapse.character==TRUE){
    dt.collapse.num <- dt[, lapply(.SD, sum), by=eval(group.by.var), .SDcols=sum.cols]
    dt.collapse.char <- dt[, lapply(.SD, function(x) paste(unique(x), collapse=";")), by=eval(group.by.var), .SDcols=paste.cols]
    dt.collapse <- merge(dt.collapse.char, dt.collapse.num)
  } else{
    dt.collapse <- dt[, lapply(.SD, sum), by=eval(group.by.var), .SDcols=sum.cols]
  }
  
  return(dt.collapse)
}
FUNCTION.filter.ann_colors <- function(ann.col.list=ann_colors, ann.dt){
  ann.vars <- names(ann.col.list)
  ann.vars <- ann.vars[ann.vars%in%colnames(ann.dt)]
  filt.list <- sapply(ann.vars, FUN=function(x){
    old.list <- ann.col.list[[x]]
    lookup.v <- if(is.factor(ann.dt[, x])){
      f <- levels(ann.dt[,x])
      f[f%in%as.character(unique(ann.dt[,x]))]
    } else{
      f <- levels(factor(ann.dt[,x], levels=names(old.list)))
      f[f%in%ann.dt[,x]]
      
    }
    new.list <- old.list[lookup.v]
    return(new.list)
  })
  return(filt.list)
}
function.Fisher.Combined.PValue <- function(P_value_list)
{
  q_this=(-2)*sum(log(P_value_list))
  df_this=2*length(P_value_list)
  P_value=pchisq(q=q_this,df=df_this,lower.tail=F)
  unlist(P_value)
}



# END FUNCTIONS -------------------------------------------------------------------
# After loading, should see the following in the global env
#ls()
#[1] "counts.dt"                 "CROSS.U01.SYNTH.METADATA"  "full.seq.info.equimolar"   "full.seq.info.ratiometric" "sequence.annot"  


# annotate number of mods for filtering later
full.seq.info.equimolar[, n.mods.total:=.N, by=sequence]
full.seq.info.equimolar[, seq.len:=nchar(sequence)]
seq.ids.equimolar.5p.only <- subset(full.seq.info.equimolar, n.mods.total==1 & modification == "5'-phosphorylation")$equimolar.seqID

#setnames(counts.dt, "pool.ID", "pool")
#setnames(CROSS.U01.METADATA.ALL.CAST, "pool.ID", "pool")
#setnames(CROSS.U01.SYNTH.METADATA, "pool.ID", "pool")
#counts.dt[, lab.libMethod.pool.replicate:=paste(lib.method.detail, Lab, pool, replicate, sep=".")]

# Add a couple more grouping variables based on lab, method, pool, replicate, etc
# Make sure they're valid names as well
#CROSS.U01.METADATA.ALL.CAST[, `:=`(lab.libMethod.pool.replicate=paste(lib.method.detail, Lab, pool, replicate, sep="."),
     #                              lab.libMethod.pool=paste(lib.method.detail, Lab, pool, sep="."),
    #                               lab.libMethod.replicate=paste(lib.method.detail, Lab , replicate, sep="."),
   #                                lab.libMethod=paste(lib.method.detail, Lab, sep="."))]
#CROSS.U01.SYNTH.METADATA[, `:=`(lab.libMethod.pool.replicate=paste(lib.method.detail, Lab, pool, replicate, sep="."),
  #                              lab.libMethod.pool=paste(lib.method.detail, Lab, pool, sep="."),
 #                               lab.libMethod.replicate=paste(lib.method.detail, Lab , replicate, sep="."),
#                                lab.libMethod=paste(lib.method.detail, Lab, sep="."))]
#counts.dt[, `:=`(lab.libMethod.pool.replicate=(paste(lib.method.detail, Lab, pool, replicate, sep=".")),
#                   lab.libMethod.pool=(paste(lib.method.detail, Lab, pool, sep=".")),
#                   lab.libMethod.replicate=(paste(lib.method.detail, Lab , replicate, sep=".")),
#                   lab.libMethod=(paste(lib.method.detail, Lab, sep=".")))]

# Collapse the metadata -- only important for the NextSeq files that came with 4 lanes per sample
CROSS.U01.METADATA.ALL.CAST.collapse <- FUNCTION.collapse.lanes(CROSS.U01.METADATA.ALL.CAST, char.override = "lane", group.by.var = "lab.libMethod.pool.replicate", paste.collapse.character = TRUE)
CROSS.U01.METADATA.ALL.CAST.full <- CROSS.U01.METADATA.ALL.CAST
CROSS.U01.METADATA.ALL.CAST <- CROSS.U01.METADATA.ALL.CAST.collapse
CROSS.U01.SYNTH.METADATA <- CROSS.U01.METADATA.ALL.CAST[pool%in%c("SynthEQ", "SynthA", "SynthB")]

# Do the same collapsing to the counts.dt
counts.dt.collapse <- FUNCTION.collapse.lanes(counts.dt, char.override="lane", group.by.var = c("lab.libMethod.pool.replicate", colnames(sequence.annot)), paste.collapse.character = FALSE)
# Now just set the collapsed one to counts.dt for compatability
counts.dt.uncollapsed <- counts.dt
counts.dt <- counts.dt.collapse # this is the one to move forward with



## Wrangle Equimolar counts
# ADD ZERO-COUNTS BACK IN
# Cast to wide matrix, and fill missing with zeros
counts.cast <- dcast.data.table(counts.dt, new.seqID+equimolar.seqID+ratio.seqID+sequence+ratio.A+ratio.B~lab.libMethod.pool.replicate, fun.aggregate="sum", fill=0, value.var="count")
# Re-merge with all input sequences, so that any missing are added back in
setkeyv(counts.cast, c("new.seqID", "equimolar.seqID", "ratio.seqID", "sequence", "ratio.A", "ratio.B"))
setkeyv(sequence.annot, c("new.seqID", "equimolar.seqID", "ratio.seqID", "sequence", "ratio.A", "ratio.B"))
counts.cast.all <- merge(sequence.annot, counts.cast, all.x=TRUE)
read.cols <- unique(CROSS.U01.SYNTH.METADATA$lab.libMethod.pool.replicate)
counts.cast.all[, (read.cols):=lapply(.SD, function(x){ ifelse(is.na(x), 0, x)}), .SDcols=read.cols] # change NA's to 0's before re-melting
# Re-melt for easier plotting in ggplots.
counts.all.dt <- melt.data.table(counts.cast.all, measure.vars=read.cols, value.name="count", id.vars=c("new.seqID", "equimolar.seqID", "ratio.seqID", "ratio.A", "ratio.B" ,"sequence"), variable.name="lab.libMethod.pool.replicate")
# Lost the sample info when we casted, so pull those back in from the file.info table
setkey(counts.all.dt, lab.libMethod.pool.replicate)
setkey(CROSS.U01.SYNTH.METADATA, lab.libMethod.pool.replicate)
# Annotate some basic sequence info
counts.all.dt[, seq.len:=nchar(sequence)]
counts.all.dt[, gc.perc:=nchar(str_replace_all(sequence, pattern="[ATU]", replacement=""))/seq.len]
counts.all.dt.sizeFilt <- subset(counts.all.dt, seq.len>=min.seqLen & seq.len<=max.seqLen) # Filter for size
counts.all.dt.add.zerocount <- merge(CROSS.U01.SYNTH.METADATA, counts.all.dt.sizeFilt)


# Supplemental Table S3a: Synthetic Pool QC Metrics ----
EXCERPT.STATS.PLASMA.AND.SYNTHETIC.collapse <- FUNCTION.collapse.lanes(EXCERPT.STATS.PLASMA.AND.SYNTHETIC, char.override = c("lane"), group.by.var = c("lab.libMethod.pool.replicate", "Stage") , paste.collapse.character = TRUE)
setnames(EXCERPT.STATS.PLASMA.AND.SYNTHETIC.collapse, c("Calibrator", "NoCalibrator"), c("Calibrator.run", "noCalibrator.run"))
EXCERPT.STATS.SYNTHETIC.collapse <- dcast.data.table(subset(EXCERPT.STATS.PLASMA.AND.SYNTHETIC.collapse, pool!="PlasmaPool", select=-noCalibrator.run), ...~Stage, value.var = "Calibrator.run")
EXCERPT.STATS.SYNTHETIC.collapse[, `:=`(tot.failing.qc=InputReads-reads_used_for_alignment-calibrator)]
EXCERPT.STATS.SYNTHETIC.collapse[, `:=`(percent.input.failing.qc=tot.failing.qc/InputReads)]
EXCERPT.STATS.SYNTHETIC.collapse[, `:=`(reads.failing.adapter.trimming=tot.failing.qc-(failed_quality_filter + failed_homopolymer_filter + UniVec_contaminants + rRNA))]
qc.fail.cols <- c("reads.failing.adapter.trimming", "failed_quality_filter", "failed_homopolymer_filter", "UniVec_contaminants", "rRNA")
perc.qc.fail.ids <- paste0("percent.fail.", qc.fail.cols)
EXCERPT.STATS.SYNTHETIC.collapse[, (perc.qc.fail.ids):=lapply(.SD, FUN=function(x){ x/tot.failing.qc }), by=1:nrow(EXCERPT.STATS.SYNTHETIC.collapse), .SDcols=qc.fail.cols]
EXCERPT.STATS.SYNTHETIC.collapse[, `:=`(percent.tot.input.aligned.calib=calibrator/input, 
                                        percent.prefilt.input.aligned=calibrator/(input-(reads.failing.adapter.trimming + failed_quality_filter + failed_homopolymer_filter)),
                                        percent.synth.unmapped.mapping.genome=genome/reads_used_for_alignment)]
# Add in library size after filtering short miRs, since this is what is used for filtering
counts.all.dt.filt.libSize <- counts.all.dt.sizeFilt[, sum(count), by=lab.libMethod.pool.replicate]
setkey(counts.all.dt.filt.libSize, lab.libMethod.pool.replicate)
setkey(EXCERPT.STATS.SYNTHETIC.collapse, lab.libMethod.pool.replicate)
EXCERPT.STATS.SYNTHETIC.collapse[counts.all.dt.filt.libSize, calibrator.sizefilt.tot:=V1]
EXCERPT.STATS.SYNTHETIC.collapse[, `:=`(percent.tot.input.aligned.calib.sizeFilt=calibrator.sizefilt.tot/input, 
                                        percent.prefilt.input.aligned.sizeFilt=calibrator.sizefilt.tot/(input-(reads.failing.adapter.trimming + failed_quality_filter + failed_homopolymer_filter)))]
synth.stats.select.outcols <- c("lab.libMethod.pool.replicate", "Lab", "lib.method.detail", "lib.method.simple", "pool", "replicate", "FileBaseName", "adapter.confidence", "lane", "Notes", "InputReads", "tot.failing.qc", "percent.input.failing.qc", perc.qc.fail.ids, "calibrator", "percent.tot.input.aligned.calib", "percent.prefilt.input.aligned", "calibrator.sizefilt.tot", "percent.tot.input.aligned.calib.sizeFilt", "percent.prefilt.input.aligned.sizeFilt", "reads_used_for_alignment", "percent.synth.unmapped.mapping.genome")
EXCERPT.STATS.SYNTHETIC.OUTPUT <- subset(EXCERPT.STATS.SYNTHETIC.collapse, select=synth.stats.select.outcols)


# Supplemental Table S3b: Plasma Pool QC Metrics ----
EXCERPT.STATS.PLASMA.collapse <- dcast.data.table(subset(EXCERPT.STATS.PLASMA.AND.SYNTHETIC.collapse, pool=="PlasmaPool", select=-Calibrator.run), ...~Stage, value.var = "noCalibrator.run")
EXCERPT.STATS.PLASMA.collapse[, `:=`(tot.failing.qc=InputReads-reads_used_for_alignment)]
EXCERPT.STATS.PLASMA.collapse[, `:=`(percent.input.failing.qc=tot.failing.qc/InputReads)]
EXCERPT.STATS.PLASMA.collapse[, `:=`(reads.failing.adapter.trimming=tot.failing.qc-(failed_quality_filter + failed_homopolymer_filter + UniVec_contaminants + rRNA))]
EXCERPT.STATS.PLASMA.collapse[, (perc.qc.fail.ids):=lapply(.SD, FUN=function(x){ x/tot.failing.qc }), by=1:nrow(EXCERPT.STATS.PLASMA.collapse), .SDcols=qc.fail.cols]
EXCERPT.STATS.PLASMA.collapse[, `:=`(percent.input.used.for.mapping=reads_used_for_alignment/input, 
                                     genome.mapped.percent.alignment.input=genome/reads_used_for_alignment,
                                     miRNA.mapped.percent.alignment.input=miRNA_sense/reads_used_for_alignment,
                                     unmapped.percent.alignment.input=not_mapped_to_genome_or_libs/reads_used_for_alignment)]
plasma.stats.select.outcols <- c("lab.libMethod.pool.replicate", "Lab", "lib.method.detail", "lib.method.simple", "pool", "replicate", "FileBaseName", "adapter.confidence", "excerpt.qc.call.genomeMapped", "lane", "Notes", "InputReads", "tot.failing.qc", "percent.input.failing.qc", perc.qc.fail.ids, "reads_used_for_alignment", "percent.input.used.for.mapping", "genome", "genome.mapped.percent.alignment.input", "miRNA_sense", "miRNA.mapped.percent.alignment.input", "not_mapped_to_genome_or_libs", "unmapped.percent.alignment.input")
EXCERPT.STATS.PLASMA.OUTPUT <- subset(EXCERPT.STATS.PLASMA.collapse, select=plasma.stats.select.outcols)
this.outfile <- paste0(outdirs["tables"], "/Table_S3_PlasmaQC_Metrics_Synth_and_Plasma.xlsx")
write.xlsx(EXCERPT.STATS.SYNTHETIC.OUTPUT, file = this.outfile, sheetName = "SYNTH_POOLS_QC_METRICS", row.names = FALSE)
write.xlsx(EXCERPT.STATS.PLASMA.OUTPUT, file = this.outfile, sheetName = "PLASMA_POOLS_QC_METRICS", row.names = FALSE, append = TRUE)



# EQUIMOLAR WRANGLING ------------------------
equimolar.counts <- subset(counts.all.dt.add.zerocount, pool=="SynthEQ" & !is.na(equimolar.seqID))
# Keep only equimolar seequences with a 5' Phosphate. Remove if there is another ID with the same sequence, but a different modification (see above)
equimolar.counts.5p.only.SizeFilt <- subset(equimolar.counts, equimolar.seqID%in%seq.ids.equimolar.5p.only)
n.unique.seqs <- equimolar.counts.5p.only.SizeFilt[, length(unique(equimolar.seqID))]
expected.cpm <- 10^6/n.unique.seqs
equimolar.counts.5p.only.SizeFilt[, count.plus1:=count+1]
equimolar.counts.5p.only.SizeFilt[, `:=`(count.total.by.sample.filt.lengths=sum(count), count.plus1.total.by.sample.filt.lengths=sum(count.plus1)), by=.(lab.libMethod.pool.replicate)]
equimolar.counts.5p.only.SizeFilt[, `:=`(
  cpm.filt.lengths=(count*10^6)/count.total.by.sample.filt.lengths,
  pseudo.cpm.filt.lengths=(count.plus1*10^6)/count.plus1.total.by.sample.filt.lengths,
  expected.count.filt.lengths=count.total.by.sample.filt.lengths/n.unique.seqs,
  expected.countplus1.filt.lengths=count.plus1.total.by.sample.filt.lengths/n.unique.seqs)]
equimolar.counts.5p.only.SizeFilt[, `:=`(logratio.count.filt.vs.expected=log2(count/expected.count.filt.lengths),
                                         logratio.countplus1.filt.vs.expected=log2(count.plus1/expected.countplus1.filt.lengths))]

equimolar.counts.5p.only.SizeFilt.cast.count <- dcast.data.table(equimolar.counts.5p.only.SizeFilt, equimolar.seqID+sequence+seq.len+gc.perc~lab.libMethod.replicate, value.var = "count")

equimolar.counts.5p.only.SizeFilt.cast.pseudo.cpm <- dcast.data.table(equimolar.counts.5p.only.SizeFilt, equimolar.seqID+sequence+seq.len+gc.perc~lab.libMethod.replicate, value.var = "pseudo.cpm.filt.lengths")

uniq.lab.lib.method.replicates <- unique(equimolar.counts.5p.only.SizeFilt$lab.libMethod.replicate)
equimolar.counts.5p.only.SizeFilt.cast.df.counts <- data.frame(subset(equimolar.counts.5p.only.SizeFilt.cast.count, select=c("equimolar.seqID", uniq.lab.lib.method.replicates)), row.names=1)
equimolar.counts.5p.only.SizeFilt.cast.df.pseudo.cpm <- data.frame(subset(equimolar.counts.5p.only.SizeFilt.cast.pseudo.cpm, select=c("equimolar.seqID", uniq.lab.lib.method.replicates)), row.names=1)
gene.info.eq <- data.frame(subset(equimolar.counts.5p.only.SizeFilt.cast.count, select=c("equimolar.seqID", colnames(equimolar.counts.5p.only.SizeFilt.cast.count)[!colnames(equimolar.counts.5p.only.SizeFilt.cast.count)%in%c("equimolar.seqID", uniq.lab.lib.method.replicates)])), row.names=1)

FINAL.EQUIMOLAR.LONG.AVERAGE.REPS.FILT <- equimolar.counts.5p.only.SizeFilt[, .(mean.cpm=mean(cpm.filt.lengths),
                                                                                mean.pseudo.cpm=mean(pseudo.cpm.filt.lengths),
                                                                                sd.cpm=sd(cpm.filt.lengths),
                                                                                sd.pseudo.cpm=sd(pseudo.cpm.filt.lengths),
                                                                                geomean.cpm=exp(mean(log(pseudo.cpm.filt.lengths))),
                                                                                geo.sd.cpm=exp(sd(log(pseudo.cpm.filt.lengths))),
                                                                                mean.logratio.countplus1.filt.vs.expected=mean(logratio.countplus1.filt.vs.expected),
                                                                                ln.cpm.test = sqrt(exp(sd(log(pseudo.cpm.filt.lengths))^2)-1),
                                                                                sd.logratio.countplus1.filt.vs.expected=sd(logratio.countplus1.filt.vs.expected),
                                                                                n.present = sum(ifelse(count>0, 1, 0))
), 
by=.(lab.libMethod, Lab, lib.method.detail, lib.method.simple, equimolar.seqID, sequence, seq.len, gc.perc)]
FINAL.EQUIMOLAR.LONG.AVERAGE.REPS.FILT[, `:=`(cv.cpm=(sd.pseudo.cpm/mean.pseudo.cpm), geo.cv.cpm=(geo.sd.cpm/geomean.cpm))]


FINAL.EQUIMOLAR.LONG.AVERAGE.REPS.FILT[, `:=`(lab.libMethod=factor(lab.libMethod, levels=new.order$equimolar),
                                              lib.method.simple=factor(lib.method.simple, levels=lib.method.simple.ordered),
                                              lib.method.detail=factor(lib.method.detail, levels=lib.method.detail.ordered))]

FINAL.EQUIMOLAR.LONG <- equimolar.counts.5p.only.SizeFilt # convert to the name we're used to from earlier versions
FINAL.EQUIMOLAR.LONG[, `:=`(lab.libMethod=factor(lab.libMethod, levels=new.order$equimolar),
                            lib.method.simple=factor(lib.method.simple, levels=lib.method.simple.ordered),
                            lib.method.detail=factor(lib.method.detail, levels=lib.method.detail.ordered))]

sample.info.eq <- data.frame(FINAL.EQUIMOLAR.LONG[, .(tot.count=max(count.total.by.sample.filt.lengths), tot.count.plus1=max(count.total.by.sample.filt.lengths)), by=.(lab.libMethod.replicate, lab.libMethod, Lab, lib.method.detail, lib.method.simple, replicate)], row.names=1)
sample.info.eq$log.lib.size <- log2(sample.info.eq$tot.count)
sample.info.eq$log.plus1.lib.size <- log2(sample.info.eq$tot.count.plus1)

sample.info.eq$lab.libMethod <- factor(sample.info.eq$lab.libMethod, levels=new.order$equimolar)
sample.info.eq$lib.method.simple <- factor(sample.info.eq$lib.method.simple, levels=lib.method.simple.ordered)
sample.info.eq$lib.method.detail <- factor(sample.info.eq$lib.method.detail, levels=lib.method.detail.ordered)

# Remove miRs missing in all samples. These are probably just technical errors
missing.mirs.from.all <- FINAL.EQUIMOLAR.LONG[, .(count.any.gt.0=sum(ifelse(count>0, 1, 0))), by=equimolar.seqID][count.any.gt.0==0]$equimolar.seqID
FINAL.EQUIMOLAR.LONG <- subset(FINAL.EQUIMOLAR.LONG, !equimolar.seqID%in%missing.mirs.from.all)
FINAL.EQUIMOLAR.LONG.AVERAGE.REPS.FILT <- subset(FINAL.EQUIMOLAR.LONG.AVERAGE.REPS.FILT, !equimolar.seqID%in%missing.mirs.from.all)

# KEY VARIABLES FROM EQUIMOLAR POOL WRANGLING
# equimolar.counts --> full results from equimolar pool. Unfiltered. Includes zerocounts. 
# equimolar.counts.5p.only --> Full equimolar pool. Filtered for 5'p only.
# FINAL.EQUIMOLAR.LONG --> Further filtered for size contraints used in paper. All replicates included. 
# FINAL.EQUIMOLAR.LONG.AVERAGE.REPS.FILT --> Averaged multiple things over replicates
# equimolar.counts.5p.only.SizeFilt.cast.df.counts --> data frame of count matrix, including all replicates
# equimolar.counts.5p.only.SizeFilt.cast.df.pseudo.cpm --> data frame of cpm, including all replicates
# sample.info.eq --> sample info for equimolar pool in data.frame format
# gene.info.eq --> gene info for equimolar pool in data.frame format





# FIG2 -------------------------------------------------------------------

# HEATMAP -- Unnormalized cpm +1 pseudo count
equimolar.counts.5p.only.SizeFilt.cast.df.pseudo.cpm.log2 <- log2(equimolar.counts.5p.only.SizeFilt.cast.df.pseudo.cpm)
colnames(equimolar.counts.5p.only.SizeFilt.cast.df.pseudo.cpm.log2) <- sub("^X", "", colnames(equimolar.counts.5p.only.SizeFilt.cast.df.pseudo.cpm.log2))
# gene.info.eq --> gene info for equimolar pool in data.frame format

# HEATMAP -- Use EdgeR to normalize first 
dge.eq <- DGEList(equimolar.counts.5p.only.SizeFilt.cast.df.counts)
dge.eq$samples <- cbind(dge.eq$samples, sample.info.eq)
dge.eq.filt <- dge.eq[, dge.eq$samples$lab.libMethod!="4N_NEXTflex.Lab8"]
dge.eq.norm.rle <- calcNormFactors(dge.eq.filt, method = "RLE")
cpm.rle <- cpm(dge.eq.norm.rle, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
eq.heat.col.annot <- dge.eq.filt$samples[, c("log.plus1.lib.size", "lib.method.simple", "lib.method.detail")]
# FIG2A -------------------------------------------------------------------
this.outfile <- paste0(outdirs["FIG2"], "/FIG2A_EQ_Heatmap_RLE_WITH_LEGEND.pdf")
pheatmap(cpm.rle, drop_levels=TRUE, treeheight_col = 20, treeheight_row = 10,
         show_rownames = FALSE, 
         show_colnames=TRUE,
         cluster_rows = TRUE,
         annotation_col = eq.heat.col.annot,
         fontsize = 8, 
         legend_breaks = seq(-4,20,2),
         annotation_colors = FUNCTION.filter.ann_colors(ann_colors, eq.heat.col.annot),
         width=6.5,
         height=4, 
         labels_col = sub("^X", "", colnames(cpm.rle)),
         border_color=NA,
         legend=FALSE, annotation_legend=FALSE)
pheatmap(cpm.rle, drop_levels=TRUE, 
         show_rownames = FALSE, 
         show_colnames=TRUE,
         cluster_rows = TRUE,
         annotation_col = eq.heat.col.annot,
         fontsize = 10, 
         legend_breaks = seq(-4,20,2),
         annotation_colors = FUNCTION.filter.ann_colors(ann_colors, eq.heat.col.annot),
         width=8.5,
         height=8, 
         labels_col = sub("^X", "", colnames(cpm.rle)),
         border_color=NA,
         filename = this.outfile
)
this.outfile <- paste0(outdirs["FIG2"], "/FIG2A_EQ_Heatmap_RLE_NO_LEGEND.pdf")
pheatmap(cpm.rle, drop_levels=TRUE, treeheight_col = 20, treeheight_row = 10,
         show_rownames = FALSE, 
         show_colnames=TRUE,
         cluster_rows = TRUE,
         annotation_col = eq.heat.col.annot,
         fontsize = 8, 
         legend_breaks = seq(-4,20,2),
         annotation_colors = FUNCTION.filter.ann_colors(ann_colors, eq.heat.col.annot),
         width=7.5,
         height=4, 
         labels_col = sub("^X", "", colnames(cpm.rle)),
         border_color=NA,
         filename = this.outfile,
         legend=FALSE, annotation_legend=FALSE)

# FIG2B -------------------------------------------------------------------
## EQUIMOLAR BIAS VIOLIN
this.outfile <- paste0(outdirs["FIG2"], "/FIG2B_ViolinPlot_equimolar_bias.pdf")
f2b <- ggplot(
  FINAL.EQUIMOLAR.LONG.AVERAGE.REPS.FILT[lab.libMethod!="4N_NEXTflex.Lab8"],
  aes(
    x=lab.libMethod,
    y=mean.pseudo.cpm,
    fill=lib.method.simple
  )
) + 
  geom_violin(
    draw_quantiles = c(0.25,0.5,0.75),
    alpha=1
  ) + 
  geom_hline(
    yintercept=expected.cpm,
    size=0.75,
    linetype = 2,
    alpha=0.8
  ) +
  scale_y_log10(
    limits=c(10^-2,10^5),
    breaks = scales::trans_breaks("log10", n = 7,  function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) + 
  annotation_logticks(sides = "l", size = 1, color="black") +
  theme(
    axis.text.x=element_text(angle=50, hjust=1),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1, color="black"),
    axis.ticks = element_line(color="black"),
    legend.position = "none"
    #legend.position = "top",
    #legend.direction = "horizontal",
    #legend.key.size=unit(1, "lines"),
    #legend.background = element_rect(fill=NA)
  ) +
  labs(
    #title="EQUIMOLAR POOL:\nAVERAGE COUNTS-PER-MILLION",
    x=NULL,
    y="Average CPM",
    fill=NULL
  ) 
print(f2b)
ggsave(this.outfile, width=4.6, height=2.75, units = "in")

## Now use a cutoff to get miRs > x fold from expected and see what % miRs in the different methods
# FIG2C -------------------------------------------------------------------
n.seqs <- length(unique(FINAL.EQUIMOLAR.LONG$equimolar.seqID))
cutoff.from.exp = log2(10)
FINAL.EQUIMOLAR.LONG.count.simple.expected <- FINAL.EQUIMOLAR.LONG[, .(n.outside.cutoff=sum(ifelse(abs(logratio.countplus1.filt.vs.expected)>=cutoff.from.exp, 1, 0)), n.total=.N), by=.(lab.libMethod, lab.libMethod.replicate, lib.method.detail, lib.method.simple, Lab)]
FINAL.EQUIMOLAR.LONG.count.simple.expected[, percent.outside.cutoff:=n.outside.cutoff/n.total]

this.outfile <- paste0(outdirs["FIG2"], "/FIG2C_SIMPLE_10x_expected_labCombined.pdf")
ggplot(FINAL.EQUIMOLAR.LONG.count.simple.expected[lab.libMethod!="4N_NEXTflex.Lab8"],
       aes(y=percent.outside.cutoff,
           x=lib.method.detail,
           fill=lib.method.simple)
) + 
  geom_boxplot() +
  scale_y_continuous(limits = c(0, .7), breaks=seq(0,1,0.1), labels = scales::percent) +
  theme(
    axis.ticks.x=element_line(color=NA),
    axis.text.x=element_text(vjust=1, angle=50, hjust=1),
    panel.border = element_rect(size=1, color="black"),
    panel.grid.major=element_line(color=NA),
    panel.grid.minor=element_line(color=NA),
    legend.position = "none"
  ) + 
  labs(
    #title="Equimolar Pool: % > 10x Expected\nSimple Count Expected",
    y="% sequences > 10x from expected",
    fill="Library Prep Method",
    x=NULL
  ) 
ggsave(filename = this.outfile, height=2.75, width=3.7)
summary.from.expected.eq <- FINAL.EQUIMOLAR.LONG.count.simple.expected[, .(min=min(percent.outside.cutoff), max=max(percent.outside.cutoff), median=median(percent.outside.cutoff), mean=mean(percent.outside.cutoff)), .(by=lib.method.detail)]

# Table with median % sequences more than 10x from expectation for in-text values ----
median.percent.10x.expected.eq.by.method <- FINAL.EQUIMOLAR.LONG.count.simple.expected[, .(median.percent.10x.exp=median(percent.outside.cutoff)), by=lib.method.detail]
this.outfile <- paste0(outdirs["tables"], "/TextSummary_Fig2C_median_percent_10x_expected_eqPool.txt")
write.table(median.percent.10x.expected.eq.by.method, file = this.outfile, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)


this.outfile <- paste0(outdirs["FIG2"], "/FIG2C_SIMPLE_10x_expected_labSeparate.pdf")
g1 <- ggplot(FINAL.EQUIMOLAR.LONG.count.simple.expected[lab.libMethod!="4N_NEXTflex.Lab8"],
       aes(y=percent.outside.cutoff,
           x=lib.method.detail,
           pos=Lab,
           color=lib.method.simple)
) + 
  stat_summary(position = position_dodge(width=1.0), size=0.01, fun.y = "median", fun.ymax = "max", fun.ymin = "min") +
  scale_y_continuous(limits = c(0, .7), breaks=seq(0,1,0.1), labels = scales::percent) +
  theme(
    axis.text.x=element_text(vjust=1, angle=50, hjust=1),
    panel.border = element_rect(size=1, color="black"),
    panel.grid.major=element_line(color=NA),
    panel.grid.minor=element_line(color=NA),
    legend.position = "none"
  ) + 
  labs(
    y="% sequences > 10x from expected",
    color="Library Prep Method",
    x=NULL
  ) 
print(g1)
ggsave(filename = this.outfile, plot = g1, height=2.75, width=3.7, units = "in")

this.outfile <- paste0(outdirs["FIG2"], "/FIG2C_SIMPLE_10x_expected_labSeparate_QuasiRandom.pdf")
g <- ggplot(FINAL.EQUIMOLAR.LONG.count.simple.expected, 
       aes(
         x=lib.method.detail,
         y=percent.outside.cutoff, 
         color=lib.method.simple,
         pos=lab.libMethod)) + 
  stat_summary(
    #position = position_quasirandom(dodge.width = 0.5, method = "quasirandom"),
    position = position_beeswarm(cex = 2.5, priority = "descending"),
    size=0.1, 
    fun.y = "median", 
    fun.ymax = "max",
    fun.ymin = "min") +
  labs(x=NULL,
       y="% sequences > 10x from expected",
       color="Library Prep Method") +
  scale_y_continuous(breaks=seq(0, 1, 0.1), labels = scales::percent) + theme(
    axis.text.x=element_text(vjust=1, angle=50, hjust=1),
    panel.border = element_rect(color="black"),
    panel.grid.major=element_line(color=NA),
    panel.grid.minor=element_line(color=NA),
    legend.position = "none"
  ); g

ggsave(this.outfile, plot = g,  width = 4.5, height = 3, units = "in")

# FIG S3 Spearman Correlation Heatmaps for Eq Pool -------------------------------------------------------------------

pheatmap(cor(cpm.rle, use = "pair", method = "spearman")^2, 
         treeheight_row = 15, 
         treeheight_col = 15, 
         legend=FALSE,
         annotation_legend=FALSE,
         breaks=seq(0,1,length.out = 101),
         annotation_col = eq.heat.col.annot,
         annotation_row = eq.heat.col.annot,
         labels_col = sub("^X", "", colnames(cpm.rle)), 
         labels_row = sub("^X", "", colnames(cpm.rle)), 
         annotation_colors =  FUNCTION.filter.ann_colors(ann_colors, eq.heat.col.annot), 
         fontsize=8,
         width=10,
         height=8)
this.outfile <- paste0(outdirs["supplemental_figures"], "/FIGS3_EQPOOL_SPEARMAN_HEATMAP_NOLEGEND.pdf")
pheatmap(cor(cpm.rle, use = "pair", method = "spearman")^2, 
         treeheight_row = 15, 
         treeheight_col = 15, 
         legend=FALSE,
         annotation_legend=FALSE,
         breaks=seq(0,1,length.out = 101),
         annotation_col = eq.heat.col.annot,
         annotation_row = eq.heat.col.annot,
         labels_col = sub("^X", "", colnames(cpm.rle)), 
         labels_row = sub("^X", "", colnames(cpm.rle)), 
         annotation_colors =  FUNCTION.filter.ann_colors(ann_colors, eq.heat.col.annot), 
         fontsize=8,
         width=10,
         height=8,
         filename=this.outfile)
this.outfile <- paste0(outdirs["supplemental_figures"], "/FIGS3_EQPOOL_SPEARMAN_HEATMAP_WITHLEGEND.pdf")
pheatmap(cor(cpm.rle, use = "pair", method = "spearman")^2, 
         treeheight_row = 15, 
         treeheight_col = 15, 
         legend=TRUE,
         breaks=seq(0,1,length.out = 101),
         annotation_col = eq.heat.col.annot,
         annotation_row = eq.heat.col.annot,
         labels_col = sub("^X", "", colnames(cpm.rle)), 
         labels_row = sub("^X", "", colnames(cpm.rle)), 
         annotation_colors =  FUNCTION.filter.ann_colors(ann_colors, eq.heat.col.annot), 
         fontsize=8,
         width=10,
         height=8,
         filename=this.outfile)


# SUMMARY CORRELATION COEFFICIENT TABLES FOR EQUIMOLAR POOL -----------------------------
# Convert the cpm.rle matrix to a melted data.frame for data.table calculations of correlations.
cor.matrix.all.eq.square <- cor(cpm.rle, use="pair", method="spearman")
cor.matrix.all.eq <- data.table(cbind(melt(cor.matrix.all.eq.square), melt(upper.tri(cor.matrix.all.eq.square, diag=TRUE))))
setnames(cor.matrix.all.eq, c("lab.libMethod.replicate.A", "lab.libMethod.replicate.B", "rho", "row", "col", "upper.tri"))
cor.matrix.all.eq[, `:=`(id.A=ifelse(upper.tri==TRUE, as.character(lab.libMethod.replicate.A), as.character(lab.libMethod.replicate.B)), id.B=ifelse(upper.tri==TRUE, as.character(lab.libMethod.replicate.B), as.character(lab.libMethod.replicate.A)))]
cor.matrix.all.eq[, `:=`(lab.libMethod.replicate.A=id.A, lab.libMethod.replicate.B=id.B)]
cor.matrix.all.eq[, `:=`(id.A=NULL, id.B=NULL)]
short.sample.info.eq.A <- dge.eq.norm.rle$samples[, c("Lab", "lib.method.detail", "lib.method.simple", "replicate")]
short.sample.info.eq.B <- dge.eq.norm.rle$samples[, c("Lab", "lib.method.detail", "lib.method.simple", "replicate")]

colnames(short.sample.info.eq.A) <- paste0(colnames(short.sample.info.eq.A), ".A")
colnames(short.sample.info.eq.B) <- paste0(colnames(short.sample.info.eq.B), ".B")

cor.matrix.all.eq.tmp <- merge.data.frame(cor.matrix.all.eq, short.sample.info.eq.A, by.x=1, by.y=0)
cor.matrix.all.eq.sample.info <- data.table(merge.data.frame(cor.matrix.all.eq.tmp, short.sample.info.eq.B, by.x=2, by.y=0))
lib.detail.levels.with.nextFlex <- levels(cor.matrix.all.eq.sample.info$lib.method.simple.A)
lib.simple.levels <- levels(cor.matrix.all.eq.sample.info$lib.method.simple.A)
reorder.paste.fun <- function(query.vec, lookup.vec){
  matches <- match(query.vec, lookup.vec)
  comp <- paste(query.vec[order(matches)], collapse=".VS.")
  return(comp)
}

cor.matrix.all.eq.sample.info[, comparison.detail:=apply(.SD, 1, reorder.paste.fun, lookup.vec=lib.detail.levels.with.nextFlex), .SDcols=c("lib.method.detail.A", "lib.method.detail.B")]
cor.matrix.all.eq.sample.info[, comparison.simple:=apply(.SD, 1, reorder.paste.fun, lookup.vec=lib.simple.levels), .SDcols=c("lib.method.simple.A", "lib.method.simple.B")]

setkeyv(cor.matrix.all.eq.sample.info, c("lab.libMethod.replicate.B", "lab.libMethod.replicate.A"))
cor.matrix.all.eq.sample.info.unique <- unique(cor.matrix.all.eq.sample.info,  by=key(cor.matrix.all.eq.sample.info))


# form meanrho
FUNCTION.meanRhoFromCoefs <- function(r, ns, nr){
  delr <- length(r) - length(r[(r < 1) & (r > -1)])
  r <- r[(r < 1) & (r > -1)]
  rz <- 1/2 * log((1 + r)/(1 - r))
  mrz <- mean(rz)
  coeff <- (exp(2 * mrz) - 1)/(exp(2 * mrz) + 1)
  SE <- sqrt(1/(ns - 3))
  u <- coeff/SE
  mean.r <- mean(r)
  r02 <- quantile(r, 0.02)
  r98 <- quantile(r, 0.98)
  p.value <- 2 * (1 - pnorm(abs(u)))
  rval <- list(irr.name = "Rho", nlibs= nr, value = coeff, stat.name = "z", statistic = u, p.value = p.value, simple.mean.r=mean.r, r02=r02, r98=r98)
}
n.miRs <- dim(cpm.rle)[1]
cor.matrix.all.eq.sample.info.cor.summary.details <- cor.matrix.all.eq.sample.info.unique[, { nr=length(unique(c(lab.libMethod.replicate.A, lab.libMethod.replicate.B))); FUNCTION.meanRhoFromCoefs(rho, ns=n.miRs, nr=nr) }, by=comparison.detail]
cor.matrix.all.eq.sample.info.cor.summary.simple <- cor.matrix.all.eq.sample.info.unique[, { nr=length(unique(c(lab.libMethod.replicate.A, lab.libMethod.replicate.B))); FUNCTION.meanRhoFromCoefs(rho, ns=n.miRs, nr=nr) }, by=comparison.simple]
setnames(cor.matrix.all.eq.sample.info.cor.summary.details, "comparison.detail", "comparison")
setnames(cor.matrix.all.eq.sample.info.cor.summary.simple, "comparison.simple", "comparison")
cor.matrix.all.eq.sample.info.cor.summary.all <- rbind(cor.matrix.all.eq.sample.info.cor.summary.simple, cor.matrix.all.eq.sample.info.cor.summary.details)

this.outfile <- paste0(outdirs["tables"], "/Spearman_correlation_equimolar_Pool_summaries.txt")
write.table(cor.matrix.all.eq.sample.info.cor.summary.all, this.outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


# meanRho Summary Equimolar: All replicates, but no comparison across technical replicates----
cor.matrix.all.eq.sample.info.unique[, `:=`(lab.libMethod.A=sub(".[0-9]$", "", lab.libMethod.replicate.A), lab.libMethod.B=sub(".[0-9]$", "", lab.libMethod.replicate.B)) ]
cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison <- cor.matrix.all.eq.sample.info.unique[lab.libMethod.A!=lab.libMethod.B]

cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.details <- cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison[, { nr=length(unique(c(lab.libMethod.replicate.A, lab.libMethod.replicate.B))); FUNCTION.meanRhoFromCoefs(rho, ns=n.miRs, nr=nr) }, by=comparison.detail]
cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.simple <- cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison[, { nr=length(unique(c(lab.libMethod.replicate.A, lab.libMethod.replicate.B))); FUNCTION.meanRhoFromCoefs(rho, ns=n.miRs, nr=nr) }, by=comparison.simple]
setnames(cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.details, "comparison.detail", "comparison")
setnames(cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.simple, "comparison.simple", "comparison")
cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.all <- unique(rbind(cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.details, cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.simple), by="comparison")

this.outfile <- paste0(outdirs["tables"], "/Spearman_correlation_equimolar_Pool_summaries_noTechRep_comparison.txt")
write.table(cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.all, this.outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

comparison.levels.order <- unique(as.character(melt(outer(lib.detail.levels, lib.detail.levels, FUN = paste, sep=".VS."))$value))
comparison.levels.order <- comparison.levels.order[comparison.levels.order%in%cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.simple$comparison]
cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.all[, comparison:=factor(comparison, levels=comparison.levels.order)]
setorder(cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.all, comparison)
# For in-text Eq correlation summary numbers----
paste(
  cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.all[comparison%in%c("TruSeq.VS.TruSeq", "NEBNext.VS.NEBNext", "4N_B.VS.4N_B"),
  paste0(
      sub(".VS.*", "", comparison),
      ": ", 
      round(simple.mean.r, digits = 2),
      " (", 
      round(r02, digits = 2), 
      ", ", 
      round(r98, digits = 2),
      "; n=",
      nlibs,  
      ")")], 
  collapse="; ")

# inter-protocol comparisons
setkey(cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.all, comparison)
gsub(".VS.", " vs ", paste(cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.all[Cs(TruSeq.VS.NEBNext, TruSeq.VS.CleanTag, TruSeq.VS.4N_B, CleanTag.VS.NEBNext, NEBNext.VS.4N_B, CleanTag.VS.4N_B), paste0(comparison, ": ", round(simple.mean.r, digits = 2), " (", round(r02, digits = 2), ", ", round(r98, digits = 2), "; n=", nlibs,  ")")], collapse="; "))



# meanRho Summary Equimolar: Average Replicates First ----
#n.miRs <- dim(cpm.rle)[1]
#repl.groups <- sub(".[0-9]$", "", colnames(cpm.rle))
#cpm.rle.mean <- FUNCTION.aveTechReps(cpm.rle, ID = repl.groups)
#cor.cpm.rle.mean <- cor(cpm.rle.mean, use = "pair", method = "spearman")
#cor.cpm.rle.mean.dt <- data.table(cbind(melt(cor.cpm.rle.mean), melt(upper.tri(cor.cpm.rle.mean, diag=TRUE))))
#setnames(cor.cpm.rle.mean.dt, c("lab.libMethod.A", "lab.libMethod.B", "rho", "row", "col", "upper.tri"))
#cor.cpm.rle.mean.dt[, `:=`(id.A=ifelse(upper.tri==TRUE, as.character(lab.libMethod.A), as.character(lab.libMethod.B)), id.B=ifelse(upper.tri==TRUE, as.character(lab.libMethod.B), as.character(lab.libMethod.A)))]
#cor.cpm.rle.mean.dt[, `:=`(lab.libMethod.A=sub("^X", "", id.A), lab.libMethod.B=sub("^X", "", id.B))][, `:=`(id.A=NULL, id.B=NULL)]
#cor.cpm.rle.mean.dt[, c("lib.method.detail.A", "Lab.A"):=tstrsplit(lab.libMethod.A, split=".", fixed=TRUE)]
#cor.cpm.rle.mean.dt[, c("lib.method.detail.B", "Lab.B"):=tstrsplit(lab.libMethod.B, split=".", fixed=TRUE)]

#comparison.levels <- cor.cpm.rle.mean.dt[, unique(lib.method.detail.A)]
#cor.cpm.rle.mean.dt[, comparison.detail:=apply(.SD, 1, reorder.paste.fun, lookup.vec=comparison.levels), .SDcols=c("lib.method.detail.A", "lib.method.detail.B")]
#setkeyv(cor.cpm.rle.mean.dt, c("lab.libMethod.A", "lab.libMethod.B"))
#cor.cpm.rle.mean.dt.unique <- unique(cor.cpm.rle.mean.dt, by=key(cor.cpm.rle.mean.dt))

#cor.cpm.rle.mean.dt.summary.details <- cor.cpm.rle.mean.dt.unique[, { nr=length(unique(c(lab.libMethod.A, lab.libMethod.B))); FUNCTION.meanRhoFromCoefs(rho, ns=n.miRs, nr=nr) }, by=comparison.detail]

#this.outfile <- paste0(outdirs["tables"], "/Spearman_correlation_equimolar_Pool_summaries_averageTechReplicatesFirst.txt")
#write.table(cor.cpm.rle.mean.dt.summary.details, this.outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)



cor.matrix.all.eq.sample.info.cor.summary.details <- cor.matrix.all.eq.sample.info.unique[, { nr=length(unique(c(lab.libMethod.replicate.A, lab.libMethod.replicate.B))); FUNCTION.meanRhoFromCoefs(rho, ns=n.miRs, nr=nr) }, by=comparison.detail]
cor.matrix.all.eq.sample.info.cor.summary.simple <- cor.matrix.all.eq.sample.info.unique[, { nr=length(unique(c(lab.libMethod.replicate.A, lab.libMethod.replicate.B))); FUNCTION.meanRhoFromCoefs(rho, ns=n.miRs, nr=nr) }, by=comparison.simple]
setnames(cor.matrix.all.eq.sample.info.cor.summary.details, "comparison.detail", "comparison")
setnames(cor.matrix.all.eq.sample.info.cor.summary.simple, "comparison.simple", "comparison")
cor.matrix.all.eq.sample.info.cor.summary.all <- rbind(cor.matrix.all.eq.sample.info.cor.summary.simple, cor.matrix.all.eq.sample.info.cor.summary.details)

this.outfile <- paste0(outdirs["tables"], "/Spearman_correlation_equimolar_Pool_summaries.txt")
write.table(cor.matrix.all.eq.sample.info.cor.summary.all, this.outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


# FIG4A--%CV -------------------------------------------------------------------
#MEAN.CV.WI.LABS <- FINAL.EQUIMOLAR.LONG[, .(mean.pseudo.cpm=mean(pseudo.cpm.filt.lengths), sd.pseudo.cpm=sd(pseudo.cpm.filt.lengths), n=.N), by=.(lab.libMethod, lib.method.simple, lib.method.detail, equimolar.seqID)]
#MEAN.CV.WI.LABS[, cv.intra:=sd.pseudo.cpm/mean.pseudo.cpm]
MEAN.CV.WI.LABS <- FINAL.EQUIMOLAR.LONG[, .(mean.cpm=mean(pseudo.cpm.filt.lengths), sd.cpm=sd(pseudo.cpm.filt.lengths), q1.cpm=quantile(pseudo.cpm.filt.lengths, 0.25), q3.cpm=quantile(pseudo.cpm.filt.lengths, 0.75), n=.N), by=.(lab.libMethod, lib.method.simple, lib.method.detail, equimolar.seqID)]
MEAN.CV.WI.LABS[, `:=`(intralab.cv=100*sd.cpm/mean.cpm, intralab.qcd=((q3.cpm-q1.cpm)/2)/((q3.cpm+q1.cpm)/2))]


cv.range <- c(min(MEAN.CV.WI.LABS[lib.method.detail!="4N_NEXTflex"]$intralab.cv), max(MEAN.CV.WI.LABS[lib.method.detail!="4N_NEXTflex"]$intralab.cv)); cv.range
this.outfile <- paste0(outdirs["FIG4"], "/FIG4A_EQUIMOLAR_PERCENT_CV_INTRALAB.pdf")
ggplot(MEAN.CV.WI.LABS[lib.method.detail!="4N_NEXTflex"],
       aes(x = lab.libMethod,
           y = intralab.cv,
           fill=lib.method.simple))  +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_y_continuous(limits = c(0,180), breaks = seq(0, 200, 20), expand=c(0,0)) +
  theme(
    axis.text.x = element_text(angle=50, hjust=1),
    panel.grid.major = element_line(color = NA),
    panel.grid.minor = element_line(color = NA),
    panel.border = element_rect(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  ) +
  labs(#title = "EQUIMOLAR: %CV Within-labs",
    #x = "Library Prep Method",
    x=NULL,
    y = "%CV") 
ggsave(this.outfile, height=2.8, width=4.7)


# FIG4A--QCD --------------------------------------------------------------
#Intralab cv, violin
this.outfile <- paste0(outdirs["FIG4"], "/FIG4A_EQUIMOLAR_QCD_INTRALAB.pdf")
qcd.range <- c(min(MEAN.CV.WI.LABS[lib.method.detail!="4N_NEXTflex"]$intralab.qcd), max(MEAN.CV.WI.LABS[lib.method.detail!="4N_NEXTflex"]$intralab.qcd)); qcd.range
ggplot(MEAN.CV.WI.LABS[lib.method.detail!="4N_NEXTflex"],
       aes(x = lab.libMethod,
           y = intralab.qcd,
           fill=lib.method.simple))  +
  geom_violin(
    draw_quantiles = c(0.25, 0.5, 0.75)
  ) + scale_y_continuous(limits = c(0,0.8), breaks = seq(0, 1, 0.1), expand=c(0.01,0)) +
  theme(
    axis.text.x = element_text(angle=50, hjust=1),
    panel.grid.major = element_line(color = NA),
    panel.grid.minor = element_line(color = NA),
    panel.border = element_rect(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  ) +
  labs(#title = "EQUIMOLAR: QCD Within-labs",
    #x = "Library Prep Method",
    x=NULL,
    y = "QCD") 
ggsave(this.outfile, height=2.8, width=4.7)


# Intralab CV Equimolar Summary Table -------------------------------
summary.quantiles <- c(0.02, 0.25, 0.50, 0.75, 0.98)
MEAN.CV.WI.LABS.copy <- copy(MEAN.CV.WI.LABS)

MEAN.CV.WI.LABS.copy.4N <- subset(MEAN.CV.WI.LABS.copy, grepl("4N_[ABCD]", lib.method.detail))
MEAN.CV.WI.LABS.copy.4N[, lib.method.group:="4N"]
MEAN.CV.WI.LABS.copy.NEBSubset <- subset(MEAN.CV.WI.LABS.copy, lab.libMethod%in%NEB.subset)
MEAN.CV.WI.LABS.copy.NEBSubset[, lib.method.group:="NEBNext_subset"]
MEAN.CV.WI.LABS.copy[, lib.method.group:=lib.method.detail]

MEAN.CV.WI.LABS.subgroups <- rbind(MEAN.CV.WI.LABS.copy.4N, MEAN.CV.WI.LABS.copy.NEBSubset, MEAN.CV.WI.LABS.copy)
MEAN.CV.WI.LABS.subgroups[, (sub("^0.", "intralab.cv", summary.quantiles)):=lapply(summary.quantiles, USE.NAMES=FALSE, SIMPLIFY=FALSE, quantile, x=intralab.cv), by=.(lib.method.group)]
MEAN.CV.WI.LABS.subgroups[, (sub("^0.", "intralab.qcd", summary.quantiles)):=lapply(summary.quantiles, USE.NAMES=FALSE, SIMPLIFY=FALSE, quantile, x=intralab.qcd), by=.(lib.method.group)]
INTRALAB.CV.EQUIMOLAR.SUMMARY.tmp <- subset(unique(MEAN.CV.WI.LABS.subgroups[lib.method.detail!="4N_NEXTflex"], by=c("lab.libMethod", "lib.method.group")), select=c("lib.method.group", "lab.libMethod", "n", sub("^0.", "intralab.cv", summary.quantiles), sub("^0.", "intralab.qcd", summary.quantiles)))
INTRALAB.CV.EQUIMOLAR.SUMMARY.tmp[, n:=sum(n), by=lib.method.group]  
INTRALAB.CV.EQUIMOLAR.SUMMARY <- subset(unique(INTRALAB.CV.EQUIMOLAR.SUMMARY.tmp, by="lib.method.group"), select=-lab.libMethod)
this.outfile <- paste0(outdirs["tables"], "/Equimolar_intralab_CV_QCD_summary.txt")
write.table(INTRALAB.CV.EQUIMOLAR.SUMMARY, this.outfile, row.names = FALSE, col.names=TRUE, sep="\t", quote=FALSE)

# For in-text equimolar CV/QCD summary numbers----
setorder(INTRALAB.CV.EQUIMOLAR.SUMMARY, intralab.cv5)
paste(
  INTRALAB.CV.EQUIMOLAR.SUMMARY[!grepl("4N_[ACDX]", lib.method.group), 
                                                   paste0(lib.method.group,
                                                     ": ", 
                                                     round(intralab.cv5, digits = 2),
                                                     " (", 
                                                     round(intralab.cv02, digits = 2), 
                                                     ", ", 
                                                     round(intralab.cv98, digits = 2),
                                                     "; n=",
                                                     n,  
                                                     ")")], 
  collapse="; ")
# TruSeq: 6.18 (1.04, 43.36; n=32); 4N_B: 7.74 (1.59, 24.64; n=16); 4N: 9.29 (1.74, 35.83; n=28); NEBNext: 11.51 (1.36, 60.44; n=24); NEBNext_subset: 12.05 (1.23, 66; n=16); CleanTag: 23.7 (3.23, 61.29; n=4)
setorder(INTRALAB.CV.EQUIMOLAR.SUMMARY, intralab.qcd5)
paste(
  INTRALAB.CV.EQUIMOLAR.SUMMARY[!grepl("4N_[ACDX]", lib.method.group), 
                                paste0(lib.method.group,
                                       ": ", 
                                       round(intralab.qcd5, digits = 2),
                                       " (", 
                                       round(intralab.qcd02, digits = 2), 
                                       ", ", 
                                       round(intralab.qcd98, digits = 2),
                                       "; n=",
                                       n,  
                                       ")")], 
  collapse="; ")
# "TruSeq: 0.03 (0, 0.3; n=32); 4N_B: 0.04 (0.01, 0.16; n=16); 4N: 0.05 (0.01, 0.24; n=28); NEBNext: 0.05 (0.01, 0.32; n=24); NEBNext_subset: 0.06 (0.01, 0.36; n=16); CleanTag: 0.09 (0.01, 0.33; n=4)"


# Equimolar Cross-Lab %CV ----
# FIG4b--Equimolar Interlab %CV --------------------------------------------------------------
# FIG4b V1--Equimolar Interlab %CV 4N=All In-house --------------------------------------------------------------
this.outfile <- paste0(outdirs["FIG4"], "/FIG4B_EQUIMOLAR_PERCENT_CV_INTERLAB_TS_NEB_4N-Inhouse.pdf")
MEAN.CV.ACROSS.LABS <- MEAN.CV.WI.LABS[(lib.method.simple=="4N" | lib.method.simple=="TruSeq" | lib.method.simple == "NEBNext") & lib.method.detail != " 4N_NEXTflex" & lib.method.detail != "4N_Xu", .(cv.inter=sd(mean.cpm)/(sum(mean.cpm*n)/sum(n))), by=.(lib.method.simple, equimolar.seqID)]
ggplot(MEAN.CV.ACROSS.LABS,
       aes(x = lib.method.simple,
           y = cv.inter * 100))  +
  geom_boxplot(
    aes(fill = lib.method.simple),
    alpha = 0.25, 
    outlier.color = NA,
    width = 0.9
  ) + 
  scale_fill_manual(values = ann_colors$lib.method.simple) +
  geom_jitter(
    size = rel(0.001),
    alpha = 0.25,
    width = 0.4,
    height = 0
  ) + scale_y_continuous(limits = c(0, 350), breaks = seq(0, 350, 25)) +
  theme(
    panel.grid.major = element_line(color = NA),
    panel.grid.minor = element_line(color = NA),
    panel.border = element_rect(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  ) + 
  labs(#title = "EQUIMOLAR:\n%CV Across Labs",
    x=NULL,
    y = "%CV") 
ggsave(this.outfile, width = 2.0, height=2.25, units = "in")


# FIG4b V2--Equimolar Interlab %CV 4N_B --------------------------------------------------------------
MEAN.CV.ACROSS.LABS2 <- MEAN.CV.WI.LABS[(lib.method.detail=="4N_B" | lib.method.simple=="TruSeq" | lib.method.simple == "NEBNext"), .(cv.inter=sd(mean.cpm)/(sum(mean.cpm*n)/sum(n))), by=.(lib.method.simple, lib.method.detail, equimolar.seqID)]
this.outfile <- paste0(outdirs["FIG4"], "/FIG4B_EQUIMOLAR_PERCENT_CV_INTERLAB_TS_NEB_4N_B.pdf")
ggplot(MEAN.CV.ACROSS.LABS2,
       aes(x = lib.method.detail,
           y = cv.inter * 100))  +
  geom_boxplot(
    aes(fill = lib.method.simple),
    alpha = 0.25, 
    outlier.color = NA,
    width = 0.9
  ) + 
  scale_fill_manual(values = ann_colors$lib.method.simple) +
  geom_jitter(
    size = rel(0.001),
    alpha = 0.25,
    width = 0.4,
    height = 0
  ) + scale_y_continuous(limits = c(0, 275), breaks = seq(0, 350, 25)) +
  theme(
    panel.grid.major = element_line(color = NA),
    panel.grid.minor = element_line(color = NA),
    panel.border = element_rect(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  ) + 
  labs(#title = "EQUIMOLAR:\n%CV Across Labs",
    x=NULL,
    y = "%CV") 
ggsave(this.outfile, width = 2.0, height=2.25, units = "in")

# test NEB subsets
# FIG4b V3--Equimolar Interlab %CV 4N_B. NEB with labs that did 1:6 dilutions--------------------------------------------------------------
MEAN.CV.WI.LABS.drop.NEB.1to2.dilutions <- subset(MEAN.CV.WI.LABS, lab.libMethod!="NEBNext.Lab9" & lab.libMethod!="NEBNext.Lab3")
MEAN.CV.ACROSS.LABS3 <- MEAN.CV.WI.LABS.drop.NEB.1to2.dilutions[(lib.method.detail=="4N_B" | lib.method.simple=="TruSeq" | lib.method.simple == "NEBNext"), .(cv.inter=sd(mean.cpm)/(sum(mean.cpm*n)/sum(n))), by=.(lib.method.simple, lib.method.detail, equimolar.seqID)]
this.outfile <- paste0(outdirs["FIG4"], "/FIG4B_EQUIMOLAR_PERCENT_CV_INTERLAB_TS_NEB_4N_B_NEB1to2DilutionsOnly.pdf")
ggplot(MEAN.CV.ACROSS.LABS3,
       aes(x = lib.method.detail,
           y = cv.inter * 100))  +
  geom_boxplot(
    aes(fill = lib.method.simple),
    alpha = 0.25, 
    outlier.color = NA,
    width = 0.9
  ) + 
  scale_fill_manual(values = ann_colors$lib.method.simple) +
  geom_jitter(
    size = rel(0.001),
    alpha = 0.25,
    width = 0.4,
    height = 0
  ) + scale_y_continuous(limits = c(0, 275), breaks = seq(0, 350, 25)) +
  theme(
    panel.grid.major = element_line(color = NA),
    panel.grid.minor = element_line(color = NA),
    panel.border = element_rect(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  ) + 
  labs(#title = "EQUIMOLAR:\n%CV Across Labs",
    x=NULL,
    y = "%CV") 
ggsave(this.outfile, width = 2.0, height=2.25, units = "in")


#QUARTILE COEFFICIENT OF DISPERSION
#q3 <- quantile(x, 0.75)
#q1 <- quantile(x, 0.25)
#iqr <- (q3-q1)/2
#mid.hinge <- (q1+q3)/2
#iqr/mid.hinge

QCD.ACROSS.LABS <- MEAN.CV.WI.LABS[(lib.method.simple=="4N" | lib.method.simple=="TruSeq" | lib.method.simple == "NEBNext") & lib.method.detail != " 4N_NEXTflex" & lib.method.detail != "4N_Xu", .(q1=quantile(mean.cpm, 0.25), q3=quantile(mean.cpm, 0.75)), by=.(lib.method.simple, equimolar.seqID)]
QCD.ACROSS.LABS[, `:=`(iqr=(q3-q1)/2, midhinge=(q1+q3)/2)]
QCD.ACROSS.LABS[, qcd:=iqr/midhinge]

# FIG4b--QCD  4N=All In-house--------------------------------------------------------------
this.outfile <- paste0(outdirs["FIG4"], "/FIG4B_EQUIMOLAR_QCD_INTERLAB_TS_NEB_4N-Inhouse.pdf")
qcd.range <- c(min(QCD.ACROSS.LABS$qcd), max(QCD.ACROSS.LABS$qcd))
ggplot(QCD.ACROSS.LABS,
       aes(x = lib.method.simple,
           y = qcd))  +
  geom_boxplot(
    aes(fill = lib.method.simple),
    alpha = 0.25, 
    outlier.color = NA,
    width = 0.9
  ) + 
  scale_fill_manual(values = ann_colors$lib.method.simple) +
  geom_jitter(
    size = rel(0.001),
    alpha = 0.25,
    width = 0.4,
    height = 0
  ) + scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.1)) +
  theme(
    panel.grid.major = element_line(color = NA),
    panel.grid.minor = element_line(color = NA),
    panel.border = element_rect(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  ) + 
  labs(#title = "EQUIMOLAR:\n%CV Across Labs",
    x=NULL,
    y = "QCD") 
ggsave(this.outfile, width = 2.0, height=2.25, units = "in")

# FIG4b--QCD  4N_B--------------------------------------------------------------
QCD.ACROSS.LABS2 <- MEAN.CV.WI.LABS[(lib.method.detail=="4N_B" | lib.method.simple=="TruSeq" | lib.method.simple == "NEBNext"), .(q1=quantile(mean.cpm, 0.25), q3=quantile(mean.cpm, 0.75)), by=.(lib.method.simple, lib.method.detail, equimolar.seqID)]
QCD.ACROSS.LABS2[, `:=`(iqr=(q3-q1)/2, midhinge=(q1+q3)/2)]
QCD.ACROSS.LABS2[, qcd:=iqr/midhinge]
this.outfile <- paste0(outdirs["FIG4"], "/FIG4B_EQUIMOLAR_QCD_INTERLAB_TS_NEB_4N_B.pdf")
qcd.range <- c(min(QCD.ACROSS.LABS2$qcd), max(QCD.ACROSS.LABS2$qcd))
ggplot(QCD.ACROSS.LABS2,
       aes(x = lib.method.detail,
           y = qcd))  +
  geom_boxplot(
    aes(fill = lib.method.simple),
    alpha = 0.25, 
    outlier.color = NA,
    width = 0.9
  ) + 
  scale_fill_manual(values = ann_colors$lib.method.simple) +
  geom_jitter(
    size = rel(0.001),
    alpha = 0.25,
    width = 0.4,
    height = 0
  ) + scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.1)) +
  theme(
    panel.grid.major = element_line(color = NA),
    panel.grid.minor = element_line(color = NA),
    panel.border = element_rect(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  ) + 
  labs(#title = "EQUIMOLAR:\n%CV Across Labs",
    x=NULL,
    y = "QCD") 
ggsave(this.outfile, width = 2.0, height=2.25, units = "in")

# test NEB subsets
# FIG4b V3--Equimolar Interlab %CV 4N_B. NEB with labs that did 1:6 dilutions--------------------------------------------------------------
QCD.ACROSS.LABS3 <- MEAN.CV.WI.LABS.drop.NEB.1to2.dilutions[(lib.method.detail=="4N_B" | lib.method.simple=="TruSeq" | lib.method.simple == "NEBNext"), .(q1=quantile(mean.cpm, 0.25), q3=quantile(mean.cpm, 0.75)), by=.(lib.method.simple, lib.method.detail, equimolar.seqID)]
QCD.ACROSS.LABS3[, `:=`(iqr=(q3-q1)/2, midhinge=(q1+q3)/2)]
QCD.ACROSS.LABS3[, qcd:=iqr/midhinge]
this.outfile <- paste0(outdirs["FIG4"], "/FIG4B_EQUIMOLAR_PERCENT_QCD_INTERLAB_TS_NEB_4N_B_NEB1to2DilutionsOnly.pdf")
ggplot(QCD.ACROSS.LABS3,
       aes(x = lib.method.detail,
           y = qcd))  +
  geom_boxplot(
    aes(fill = lib.method.simple),
    alpha = 0.25, 
    outlier.color = NA,
    width = 0.9
  ) + 
  scale_fill_manual(values = ann_colors$lib.method.simple) +
  geom_jitter(
    size = rel(0.001),
    alpha = 0.25,
    width = 0.4,
    height = 0
  ) + scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.1)) +
  theme(
    panel.grid.major = element_line(color = NA),
    panel.grid.minor = element_line(color = NA),
    panel.border = element_rect(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  ) + 
  labs(#title = "EQUIMOLAR:\n%CV Across Labs",
    x=NULL,
    y = "QCD") 
ggsave(this.outfile, width = 2.0, height=2.25, units = "in")


# INTERLAB CV/QCD Equimolar Summary Table -------------------------------
MEAN.CV.WI.LABS.subgroups[, n:=.N, by=.(lib.method.group, equimolar.seqID)]

MEAN.CV.ACROSS.LABS.subgroups <- MEAN.CV.WI.LABS.subgroups[n>2, .(mean.cpm=mean(mean.cpm), sd.cpm=sd(mean.cpm), q1.cpm=quantile(mean.cpm, 0.25), q3.cpm=quantile(mean.cpm, 0.75), n=.N), by=.(lib.method.group, equimolar.seqID)]
MEAN.CV.ACROSS.LABS.subgroups[, `:=`(interlab.cv=100*sd.cpm/mean.cpm, interlab.qcd=((q3.cpm-q1.cpm)/2)/((q3.cpm+q1.cpm)/2))]
MEAN.CV.ACROSS.LABS.subgroups[, (sub("^0.", "interlab.cv", summary.quantiles)):=lapply(summary.quantiles, USE.NAMES=FALSE, SIMPLIFY=FALSE, quantile, x=interlab.cv), by=.(lib.method.group)]
MEAN.CV.ACROSS.LABS.subgroups[, (sub("^0.", "interlab.qcd", summary.quantiles)):=lapply(summary.quantiles, USE.NAMES=FALSE, SIMPLIFY=FALSE, quantile, x=interlab.qcd), by=.(lib.method.group)]

INTERLAB.CV.EQUIMOLAR.SUMMARY <- subset(unique(MEAN.CV.ACROSS.LABS.subgroups, by=c("lib.method.group", "n")), select=c("lib.method.group", "n", sub("^0.", "interlab.cv", summary.quantiles), sub("^0.", "interlab.qcd", summary.quantiles)))
this.outfile <- paste0(outdirs["tables"], "/Equimolar_interlab_CV_QCD_summary.txt")
write.table(INTERLAB.CV.EQUIMOLAR.SUMMARY, this.outfile, row.names = FALSE, col.names=TRUE, sep="\t", quote=FALSE)

# For in-text INTER-lab equimolar CV/QCD summary numbers----
setorder(INTERLAB.CV.EQUIMOLAR.SUMMARY, interlab.cv5)
paste(
  INTERLAB.CV.EQUIMOLAR.SUMMARY[, 
                                paste0(lib.method.group,
                                       ": ", 
                                       round(interlab.cv5, digits = 2),
                                       " (", 
                                       round(interlab.cv02, digits = 2), 
                                       ", ", 
                                       round(interlab.cv98, digits = 2),
                                       "; n=",
                                       n,  
                                       ")")], 
  collapse="; ")
# 4N_B: 30.42 (6.07, 74.24; n=4); TruSeq: 35.09 (16.32, 142.72; n=8); NEBNext_subset: 35.26 (9.75, 69.53; n=4); 4N: 54.1 (20.84, 108.59; n=7); NEBNext: 83.05 (19.3, 177.42; n=6)

setorder(INTERLAB.CV.EQUIMOLAR.SUMMARY, interlab.qcd5)
paste(
  INTERLAB.CV.EQUIMOLAR.SUMMARY[, 
                                paste0(lib.method.group,
                                       ": ", 
                                       round(interlab.qcd5, digits = 2),
                                       " (", 
                                       round(interlab.qcd02, digits = 2), 
                                       ", ", 
                                       round(interlab.qcd98, digits = 2),
                                       "; n=",
                                       n,  
                                       ")")], 
  collapse="; ")
# 4"4N_B: 0.13 (0.03, 0.42; n=4); NEBNext_subset: 0.18 (0.04, 0.46; n=4); TruSeq: 0.18 (0.06, 0.45; n=8); NEBNext: 0.25 (0.06, 0.56; n=6); 4N: 0.3 (0.07, 0.82; n=7)"

# missing miRNAs
total.eq.seqs <- length(unique(FINAL.EQUIMOLAR.LONG$equimolar.seqID))
FINAL.EQUIMOLAR.LONG[, `:=`(n.replicates.detected=sum(ifelse(count>0, 1, 0)), total.replicates=.N, mean.lib.size=exp(mean(log( count.total.by.sample.filt.lengths)))), by=.(lab.libMethod, Lab, lib.method.detail, lib.method.simple, equimolar.seqID, sequence, seq.len, gc.perc)]
FINAL.EQUIMOLAR.LONG[, mean.lib.size:=ifelse(n.replicates.detected==total.replicates, mean.lib.size, 0)]
FINAL.EQUIMOLAR.LONG[, mean.lib.size:=max(mean.lib.size), by=.(lab.libMethod)]
EQUIMOLAR.MISSING.MIR.SUMMARY.BY.MIR <- FINAL.EQUIMOLAR.LONG[, .(mean.lib.size=max(mean.lib.size)), by=.(lab.libMethod, Lab, lib.method.detail, lib.method.simple, equimolar.seqID, sequence, seq.len, gc.perc, n.replicates.detected, total.replicates)]
EQUIMOLAR.MISSING.MIR.SUMMARY <- EQUIMOLAR.MISSING.MIR.SUMMARY.BY.MIR[, .(n.miRs=.N), by=.(lab.libMethod, Lab, lib.method.detail, lib.method.simple, mean.lib.size, n.replicates.detected, total.replicates)]
EQUIMOLAR.MISSING.MIR.SUMMARY.drop.nonMissing <- EQUIMOLAR.MISSING.MIR.SUMMARY[n.miRs==total.eq.seqs | n.replicates.detected<total.replicates]
EQUIMOLAR.MISSING.MIR.SUMMARY.total.missing.in.any <- EQUIMOLAR.MISSING.MIR.SUMMARY[, .(n.miRs.missing=sum(ifelse(n.replicates.detected<total.replicates, n.miRs, 0))), by=.(lab.libMethod, Lab, lib.method.detail, lib.method.simple, mean.lib.size, total.replicates)]

# SUP Table 4: Number of undetected sequences in the equimolar pool -----
this.outfile <- paste0(outdirs["tables"], "/TABLES4_EQUIMOLAR_POOL_MISSING_MIR_COUNTS.xlsx")
write.xlsx(EQUIMOLAR.MISSING.MIR.SUMMARY.total.missing.in.any, this.outfile, row.names=FALSE, col.names=TRUE)


# Top and bottom miRs by protocol
# Drop the low-coverage truseq and 4N protocols
dge.eq.filt <- dge.eq[rowSums(dge.eq$counts)>0,]
dge.eq.filt.lrt.norm <- calcNormFactors(dge.eq.filt, "RLE")
cpm.lrt <- cpm(dge.eq.filt.lrt.norm, normalized.lib.sizes = TRUE, prior.count = 1)
drop.eq.samples <- row.names(dge.eq$samples[dge.eq$samples$lib.method.detail%in%c("4N_NEXTflex", "4N_Xu") | dge.eq$samples$lab.libMethod=="TruSeq.Lab8" ,])
cpm.lrt.filt <- cpm.lrt[, !colnames(cpm.lrt)%in%drop.eq.samples]
cpm.lrt.eq.sample.info <- dge.eq$samples[!row.names(dge.eq$samples)%in%drop.eq.samples,]
rank.matrix.all <- apply(-1*cpm.lrt.filt, 2, rank, ties.method = "min")
colnames(rank.matrix.all) <- sub("^X", "", colnames(rank.matrix.all))
in.top.ten <- ifelse(rank.matrix.all<=10, 1, 0)

n.samples <- ifelse(rank.matrix.all<=10, 1, 1)
in.top.ten.agreement <- sumTechReps(in.top.ten, ID = sample.info.eq[colnames(in.top.ten),]$lib.method.simple)
n.samples.summary <- sumTechReps(n.samples, ID = sample.info.eq[colnames(n.samples),]$lib.method.simple)
in.top.ten.agreement.perc <- in.top.ten.agreement/n.samples.summary
mirs.with.75.perc.agreement.in.top.10 <- names(which(apply(in.top.ten.agreement.perc, 1, max)>=0.75))

FINAL.EQUIMOLAR.LONG.25nt.disp.replicates.top.10 <- subset(FINAL.EQUIMOLAR.LONG.AVERAGE.REPS.FILT,  equimolar.seqID%in%mirs.with.75.perc.agreement.in.top.10)
FINAL.EQUIMOLAR.LONG.25nt.disp.replicates.top.10.mean <- FINAL.EQUIMOLAR.LONG.25nt.disp.replicates.top.10[, .(lib.meanCPM=mean(mean.pseudo.cpm)), by=.(lib.method.simple, equimolar.seqID)]

# now, instead, divide by the expected value to get a difference from expected
FINAL.EQUIMOLAR.LONG.25nt.disp.replicates.top.10.mean[, `:=`(meanCPM.from.expected=(lib.meanCPM-expected.cpm)/expected.cpm, log.CPM.from.expected=log2(lib.meanCPM/expected.cpm))]


# FIG S1a: Over-represented miRs--------------------------------------------------------------
# log CPM from expected
FINAL.EQUIMOLAR.LONG.25nt.disp.replicates.top.10.mean[, `:=`(meanCPM.from.expected=(lib.meanCPM-expected.cpm)/expected.cpm, log.CPM.from.expected=log2(lib.meanCPM/expected.cpm))]
this.outfile <- paste0(outdirs["supplemental_figures"], "/SUPFIG1_MANUSCRIPT_Top_equimolar_miRs_expressed_by_protocol_log_fold_diff_from_expected.pdf")
fold.min.max <- c(min(FINAL.EQUIMOLAR.LONG.25nt.disp.replicates.top.10.mean$meanCPM.from.expected), max(FINAL.EQUIMOLAR.LONG.25nt.disp.replicates.top.10.mean$meanCPM.from.expected))
g <- ggplot(FINAL.EQUIMOLAR.LONG.25nt.disp.replicates.top.10.mean, 
            aes(x=sub("MIR", "miR", equimolar.seqID), y=log.CPM.from.expected, fill= lib.method.simple)) + 
  geom_bar(stat="identity", color="black", pos="dodge", width=0.6) +
  scale_y_continuous(breaks = seq(-5, 10, 0.5)) +
  theme_bw() + 
  theme(
    axis.text=element_text(size=8, color="black"),
    axis.ticks.x=element_line(color=NA),
    axis.ticks.y=element_line(color="black", size=1.2),
    axis.line.x=element_line(color="black", size=1.2),
    axis.line.y=element_line(color="black", size=1.2),
    axis.text.x=element_text(vjust=1, angle=50, hjust=1),
    axis.title=element_text(size=8),
    plot.title=element_text(size=8),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=8),
    legend.title=element_text(size=8),
    legend.key.size=unit(1.2, "lines")
  ) + labs(
    x=NULL,
    #title="Equimolar Pool:\nOver-represented Sequences",
    y="mean CPM observed / expected (log2)",
    fill="Library Prep Method"
  ) + geom_hline(yintercept = 0, size=1.5)
ggsave(plot = g, filename = this.outfile, width=5, height=4, units="in")


# Bottom miRs by protocol
rank.matrix.all.down <- apply(cpm.lrt.filt, 2, rank, ties.method = "min")
in.bottom.ten <- ifelse(rank.matrix.all.down<=10, 1, 0)
n.samples <- ifelse(rank.matrix.all.down<=10, 1, 1)
in.bottom.ten.agreement <- sumTechReps(in.bottom.ten, ID = sample.info.eq[colnames(in.top.ten),]$lib.method.simple)
n.samples.summary <- sumTechReps(n.samples, ID = sample.info.eq[colnames(in.top.ten),]$lib.method.simple)
in.bottom.ten.agreement.perc <- in.bottom.ten.agreement/n.samples.summary
mirs.with.75.perc.agreement.in.bottom.10 <- names(which(apply(in.bottom.ten.agreement.perc, 1, max)>=0.75))

FINAL.EQUIMOLAR.LONG.25nt.disp.replicates.bottom.10 <- subset(FINAL.EQUIMOLAR.LONG.AVERAGE.REPS.FILT,  equimolar.seqID%in%mirs.with.75.perc.agreement.in.bottom.10)
FINAL.EQUIMOLAR.LONG.25nt.disp.replicates.bottom.10.mean <- FINAL.EQUIMOLAR.LONG.25nt.disp.replicates.bottom.10[, .(lib.meanCPM=mean(mean.pseudo.cpm)), by=.(lib.method.simple, equimolar.seqID)]
max.cpm.bottom <- max(FINAL.EQUIMOLAR.LONG.25nt.disp.replicates.bottom.10.mean$lib.meanCPM)
FINAL.EQUIMOLAR.LONG.25nt.disp.replicates.bottom.10.mean[, `:=`(percent.of.expected.value=(lib.meanCPM)/expected.cpm, log.CPM.from.expected=log2(lib.meanCPM/expected.cpm))]

# FIG S1b: Under-represented miRs--------------------------------------------------------------
this.outfile <- paste0(outdirs["supplemental_figures"], "/SUP1B_Bottom_equimolar_miRs_expressed_by_protocol_log_fold_diff_from_expected.pdf")
  g <- ggplot(FINAL.EQUIMOLAR.LONG.25nt.disp.replicates.bottom.10.mean, 
              aes(x=sub("MIR", "miR", equimolar.seqID), y=log.CPM.from.expected, fill= lib.method.simple)) + 
    geom_bar(stat="identity", color="black", pos="dodge", width=0.6) +
    scale_y_continuous(breaks = seq(0,-12, -2)) +
    theme_bw() + 
    theme(
      axis.text=element_text(size=8, color="black"),
      axis.ticks.x=element_line(color=NA),
      axis.ticks.y=element_line(color="black"),
      axis.line.x=element_line(color="black"),
      axis.line.y=element_line(color="black"),
      axis.text.x=element_text(vjust=1, angle=50, hjust=1),
      axis.title=element_text(size=8),
      plot.title=element_text(size=8),
      panel.grid.major = element_line(color = NA),
      panel.grid.minor = element_line(color = NA),
      legend.position = "none"
    ) + labs(
      x=NULL,
      #title="Equimolar Pool:\nOver-represented Sequences",
      y="mean CPM observed / expected (log2)",
      fill="Library Prep Method"
    ) + geom_hline(yintercept = 0, size=1.5)
  g
  ggsave(plot = g, filename = this.outfile, width=6, height=4, units="in")
  
# RATIOMETRIC WRANGLING --------------------------

ratiometric.counts <- subset(counts.all.dt.add.zerocount, pool!="SynthEQ" & !is.na(ratio.seqID))
ratiometric.counts[, ratioGroup:=ifelse(ratio.A=="1x", ratio.B, ratio.A)]
ratio.missing.in.all <- ratiometric.counts[, sum(ifelse(count>0,1,0)), by=ratio.seqID][V1==0]$ratio.seqID

# Skipping a few steps, since it was taken care of up-front. Name is kept to simplify downstream 
ratiometric.counts.addinfo.sizeFilt <- subset(ratiometric.counts, !ratio.seqID%in%ratio.missing.in.all)

ratiometric.counts.addinfo.sizeFilt.cast <- dcast.data.table(ratiometric.counts.addinfo.sizeFilt, ratio.seqID+ratio.A+ratio.B+ratioGroup+sequence+seq.len+gc.perc~lab.libMethod.pool.replicate, fun.aggregate="sum", fill=0, value.var="count")

FINAL.RATIOMETRIC.LONG <- ratiometric.counts.addinfo.sizeFilt
FINAL.RATIOMETRIC.LONG[, count.plus1:=count+1]
FINAL.RATIOMETRIC.LONG[, `:=`(count.total.by.sample.filt.lengths=sum(count),
                              count.plus1.total.by.sample.all.lengths=sum(count.plus1)
), by=.(lab.libMethod.pool.replicate, lab.libMethod.pool, lab.libMethod.replicate, lab.libMethod, Lab, lib.method.detail, lib.method.simple, replicate)]
FINAL.RATIOMETRIC.LONG[, `:=`(cpm.filt.lengths=(count*10^6)/count.total.by.sample.filt.lengths,
                              pseudo.cpm.filt.lengths=(count.plus1*10^6)/count.plus1.total.by.sample.all.lengths)]

# data frame / matrix for downstream diff expr and heatmaps
# also row (gene) and column (sample) annotation data frames
sample.info.ratio <- data.frame(FINAL.RATIOMETRIC.LONG[, .(tot.count=max(count.total.by.sample.filt.lengths), tot.count.plus1=max(count.plus1.total.by.sample.all.lengths)), by=.(lab.libMethod.pool.replicate, lab.libMethod.pool, lab.libMethod.replicate, lab.libMethod, Lab, lib.method.detail, lib.method.simple, pool, replicate)], row.names=1)
ratiometric.counts.addinfo.sizeFilt.cast.df <- data.frame(subset(ratiometric.counts.addinfo.sizeFilt.cast, select=c("ratio.seqID", row.names(sample.info.ratio))), row.names=1)
gene.info.ratio <- data.frame(subset(ratiometric.counts.addinfo.sizeFilt.cast, select=colnames(ratiometric.counts.addinfo.sizeFilt.cast)[!colnames(ratiometric.counts.addinfo.sizeFilt.cast)%in%row.names(sample.info.ratio)]), row.names=1)

FINAL.RATIOMETRIC.LONG[, `:=`(lab.libMethod=factor(lab.libMethod, levels=new.order$ratiometric),
                              lib.method.simple=factor(lib.method.simple, levels=lib.method.simple.ordered),
                              lib.method.detail=factor(lib.method.detail, levels=lib.method.detail.ordered))]

sample.info.ratio$log.lib.size <- log2(sample.info.ratio$tot.count)
sample.info.ratio$log.plus1.lib.size <- log2(sample.info.ratio$tot.count.plus1)

sample.info.ratio$lab.libMethod <- factor(sample.info.ratio$lab.libMethod, levels=new.order$ratiometric)
sample.info.ratio$lib.method.simple <- factor(sample.info.ratio$lib.method.simple, levels=lib.method.simple.ordered)
sample.info.ratio$lib.method.detail <- factor(sample.info.ratio$lib.method.detail, levels=lib.method.detail.ordered)


ratio.ns <- c(10,8,5,4,3,2,1.5)
ratio.all <- c(1/ratio.ns, 1, rev(ratio.ns))
ratio.order <- c("1x:10x", "1x:8x", "1x:5x", "1x:4x", "1x:3x", "1x:2x", "1x:1.5x", "1x:1x", "1.5x:1x", "2x:1x", "3x:1x", "4x:1x", "5x:1x", "8x:1x", "10x:1x")
ratio.order.fract <- c("1/10", "1/8", "1/5", "1/4", "1/3", "1/2", "1/1.5", "1", "1.5", "2", "3", "4", "5", "8", "10")
FINAL.RATIOMETRIC.LONG.AVERAGE.REPS <- FINAL.RATIOMETRIC.LONG[, .(mean.cpm.filt.lengths=mean(cpm.filt.lengths),
                                                                  mean.pseudo.cpm.filt.lengths=mean(pseudo.cpm.filt.lengths),
                                                                  mean.lib.size=mean(count.total.by.sample.filt.lengths),
                                                                  mean.lib.size.plus1=mean(count.plus1.total.by.sample.all.lengths)),
                                                              by=.(ratio.seqID, sequence, lab.libMethod.pool, lab.libMethod, Lab, lib.method.detail, lib.method.simple, pool, ratio.A, ratio.B, seq.len, gc.perc, ratioGroup)]
FINAL.RATIOMETRIC.LONG[, log.pseudo.cpm:=log2(pseudo.cpm.filt.lengths)]
FINAL.RATIOMETRIC.LONG.AVERAGE.REPS.CAST.AvsB <- dcast.data.table(FINAL.RATIOMETRIC.LONG, ratio.seqID+ratio.A+ratio.B+sequence+seq.len+gc.perc+ratioGroup+lab.libMethod+Lab+lib.method.detail+lib.method.simple~pool, value.var=c("cpm.filt.lengths"), fun.aggregate = mean, fill=NA)
setnames(FINAL.RATIOMETRIC.LONG.AVERAGE.REPS.CAST.AvsB, c("SynthA", "SynthB"), c("mean.cpm_SynthA", "mean.cpm_SynthB"))
FINAL.RATIOMETRIC.LONG.AVERAGE.REPS.CAST.AvsB[, ratio.AvsB.Measured:=ifelse(mean.cpm_SynthA==0 | mean.cpm_SynthB==0, NA, (mean.cpm_SynthA)/(mean.cpm_SynthB))]
FINAL.RATIOMETRIC.LONG.AVERAGE.REPS.CAST.AvsB[, RATIO.AvsB.Expected:=factor(paste(ratio.A, ratio.B, sep=":"), levels=ratio.order)]

# FIG3A -------------------------------------------------------------------
# Specify new order to fit into grid in a semi-resonable way
group.split <- split(new.order$ratiometric, strsplit2(new.order$ratiometric, split=".", fixed=TRUE)[,1])
new.order.ratio.2 <- c(group.split[["TruSeq"]][1:2],
                       group.split[["NEBNext"]][1:3],
                       group.split[["TruSeq"]][3:4],
                       group.split[["NEBNext"]][4:6],
                       group.split[["TruSeq"]][5:6],
                       group.split[["4N_B"]][1:2],
                       group.split[["4N_NEXTflex"]][1],
                       group.split[["TruSeq"]][7:8],
                       group.split[["4N_B"]][3:4],
                       group.split[["4N_NEXTflex"]][2],
                       group.split[["CleanTag"]][1],
                       group.split[["4N_A"]][1],
                       group.split[["4N_C"]][1],
                       group.split[["4N_D"]][1],
                       group.split[["4N_Xu"]][1])
new.order.ratio.2b  <- c(group.split[["TruSeq"]][1:2],
                       group.split[["NEBNext"]][1:3],
                       group.split[["TruSeq"]][3:4],
                       group.split[["NEBNext"]][4:6],
                       group.split[["TruSeq"]][5:6],
                       group.split[["4N_B"]][1:2],
                       group.split[["4N_NEXTflex"]][1],
                       group.split[["TruSeq"]][7:8],
                       group.split[["4N_B"]][3:4],
                       group.split[["4N_NEXTflex"]][2],
                       group.split[["CleanTag"]][1],
                       group.split[["4N_A"]][1],
                       group.split[["4N_C"]][1],
                       group.split[["4N_D"]][1],
                       group.split[["4N_Xu"]][1])

new.order.ratio.3 <- c(group.split[["TruSeq"]][1:4],
                       group.split[["CleanTag"]][1],
                       group.split[["TruSeq"]][5:8],
                       group.split[["4N_A"]][1],
                       group.split[["NEBNext"]][1:3],
                       group.split[["4N_NEXTflex"]][1],
                       group.split[["4N_C"]][1],
                       group.split[["NEBNext"]][4:6],
                       group.split[["4N_NEXTflex"]][2],
                       group.split[["4N_D"]][1],
                       group.split[["4N_B"]][1:4],
                       group.split[["4N_Xu"]][1])

FINAL.RATIOMETRIC.LONG.AVERAGE.REPS.CAST.AvsB[, lab.libMethod2:=factor(lab.libMethod, levels=new.order.ratio.2)]
ratio.labels <- ratio.order.fract
names(ratio.labels) <- ratio.order
undetected.mirs <- FINAL.RATIOMETRIC.LONG.AVERAGE.REPS.CAST.AvsB[, .(n.missing=sum(ifelse(is.na(log2(ratio.AvsB.Measured)), 1, 0))), by=.(lab.libMethod2)]
g <- ggplot(FINAL.RATIOMETRIC.LONG.AVERAGE.REPS.CAST.AvsB, 
              aes(
                x=RATIO.AvsB.Expected,
                y=ratio.AvsB.Measured,
                fill=lib.method.simple
              )) + scale_x_discrete(expand=c(0,1.2)) +
    geom_boxplot(show.legend=FALSE, outlier.size = 0.75, outlier.stroke=0, outlier.fill="black") + 
    theme(
      axis.text.x=element_text(angle=90, hjust=1),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.background=element_blank(),
      strip.text=element_text(size=8),
      legend.position = "none"
    ) + scale_y_continuous(trans=log2_trans(),
                           #limits=c(1/32, 32),
                           breaks = trans_breaks("log2", n = 9, function(x) 2^x),
                           labels = function(x){ fractions(x) }
    ) + labs(x=NULL, y=NULL) + 
    annotation_logticks(sides = "l", short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm")) + facet_wrap(~lab.libMethod2, dir = "h", ncol=5); g
this.plot <- paste0(outdirs["FIG3"], "/FIG3A_RatioVSExpected.pdf")
ggsave(this.plot, plot = g, width = 7.5, height=6.5, units = "in")


# Differential Expression Ratio Pool -------------------------------------------------------------------
dge.ratio <- DGEList(ratiometric.counts.addinfo.sizeFilt.cast.df, genes=gene.info.ratio, remove.zeros = TRUE)
dge.ratio$samples <- cbind(dge.ratio$samples, sample.info.ratio)

#dge.ratio$samples$lab.libMethod <- droplevels(dge.ratio$samples$lab.libMethod)
dge.ratio$samples$pool <- factor(dge.ratio$samples$pool, c("SynthA", "SynthB"))

ratio.groups <- unique(dge.ratio$genes$ratioGroup)
names(ratio.groups) <- ratio.groups
ratio.groups.list <- sapply(ratio.groups, USE.NAMES=FALSE, simplify = FALSE, FUN=function(x){
  if(x=="1x"){
    Up.in.A <- !(dge.ratio$genes$ratio.A == x & dge.ratio$genes$ratio.B == x)
    Up.in.B <- dge.ratio$genes$ratio.A == x & dge.ratio$genes$ratio.B == x
  } else{
    Up.in.A <- dge.ratio$genes$ratio.A == x
    Up.in.B <- dge.ratio$genes$ratio.B == x
  }
  
  return(list(Up.in.A=Up.in.A,
              up.in.B=Up.in.B))
})

fdr.cutoff=0.01
labs.all <- as.character(unique(dge.ratio$samples$lab.libMethod))

all.ratio.topTables <- do.call(rbind, sapply(labs.all, USE.NAMES=FALSE, simplify = FALSE, FUN=function(this.lab){
  dge.subset <- dge.ratio[, dge.ratio$samples$lab.libMethod==this.lab]  
  subset.design <- model.matrix(~pool, dge.subset$samples)
  colnames(subset.design) <- sub("pool", "", colnames(subset.design))
  dge.subset <- calcNormFactors(dge.subset, method="RLE")
  dge.subset <- estimateDisp(dge.subset, subset.design, robust = TRUE)  
  fit.subset <- glmFit(dge.subset, subset.design)
  lrt.subset <- glmLRT(fit.subset)  
  tt.subset <- data.table(data.frame(ratio.seqID=row.names(dge.subset), topTags(lrt.subset, n = Inf, sort.by = "none")))
  tt.subset[, lab.libMethod:=this.lab]
}))

all.ratio.topTables.limma <- do.call(rbind, sapply(labs.all, USE.NAMES=FALSE, simplify = FALSE, FUN=function(this.lab){
  dge.subset <- dge.ratio[, dge.ratio$samples$lab.libMethod==this.lab]
  subset.design <- model.matrix(~pool, dge.subset$samples)
  colnames(subset.design) <- sub("pool", "", colnames(subset.design))
  dge.subset <- calcNormFactors(dge.subset, method="RLE")
  
  #logCPM <- cpm(dge.subset, log=TRUE, prior.count=1)
  #fit <- lmFit(logCPM, subset.design)
  #fit <- eBayes(fit, trend=TRUE)
  v <- voom(dge.subset, subset.design)
  fit <- lmFit(v)
  fit <- eBayes(fit, trend = TRUE, robust=TRUE)
  tt.subset <- data.table(data.frame(ratio.seqID=row.names(dge.subset), topTable(fit, n = Inf, confint=0.95, sort.by = "none", genelist = dge.subset$genes)))
  tt.subset[, lab.libMethod:=this.lab]
  tt.subset[, `:=`(ratio.A.num=as.numeric(sub("x", "", ratio.A)),
                   ratio.B.num=as.numeric(sub("x", "", ratio.B)))]
  tt.subset[, logratio.expected:=log2(ratio.B.num/ratio.A.num)]
  # barcodeplot(tt.subset$logFC, ratio.groups.list[["1x"]][[1]], ratio.groups.list[["1x"]][[2]])
}))
all.ratio.topTables.deseq2 <- do.call(rbind, sapply(labs.all, USE.NAMES=FALSE, simplify = FALSE, FUN=function(this.lab){
  dge.subset <- dge.ratio[, dge.ratio$samples$lab.libMethod==this.lab]
  dds <- DESeqDataSetFromMatrix(dge.subset$counts, design = ~ pool, colData = data.frame(dge.subset$samples))
  mcols(dds) <- DataFrame(dge.subset$genes, check.names = TRUE, row.names = row.names(dge.subset$genes))
  dds <- DESeq(dds)
  res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)
  res$padj <- ifelse(is.na(res$padj), 1, res$padj)
  tt.subset <- data.table(data.frame(merge(data.frame(res), mcols(dds, use.names = TRUE), by=0)))
  setnames(tt.subset, "Row.names", "ratio.seqID")
  tt.subset[, lab.libMethod:=this.lab]
  tt.subset[, `:=`(ratio.A.num=as.numeric(sub("x", "", ratio.A)),
                   ratio.B.num=as.numeric(sub("x", "", ratio.B)))]
  tt.subset[, logratio.expected:=log2(ratio.B.num/ratio.A.num)]
  # barcodeplot(tt.subset$logFC, ratio.groups.list[["1x"]][[1]], ratio.groups.list[["1x"]][[2]])
}))

all.ratio.topTables.summary <- all.ratio.topTables[, .(n.signif=sum(ifelse(FDR<fdr.cutoff, 1, 0)), n.seqs=.N), by=.(ratioGroup, lab.libMethod)]
all.ratio.topTables.summary[, `:=`(percent.signif=n.signif/n.seqs)]
all.ratio.topTables.summary[, c("lib.method.detail", "Lab"):=tstrsplit(lab.libMethod, split=".", fixed=TRUE)]
all.ratio.topTables.summary[, lib.method.simple:=tstrsplit(lib.method.detail, split="_", fixed=TRUE)[1]]
all.ratio.topTables.summary[, `:=`(lab.libMethod=factor(lab.libMethod, levels=new.order$ratiometric),
                                   lib.method.detail=factor(lib.method.detail, levels=lib.detail.levels),
                                   lib.method.simple=factor(lib.method.simple, levels=lib.simple.levels),
                                   ratioGroup=factor(ratioGroup, levels=c("1x", "1.5x", "2x", "3x", "4x", "5x", "8x", "10x")))]
all.ratio.topTables.summary[, algorithm:="edgeR"]
all.ratio.topTables.limma.summary <- all.ratio.topTables.limma[, .(n.signif=sum(ifelse(adj.P.Val<fdr.cutoff, 1, 0)), n.seqs=.N), by=.(ratioGroup, lab.libMethod)]
all.ratio.topTables.limma.summary[, `:=`(percent.signif=n.signif/n.seqs)]
all.ratio.topTables.limma.summary[, c("lib.method.detail", "Lab"):=tstrsplit(lab.libMethod, split=".", fixed=TRUE)]
all.ratio.topTables.limma.summary[, lib.method.simple:=tstrsplit(lib.method.detail, split="_", fixed=TRUE)[1]]
all.ratio.topTables.limma.summary[, `:=`(lab.libMethod=factor(lab.libMethod, levels=new.order$ratiometric),
                                         lib.method.detail=factor(lib.method.detail, levels=lib.detail.levels),
                                         lib.method.simple=factor(lib.method.simple, levels=lib.simple.levels),
                                         ratioGroup=factor(ratioGroup, levels=c("1x", "1.5x", "2x", "3x", "4x", "5x", "8x", "10x")))]
all.ratio.topTables.limma.summary[, algorithm:="limma"]
all.ratio.topTables.deseq2.summary <- all.ratio.topTables.deseq2[, .(n.signif=sum(ifelse(padj<fdr.cutoff, 1, 0)), n.seqs=.N), by=.(ratioGroup, lab.libMethod)]
all.ratio.topTables.deseq2.summary[, `:=`(percent.signif=n.signif/n.seqs)]
all.ratio.topTables.deseq2.summary[, c("lib.method.detail", "Lab"):=tstrsplit(lab.libMethod, split=".", fixed=TRUE)]
all.ratio.topTables.deseq2.summary[, lib.method.simple:=tstrsplit(lib.method.detail, split="_", fixed=TRUE)[1]]
all.ratio.topTables.deseq2.summary[, `:=`(lab.libMethod=factor(lab.libMethod, levels=new.order$ratiometric),
                                          lib.method.detail=factor(lib.method.detail, levels=lib.detail.levels),
                                          lib.method.simple=factor(lib.method.simple, levels=lib.simple.levels),
                                          ratioGroup=factor(ratioGroup, levels=c("1x", "1.5x", "2x", "3x", "4x", "5x", "8x", "10x")))]
all.ratio.topTables.deseq2.summary[, algorithm:="deseq2"]
all.ratio.topTables.both <- rbind(all.ratio.topTables.limma.summary, all.ratio.topTables.summary, all.ratio.topTables.deseq2.summary)
all.ratio.topTables.both[, lab.libMethod2:=factor(lab.libMethod, levels=new.order.ratio.2)]
# FIG S4: Ratiometric Differential Expression Barplot -------------------------------------------------------------------

this.outfile <- paste0(outdirs["supplemental_figures"], "/SUPP_FIGS4_DIFFExprs_Results_RATIOPOOL_barplot.pdf")
ggplot(all.ratio.topTables.both, 
       aes(x=ratioGroup, y=percent.signif, fill=algorithm)) + 
  geom_bar(stat="identity", position="dodge", width=0.75) +
  facet_wrap(~lab.libMethod2, ncol=5, dir = "v")  + theme(
    axis.text.x=element_text(angle=90, hjust=1, color="black"),
    axis.text.y=element_text(color="black"),
    panel.grid.major=element_line(color=NA),
    panel.grid.minor=element_line(color=NA),
    axis.title=element_text(size=10),
    plot.title=element_text(size=10),
    panel.background=element_blank(),
    strip.text=element_text(size=8),
    legend.position = "top"
  ) + scale_y_continuous(labels=scales::percent_format(), expand=c(0,0), limits=c(0,1.05))
ggsave(this.outfile, height=6, width=8)


# Maria asked for a table to show the fold changes/cutoffs.
all.ratio.topTables.limma.voom.qual <- do.call(rbind, sapply(labs.all, USE.NAMES=FALSE, simplify = FALSE, FUN=function(this.lab){
  dge.subset <- dge.ratio[, dge.ratio$samples$lab.libMethod==this.lab]  
  subset.design <- model.matrix(~pool, dge.subset$samples)
  colnames(subset.design) <- unique(dge.subset$samples$pool)
  dge.subset <- calcNormFactors(dge.subset, method="RLE")
  v.subset <- voom(dge.subset, subset.design, normalize.method = "none")
  fit <- lmFit(v.subset, subset.design)
  fit <- eBayes(fit, trend=TRUE)
  tt.subset <- data.table(data.frame(ratio.seqID=row.names(dge.subset), topTable(fit, 2, n = Inf, confint=0.95, sort.by = "none", genelist = dge.subset$genes)))
  tt.subset[, lab.libMethod:=this.lab]
  tt.subset[, `:=`(ratio.A.num=as.numeric(sub("x", "", ratio.A)),
                   ratio.B.num=as.numeric(sub("x", "", ratio.B)))]
  tt.subset[, logratio.expected:=log2(ratio.B.num/ratio.A.num)]
  # barcodeplot(tt.subset$logFC, ratio.groups.list[["1x"]][[1]], ratio.groups.list[["1x"]][[2]])
}))
all.ratio.topTables.limma[, `:=`(logFC.vs.exp=logFC-logratio.expected, CI.L.vs.exp=CI.L-logratio.expected, CI.R.vs.exp=CI.R-logratio.expected)]
all.ratio.topTables.limma[, ratio.subpool:=factor(paste(ratio.A, ratio.B, sep=":"), levels=ratio.order)]
all.ratio.topTables.limma[, `:=`(expected.within.95.CI=ifelse(logratio.expected>=CI.L & logratio.expected<=CI.R, 1, 0), span.95.CI=CI.R-CI.L)]
this.outfile <- paste0(outdirs["tables"], "/RatioPool_limma_voom_calls_for_sup_table.txt")
write.table(all.ratio.topTables.limma, this.outfile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

# Table S5 Ratiometric Pool Observed VS expected ratios
all.ratio.topTables.limma[, abs.logFC.vs.exp:=abs(logFC.vs.exp)]
all.ratio.topTables.limma.summary <- all.ratio.topTables.limma[, .(n.WI.1.5x=sum(ifelse(abs.logFC.vs.exp<=log2(1.5), 1, 0)),
                                                                   n.WI.2.0x=sum(ifelse(abs.logFC.vs.exp<=log2(2.0), 1, 0)),
                                                                   n.WI.2.5x=sum(ifelse(abs.logFC.vs.exp<=log2(2.5), 1, 0)),
                                                                   n.tot=.N), by=.(lab.libMethod, ratio.A, ratio.B)]
all.ratio.topTables.limma.summary[, c("Protocol", "Lab"):=tstrsplit(lab.libMethod, split=".", fixed=TRUE)]
this.outfile <- paste0(outdirs["tables"], "/TABLE_S5_RatioPool_limma_voom_calls_obs_vs_exp.xlsx")
write.xlsx(all.ratio.topTables.limma.summary, this.outfile, row.names=FALSE, col.names=TRUE)

# COR HEATMAPS FOR MAN
#dge.ratio.filt <- dge.ratio[rowSums(dge.ratio$counts>0)>0,]

dge.ratio <- calcNormFactors(dge.ratio, method="RLE")
cpm.ratio.rle <- cpm(dge.ratio, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)


# Now average over replicates for CPM as well
cpm.ratio.sampleA <- cpm.ratio.rle[, grep("SynthA", colnames(cpm.ratio.rle), value = TRUE)]
sample.info.ratioA <- dge.ratio$samples[dge.ratio$samples$pool=="SynthA",]
cpm.ratio.sampleB <- cpm.ratio.rle[, grep("SynthB", colnames(cpm.ratio.rle), value = TRUE)]
sample.info.ratioB <- dge.ratio$samples[dge.ratio$samples$pool=="SynthB",]
sample.groups.A <- gsub(pattern = ".Synth[AB].[1-4]", replacement = "", x = colnames(cpm.ratio.sampleA))
sample.groups.B <- gsub(pattern = ".Synth[AB].[1-4]", replacement = "", x = colnames(cpm.ratio.sampleB))
n.reps.A <- as.numeric(table(sample.groups.A))
n.reps.B <- as.numeric(table(sample.groups.B))
sum.cpm.ratio.sampleA <- sumTechReps(cpm.ratio.sampleA, ID = sample.groups.A)
sum.cpm.ratio.sampleB <- sumTechReps(cpm.ratio.sampleB, ID = sample.groups.B)
n.reps.A <- matrix(as.numeric(table(sample.groups.A)), nrow = dim(sum.cpm.ratio.sampleA)[1], ncol=dim(sum.cpm.ratio.sampleA)[2], byrow = TRUE)
n.reps.B <- matrix(as.numeric(table(sample.groups.B)), nrow = dim(sum.cpm.ratio.sampleB)[1], ncol=dim(sum.cpm.ratio.sampleB)[2], byrow = TRUE)
ave.cpm.ratio.sampleA <- (sum.cpm.ratio.sampleA/n.reps.A)
ave.cpm.ratio.sampleB <- (sum.cpm.ratio.sampleB/n.reps.B)
colnames(ave.cpm.ratio.sampleA) <- sub("^X", "", colnames(ave.cpm.ratio.sampleA))
colnames(ave.cpm.ratio.sampleB) <- sub("^X", "", colnames(ave.cpm.ratio.sampleB))
ave.cpm.ratio.sampleAB.ratio <- ave.cpm.ratio.sampleA-ave.cpm.ratio.sampleB

if(!all.equal(colnames(sum.cpm.ratio.sampleB), colnames(sum.cpm.ratio.sampleA)) | !all.equal(row.names(sum.cpm.ratio.sampleB), row.names(sum.cpm.ratio.sampleA))){
  warning("NOT ALL RATIO MATRIX GROUPS ARE EQUAL")
}

full.spearman.ave.CPM.A <- cor(ave.cpm.ratio.sampleA, use = "pair", method = "spearman")^2
full.spearman.ave.CPM.B <- cor(ave.cpm.ratio.sampleB, use = "pair", method = "spearman")^2
full.spearman.ratio.sampleAB.ratio <- cor(ave.cpm.ratio.sampleAB.ratio, use = "pair", method = "spearman")^2


sample.info.logfc <- data.frame(row.names=colnames(full.spearman.ave.CPM.A), 
                                lib.method.simple=strsplit2(strsplit2(colnames(full.spearman.ave.CPM.A), split=".", fixed=TRUE)[,1], split="_")[,1],
                                lib.method.detail=strsplit2(colnames(full.spearman.ave.CPM.A), split=".", fixed=TRUE)[,1])
sample.info.logfc$lib.method.simple <- factor(sample.info.logfc$lib.method.simple, levels = lib.simple.levels)
sample.info.logfc$lib.method.detail <- factor(sample.info.logfc$lib.method.detail, levels = lib.detail.levels)

pheatmap(full.spearman.ave.CPM.A, fontsize = 8, annotation_row = sample.info.logfc, annotation_col=sample.info.logfc, breaks=seq(0,1, length.out = 101), annotation_colors=FUNCTION.filter.ann_colors(ann_colors, sample.info.logfc), main = "RATIOMETRIC Pool Correlation\n Log CPM Sample A Spearman rho^2", labels_row = sub("^X", "", row.names(sample.info.logfc)), labels_col = sub("^X", "", row.names(sample.info.logfc)) )
pheatmap(full.spearman.ave.CPM.B, fontsize = 8, annotation_row = sample.info.logfc, annotation_col=sample.info.logfc, annotation_colors=FUNCTION.filter.ann_colors(ann_colors, sample.info.logfc), breaks=seq(0,1, length.out = 101), main = "RATIOMETRIC Pool Correlation\n Log CPM Sample B Spearman rho^2", labels_row = sub("^X", "", row.names(sample.info.logfc)), labels_col = sub("^X", "", row.names(sample.info.logfc)) )
pheatmap(full.spearman.ratio.sampleAB.ratio, fontsize = 8, annotation_row = sample.info.logfc, annotation_col=sample.info.logfc,  annotation_colors=FUNCTION.filter.ann_colors(ann_colors, sample.info.logfc), main = "RATIOMETRIC Pool Correlation\n Log2 Ratio SampleB:A Spearman rho^2", labels_row = sub("^X", "", row.names(sample.info.logfc)), breaks=seq(0,1, length.out = 101), labels_col = sub("^X", "", row.names(sample.info.logfc)))

# FIG3B -------------------------------------------------------------------
# For manuscript
library(dendsort)
dA <- dist(full.spearman.ave.CPM.A)
dB <- dist(full.spearman.ave.CPM.B)
d.mean <- as.dist((as.matrix(dA)+as.matrix(dB))/2)
this.outfile <- paste0(outdirs["FIG3"], "/FIG3B_SPEARMAN_COR_HEATMAP_CPM.pdf")
pA <- pheatmap(full.spearman.ave.CPM.A, fontsize = 8,
         annotation_legend = FALSE,
         legend = FALSE, 
         annotation_row = sample.info.logfc,
         annotation_colors=ann_colors, 
         annotation_col=sample.info.logfc,
         breaks=seq(0,1, length.out = 101),
         treeheight_row=10,
         treeheight_col=10,
         width=3,
         height=3,
         clustering_distance_rows = d.mean,
         clustering_distance_cols = d.mean)$gtable

pB <-pheatmap(full.spearman.ave.CPM.B, fontsize = 8,
         annotation_legend = FALSE,
         legend = FALSE, 
         annotation_row = sample.info.logfc,
         annotation_colors=ann_colors, 
         annotation_col=sample.info.logfc,
         breaks=seq(0,1, length.out = 101),
         treeheight_row=10,
         treeheight_col=10,
         width=3,
         height=3,
         clustering_distance_rows = d.mean,
         clustering_distance_cols = d.mean)$gtable
pRatio <- pheatmap(full.spearman.ratio.sampleAB.ratio, fontsize = 8,
         annotation_legend = FALSE,
         legend = FALSE, 
         annotation_row = sample.info.logfc,
         annotation_colors=ann_colors, 
         annotation_col=sample.info.logfc,
         breaks=seq(0,1, length.out = 101),
         treeheight_row=10,
         treeheight_col=10,
         width=3,
         height=3,
         clustering_distance_rows = d.mean,
         clustering_distance_cols = d.mean)$gtable
ratio.cor.heatmap <- plot_grid(plotlist = list(pA, pB, pRatio), nrow=1)
save_plot(this.outfile, ratio.cor.heatmap, ncol = 3, nrow = 1, base_height = 4)

# FIG3B Alternative using unnormalized ratios-------------------------------------------------------------------
# For manuscript
library(dendsort)
ave.cpm.ratio.sampleA2 <- (data.frame(dcast.data.table(FINAL.RATIOMETRIC.LONG.AVERAGE.REPS.CAST.AvsB, ratio.seqID~lab.libMethod, value.var="mean.cpm_SynthA", fill=0), row.names=1))
ave.cpm.ratio.sampleB2 <- (data.frame(dcast.data.table(FINAL.RATIOMETRIC.LONG.AVERAGE.REPS.CAST.AvsB, ratio.seqID~lab.libMethod, value.var="mean.cpm_SynthB", fill=0), row.names=1))
full.spearman.ratio.sampleAB.ratio2 <- (data.frame(dcast.data.table(FINAL.RATIOMETRIC.LONG.AVERAGE.REPS.CAST.AvsB, ratio.seqID~lab.libMethod, value.var="ratio.AvsB.Measured", fill=0), row.names=1))
full.spearman.ave.CPM.A2 <- cor(ave.cpm.ratio.sampleA2, use = "pair", method = "spearman")^2
full.spearman.ave.CPM.B2 <- cor(ave.cpm.ratio.sampleB2, use = "pair", method = "spearman")^2
full.spearman.ratio.sampleAB.ratio2 <- cor(full.spearman.ratio.sampleAB.ratio2, use = "pair", method = "spearman")^2
dimnames(full.spearman.ave.CPM.A2) <- lapply(dimnames(full.spearman.ave.CPM.A2), sub, pattern = "^X", replacement =  "")
dimnames(full.spearman.ave.CPM.B2) <- lapply(dimnames(full.spearman.ave.CPM.B2), sub, pattern = "^X", replacement =  "")
dimnames(full.spearman.ratio.sampleAB.ratio2) <- lapply(dimnames(full.spearman.ratio.sampleAB.ratio2), sub, pattern = "^X", replacement =  "")

sample.info.logfc <- data.frame(row.names=colnames(full.spearman.ave.CPM.A2), 
                                lib.method.simple=strsplit2(strsplit2(colnames(full.spearman.ave.CPM.A2), split=".", fixed=TRUE)[,1], split="_")[,1],
                                lib.method.detail=strsplit2(colnames(full.spearman.ave.CPM.A2), split=".", fixed=TRUE)[,1])
sample.info.logfc$lib.method.simple <- factor(sample.info.logfc$lib.method.simple, levels = lib.simple.levels)
sample.info.logfc$lib.method.detail <- factor(sample.info.logfc$lib.method.detail, levels = lib.detail.levels)


dA <- dist(full.spearman.ave.CPM.A2)
dB <- dist(full.spearman.ave.CPM.B2)
d.mean <- as.dist((as.matrix(dA)+as.matrix(dB))/2)
this.outfile <- paste0(outdirs["FIG3"], "/FIG3B_SPEARMAN_COR_HEATMAP_CPM_unnormalized.pdf")
pA <- pheatmap(full.spearman.ave.CPM.A2, fontsize = 8,annotation_names_row = FALSE, annotation_names_col = FALSE,
               annotation_legend = FALSE, 
               labels_row = substring(colnames(full.spearman.ave.CPM.A2), nchar(colnames(full.spearman.ave.CPM.A2))), 
               labels_col = substring(row.names(full.spearman.ave.CPM.A2), nchar(row.names(full.spearman.ave.CPM.A2))), 
               legend = FALSE, 
               annotation_row = sample.info.logfc,
               annotation_colors=ann_colors, 
               annotation_col=sample.info.logfc,
               breaks=seq(0,1, length.out = 101),
               treeheight_row=10,
               treeheight_col=10,
               width=3,
               height=3,
               clustering_distance_rows = d.mean,
               clustering_distance_cols = d.mean)$gtable

pB <-pheatmap(full.spearman.ave.CPM.B2, fontsize = 8,annotation_names_row = FALSE, annotation_names_col = FALSE,
              annotation_legend = FALSE, 
              labels_row = substring(colnames(full.spearman.ave.CPM.A2), nchar(colnames(full.spearman.ave.CPM.A2))), 
              labels_col = substring(row.names(full.spearman.ave.CPM.A2), nchar(row.names(full.spearman.ave.CPM.A2))), 
              legend = FALSE, 
              annotation_row = sample.info.logfc,
              annotation_colors=ann_colors, 
              annotation_col=sample.info.logfc,
              breaks=seq(0,1, length.out = 101),
              treeheight_row=10,
              treeheight_col=10,
              width=3,
              height=3,
              clustering_distance_rows = d.mean,
              clustering_distance_cols = d.mean)$gtable
pRatio <- pheatmap(full.spearman.ratio.sampleAB.ratio2, fontsize = 8, annotation_names_row = FALSE, annotation_names_col = FALSE,
                   annotation_legend = FALSE, 
                   labels_row = substring(colnames(full.spearman.ave.CPM.A2), nchar(colnames(full.spearman.ave.CPM.A2))), 
                   labels_col = substring(row.names(full.spearman.ave.CPM.A2), nchar(row.names(full.spearman.ave.CPM.A2))), 
                   legend = FALSE, 
                   annotation_row = sample.info.logfc,
                   annotation_colors=ann_colors, 
                   annotation_col=sample.info.logfc,
                   #breaks=seq(0,1, length.out = 101),
                   treeheight_row=10,
                   treeheight_col=10,
                   width=3,
                   height=3,
                   clustering_distance_rows = d.mean,
                   clustering_distance_cols = d.mean)$gtable
ratio.cor.heatmap <- plot_grid(plotlist = list(pA, pB, pRatio), nrow=1)
save_plot(this.outfile, ratio.cor.heatmap, ncol = 3, nrow = 1, base_height = 3, base_width = 2.5)
# Contains all the annotation/label info for adding later
this.outfile <- paste0(outdirs["FIG3"], "/FULL_ANNOT_FIG3B_SPEARMAN_COR_HEATMAP_CPM_unnormalized.pdf")
pheatmap(full.spearman.ave.CPM.A2, fontsize = 8,
               legend = FALSE, 
          annotation_legend=FALSE,
         labels_col = sub(".Lab[0-9]", "", colnames(full.spearman.ave.CPM.A2)),
               annotation_row = sample.info.logfc,
               annotation_colors=ann_colors, 
               annotation_col=sample.info.logfc,
               breaks=seq(0,1, length.out = 101),
               treeheight_row=10,
               treeheight_col=10,
               width=3,
               height=3,
               clustering_distance_rows = d.mean,
               clustering_distance_cols = d.mean, 
              filename=this.outfile)




# SUMMARY CORRELATION COEFFICIENT TABLES FOR RATIOMETRIC POOL With Replicates -------------
# WITH REPLICATES
dge.ratio <- DGEList(ratiometric.counts.addinfo.sizeFilt.cast.df, genes=gene.info.ratio, remove.zeros = TRUE)
dge.ratio$samples <- cbind(dge.ratio$samples, sample.info.ratio)
dge.ratio$samples$pool <- factor(dge.ratio$samples$pool, c("SynthA", "SynthB"))
cpm.ratio <- cpm(dge.ratio, normalized.lib.sizes = FALSE)

cor.matrix.all.ratio.square <- cor(cpm.ratio.rle, use="pair", method="spearman")
cor.matrix.all.ratio <- data.table(cbind(melt(cor.matrix.all.ratio.square), melt(upper.tri(cor.matrix.all.ratio.square, diag=TRUE))))
setnames(cor.matrix.all.ratio, c("lab.libMethod.replicate.A", "lab.libMethod.replicate.B", "rho", "row", "col", "upper.tri"))
cor.matrix.all.ratio[, `:=`(id.A=ifelse(upper.tri==TRUE, as.character(lab.libMethod.replicate.A), as.character(lab.libMethod.replicate.B)), id.B=ifelse(upper.tri==TRUE, as.character(lab.libMethod.replicate.B), as.character(lab.libMethod.replicate.A)))]
cor.matrix.all.ratio[, `:=`(lab.libMethod.replicate.A=id.A, lab.libMethod.replicate.B=id.B)]
cor.matrix.all.ratio[, `:=`(id.A=NULL, id.B=NULL)]
dge.ratio$samples$lib.detail.pool <- 
  short.sample.info.ratio.A <- dge.ratio$samples[, c("Lab", "lab.libMethod.pool", "pool", "lib.method.detail", "lib.method.simple", "replicate")]
short.sample.info.ratio.B <- dge.ratio$samples[, c("Lab", "lab.libMethod.pool", "pool", "lib.method.detail", "lib.method.simple", "replicate")]

colnames(short.sample.info.ratio.A) <- paste0(colnames(short.sample.info.ratio.A), ".A")
colnames(short.sample.info.ratio.B) <- paste0(colnames(short.sample.info.ratio.B), ".B")

cor.matrix.all.ratio.tmp <- merge.data.frame(cor.matrix.all.ratio, short.sample.info.ratio.A, by.x=1, by.y=0)
cor.matrix.all.ratio.sample.info <- data.table(merge.data.frame(cor.matrix.all.ratio.tmp, short.sample.info.ratio.B, by.x=2, by.y=0))
lib.detail.levels <- levels(cor.matrix.all.ratio.sample.info$lib.method.detail.A)
reorder.paste.fun <- function(query.vec, lookup.vec){
  matches <- match(query.vec, lookup.vec)
  comp <- paste(query.vec[order(matches)], collapse=".VS.")
  return(comp)
}

cor.matrix.all.ratio.sample.info[, comparison.detail:=apply(.SD, 1, reorder.paste.fun, lookup.vec=lib.detail.levels.with.nextFlex), .SDcols=c("lib.method.detail.A", "lib.method.detail.B")]
cor.matrix.all.ratio.sample.info[, comparison.simple:=apply(.SD, 1, reorder.paste.fun, lookup.vec=lib.simple.levels), .SDcols=c("lib.method.simple.A", "lib.method.simple.B")]


setkeyv(cor.matrix.all.ratio.sample.info, c("lab.libMethod.replicate.B", "lab.libMethod.replicate.A"))
cor.matrix.all.ratio.sample.info.unique <- unique(cor.matrix.all.ratio.sample.info)
cor.matrix.all.ratio.sample.info.unique[, comparison.detail.pool:=ifelse(pool.A==pool.B, paste(comparison.detail, pool.A, sep="."), paste0(comparison.detail, ".AvsB")), by=1:nrow(cor.matrix.all.ratio.sample.info.unique)]
cor.matrix.all.ratio.sample.info.unique[, comparison.simple.pool:=ifelse(pool.A==pool.B, paste(comparison.simple, pool.A, sep="."), paste0(comparison.simple, ".AvsB")), by=1:nrow(cor.matrix.all.ratio.sample.info.unique)]

n.miRs <- dim(cpm.ratio.rle)[1]
cor.matrix.all.ratio.sample.info.cor.summary.details <- cor.matrix.all.ratio.sample.info.unique[, { nr=length(unique(c(lab.libMethod.replicate.A, lab.libMethod.replicate.B))); FUNCTION.meanRhoFromCoefs(rho, ns=n.miRs, nr=nr) }, by=comparison.detail.pool]
cor.matrix.all.ratio.sample.info.cor.summary.simple <- cor.matrix.all.ratio.sample.info.unique[, { nr=length(unique(c(lab.libMethod.replicate.A, lab.libMethod.replicate.B))); FUNCTION.meanRhoFromCoefs(rho, ns=n.miRs, nr=nr) }, by=comparison.simple.pool]
setnames(cor.matrix.all.ratio.sample.info.cor.summary.details, "comparison.detail.pool", "comparison")
setnames(cor.matrix.all.ratio.sample.info.cor.summary.simple, "comparison.simple.pool", "comparison")
cor.matrix.all.ratio.sample.info.cor.summary.all <- rbind(cor.matrix.all.ratio.sample.info.cor.summary.simple, cor.matrix.all.ratio.sample.info.cor.summary.details)
cor.matrix.all.ratio.sample.info.cor.summary.all[, c("MethodA", "MethodB", "Type"):=tstrsplit(comparison, split=".", fixed=TRUE)[c(1,3,4)]]

this.outfile <- paste0(outdirs["tables"], "/Spearman_correlation_ratiometric_Pool_summaries_with_replicates.txt")
write.table(cor.matrix.all.ratio.sample.info.cor.summary.all, this.outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


# SUP FigS4 Intralab %CV Ratio -------------------------------
INTRALAB.CV.RATIOMETRIC <- FINAL.RATIOMETRIC.LONG[, .(mean.cpm=mean(pseudo.cpm.filt.lengths), sd.cpm=sd(pseudo.cpm.filt.lengths), q1.cpm=quantile(pseudo.cpm.filt.lengths, 0.25), q3.cpm=quantile(pseudo.cpm.filt.lengths, 0.75)), by=.(lab.libMethod,  Lab, lib.method.detail, lib.method.simple, ratio.seqID, pool)]
INTRALAB.CV.RATIOMETRIC[, `:=`(intralab.cv=100*sd.cpm/mean.cpm, intralab.qcd=((q3.cpm-q1.cpm)/2)/((q3.cpm+q1.cpm)/2))]

this.outfile <- paste0(outdirs["supplemental_figures"], "/FIG_S4_RATIOMETRIC_CV_INTRALAB_Violin.pdf"); this.outfile
ggplot(INTRALAB.CV.RATIOMETRIC,
       aes(x = lab.libMethod,
           y = intralab.cv,
           alpha = pool,
           pos = pool,
           fill=lib.method.simple))  +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + scale_alpha_discrete(range=c(0.25, 1)) +
  scale_y_continuous(limits = c(0,190), breaks = seq(0, 180, 20), expand=c(0,0)) +
  theme(
    axis.text.x = element_text(angle=50, hjust=1),
    panel.grid.major = element_line(color = NA),
    panel.grid.minor = element_line(color = NA),
    panel.border = element_rect(size = 1, color = "black"),
    axis.ticks = element_line(color = "black", size = 1),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key.size=unit(1.0, "lines")
  ) +
  labs(title = "RATIOMETRIC: %CV Within-labs",
       x = NULL,
       y = "%CV") 
ggsave(this.outfile, width = 7.5, height=5.0)


# SUP FigS4 Intralab QCD Ratio -------------------------------
this.outfile <- paste0(outdirs["supplemental_figures"], "/FIGS4_RATIOMETRIC_QCD_INTRALAB.pdf"); this.outfile

#Intralab qcd violin
ggplot(INTRALAB.CV.RATIOMETRIC,
       aes(x = lab.libMethod,
           y = intralab.qcd,
           pos = pool,
           alpha = pool,
           fill=lib.method.simple))  +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_y_continuous(limits = c(0,0.9), breaks = seq(0, 1, 0.1), expand=c(0,0)) +
  scale_alpha_discrete(range=c(0.25, 1)) +
  theme(
    axis.text.x = element_text(angle=50, hjust=1),
    panel.grid.major = element_line(color = NA),
    panel.grid.minor = element_line(color = NA),
    panel.border = element_rect(size = 1, color = "black"),
    axis.ticks = element_line(color = "black", size = 1),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key.size=unit(1.0, "lines")
  ) +
  labs(title = "RATIOMETRIC: QCD Within-labs",
       x = NULL,
       y = "QCD") 
ggsave(this.outfile, width = 7.5, height=5.0)


# Intralab CV Ratio Summary Table -------------------------------
summary.quantiles <- c(0.02, 0.25, 0.50, 0.75, 0.98)
INTRALAB.CV.RATIOMETRIC[lib.method.detail!="4N_NEXTflex", (sub("^0.", "intralab.cv", summary.quantiles)):=lapply(summary.quantiles, USE.NAMES=FALSE, SIMPLIFY=FALSE, quantile, x=intralab.cv), by=.(lib.method.detail, pool)]
INTRALAB.CV.RATIOMETRIC[lib.method.detail!="4N_NEXTflex", (sub("^0.", "intralab.qcd", summary.quantiles)):=lapply(summary.quantiles, USE.NAMES=FALSE, SIMPLIFY=FALSE, quantile, x=intralab.qcd), by=.(lib.method.detail, pool)]
INTRALAB.CV.RATIOMETRIC.SUMMARY <- INTRALAB.CV.RATIOMETRIC[lib.method.detail!="4N_NEXTflex", lapply(.SD, max), by=.(lib.method.detail, pool), .SDcols=c(sub("^0.", "intralab.cv", summary.quantiles), sub("^0.", "intralab.qcd", summary.quantiles))]
this.outfile <- paste0(outdirs["tables"], "/Ratiometric_intralab_CV_QCD_summary.txt")
write.table(INTRALAB.CV.RATIOMETRIC.SUMMARY, this.outfile, row.names = FALSE, col.names=TRUE, sep="\t", quote=FALSE)



# SUMMARY CORRELATION COEFFICIENT TABLES FOR RATIOMETRIC POOL -- USING AVERAGE COUNTS AND FOLD CHANGES --------------------------------------------------------------
# AVERAGED OVER REPLICATES 
# Convert to matrix format, casting samples and measures
MEAN.CPM.RATIO.FOR.CPM.FOR.RATIO.ALL <- FINAL.RATIOMETRIC.LONG[, .(mean.cpm.over.pool.repl=mean(pseudo.cpm.filt.lengths)), by=.(lab.libMethod,  Lab, lib.method.detail, lib.method.simple, ratio.seqID, pool)]

MEAN.CPM.RATIO.FOR.CPM.FOR.RATIO.castAB.ALL <- dcast.data.table(MEAN.CPM.RATIO.FOR.CPM.FOR.RATIO.ALL, lab.libMethod+Lab+lib.method.detail+lib.method.simple+ratio.seqID~pool)
MEAN.CPM.RATIO.FOR.CPM.FOR.RATIO.castAB.ALL[, ratioAB:=SynthA/SynthB]
MEAN.CPM.RATIO.FOR.CPM.FOR.RATIO.cast.lib.and.measure <- as.matrix(data.frame(dcast.data.table(MEAN.CPM.RATIO.FOR.CPM.FOR.RATIO.castAB.ALL, ratio.seqID ~ lab.libMethod, sep="_log2_", fun.aggregate = log2, fill = 0, value.var = c("SynthA", "SynthB", "ratioAB")),  row.names=1))


cor.matrix.all.ratio.square <- cor(MEAN.CPM.RATIO.FOR.CPM.FOR.RATIO.cast.lib.and.measure, use="pair", method="spearman")
cor.matrix.all.ratio <- data.table(cbind(melt(cor.matrix.all.ratio.square), melt(upper.tri(cor.matrix.all.ratio.square, diag=TRUE))))
setnames(cor.matrix.all.ratio, c("measure.lab.libMethod.replicate.A", "measure.lab.libMethod.replicate.B", "rho", "row", "col", "upper.tri"))
cor.matrix.all.ratio[, `:=`(id.A=ifelse(upper.tri==TRUE, as.character(measure.lab.libMethod.replicate.A), as.character(measure.lab.libMethod.replicate.B)), id.B=ifelse(upper.tri==TRUE, as.character(measure.lab.libMethod.replicate.B), as.character(measure.lab.libMethod.replicate.A)))]
cor.matrix.all.ratio[, `:=`(measure.lab.libMethod.replicate.A=id.A, measure.lab.libMethod.replicate.B=id.B)]
cor.matrix.all.ratio[, `:=`(id.A=NULL, id.B=NULL)]
cor.matrix.all.ratio[, c("measure.A", "lab.libMethod.replicate.A"):=tstrsplit(measure.lab.libMethod.replicate.A, split="_log2_")]
cor.matrix.all.ratio[, c("measure.B", "lab.libMethod.replicate.B"):=tstrsplit(measure.lab.libMethod.replicate.B, split="_log2_")]


short.sample.info.ratio.A <- data.frame(unique(subset(MEAN.CPM.RATIO.FOR.CPM.FOR.RATIO.castAB.ALL, select=c("Lab", "lab.libMethod", "lib.method.detail", "lib.method.simple"))), row.names=2)
short.sample.info.ratio.B <- data.frame(unique(subset(MEAN.CPM.RATIO.FOR.CPM.FOR.RATIO.castAB.ALL, select=c("Lab", "lab.libMethod", "lib.method.detail", "lib.method.simple"))), row.names=2)

colnames(short.sample.info.ratio.A) <- paste0(colnames(short.sample.info.ratio.A), ".A")
colnames(short.sample.info.ratio.B) <- paste0(colnames(short.sample.info.ratio.B), ".B")

x.col.1 <- which(colnames(cor.matrix.all.ratio)=="lab.libMethod.replicate.A")
x.col.2 <- which(colnames(cor.matrix.all.ratio)=="lab.libMethod.replicate.B")
cor.matrix.all.ratio.tmp <- merge.data.frame(cor.matrix.all.ratio, short.sample.info.ratio.A, by.x=x.col.1, by.y=0)
cor.matrix.all.ratio.sample.info <- data.table(merge.data.frame(cor.matrix.all.ratio.tmp, short.sample.info.ratio.B, by.x=x.col.2, by.y=0))
reorder.paste.fun <- function(query.vec, lookup.vec){
  matches <- match(query.vec, lookup.vec)
  comp <- paste(query.vec[order(matches)], collapse=".VS.")
  return(comp)
}

cor.matrix.all.ratio.sample.info[, comparison.detail:=apply(.SD, 1, reorder.paste.fun, lookup.vec=lib.detail.levels.with.nextFlex), .SDcols=c("lib.method.detail.A", "lib.method.detail.B")]
cor.matrix.all.ratio.sample.info[, comparison.simple:=apply(.SD, 1, reorder.paste.fun, lookup.vec=lib.simple.levels), .SDcols=c("lib.method.simple.A", "lib.method.simple.B")]

cor.matrix.all.ratio.sample.info[, comparison.detail:=apply(.SD, 1, paste, collapse=".VS."), .SDcols=c("lib.method.detail.A", "lib.method.detail.B")]
cor.matrix.all.ratio.sample.info[, comparison.simple:=apply(.SD, 1, paste, collapse=".VS."), .SDcols=c("lib.method.simple.A", "lib.method.simple.B")]

setkeyv(cor.matrix.all.ratio.sample.info, c("measure.lab.libMethod.replicate.B", "measure.lab.libMethod.replicate.A"))
cor.matrix.all.ratio.sample.info.unique.selfCompare <- subset(unique(cor.matrix.all.ratio.sample.info), measure.A==measure.B)
cor.matrix.all.ratio.sample.info.unique.selfCompare[, comparison.detail.measure:=paste(comparison.detail, measure.A, sep="."), by=1:nrow(cor.matrix.all.ratio.sample.info.unique.selfCompare)]
cor.matrix.all.ratio.sample.info.unique.selfCompare[, comparison.simple.measure:=paste(comparison.simple, measure.A, sep="."), by=1:nrow(cor.matrix.all.ratio.sample.info.unique.selfCompare)]

n.miRs <- dim(cpm.ratio.rle)[1]
cor.matrix.all.ratio.sample.info.cor.summary.details <- cor.matrix.all.ratio.sample.info.unique.selfCompare[, { nr=length(unique(c(lab.libMethod.replicate.A, lab.libMethod.replicate.B))); FUNCTION.meanRhoFromCoefs(rho, ns=n.miRs, nr=nr) }, by=comparison.detail.measure]
cor.matrix.all.ratio.sample.info.cor.summary.simple <- cor.matrix.all.ratio.sample.info.unique.selfCompare[, { nr=length(unique(c(lab.libMethod.replicate.A, lab.libMethod.replicate.B))); FUNCTION.meanRhoFromCoefs(rho, ns=n.miRs, nr=nr) }, by=comparison.simple.measure]
setnames(cor.matrix.all.ratio.sample.info.cor.summary.details, "comparison.detail.measure", "comparison")
setnames(cor.matrix.all.ratio.sample.info.cor.summary.simple, "comparison.simple.measure", "comparison")
cor.matrix.all.ratio.sample.info.cor.summary.all <- unique(rbind(cor.matrix.all.ratio.sample.info.cor.summary.simple, cor.matrix.all.ratio.sample.info.cor.summary.details), by="comparison")
cor.matrix.all.ratio.sample.info.cor.summary.all[, c("MethodA", "MethodB", "Type"):=tstrsplit(comparison, split=".", fixed=TRUE)[c(1,3,4)]]

this.outfile <- paste0(outdirs["tables"], "/Spearman_correlation_ratiometric_Pool_summaries_mean_cpm_and_ratios.txt")
write.table(cor.matrix.all.ratio.sample.info.cor.summary.all, this.outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# Ratio correlatoin summary for in-text summary ----
setkey(cor.matrix.all.ratio.sample.info.cor.summary.all, comparison)
comparison.levels.order.ratio <- paste0(unique(as.character(melt(outer(lib.detail.levels, lib.detail.levels, FUN = paste, sep=".VS."))$value)), ".ratioAB")
comparison.levels.order.ratio <- comparison.levels.order[comparison.levels.order%in%cor.matrix.all.ratio.sample.info.cor.summary.all$comparison]



# For in-text ratiometric correlation summary numbers----
paste(
  cor.matrix.all.ratio.sample.info.cor.summary.all[paste0(c("TruSeq.VS.TruSeq", "NEBNext.VS.NEBNext", "4N_B.VS.4N_B"), ".ratioAB"),
                                                                              paste0(
                                                                                sub(".VS.*", "", comparison),
                                                                                ": ", 
                                                                                round(simple.mean.r, digits = 2),
                                                                                " (", 
                                                                                round(r02, digits = 2), 
                                                                                ", ", 
                                                                                round(r98, digits = 2),
                                                                                "; n=",
                                                                                nlibs,  
                                                                                ")")], 
  collapse="; ")
# TruSeq: 0.99 (0.98, 1; n=8); NEBNext: 0.98 (0.97, 0.99; n=6); 4N_B: 0.99 (0.98, 0.99; n=4); 4N: 0.98 (0.95, 0.99; n=10)

# inter-protocol comparisons
gsub(pattern = ".VS.",
     replacement = " vs ", 
     x=paste(
        cor.matrix.all.ratio.sample.info.cor.summary.all[paste0(Cs(TruSeq.VS.NEBNext, TruSeq.VS.CleanTag, TruSeq.VS.4N_B, CleanTag.VS.NEBNext, NEBNext.VS.4N_B, CleanTag.VS.4N_B), 
                                                                ".ratioAB"),
                                                         paste0(
                                                           sub(".ratioAB", "", comparison),
                                                           ": ",
                                                                round(simple.mean.r, digits = 2), 
                                                                " (",
                                                                round(r02, digits = 2),
                                                                ", ",
                                                                round(r98, digits = 2),
                                                                "; n=", 
                                                                nlibs,
                                                                ")")],
        collapse="; "))
# TruSeq vs NEBNext: 0.96 (0.94, 0.98; n=14); TruSeq vs CleanTag: 0.95 (0.94, 0.96; n=9); TruSeq vs 4N_B: 0.97 (0.96, 0.98; n=12); CleanTag vs NEBNext: 0.95 (0.93, 0.96; n=7); NEBNext vs 4N_B: 0.96 (0.94, 0.98; n=10); CleanTag vs 4N_B: 0.96 (0.96, 0.97; n=5)
cor.matrix.all.ratio.sample.info.cor.summary.all.switch <- copy(cor.matrix.all.ratio.sample.info.cor.summary.all)
cor.matrix.all.ratio.sample.info.cor.summary.all.switch[, comparison:=paste0(MethodB, ".VS.", MethodA, ".", Type)]
setnames(cor.matrix.all.ratio.sample.info.cor.summary.all.switch, c("MethodA", "MethodB"), c("MethodB", "MethodA"))
cor.matrix.all.ratio.sample.info.cor.summary.all.dupl <- unique(rbind(cor.matrix.all.ratio.sample.info.cor.summary.all, cor.matrix.all.ratio.sample.info.cor.summary.all.switch), by=c("comparison"))

levels.use <- lib.detail.levels
comparison.levels.order.ratio.show.matr <- expand.grid(levels.use, levels.use)
cor.matrix.all.ratio.sample.info.cor.summary.all.dupl.filt <- subset(cor.matrix.all.ratio.sample.info.cor.summary.all.dupl, MethodA%in%levels.use & MethodB%in%levels.use)
cor.matrix.all.ratio.sample.info.cor.summary.all.dupl.filt[, `:=`(MethodA=factor(MethodA, levels=levels.use), MethodB=factor(MethodB, levels=levels.use))]
# Not shown, but makes grid showing correlation summaries for hte ratio and the A and B separate
g <- ggplot(cor.matrix.all.ratio.sample.info.cor.summary.all.dupl.filt, aes(x=Type)) +
geom_pointrange(aes(y=simple.mean.r, ymin=r02, ymax=r98)) + facet_grid(MethodA~MethodB); g


# OVerall metric idea for correlation summary 
# 
cor.matrix.all.ratio.sample.info.unique.selfCompare[, overall.measure.type:=paste0(measure.A, ifelse(lib.method.detail.A==lib.method.detail.B, ".WithinMethod", ".AcrossMethods"))]

cor.matrix.all.ratio.sample.info.cor.summary.details.overall <- cor.matrix.all.ratio.sample.info.unique.selfCompare[, { nr=length(unique(c(lab.libMethod.replicate.A, lab.libMethod.replicate.B))); FUNCTION.meanRhoFromCoefs(rho, ns=n.miRs, nr=nr) }, by=overall.measure.type]



# DOWNSAMPLING: Equimolar ------------------------------
dge.eq.counts <- dge.eq$counts
# SO REPEAT SAME PARAMS FOR EQUIMOLAR
downsample.to.vec <- c(Inf, round(10^seq(4, 6.5, 0.5), 0))
drr.list <- sapply(downsample.to.vec, simplify=FALSE, function(x) t(drarefy(t(dge.eq.counts), x)))

names(drr.list) <- downsample.to.vec
drr.all <- data.table(rbindlist(lapply(drr.list, melt), idcol=TRUE, use.names=TRUE), keep.rownames=TRUE)
setnames(drr.all, c("downsample.to", "miR.id", "lab.libMethod.replicate", "est.prob.detected"))
setkey(drr.all, lab.libMethod.replicate)
sample.info.eq.dt <- data.table(dge.eq$samples, keep.rownames=TRUE)
setnames(sample.info.eq.dt, "rn", "lab.libMethod.replicate" )

drr.all.merge <- merge(drr.all, sample.info.eq.dt)
drr.all.merge[, downsample.to:=as.numeric(downsample.to)]
drr.all.filt <- subset(drr.all.merge, lib.size>=downsample.to|lib.size==Inf)

# plot cutoff detected
min.prob.detected=0.90
#drr.all.filt.summarize.repl <- drr.all.filt[, .(est.prob.detected.all=prod(est.prob.detected)), by=.(lab.libMethod, downsample.to, Lab, lib.method.detail, lib.method.simple,  miR.id)]
#est.prob.detected.n.miRs <- drr.all.filt.summarize.repl[, .(n.detected=sum(ifelse(est.prob.detected.all>=min.prob.detected, 1, 0)), n.not.detected=sum(ifelse(est.prob.detected.all>=min.prob.detected, 0, 1)), n.total=.N), by=.(lab.libMethod, downsample.to, Lab, lib.method.detail, lib.method.simple)]


est.prob.detected.n.miRs <- drr.all.filt[, .(n.detected=sum(ifelse(est.prob.detected>=min.prob.detected, 1, 0)), n.not.detected=sum(ifelse(est.prob.detected>=min.prob.detected, 0, 1)), n.total=.N), by=.(lab.libMethod.replicate, lab.libMethod, downsample.to, Lab, lib.method.detail, lib.method.simple)]

est.prob.detected.n.miRs[, percent.detected:=n.detected/max(n.total)]
est.prob.detected.n.miRs[, downsampled.to.group:=factor(round(log10(downsample.to), digits=1))]
est.prob.detected.n.miRs[, lib.method.simple:=factor(lib.method.simple, c("TruSeq", "CleanTag", "NEBNext", "4N"))]
est.prob.detected.n.miRs[, n.not.detected:=n.total-n.detected]

est.prob.detected.n.miRs[, downsampled.to.group:=factor(downsampled.to.group, levels=rev(levels(downsampled.to.group)))]
group.levels <- rev(downsample.to.vec[-1])

g <- ggplot(est.prob.detected.n.miRs,
            aes(
              x=lab.libMethod,
              y=n.not.detected, fill=lib.method.simple)) + 
  #geom_bar(stat="identity", color="black", size=1) + 
  stat_summary(fun.y = median,fun.ymax = max, fun.ymin = min) +
  theme_bw() +
  scale_x_discrete(drop=FALSE) +
  facet_wrap(~downsampled.to.group, ncol=2) +
  labs(x=NULL,
       y="# miRNAs with < 90% Probability of Detection") +
  #scale_y_continuous(expand=c(0, 0)) +
  theme(panel.grid.minor = element_blank(),
        axis.text=element_text(size=8, color="black"),
        axis.text.x = element_text(angle=50, hjust=1, face="bold", color="black"),
        axis.line=element_line(color="black"),
        axis.title=element_text(size=8),
        plot.title=element_text(size=8),
        panel.grid.major = element_line(color=NA),
        #panel.border = element_blank(),
        #panel.spacing = unit(0, "lines"),
        #strip.background = element_blank(),
        strip.text=element_text(size=8, color="black"),
        axis.ticks = element_line(color="black", size=1.5),
        legend.position = "none"
        ) 
plot_grid(plotlist = list(g %+% est.prob.detected.n.miRs[downsample.to%in%group.levels[1:2],], 
g %+% est.prob.detected.n.miRs[downsample.to%in%group.levels[3:4],],
g %+% est.prob.detected.n.miRs[downsample.to%in%group.levels[5:6],]), nrow = 3)

#ggsave("SUPPLEMENTAL_FIGS2_downsampling_not_detected_1.svg", g, width=8, height=4, units="in")
#ggsave("SUPPLEMENTAL_FIGS2_downsampling_not_detected_2.svg", g %+% est.prob.detected.n.miRs[downsample.to%in%group.levels[3:4],], width=8, height=4, units="in")
#ggsave("SUPPLEMENTAL_FIGS2_downsampling_not_detected_3.svg", g %+% est.prob.detected.n.miRs[downsample.to%in%group.levels[5:6],], width=8, height=4, units="in")
# SUP FIGS2: Equimolar pool downsample percent detected ------------------------------
this.outfile <- paste0(outdirs["supplemental_figures"], "/SUPPLEMENTAL_FIGS2_equimolar_percent_detected_downsample.pdf")
g <- ggplot(est.prob.detected.n.miRs,
            aes(
              x=lab.libMethod,
              y=percent.detected, color=lib.method.simple)) + 
  stat_summary(fun.y = median,fun.ymax = max, fun.ymin = min) +
  scale_x_discrete(drop=FALSE) +
  facet_wrap(~downsampled.to.group, ncol=2) +
  labs(x=NULL,
       y="# miRNAs with > 90% Probability of Detection") +
  scale_y_continuous(labels = percent_format(), breaks = seq(0,1,0.1), limits=c(0,1)) +
  theme(panel.grid.minor = element_blank(),
        axis.text=element_text(size=8, color="black"),
        axis.text.x = element_text(angle=50, hjust=1, color="black"),
        axis.line=element_line(color="black"),
        axis.title=element_text(size=8),
        plot.title=element_text(size=8),
        panel.grid.major = element_line(color=NA),
        #panel.border = element_blank(),
        #panel.spacing = unit(0, "lines"),
        #strip.background = element_blank(),
        strip.text=element_text(size=8, color="black"),
        axis.ticks = element_line(color="black"),
        legend.position = "none"
  ); g
ggsave(this.outfile, g, width=7.5, height=7)

# Plasma Pool Wrangling ------------------------------
# Plasma sample filtering by counts ------------------
total.count.cutoff = 100000

mirs.countsummary <- counts.dt.endog.miRs.all[grepl("MIMAT", ReferenceID), .(total.mir.adjusted.counts=sum(multimapAdjustedReadCount)), by=.(lab.libMethod.pool.replicate, Lab, lib.method.detail, pool)]

mirs.countsummary[, n.replicates.below.cutoff:=sum(ifelse(total.mir.adjusted.counts<total.count.cutoff, 1, 0)), by=.(Lab, lib.method.detail, pool)]
# If only 1 replicate has low counts, return just that sample replicate ID. If more than 1 replicate is below the cutoff for total counts, return all the samples to filter
sample.filt <- as.character(mirs.countsummary[(n.replicates.below.cutoff==1 & total.mir.adjusted.counts<total.count.cutoff) | n.replicates.below.cutoff>1]$lab.libMethod)


# Casts the endogenous miR counts from both synthetic and plasma sources, summing over lanes where needed
# Filters out counts to stem loop sequences (only mature mIRs)
this.value.var <- "multimapAdjustedReadCount" # value to use in in the casted table
plasma.use.pools <- c("PlasmaPool", "SynthEQ")
mature.counts.dt.cast <- dcast.data.table(counts.dt.endog.miRs.all, 
                                          subset = .(grepl("MIMAT", ReferenceID) & pool%in%plasma.use.pools & !lab.libMethod.pool.replicate%in%sample.filt),
                                          miR.ID~lab.libMethod.pool.replicate,
                                          value.var = this.value.var,
                                          fun.aggregate = sum, 
                                          fill=0)


# Re-melt to add zerocounts back in
mature.counts.dt.zerocounts.tmp <- melt(mature.counts.dt.cast, value.name="multimapAdjustedReadCount", variable.name="lab.libMethod.pool.replicate", id.vars=c("miR.ID"))

setkey(mature.counts.dt.zerocounts.tmp, lab.libMethod.pool.replicate)
setkey(CROSS.U01.METADATA.ALL.CAST, lab.libMethod.pool.replicate)
ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE <- merge(mature.counts.dt.zerocounts.tmp, CROSS.U01.METADATA.ALL.CAST)
ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE[, total.mir.adjusted.counts:=sum(multimapAdjustedReadCount), by=.(lab.libMethod.pool.replicate)]
ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE[, cpm.multimap.adjust:=(multimapAdjustedReadCount*10^6)/total.mir.adjusted.counts]
ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE[, `:=`(lab.libMethod=factor(lab.libMethod, levels=new.order$all),
                                   lib.method.detail=factor(lib.method.detail, levels=lib.detail.levels),
                                   lib.method.simple=factor(lib.method.simple, levels=lib.simple.levels))]

# Filter to keep only plasma pool
mature.counts.dt.zerocounts <- subset(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE, pool=="PlasmaPool")
mature.counts.dt.zerocounts <- mature.counts.dt.zerocounts[, lab.libMethod:=factor(lab.libMethod, levels=new.order$plasma)]
mature.counts.dt.zerocounts[, count.rank.in.sample:=rank(-multimapAdjustedReadCount, ties.method = "min"), by=.(lab.libMethod.replicate)]
plasma.pool.top10.miR.summary <- mature.counts.dt.zerocounts[count.rank.in.sample<=10, .(sum.multimapAdjustedReadCount.top10.miRs=sum(multimapAdjustedReadCount)), by=.(lab.libMethod, lib.method.detail, lib.method.simple, replicate, total.mir.adjusted.counts)]
  
plasma.pool.top10.miR.summary[, percent.reads.top10.miRs:=sum.multimapAdjustedReadCount.top10.miRs/total.mir.adjusted.counts]
plasma.pool.top10.miR.summary[, log.total.adjusted.counts:=log10(total.mir.adjusted.counts)]

# SUP FIG S6 % contribution of top-10 miRs -- 4 separate bars -------------
this.outfile <- paste0(outdirs["supplemental_figures"], "/FIGS6_Percent_plasma_reads_from_top10_expressed_mirs_barchart.pdf")
ggplot(plasma.pool.top10.miR.summary, aes(x=lab.libMethod, y=percent.reads.top10.miRs, pos=factor(replicate, levels=seq(1,4)), alpha=1/as.numeric(replicate), fill=lib.method.simple)) + 
  geom_bar(stat="identity", color="black", pos="dodge", width=0.8) +
  scale_y_continuous(
    limits=c(0, 1),
    expand=c(0, 0),
    breaks = seq(0,1,0.1),
    labels=scales::percent
  ) +
  theme(
    axis.text=element_text(size=10, color="black"),
    axis.text.x=element_text(angle=50, hjust=1),
    axis.title=element_text(size=10),
    plot.title=element_text(size=10),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.0, color="black"),
    axis.ticks = element_line(color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=10),
    legend.title=element_text(size=10)
    #legend.key.size=unit(1.2, "lines")
  ) + labs(x=NULL, y="%Reads from Top-10 Expressed miRNAs")
ggsave(this.outfile, width=8, height=6)

# COMPARE SYNTHETIC ENDOGENOUS COUNT MAPPING TO PLASMA COUNTS 4N VS TRUSEQ -----------
drop.libs <- c("4N_Xu", "4N_NEXTflex")
ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.filt <- subset(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE, !lib.method.detail%in%drop.libs)
ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.counts <- data.frame(dcast.data.table(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.filt, miR.ID~lab.libMethod.pool.replicate, fill=0, fun.aggregate = round, digits=0, value.var=this.value.var), row.names=1)
ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.sample.info <- data.frame(subset(CROSS.U01.METADATA.ALL.CAST, lab.libMethod.pool.replicate%in%ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.filt[, unique(lab.libMethod.pool.replicate)]))
row.names(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.sample.info) <- make.names(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.sample.info$lab.libMethod.pool.replicate)
ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.sample.info <- ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.sample.info[colnames(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.counts),]
ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE <- DGEList(counts = ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.counts, samples = ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.sample.info, remove.zeros = TRUE)

ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE <- ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE[ , ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE$samples$pool=="PlasmaPool"]
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE[rowSums(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE$counts>0)>0,]
# USED THIS CUTOFF
cpm.cutoff=100

## 1-8-17 Think we should use a CPM>100 cutoff. The main concern is that it makes a messy plot when we use our standard cutoff of cpm>100 in all libraries of >= 1 protocol. Therefore, 
# let's try a little different one. On the downsampling cutoff, we chose to plot only commonly-detected miRs for the main text and supplement a heatmap with those detected accross at least one protocol.
# So lets to a cutoff of >= 90% libraries detecting the miR (at >=100cpm), and then plot the others (>100CPM in at least one method) as a supplement
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.norm <- calcNormFactors(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA, method="RLE")
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm <- cpm(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.norm, normalized.lib.sizes = TRUE, prior.count = 1, log = FALSE)
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.logcpm <- cpm(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.norm, normalized.lib.sizes = TRUE, prior.count = 1, log = TRUE)

# Now get the protocol-specific ones
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.gt100 <- ifelse(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm>=100, 1, 0)
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.all <- ifelse(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm>=-Inf, 1, 0)
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple <- factor(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple, lib.simple.levels)
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.n.libs.gt100 <- sumTechReps(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.gt100, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple)
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.n.libs <- sumTechReps(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.all, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple)

ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100 <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.n.libs.gt100/ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.n.libs
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.filt <- rowMaxs(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100)>=0.9

LOG.CPM.mirs.passing.100cpm.by.method <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.logcpm[ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.filt,]
#LOG.CPM.mirs.passing.100cpm.by.method.drop.common.miRs <- LOG.CPM.mirs.passing.100cpm.by.method[!row.names(LOG.CPM.mirs.passing.100cpm.by.method)%in%row.names(LOG.CPM.mirs.passing.90perc.all),]
# FIG5A -------------------------------------------------------------------
# Note: The "cutree_rows" parameter was used to split the dendrogram rows into two chunks. In illustrator, these two groups were 
# switched in position to put the group with the higher expression at the top--just for visualization purposes.
LOG.CPM.mirs.passing.100cpm.by.method.row.max <- rowMedians(LOG.CPM.mirs.passing.100cpm.by.method)

callback.mirs <- function(tree_row, mat){
  mir.wts <- LOG.CPM.mirs.passing.100cpm.by.method.row.max[tree_row$labels]
  dend <- reorder(as.dendrogram(tree_row), wts=mir.wts)
  as.hclust(dend)
}

plasma.col.info <- data.frame(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples[, c("lib.size", "Lab", "lib.method.detail", "lib.method.simple")], row.names=row.names(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples))
plasma.col.info$lib.method.simple <- factor(plasma.col.info$lib.method.simple, levels=lib.simple.levels)
plasma.col.info$lib.method.detail <- factor(plasma.col.info$lib.method.detail, levels=lib.detail.levels)
plasma.col.info$log.lib.size <- log2(plasma.col.info$lib.size)


pheatmap(LOG.CPM.mirs.passing.100cpm.by.method, 
         show_rownames = FALSE,
         breaks=seq(min(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.logcpm), max(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.logcpm), length.out=101),
         show_colnames=TRUE,
         cluster_rows = TRUE, 
         legend_breaks = seq(0, 20, 2),
         annotation_colors=FUNCTION.filter.ann_colors(ann_colors, plasma.col.info[, c("log.lib.size", "lib.method.detail", "lib.method.simple")]),
         labels_col = sub("PlasmaPool.", "", sub("^X", "", colnames(LOG.CPM.mirs.passing.100cpm.by.method))),
         annotation_col = plasma.col.info[, c("log.lib.size", "lib.method.detail", "lib.method.simple")],
         fontsize = 8,
         treeheight_row = 15,
         treeheight_col = 15,
         width=14, 
         height=14,
         clustering_callback = callback.mirs)
this.outfile <- paste0(outdirs["FIG5"], "/FIG5A_Plasma_Heatmap_RLE_WITH_LEGEND.pdf")
pheatmap(LOG.CPM.mirs.passing.100cpm.by.method, 
         show_rownames = FALSE,
         breaks=seq(min(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.logcpm), max(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.logcpm), length.out=101),
         show_colnames=TRUE,
         cluster_rows = TRUE, 
         legend_breaks = seq(0, 20, 2),
         annotation_colors=FUNCTION.filter.ann_colors(ann_colors, plasma.col.info[, c("log.lib.size", "lib.method.detail", "lib.method.simple")]),
         labels_col = sub("PlasmaPool.", "", sub("^X", "", colnames(LOG.CPM.mirs.passing.100cpm.by.method))),
         annotation_col = plasma.col.info[, c("log.lib.size", "lib.method.detail", "lib.method.simple")],
         fontsize = 8,
         treeheight_row = 15,
         treeheight_col = 15,
         width=14, 
         height=14, 
         cutree_rows = 2,
         clustering_callback = callback.mirs,
         filename=this.outfile)
this.outfile <- paste0(outdirs["FIG5"], "/FIG5A_Plasma_Heatmap_RLE_NO_LEGEND.pdf")
pheatmap(LOG.CPM.mirs.passing.100cpm.by.method, annotation_legend = FALSE, legend = FALSE,
         show_rownames = FALSE,
         breaks=seq(min(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.logcpm), max(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.logcpm), length.out=101),
         show_colnames=TRUE,
         cluster_rows = TRUE, 
         legend_breaks = seq(0, 20, 2),
         annotation_colors=FUNCTION.filter.ann_colors(ann_colors, plasma.col.info[, c("log.lib.size", "lib.method.detail", "lib.method.simple")]),
         labels_col = sub("PlasmaPool.", "", sub("^X", "", colnames(LOG.CPM.mirs.passing.100cpm.by.method))),
         annotation_col = plasma.col.info[, c("log.lib.size", "lib.method.detail", "lib.method.simple")],
         fontsize = 8,
         treeheight_row = 15,
         treeheight_col = 15,
         width=7.0, 
         height=4.0,
         clustering_callback = callback.mirs,
         filename=this.outfile)


# Plasma Intra- and Inter-lab variation
plasma.libs.used <- sub("^X", "", colnames(LOG.CPM.mirs.passing.100cpm.by.method))
ENDOG.MIR.COUNTS.PLASMA.MATURE <- subset(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE, lab.libMethod.pool.replicate %in%plasma.libs.used )
# In initial submission, did CV calculations for miRs with cutoff > 100 CPM in all samples for the method specified.
ENDOG.MIR.COUNTS.PLASMA.MATURE[, passes.cutoff.within.lab:=(sum(ifelse(cpm.multimap.adjust>=cpm.cutoff, 1, 0))/.N)==1.0, by=.(lab.libMethod, miR.ID)]
ENDOG.MIR.COUNTS.PLASMA.MATURE[, all.gt0.within.lab:=(sum(ifelse(cpm.multimap.adjust>0, 1, 0))/.N)==1.0, by=.(lab.libMethod, miR.ID)]
ENDOG.MIR.COUNTS.PLASMA.MATURE[, passes.cutoff.between.lab.simple:=(sum(ifelse(cpm.multimap.adjust>=cpm.cutoff, 1, 0))/sum(ifelse(cpm.multimap.adjust>=cpm.cutoff, 1, 1)))>=0.9, by=.(lib.method.simple, miR.ID)]
ENDOG.MIR.COUNTS.PLASMA.MATURE[, passes.cutoff.between.lab.detail:=(sum(ifelse(cpm.multimap.adjust>=cpm.cutoff, 1, 0))/sum(ifelse(cpm.multimap.adjust>=cpm.cutoff, 1, 1)))>=0.9, by=.(lib.method.detail, miR.ID)]
MEAN.CV.WI.LABS.plasma <- ENDOG.MIR.COUNTS.PLASMA.MATURE[, .(mean.cpm=mean(cpm.multimap.adjust), sd.cpm=sd(cpm.multimap.adjust), q1.cpm=quantile(cpm.multimap.adjust, 0.25), q3.cpm=quantile(cpm.multimap.adjust, 0.75), n=.N), by=.(lab.libMethod, lib.method.simple, lib.method.detail, miR.ID, all.gt0.within.lab, passes.cutoff.within.lab, passes.cutoff.between.lab.simple, passes.cutoff.between.lab.detail)]
MEAN.CV.WI.LABS.plasma[, `:=`(intralab.cv=100*sd.cpm/mean.cpm, intralab.qcd=((q3.cpm-q1.cpm)/2)/((q3.cpm+q1.cpm)/2))]


# Intralab CV PLASMA Summary Table -------------------------------

MEAN.CV.WI.LABS.plasma.copy <- copy(MEAN.CV.WI.LABS.plasma)

MEAN.CV.WI.LABS.plasma.copy.4N <- subset(MEAN.CV.WI.LABS.plasma.copy, grepl("4N_[ABCD]", lib.method.detail) &  passes.cutoff.between.lab.simple == TRUE)
MEAN.CV.WI.LABS.plasma.copy.4N[, lib.method.group:="4N"]
MEAN.CV.WI.LABS.plasma.copy.filt <- subset(MEAN.CV.WI.LABS.plasma.copy, passes.cutoff.between.lab.detail==TRUE)
MEAN.CV.WI.LABS.plasma.copy.filt[, lib.method.group:=lib.method.detail]

MEAN.CV.WI.LABS.plasma.subgroups <- rbind(MEAN.CV.WI.LABS.plasma.copy.4N, MEAN.CV.WI.LABS.plasma.copy.filt)
MEAN.CV.WI.LABS.plasma.subgroups[, (sub("^0.", "intralab.cv", summary.quantiles)):=lapply(summary.quantiles, USE.NAMES=FALSE, SIMPLIFY=FALSE, quantile, x=intralab.cv), by=.(lib.method.group)]
MEAN.CV.WI.LABS.plasma.subgroups[, (sub("^0.", "intralab.qcd", summary.quantiles)):=lapply(summary.quantiles, USE.NAMES=FALSE, SIMPLIFY=FALSE, quantile, x=intralab.qcd), by=.(lib.method.group)]
MEAN.CV.WI.LABS.plasma.subgroups.tmp <- subset(unique(MEAN.CV.WI.LABS.plasma.subgroups[lib.method.detail!="4N_NEXTflex"], by=c("lab.libMethod", "lib.method.group")), select=c("lib.method.group", "lab.libMethod", "n", sub("^0.", "intralab.cv", summary.quantiles), sub("^0.", "intralab.qcd", summary.quantiles)))
MEAN.CV.WI.LABS.plasma.subgroups.tmp[, n:=sum(n), by=lib.method.group]  
INTRALAB.CV.PLASMA.SUMMARY <- subset(unique(MEAN.CV.WI.LABS.plasma.subgroups.tmp, by="lib.method.group"), select=-lab.libMethod)
this.outfile <- paste0(outdirs["tables"], "/Plasma_intralab_CV_QCD_summary.txt")
write.table(INTRALAB.CV.PLASMA.SUMMARY, this.outfile, row.names = FALSE, col.names=TRUE, sep="\t", quote=FALSE)

# For in-text PLASMA Intralab CV/QCD summary numbers----
setorder(INTRALAB.CV.PLASMA.SUMMARY, intralab.cv5)
paste(
  INTRALAB.CV.PLASMA.SUMMARY[!grepl("4N_[ACDX]", lib.method.group), 
                                paste0(lib.method.group,
                                       ": ", 
                                       round(intralab.cv5, digits = 2),
                                       " (", 
                                       round(intralab.cv02, digits = 2), 
                                       ", ", 
                                       round(intralab.cv98, digits = 2),
                                       "; n=",
                                       n,  
                                       ")")], 
  collapse="; ")
# TruSeq: 7.7 (1.55, 29.3; n=23); 4N_B: 9.28 (2.16, 30.74; n=16); 4N: 9.49 (1.77, 34.54; n=28); NEBNext: 10.86 (2.44, 37.36; n=16); CleanTag: 24.9 (5.47, 52.39; n=4)
setorder(INTRALAB.CV.PLASMA.SUMMARY, intralab.qcd5)
paste(
  INTRALAB.CV.PLASMA.SUMMARY[!grepl("4N_[ACDX]", lib.method.group), 
                                paste0(lib.method.group,
                                       ": ", 
                                       round(intralab.qcd5, digits = 2),
                                       " (", 
                                       round(intralab.qcd02, digits = 2), 
                                       ", ", 
                                       round(intralab.qcd98, digits = 2),
                                       "; n=",
                                       n,  
                                       ")")], 
  collapse="; ")
# TruSeq: 0.04 (0.01, 0.14; n=23); 4N_B: 0.05 (0.01, 0.17; n=16); 4N: 0.05 (0.01, 0.21; n=28); NEBNext: 0.05 (0.01, 0.24; n=16); CleanTag: 0.1 (0.02, 0.24; n=4)



# FIG5A--%CV -------------------------------------------------------------------
cv.range <- c(min(MEAN.CV.WI.LABS.plasma[lib.method.detail!="4N_NEXTflex" & passes.cutoff.between.lab.simple==TRUE]$intralab.cv, na.rm = TRUE), max(MEAN.CV.WI.LABS.plasma[lib.method.detail!="4N_NEXTflex" & passes.cutoff.between.lab.simple==TRUE]$intralab.cv, na.rm=TRUE))
# Use the cutoff from 5A--that is, find those where >90% of samples have >=100CPM and plot those. Just subset by method.
this.outfile <- paste0(outdirs["FIG5"], "/FIG5B_PLASMA_PERCENT_CV_INTRALAB_USE5ACutoff.pdf")
ggplot(MEAN.CV.WI.LABS.plasma[lib.method.detail!="4N_NEXTflex" & passes.cutoff.between.lab.simple==TRUE],
       aes(x = lab.libMethod,
           y = intralab.cv,
           fill=lib.method.simple))  +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0, 200, 20), expand=c(0,0)) +
  theme(
    axis.text.x = element_text(angle=50, hjust=1),
    panel.grid.major = element_line(color = NA),
    panel.grid.minor = element_line(color = NA),
    panel.border = element_rect(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  ) +
  labs(
    x=NULL,
    y = "%CV") 
ggsave(this.outfile, height=2.6, width=4.5)


# FIG5A--QCD --------------------------------------------------------------
#Intralab cv, violin
this.outfile <- paste0(outdirs["FIG5"], "/FIG5B_PLASMA_QCD_INTRALAB_USE5ACutoff.pdf")
QCD.WI.LABS.plasma <- ENDOG.MIR.COUNTS.PLASMA.MATURE[, .(q1=quantile(cpm.multimap.adjust, 0.25), q3=quantile(cpm.multimap.adjust, 0.75)), by=.(lab.libMethod, lib.method.simple, lib.method.detail, miR.ID, all.gt0.within.lab, passes.cutoff.within.lab, passes.cutoff.between.lab.simple, passes.cutoff.between.lab.detail)]
QCD.WI.LABS.plasma[, `:=`(iqr=(q3-q1)/2, midhinge=(q1+q3)/2)]
QCD.WI.LABS.plasma[, intralab.qcd:=iqr/midhinge]
ggplot(QCD.WI.LABS.plasma[lib.method.detail!="4N_NEXTflex" & passes.cutoff.between.lab.simple==TRUE],
       aes(x = lab.libMethod,
           y = intralab.qcd,
           fill=lib.method.simple))  +
  geom_violin(
    draw_quantiles = c(0.25, 0.5, 0.75)
  ) + scale_y_continuous(limits = c(0,1.0), breaks = seq(0, 1, 0.1), expand=c(0.01,0)) +
  theme(
    axis.text.x = element_text(angle=50, hjust=1),
    panel.grid.major = element_line(color = NA),
    panel.grid.minor = element_line(color = NA),
    panel.border = element_rect(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  ) +
  labs(#title = "EQUIMOLAR: QCD Within-labs",
    #x = "Library Prep Method",
    x=NULL,
    y = "QCD") 
ggsave(this.outfile, height=2.6, width=4.5)

# Cross-Lab %CV

# FIG5C--Plasma Pool Interlab %CV --------------------------------------------------------------

MEAN.CV.ACROSS.LABS2.plasma <- MEAN.CV.WI.LABS.plasma[(lib.method.detail=="4N_B" | lib.method.simple=="TruSeq" | lib.method.simple == "NEBNext"), .(cv.inter=sd(mean.cpm)/(sum(mean.cpm*n)/sum(n))), by=.(lib.method.simple, lib.method.detail, miR.ID, passes.cutoff.between.lab.detail)]
this.outfile <- paste0(outdirs["FIG5"], "/FIG5C_PLASMA_PERCENT_CV_INTERLAB_TS_NEB_4N_B.pdf")
ggplot(MEAN.CV.ACROSS.LABS2.plasma[passes.cutoff.between.lab.detail==TRUE],
       aes(x = lib.method.detail,
           y = cv.inter * 100))  +
  geom_boxplot(
    aes(fill = lib.method.simple),
    alpha = 0.25, 
    outlier.color = NA,
    width = 0.9
  ) + 
  scale_fill_manual(values = ann_colors$lib.method.simple) +
  geom_jitter(
    size = rel(0.001),
    alpha = 0.25,
    width = 0.4,
    height = 0
  ) + scale_y_continuous(limits = c(0, 150), breaks = seq(0, 200, 25)) +
  theme(
    panel.grid.major = element_line(color = NA),
    panel.grid.minor = element_line(color = NA),
    panel.border = element_rect(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  ) + 
  labs(#title = "EQUIMOLAR:\n%CV Across Labs",
    x=NULL,
    y = "%CV") 
ggsave(this.outfile, width = 3.0, height=2.25, units = "in")


# FIG5C--Plasma QCD  4N_B--------------------------------------------------------------
QCD.ACROSS.LABS2.plasma <- MEAN.CV.WI.LABS.plasma[(lib.method.detail=="4N_B" | lib.method.simple=="TruSeq" | lib.method.simple == "NEBNext"), .(q1=quantile(mean.cpm, 0.25), q3=quantile(mean.cpm, 0.75)), by=.(lib.method.simple, lib.method.detail, miR.ID, passes.cutoff.between.lab.detail)]
QCD.ACROSS.LABS2.plasma[, `:=`(iqr=(q3-q1)/2, midhinge=(q1+q3)/2)]
QCD.ACROSS.LABS2.plasma[, qcd:=iqr/midhinge]
this.outfile <- paste0(outdirs["FIG5"], "/FIG5C_PLASMA_QCD_INTERLAB_TS_NEB_4N_B.pdf")

ggplot(QCD.ACROSS.LABS2.plasma[passes.cutoff.between.lab.detail==TRUE],
       aes(x = lib.method.detail,
           y = qcd))  +
  geom_boxplot(
    aes(fill = lib.method.simple),
    alpha = 0.25, 
    outlier.color = NA,
    width = 0.9
  ) + 
  scale_fill_manual(values = ann_colors$lib.method.simple) +
  geom_jitter(
    size = rel(0.001),
    alpha = 0.25,
    width = 0.4,
    height = 0
  ) + scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.1)) +
  theme(
    panel.grid.major = element_line(color = NA),
    panel.grid.minor = element_line(color = NA),
    panel.border = element_rect(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  ) + 
  labs(#title = "EQUIMOLAR:\n%CV Across Labs",
    x=NULL,
    y = "QCD") 
ggsave(this.outfile, width = 3.0, height=2.25, units = "in")

# INTERLAB CV/QCD Plasma Summary Table -------------------------------
MEAN.CV.ACROSS.LABS.plasma.subgroups <- MEAN.CV.WI.LABS.plasma.subgroups[n>2, .(mean.cpm=mean(mean.cpm), sd.cpm=sd(mean.cpm), q1.cpm=quantile(mean.cpm, 0.25), q3.cpm=quantile(mean.cpm, 0.75), n=.N), by=.(lib.method.group, miR.ID)]
MEAN.CV.ACROSS.LABS.plasma.subgroups[, `:=`(interlab.cv=100*sd.cpm/mean.cpm, interlab.qcd=((q3.cpm-q1.cpm)/2)/((q3.cpm+q1.cpm)/2))]
MEAN.CV.ACROSS.LABS.plasma.subgroups[, (sub("^0.", "interlab.cv", summary.quantiles)):=lapply(summary.quantiles, USE.NAMES=FALSE, SIMPLIFY=FALSE, quantile, x=interlab.cv), by=.(lib.method.group)]
MEAN.CV.ACROSS.LABS.plasma.subgroups[, (sub("^0.", "interlab.qcd", summary.quantiles)):=lapply(summary.quantiles, USE.NAMES=FALSE, SIMPLIFY=FALSE, quantile, x=interlab.qcd), by=.(lib.method.group)]

INTERLAB.CV.PLASMA.SUMMARY <- subset(unique(MEAN.CV.ACROSS.LABS.plasma.subgroups, by=c("lib.method.group", "n")), select=c("lib.method.group", "n", sub("^0.", "interlab.cv", summary.quantiles), sub("^0.", "interlab.qcd", summary.quantiles)))
this.outfile <- paste0(outdirs["tables"], "/Plasma_interlab_CV_QCD_summary.txt")
write.table(INTERLAB.CV.PLASMA.SUMMARY, this.outfile, row.names = FALSE, col.names=TRUE, sep="\t", quote=FALSE)

# IN-TEXT Plasma INTER-lab CV/QCD summary numbers----
setorder(INTERLAB.CV.PLASMA.SUMMARY, interlab.cv5)
paste(
  INTERLAB.CV.PLASMA.SUMMARY[, 
                                paste0(lib.method.group,
                                       ": ", 
                                       round(interlab.cv5, digits = 2),
                                       " (", 
                                       round(interlab.cv02, digits = 2), 
                                       ", ", 
                                       round(interlab.cv98, digits = 2),
                                       "; n=",
                                       n,  
                                       ")")], 
  collapse="; ")
# 4N_B: 25.66 (6.72, 68.17; n=4); NEBNext: 29.17 (9.71, 66.35; n=4); TruSeq: 32.87 (12.03, 110.34; n=6); 4N: 47.46 (21.29, 95.67; n=7)

setorder(INTERLAB.CV.PLASMA.SUMMARY, interlab.qcd5)
paste(
  INTERLAB.CV.PLASMA.SUMMARY[, 
                                paste0(lib.method.group,
                                       ": ", 
                                       round(interlab.qcd5, digits = 2),
                                       " (", 
                                       round(interlab.qcd02, digits = 2), 
                                       ", ", 
                                       round(interlab.qcd98, digits = 2),
                                       "; n=",
                                       n,  
                                       ")")], 
  collapse="; ")
# 4N_B: 0.11 (0.03, 0.3; n=4); NEBNext: 0.15 (0.04, 0.44; n=4); TruSeq: 0.19 (0.04, 0.61; n=6); 4N: 0.26 (0.06, 0.75; n=7)

# DOWNSAMPLING: Plasma pool ------------------------------
dge.plasma <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, colnames(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA)%in%colnames(LOG.CPM.mirs.passing.100cpm.by.method)]
dge.plasma.counts <- dge.plasma$counts
# SO REPEAT SAME PARAMS FOR EQUIMOLAR
downsample.to.vec <- c(Inf, round(10^seq(4, 5.5, 0.5), 0))
drr.list.plasma <- sapply(downsample.to.vec, simplify=FALSE, function(x) t(drarefy(t(dge.plasma.counts), x)))

names(drr.list.plasma) <- downsample.to.vec
drr.all.plasma <- data.table(rbindlist(lapply(drr.list.plasma, melt), idcol=TRUE, use.names=TRUE), keep.rownames=TRUE)
setnames(drr.all.plasma, c("downsample.to", "miR.id", "lab.libMethod.replicate", "est.prob.detected"))
drr.all.plasma[, lab.libMethod.replicate:=sub("PlasmaPool.", "", sub("^X", "", as.character(lab.libMethod.replicate)))]
setkey(drr.all.plasma, lab.libMethod.replicate)
sample.info.plasma.dt <- data.table(dge.plasma$samples, keep.rownames=FALSE)
setkey(sample.info.plasma.dt, lab.libMethod.replicate)
drr.all.plasma.merge <- merge(drr.all.plasma, sample.info.plasma.dt)
drr.all.plasma.merge[, downsample.to:=as.numeric(downsample.to)]
drr.all.plasma.filt <- subset(drr.all.plasma.merge, lib.size>=downsample.to|lib.size==Inf)

# plot cutoff detected
min.prob.detected=0.90
#drr.all.filt.summarize.repl <- drr.all.filt[, .(est.prob.detected.all=prod(est.prob.detected)), by=.(lab.libMethod, downsample.to, Lab, lib.method.detail, lib.method.simple,  miR.id)]
#est.prob.detected.n.miRs <- drr.all.filt.summarize.repl[, .(n.detected=sum(ifelse(est.prob.detected.all>=min.prob.detected, 1, 0)), n.not.detected=sum(ifelse(est.prob.detected.all>=min.prob.detected, 0, 1)), n.total=.N), by=.(lab.libMethod, downsample.to, Lab, lib.method.detail, lib.method.simple)]
est.prob.detected.n.miRs <- drr.all.plasma.filt[, .(n.detected=sum(ifelse(est.prob.detected>=min.prob.detected, 1, 0)), n.not.detected=sum(ifelse(est.prob.detected>=min.prob.detected, 0, 1)), n.total=.N), by=.(lab.libMethod.replicate, lab.libMethod, downsample.to, Lab, lib.method.detail, lib.method.simple)]

this.outfile <- paste0(outdirs["FIG5"], "/FIG5D_miRs_detected_plasmaPool.pdf")
g <- ggplot(est.prob.detected.n.miRs,
       aes(
         x=lib.method.simple,
         y=n.detected, fill=lib.method.simple)) + 
  geom_boxplot() + 
  facet_wrap(~downsample.to, strip.position = "top",  nrow = 1) +
  labs(x=NULL,
    y="# miRNAs detected in all samples",
    strip="log10 Sequencing Depth (total miRNA-mapping reads)") +
  scale_y_continuous(breaks = seq(0,1000,100)) +
  theme(panel.grid.minor = element_blank(),
        axis.text=element_text(color="black"),
        axis.text.x = element_text(angle=50, hjust=1, color="black"),
        axis.line=element_line(color="black"),
        panel.grid.major = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text=element_text(size=8, color="black"),
        axis.ticks = element_line(color="black"),
        legend.position = "none") ; g
ggsave(this.outfile, g, width = 7, height = 2); g

# Scaling factors ----
# Redo 9-14-17 ----
# Remove the two NEBNext labs that didn't do the same adapter dilutions, as we did for the inter-lab QCD
ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.drop.NEB.1to2.dilutions <- ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE[, ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE$samples$lab.libMethod!="NEBNext.Lab9" & ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE$samples$lab.libMethod!="NEBNext.Lab3"]

# Subset matrices for plasma and equimolar
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE <- ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.drop.NEB.1to2.dilutions[ , ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.drop.NEB.1to2.dilutions$samples$pool=="PlasmaPool"]
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE[rowSums(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE$counts>0)>0,]
ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR <- ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.drop.NEB.1to2.dilutions[ , ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.drop.NEB.1to2.dilutions$samples$pool=="SynthEQ"]
ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[rowSums(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$counts>0)>0,]
lib.detail.for.comp <- levels(factor(unique(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.detail), levels=lib.detail.levels[!lib.detail.levels%in%c("4N_Xu", "4N_NEXTflex")]))
lib.detail.matr <- matrix(data = 1, nrow = length(lib.detail.for.comp), ncol= length(lib.detail.for.comp))
dimnames(lib.detail.matr) <- list(lib.detail.for.comp, lib.detail.for.comp)
comparisons.detail <- data.frame(melt(lib.detail.matr)[melt(upper.tri(lib.detail.matr, diag = FALSE))[,3],1:2], stringsAsFactors = FALSE)
row.names(comparisons.detail) <- NULL

# EQ Pool: Calculate scaling factors for all pairwise comparisons-----
ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG <- do.call(rbind, 
                                                             sapply(
                                                               1:dim(comparisons.detail)[1], 
                                                               USE.NAMES=FALSE, 
                                                               simplify = FALSE,
                                                               FUN=function(x.row){
                                                               x.ids <- c(as.character(comparisons.detail[x.row,]$Var1), as.character(comparisons.detail[x.row,]$Var2))
                                                               print(x.ids)
                                                               x.id.list <- paste0(x.ids[1], ".V.", x.ids[2])
                                                               dge.comp <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[, ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.detail%in%x.ids]
                                                               design.comp <- model.matrix(~0+lab.libMethod, dge.comp$samples)
                                                               if(length(unique(colnames(design.comp)))<2){
                                                                 warning(paste0("skipping ", x.ids[1], " ", x.ids[2]))
                                                                 return(NULL)
                                                               }
                                                               colnames(design.comp) <- sub("lab.libMethod", "", colnames(design.comp))
                                                               lib.methods.comp <- strsplit2(colnames(design.comp), split=".", fixed=TRUE)[,1]
                                                               
                                                               comp.contrasts.fodder <- sapply(unique(lib.methods.comp), USE.NAMES=FALSE, FUN=function(lib.m){
                                                                 mat <- design.comp
                                                                 y.sum <- rowSums(apply(mat, 1, FUN=function(x){
                                                                   x*ifelse(lib.methods.comp==lib.m, 1, 0)
                                                                 }))
                                                                 y.sum/sum(y.sum)
                                                               })
                                                               colnames(comp.contrasts.fodder) <- unique(lib.methods.comp)
                                                               out.contrast <- as.numeric(comp.contrasts.fodder[, x.ids[1]]-comp.contrasts.fodder[, x.ids[2]])
                                                               # Now filter miRs
                                                               dge.comp.filt <- dge.comp[rowSums(dge.comp$counts>0)==dim(dge.comp)[2],]
                                                               # v.comp <- voomWithQualityWeights(dge.comp.filt, design.comp, normalize.method = "scale")
                                                               # Drop the quality weighting, just for simplicity
                                                               dge.comp.filt <- calcNormFactors(dge.comp.filt, method = "RLE")
                                                               v.comp <- voom(dge.comp.filt, design.comp, normalize.method = "scale")
                                                               fit.comp <- lmFit(v.comp, design.comp)
                                                               fit.comp <- eBayes(fit.comp, robust = TRUE, trend = TRUE)
                                                               fit2.comp <- contrasts.fit(fit.comp, out.contrast)
                                                               fit2.comp <- eBayes(fit2.comp, trend = TRUE, robust = TRUE)
                                                               tt.comp <- data.frame(topTable(fit2.comp, n = Inf, confint=0.95, sort.by = "none"))
                                                               tt.comp.dt <- data.table(data.frame(miR.id.full=row.names(tt.comp), tt.comp))
                                                               tt.comp.dt[, comparison:=x.id.list]
                                                               return(tt.comp.dt)
                                                             })
)



# save the data from the original correction factor calculation input so that we can try alternatives
#save(file = "20170810_CrossU01_ScalingFactor_Calculation_Input_Data.RData",
#     list = Hmisc::Cs(comparisons.detail,
 #                     ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR,
  #                    downbias.labLibMethod,
   #                   ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA,
    #                  CROSS.U01.PLASM.METADATA,
     #                 CROSS.U01.SYNTH.METADATA,
      #                CROSS.U01.SYNTH.METADATA.WITH.MIR.SENSE
    # ))

ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG[, c("from.method", "to.method"):=tstrsplit(comparison, split=".V.")]
ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset <- subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG, select=c("miR.id.full", "comparison", "from.method", "to.method", "logFC", "CI.L", "CI.R"))
setnames(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, c("CI.L", "CI.R", "comparison"), c("CI.95L", "CI.95R", "correction"))
ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset[, from.method:=sub("^X4N", "4N", from.method)]
ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset[, to.method:=sub("^X4N", "4N", to.method)]
ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset[, correction:=sub("^X4N", "4N", correction)]
this.outfile <- paste0(outdirs["tables"],"/20170914_INTER-METHOD_Correction_Factors_equimolar_pool_voom_limma.txt")
write.table(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, this.outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# Now remove clean tag & 4N sub-methods from plots to show
plot.methods <- c("TruSeq", "NEBNext", "4N_B")
ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset.FOR.PLOTTING <- subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, from.method%in%plot.methods & (to.method%in%plot.methods | to.method=="DOWNBIAS"))
lib.comparison.groups.for.plot <- unique(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset.FOR.PLOTTING$correction)

# Apply Correction factors to plasma data----
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.QUANTILENORM <- do.call(rbind,  sapply(lib.comparison.groups.for.plot, USE.NAMES=FALSE, simplify = FALSE, FUN=function(x.group){
  
  correct.dt <- data.frame(subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset.FOR.PLOTTING, correction==x.group), row.names=1)
  x.ids <- strsplit2(x.group, split=".V.")
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.detail%in%x.ids]
  keep.mirs <- (rowSums(dge.test$counts>=5)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  v.test.norm <- voom(dge.test.filt, design.test, normalize.method="scale")
  cpm.test.norm <- v.test.norm$E
  
  merged.miR.ids <- merge(cpm.test.norm, correct.dt, by=0)$Row.names
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  cpm.test.norm.filt.md1 <- cpm.test.norm.filt[, grep(x.ids[1], colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)]
  tt.comp.scl.factor.filt <- correct.dt[merged.miR.ids,]
  cpm.test.norm.filt.md.1.scaled <-
    apply(
      cpm.test.norm.filt.md1,
      2,
      FUN = function(x) {
        x.log.scl <- x - tt.comp.scl.factor.filt$logFC
      })
  colnames(cpm.test.norm.filt) <-
    paste0(colnames(cpm.test.norm.filt), ".orig")
  colnames(cpm.test.norm.filt.md.1.scaled) <-
    paste0(colnames(cpm.test.norm.filt.md.1.scaled), ".transf")
  cpm.test.orig.and.scaled <-
    cbind(cpm.test.norm.filt, cpm.test.norm.filt.md.1.scaled)
  col.names.combined <- colnames(cpm.test.orig.and.scaled)
  group.des.df <- ifelse(grepl("transf$", col.names.combined), paste0(x.ids[1], ".to.", x.ids[2]), ifelse(grepl(x.ids[1], col.names.combined), paste0(x.ids[1], ".orig"), paste0(x.ids[2], ".orig")))
  group.design <- model.matrix(~0+group.des.df)
  colnames(group.design) <- sub("group.des.df", "", colnames(group.design))
  cpm.test.orig.and.scaled.rescaled <-
    voom(2^cpm.test.orig.and.scaled, group.design, normalize.method=ifelse(x.ids[2]=="DOWNBIAS", "scale", "quantile"))$E
  ids.reference <- if(x.ids[2]=="DOWNBIAS"){
    paste0(row.names(dge.test.filt$samples[dge.test.filt$samples$GROUND.TRUTH.SAMPLE=="DOWNBIAS",]), ".orig")
  } else{
    grep(
      x.ids[1],
      colnames(cpm.test.orig.and.scaled.rescaled),
      value = TRUE,
      ignore.case = TRUE,
      invert = TRUE
    )
  }
  
  ids.test.orig <-
    grep(x.ids[1],
         colnames(cpm.test.norm.filt),
         value = TRUE,
         ignore.case = TRUE)
  ids.test.scaled <- colnames(cpm.test.norm.filt.md.1.scaled)
  # Return this to get a matrix of transformed and original values
  # return(cpm.test.orig.and.scaled.rescaled)
  # Now get an average of the values only did intra-lab scaling
  cpm.rescaled.test.orig <-
    cpm.test.orig.and.scaled.rescaled[, ids.test.orig]
  cpm.rescaled.test.scaled <-
    cpm.test.orig.and.scaled.rescaled[, ids.test.scaled]
  cpm.rescaled.ref.orig <-
    cpm.test.orig.and.scaled.rescaled[, ids.reference]
  cpm.mean.reference <-
    rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.reference])
  
  # Check scale norm too. Didnt look quite as good, so removed.
  cpm.test.orig.and.scaled.rescaled.med <-
    voom(2^cpm.test.orig.and.scaled, design = group.design, normalize.method = "scale")$E
  cpm.rescaled.test.orig.med <-
    cpm.test.orig.and.scaled.rescaled.med[, ids.test.orig]
  cpm.rescaled.test.scaled.med <-
    cpm.test.orig.and.scaled.rescaled.med[, ids.test.scaled]
  cpm.rescaled.ref.orig.med <-
    cpm.test.orig.and.scaled.rescaled.med[, ids.reference] 
  cpm.mean.reference.med <-
    rowMeans(cpm.test.orig.and.scaled.rescaled.med[, ids.reference])
  
  
  cpm.test.orig.and.scaled.rescaled.mean <-
    data.table(
      data.frame(
        miRNA.full.ID = row.names(cpm.test.orig.and.scaled.rescaled),
        scaling.factor.logFC =
          tt.comp.scl.factor.filt$logFC,
        cpm.mean.reference = cpm.mean.reference,
        cpm.mean.test.orig = rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.test.orig]),
        cpm.mean.test.scaled = rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.test.scaled]),
        abs.dif.orig = rowMeans(apply(
          cpm.rescaled.test.orig,
          2,
          FUN = function(x) {
            abs(x - cpm.mean.reference)
          }
        )),
        abs.dif.transf = rowMeans(apply(
          cpm.rescaled.test.scaled,
          2,
          FUN = function(x) {
            abs(x - cpm.mean.reference)
          }
        ))
      )
    )
  cpm.test.orig.and.scaled.rescaled.mean[, comparison := x.group]
  return(cpm.test.orig.and.scaled.rescaled.mean)
  #return(cpm.test.orig.and.scaled.rescaled)
}))
#ggplot(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.LIST[[3]]) + geom_density(aes(color="VS.ORIG", x=abs.dif.orig)) + geom_density(aes(color="VS.TRANSF", x=abs.dif.transf))

ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif <- melt(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.QUANTILENORM, measure.vars=c("abs.dif.orig", "abs.dif.transf"), variable.name = "ORIG.or.TRANSF.and.SCL", value.name = "mean.abs.diff")
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif[, SCALE.METHOD:="VOOM.NO.INTERLAB"]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif[, ORIG.or.TRANSF:=sub(paste0(".", SCALE.METHOD), "", ORIG.or.TRANSF.and.SCL), by=1:nrow(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif)]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif[, `:=`(SCALE.METHOD=toupper(SCALE.METHOD), ORIG.or.TRANSF=ifelse(ORIG.or.TRANSF=="abs.dif.orig", "ORIGINAL", "TRANSFORMED"))]
#ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile <- subset(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif, SCALE.METHOD=="QUANTILE")
ggplot(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif, aes(x=mean.abs.diff, color=ORIG.or.TRANSF)) + geom_density(size=1.5) + facet_wrap(~comparison) + theme_bw() + 
  scale_y_continuous(expand=c(0.0, 0.0), breaks=seq(0.0, 1.5, 0.1), limits=c(0.0, 1.4)) + 
  scale_x_continuous(expand=c(0.01,0), limits=c(-0.1, 9.0), breaks=seq(0,9, 1.0)) + 
  theme(
    axis.text=element_text(size=14, face="bold", color="black"),
    axis.title=element_text(face="bold", size=16),
    plot.title=element_text(face="bold", size=18),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=16, face="bold", color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=14, face="bold"),
    legend.title=element_text(size=15, face="bold"),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(y="probability density",
       x="|Difference VS Reference|(log2)")

# Get the matrices for plotting
ALL.COMPARISON.SCALING.FACTOR.CORRECTED.MATR.LIST.UNRESCALED <- sapply(lib.comparison.groups.for.plot, USE.NAMES=FALSE, simplify = FALSE, FUN=function(x.group){
  
  correct.dt <- data.frame(subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset.FOR.PLOTTING, correction==x.group), row.names=1)
  x.ids <- strsplit2(x.group, split=".V.")
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.detail%in%x.ids]
  keep.mirs <- (rowSums(dge.test$counts>=1)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  v.test.norm <- voom(dge.test.filt, design.test, normalize.method="scale")
  cpm.test.norm <- v.test.norm$E
  
  merged.miR.ids <- merge(cpm.test.norm, correct.dt, by=0)$Row.names
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  cpm.test.norm.filt.md1 <- cpm.test.norm.filt[, grep(x.ids[1], colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)]
  tt.comp.scl.factor.filt <- correct.dt[merged.miR.ids,]
  cpm.test.norm.filt.md.1.scaled <-
    apply(
      cpm.test.norm.filt.md1,
      2,
      FUN = function(x) {
        x.log.scl <- x - tt.comp.scl.factor.filt$logFC
      })
  colnames(cpm.test.norm.filt) <-
    paste0(colnames(cpm.test.norm.filt), ".orig")
  colnames(cpm.test.norm.filt.md.1.scaled) <-
    paste0(colnames(cpm.test.norm.filt.md.1.scaled), ".transf")
  cpm.test.orig.and.scaled <-
    cbind(cpm.test.norm.filt, cpm.test.norm.filt.md.1.scaled)
})
names(ALL.COMPARISON.SCALING.FACTOR.CORRECTED.MATR.LIST.UNRESCALED) <- lib.comparison.groups.for.plot
library(scales)
ann.colors.plasm = list(
  from.method = ann_colors$lib.method.detail,
  to.method = ann_colors$lib.method.detail,
  is.corrected = grey_pal()(2)
)

# Plot heatmaps Plasma Pool Corrected VS Uncorrected ----
sapply(lib.comparison.groups.for.plot, simplify = FALSE, FUN=function(x){
  this.comparison <- x
  from.method.id <-  strsplit2(this.comparison, split=".V.", fixed=TRUE)[,1]
  to.method.id <- strsplit2(this.comparison, split=".V.", fixed=TRUE)[,2]
  this.matr <- ALL.COMPARISON.SCALING.FACTOR.CORRECTED.MATR.LIST.UNRESCALED[[x]]
  this.sample.info <- data.frame(strsplit2(colnames(this.matr), split=".", fixed=TRUE), stringsAsFactors = FALSE, row.names=colnames(this.matr))
  colnames(this.sample.info) <- c("lib.method.detail", "Lab", "pool", "replicate", "orig.or.transf")
  this.sample.info$lib.method.detail <- sub("^X", "", this.sample.info$lib.method.detail)
  this.sample.info$lab.libMethod.pool.replicate <- apply(this.sample.info[, c("lib.method.detail", "Lab", "pool", "replicate")], 1, paste, collapse=".")
  this.sample.info$lab.libMethod.pool.replicate.transf <- row.names(this.sample.info)
  this.sample.info <- data.table(this.sample.info)
  this.sample.info.full1 <- data.table(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples)
  merge.cols <- intersect(colnames(this.sample.info.full1), colnames(this.sample.info))

  setkeyv(this.sample.info.full1, merge.cols)
  setkeyv(this.sample.info, merge.cols)
 
  this.sample.info.full <-  merge(this.sample.info, this.sample.info.full1)
  this.sample.info.full[, is.corrected:=ifelse(orig.or.transf=="orig", "ORIGINAL", "CORRECTED")]
  this.sample.info.full[, transf.LibMethod.id:=sub("X4N", "4N", sub(".PlasmaPool", "", make.names(paste(lab.libMethod.pool.replicate, substr(is.corrected, 1, 4), sep="."), unique = TRUE)))]
  this.sample.info.full[, `:=`(from.method=lib.method.detail, to.method=ifelse(is.corrected=="ORIGINAL", as.character(lib.method.detail), to.method.id))]
  this.sample.info.full.df <- data.frame(this.sample.info.full)
  row.names(this.sample.info.full.df) <- make.names(this.sample.info.full$lab.libMethod.pool.replicate.transf)
  this.sample.info.full.df <- this.sample.info.full.df[colnames(this.matr), ]
  this.sample.info.full.df.for.plot <- this.sample.info.full.df[, c("to.method", "from.method", "transf.LibMethod.id")]
  this.sample.info.full.df.for.plot$from.method <- factor(this.sample.info.full.df.for.plot$from.method, levels=lib.detail.levels)
  this.sample.info.full.df.for.plot$to.method <- factor(this.sample.info.full.df.for.plot$to.method, levels=lib.detail.levels)
  groups <- apply(this.sample.info.full.df.for.plot[,1:2], 1, FUN=function(x) ifelse(x[1]==x[2], paste0(x[1], ".orig"), paste0(x[1], ".to.", x[2])))
  design.scl <- model.matrix(~0+groups)
  this.matr.RESCALE <- voom(2^this.matr, design.scl, normalize.method="quantile", plot=FALSE)$E
  
  row.mean.abund.df <- data.frame(FUNCTION.aveTechReps(this.matr.RESCALE, groups))
  this.matr.RESCALE.sort <- this.matr.RESCALE[order(-row.mean.abund.df[, paste0(from.method.id, ".orig")]),]
  
  this.outfile <- paste0(outdirs["supplemental_figures"], "/FIG_S8_", this.comparison, "_Correction_Heatmap_PLASMAPOOL.pdf")
  pheatmap(this.matr.RESCALE.sort, height = 7.0, width=8.5,
           cluster_rows = FALSE,
           show_rownames = FALSE,
           fontsize = 8,
           fontsize_col = 8,
           treeheight_col = 15,
           labels_col = this.sample.info.full.df.for.plot$transf.LibMethod.id,
           annotation_col = this.sample.info.full.df.for.plot[, 1:2], 
           annotation_row = row.mean.abund.df,
           annotation_colors=ann.colors.plasm, 
           filename = this.outfile)
  
})
                                                                                
# REPEAT FOR EQUIMOLAR POOL 12-11-16----

ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.QUANTILENORM.EQ <- do.call(rbind,  sapply(lib.comparison.groups.for.plot, USE.NAMES=FALSE, simplify = FALSE, FUN=function(x.group){
  correct.dt <- data.frame(subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset.FOR.PLOTTING, correction==x.group), row.names=1)
  x.ids <- strsplit2(x.group, split=".V.")
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[, (ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.detail%in%x.ids)]
  keep.mirs <- (rowSums(dge.test$counts>=1)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  v.test.norm <- voom(dge.test.filt, design.test, normalize.method="scale")
  cpm.test.norm <- v.test.norm$E
  
  merged.miR.ids <- merge(cpm.test.norm, correct.dt, by=0)$Row.names
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  cpm.test.norm.filt.md1 <- cpm.test.norm.filt[, grep(x.ids[1], colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)]
  tt.comp.scl.factor.filt <- correct.dt[merged.miR.ids,]
  cpm.test.norm.filt.md.1.scaled <-
    apply(
      cpm.test.norm.filt.md1,
      2,
      FUN = function(x) {
        x.log.scl <- x - tt.comp.scl.factor.filt$logFC
      })
  colnames(cpm.test.norm.filt) <-
    paste0(colnames(cpm.test.norm.filt), ".orig")
  colnames(cpm.test.norm.filt.md.1.scaled) <-
    paste0(colnames(cpm.test.norm.filt.md.1.scaled), ".transf")
  cpm.test.orig.and.scaled <-
    cbind(cpm.test.norm.filt, cpm.test.norm.filt.md.1.scaled)
  col.names.combined <- colnames(cpm.test.orig.and.scaled)
  group.des.df <- ifelse(grepl("transf$", col.names.combined), paste0(x.ids[1], ".to.", x.ids[2]), ifelse(grepl(x.ids[1], col.names.combined), paste0(x.ids[1], ".orig"), paste0(x.ids[2], ".orig")))
  group.design <- model.matrix(~0+group.des.df)
  colnames(group.design) <- sub("group.des.df", "", colnames(group.design))
  cpm.test.orig.and.scaled.rescaled <-
    voom(2^cpm.test.orig.and.scaled, group.design, normalize.method=ifelse(x.ids[2]=="DOWNBIAS", "scale", "quantile"))$E
  ids.reference <- if(x.ids[2]=="DOWNBIAS"){
    grep(
      downbias.labLibMethod,
      colnames(cpm.test.orig.and.scaled.rescaled),
      value=TRUE,
      ignore.case=TRUE
    )
  } else{
    grep(
      x.ids[1],
      colnames(cpm.test.orig.and.scaled.rescaled),
      value = TRUE,
      ignore.case = TRUE,
      invert = TRUE
    )
  }
  
  ids.test.orig <-
    grep(x.ids[1],
         colnames(cpm.test.norm.filt),
         value = TRUE,
         ignore.case = TRUE)
  ids.test.scaled <- colnames(cpm.test.norm.filt.md.1.scaled)
  # Return this to get a matrix of transformed and original values
  # return(cpm.test.orig.and.scaled.rescaled)
  # Now get an average of the values only did intra-lab scaling
  cpm.rescaled.test.orig <-
    cpm.test.orig.and.scaled.rescaled[, ids.test.orig]
  cpm.rescaled.test.scaled <-
    cpm.test.orig.and.scaled.rescaled[, ids.test.scaled]
  cpm.rescaled.ref.orig <-
    cpm.test.orig.and.scaled.rescaled[, ids.reference]
  cpm.mean.reference <-
    rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.reference])
  
  # Check scale norm too. Didnt look quite as good, so removed.
  cpm.test.orig.and.scaled.rescaled.med <-
    voom(2^cpm.test.orig.and.scaled, design = group.design, normalize.method = "scale")$E
  cpm.rescaled.test.orig.med <-
    cpm.test.orig.and.scaled.rescaled.med[, ids.test.orig]
  cpm.rescaled.test.scaled.med <-
    cpm.test.orig.and.scaled.rescaled.med[, ids.test.scaled]
  cpm.rescaled.ref.orig.med <-
    cpm.test.orig.and.scaled.rescaled.med[, ids.reference] 
  cpm.mean.reference.med <-
    rowMeans(cpm.test.orig.and.scaled.rescaled.med[, ids.reference])
  
  
  cpm.test.orig.and.scaled.rescaled.mean <-
    data.table(
      data.frame(
        miRNA.full.ID = row.names(cpm.test.orig.and.scaled.rescaled),
        scaling.factor.logFC =
          tt.comp.scl.factor.filt$logFC,
        cpm.mean.reference = cpm.mean.reference,
        cpm.mean.test.orig = rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.test.orig]),
        cpm.mean.test.scaled = rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.test.scaled]),
        abs.dif.orig = rowMeans(apply(
          cpm.rescaled.test.orig,
          2,
          FUN = function(x) {
            abs(x - cpm.mean.reference)
          }
        )),
        abs.dif.transf = rowMeans(apply(
          cpm.rescaled.test.scaled,
          2,
          FUN = function(x) {
            abs(x - cpm.mean.reference)
          }
        ))
      )
    )
  cpm.test.orig.and.scaled.rescaled.mean[, comparison := x.group]
  return(cpm.test.orig.and.scaled.rescaled.mean)
  #return(cpm.test.orig.and.scaled.rescaled)
}))

ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQ <- melt(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.QUANTILENORM.EQ, measure.vars=c("abs.dif.orig", "abs.dif.transf"), variable.name = "ORIG.or.TRANSF.and.SCL", value.name = "mean.abs.diff")
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQ[, SCALE.METHOD:="QUANTILE"]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQ[, ORIG.or.TRANSF:=sub(paste0(".", SCALE.METHOD), "", ORIG.or.TRANSF.and.SCL), by=1:nrow(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQ)]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQ[, `:=`(SCALE.METHOD=toupper(SCALE.METHOD), ORIG.or.TRANSF=ifelse(ORIG.or.TRANSF=="abs.dif.orig", "ORIGINAL", "TRANSFORMED"))]
#ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile <- subset(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif, SCALE.METHOD=="QUANTILE")
ggplot(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQ, aes(x=mean.abs.diff, color=ORIG.or.TRANSF)) + geom_density(size=1.5) + facet_wrap(~comparison) + theme_bw() + 
  scale_y_continuous(expand=c(0.0, 0.0), breaks=seq(0.0, 2.0, 0.2), limits=c(0.0, 1.6)) + 
  scale_x_continuous(expand=c(0.01,0), limits=c(-0.1, 8.0), breaks=seq(0,9, 1.0)) + 
  theme(
    axis.text=element_text(color="black"),
    axis.title=element_text(size=10),
    plot.title=element_text(size=10),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(color="black"),
    axis.ticks = element_line(color="black"),
    strip.text = element_text(size=10, color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=10),
    legend.title=element_text(size=10),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(y="probability density",
       x="|Difference VS Reference|(log2)")
# Get the matrices for plotting EQUIMOLAR
ALL.COMPARISON.SCALING.FACTOR.CORRECTED.MATR.LIST.UNRESCALED.EQ <- sapply(lib.comparison.groups.for.plot, USE.NAMES=FALSE, simplify = FALSE, FUN=function(x.group){
  
  correct.dt <- data.frame(subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset.FOR.PLOTTING, correction==x.group), row.names=1)
  x.ids <- strsplit2(x.group, split=".V.")
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[, ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.detail%in%x.ids ]
  keep.mirs <- (rowSums(dge.test$counts>=1)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  v.test.norm <- voom(dge.test.filt, design.test, normalize.method="scale")
  cpm.test.norm <- v.test.norm$E
  
  merged.miR.ids <- merge(cpm.test.norm, correct.dt, by=0)$Row.names
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  cpm.test.norm.filt.md1 <- cpm.test.norm.filt[, grep(x.ids[1], colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)]
  tt.comp.scl.factor.filt <- correct.dt[merged.miR.ids,]
  cpm.test.norm.filt.md.1.scaled <-
    apply(
      cpm.test.norm.filt.md1,
      2,
      FUN = function(x) {
        x.log.scl <- x - tt.comp.scl.factor.filt$logFC
      })
  colnames(cpm.test.norm.filt) <-
    paste0(colnames(cpm.test.norm.filt), ".orig")
  colnames(cpm.test.norm.filt.md.1.scaled) <-
    paste0(colnames(cpm.test.norm.filt.md.1.scaled), ".transf")
  cpm.test.orig.and.scaled <-
    cbind(cpm.test.norm.filt, cpm.test.norm.filt.md.1.scaled)
})
names(ALL.COMPARISON.SCALING.FACTOR.CORRECTED.MATR.LIST.UNRESCALED.EQ) <- lib.comparison.groups.for.plot
# Plot heatmaps Equimolar Pool Corrected VS Uncorrected ----
ALL.EQ.Correction.summary <- do.call(rbind, sapply(lib.comparison.groups.for.plot, simplify = FALSE, FUN=function(x){
  this.comparison <- x
  from.method.id <-  strsplit2(this.comparison, split=".V.", fixed=TRUE)[,1]
  to.method.id <- strsplit2(this.comparison, split=".V.", fixed=TRUE)[,2]
  this.matr <- ALL.COMPARISON.SCALING.FACTOR.CORRECTED.MATR.LIST.UNRESCALED.EQ[[x]]
  this.sample.info <- data.frame(strsplit2(colnames(this.matr), split=".", fixed=TRUE), stringsAsFactors = FALSE, row.names=colnames(this.matr))
  colnames(this.sample.info) <- c("lib.method.detail", "Lab", "pool", "replicate", "orig.or.transf")
  this.sample.info$lib.method.detail <- sub("^X", "", this.sample.info$lib.method.detail)
  this.sample.info$lab.libMethod.pool.replicate <- apply(this.sample.info[, c("lib.method.detail", "Lab", "pool", "replicate")], 1, paste, collapse=".")
  this.sample.info$lab.libMethod.pool.replicate.transf <- row.names(this.sample.info)
  this.sample.info <- data.table(this.sample.info)
  this.sample.info.full1 <- data.table(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples)
  merge.cols <- intersect(colnames(this.sample.info.full1), colnames(this.sample.info))
  
  setkeyv(this.sample.info.full1, merge.cols)
  setkeyv(this.sample.info, merge.cols)
  
  this.sample.info.full <-  merge(this.sample.info, this.sample.info.full1)
  this.sample.info.full[, is.corrected:=ifelse(orig.or.transf=="orig", "ORIGINAL", "CORRECTED")]
  this.sample.info.full[, transf.LibMethod.id:=sub("X4N", "4N", sub(".PlasmaPool", "", make.names(paste(lab.libMethod.pool.replicate, substr(is.corrected, 1, 4), sep="."), unique = TRUE)))]
  this.sample.info.full[, `:=`(from.method=lib.method.detail, to.method=ifelse(is.corrected=="ORIGINAL", as.character(lib.method.detail), to.method.id))]
  this.sample.info.full.df <- data.frame(this.sample.info.full)
  row.names(this.sample.info.full.df) <- make.names(this.sample.info.full$lab.libMethod.pool.replicate.transf)
  this.sample.info.full.df <- this.sample.info.full.df[colnames(this.matr), ]
  this.sample.info.full.df.for.plot <- this.sample.info.full.df[, c("to.method", "from.method", "transf.LibMethod.id")]
  this.sample.info.full.df.for.plot$from.method <- factor(this.sample.info.full.df.for.plot$from.method, levels=lib.detail.levels)
  this.sample.info.full.df.for.plot$to.method <- factor(this.sample.info.full.df.for.plot$to.method, levels=lib.detail.levels)
  groups <- apply(this.sample.info.full.df.for.plot[,1:2], 1, FUN=function(x) ifelse(x[1]==x[2], paste0(x[1], ".orig"), paste0(x[1], ".to.", x[2])))
  design.scl <- model.matrix(~0+groups)
  #this.matr.RESCALE <- voom(2^this.matr, design.scl, normalize.method="none", plot=FALSE)$E
  this.matr.RESCALE <- this.matr
  row.mean.abund.df <- data.frame(FUNCTION.aveTechReps(this.matr.RESCALE, groups))
  this.matr.RESCALE.sort <- this.matr.RESCALE[order(-row.mean.abund.df[, paste0(from.method.id, ".orig")]),]
  
  this.outfile <- paste0(outdirs["supplemental_figures"], "/FIG_S7_", this.comparison, "_Correction_Heatmap_EQUIMOLAR_POOL_WITH_LEGEND.pdf")
  pheatmap(this.matr.RESCALE.sort,
           cluster_rows = FALSE,
           show_rownames = FALSE,
           fontsize = 8,
           fontsize_col = 8,
           treeheight_col = 15,
           labels_col = this.sample.info.full.df.for.plot$transf.LibMethod.id,
           annotation_col = this.sample.info.full.df.for.plot[, 1:2], 
           annotation_row = row.mean.abund.df,
           annotation_colors=ann.colors.plasm, 
           filename = this.outfile)
  this.outfile <- paste0(outdirs["supplemental_figures"], "/FIG_S7_", this.comparison, "_Correction_Heatmap_EQUIMOLAR_POOL_NOLEGEND.pdf")
  pheatmap(this.matr.RESCALE.sort, legend = FALSE, annotation_legend = FALSE, height = 4, width = 7.5,
           cluster_rows = FALSE,
           show_rownames = FALSE,
           fontsize = 8,
           fontsize_col = 8,
           treeheight_col = 15,
           labels_col = sub(".SynthEQ", "", this.sample.info.full.df.for.plot$transf.LibMethod.id),
           annotation_col = this.sample.info.full.df.for.plot[, 1:2], 
           annotation_row = row.mean.abund.df,
           annotation_colors=ann.colors.plasm, 
           filename = this.outfile)
  
  if(to.method.id=="4N_B"){
    this.matr.RESCALE.sort.melt <- data.table(melt(this.matr.RESCALE.sort))
    setnames(this.matr.RESCALE.sort.melt,  c("miR.id.full", "sample.ID", "CPM"))
    this.matr.RESCALE.sort.melt[, c("lib.method.detail", "Lab", "Pool", "replicate", "correction"):=tstrsplit(sample.ID, split=".", fixed=TRUE)]
    this.matr.RESCALE.sort.melt[, lib.method.detail:=sub("^X", "", lib.method.detail)]
    this.matr.RESCALE.sort.melt[, lib.method.simple:=tstrsplit(lib.method.detail, split="_")[1]]
    this.matr.RESCALE.sort.melt[, `:=`(lib.method.detail=factor(lib.method.detail, lib.detail.levels), lib.method.simple=factor(lib.method.simple, lib.simple.levels))]
    this.matr.RESCALE.sort.melt.summary <-  this.matr.RESCALE.sort.melt[, .(ave.cpm=mean(CPM)), by=.(miR.id.full, lib.method.detail, lib.method.simple, Lab, Pool, correction)]
    this.matr.RESCALE.sort.melt.summary[, lab.libMethod:=paste(lib.method.detail, Lab, sep="."), by=1:nrow(this.matr.RESCALE.sort.melt.summary)]
    this.matr.RESCALE.sort.melt.summary[, lab.libMethod.correction:=paste(lib.method.detail, Lab, correction, sep="."), by=1:nrow(this.matr.RESCALE.sort.melt.summary)]
    lab.lib.levels <- unique(this.matr.RESCALE.sort.melt.summary$lab.libMethod.correction)
    lab.lib.levels.ordered <- c(grep("orig", grep(from.method.id , lab.lib.levels, value = TRUE), value=TRUE),
                                grep("transf", grep(from.method.id , lab.lib.levels, value = TRUE), value=TRUE),
                                grep(to.method.id , lab.lib.levels, value = TRUE))
    this.matr.RESCALE.sort.melt.summary[, lab.libMethod.correction:=factor(lab.libMethod.correction, levels=lab.lib.levels.ordered)]
    this.matr.RESCALE.sort.melt.summary[, correction:=factor(correction, levels=c("orig", "transf"))]
    this.matr.RESCALE.sort.melt.summary[, comparison:=this.comparison]
    return(this.matr.RESCALE.sort.melt.summary)
  } else{
    return(NULL)
  }
  
  }))

# FIG S7 Bias_reduction violin plot equimolar----
g1 <- ggplot(ALL.EQ.Correction.summary[comparison=="TruSeq.V.4N_B"],
             aes(x=lab.libMethod.correction, alpha=correction, fill=lib.method.detail, y=ave.cpm)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_alpha_discrete(range=c(1,0.25)) + 
  scale_y_continuous(breaks = seq(-4,30,2)) +
  scale_fill_manual(values = ann_colors$lib.method.detail) + 
  theme(axis.text.x=element_text(angle=50, hjust=1),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(color="black"),
  axis.ticks = element_line(color="black"),
  legend.position = "top",
  legend.direction = "horizontal"
) + labs(y="Mean CPM (log2)", x=NULL); g1

neb.vs.4n.correction.summary <- ALL.EQ.Correction.summary[comparison=="NEBNext.V.4N_B"]
unique.neb.levels <- neb.vs.4n.correction.summary[, unique(as.character(lab.libMethod.correction))]
unique.neb.levels.order <- c(grep("NEBNext.*.orig", unique.neb.levels, value = TRUE),
                             grep("NEBNext.*.transf", unique.neb.levels, value = TRUE),
                             grep("4N_B.*.orig", unique.neb.levels, value = TRUE))
neb.vs.4n.correction.summary[, lab.libMethod.correction:=factor(lab.libMethod.correction, levels=unique.neb.levels.order)]
g2 <- ggplot(neb.vs.4n.correction.summary, aes(x=lab.libMethod.correction, alpha=correction, fill=lib.method.detail, y=ave.cpm)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + scale_alpha_discrete(range=c(1,0.5)) + 
  scale_fill_manual(values = ann_colors$lib.method.detail) + 
  scale_y_continuous(breaks = seq(-4,30,2)) +
  theme(axis.text.x=element_text(angle=50, hjust=1),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(color="black"),
  axis.ticks = element_line(color="black"),
  legend.position = "top", 
  legend.direction = "horizontal"
) + labs(y="Mean CPM (log2)", x=NULL); g2
g.both <- plot_grid(plotlist = list(g1, g2), nrow=1, ncol=2, rel_widths = c(5, 3)); g.both
this.outfile <- paste0(outdirs[["supplemental_figures"]], "/FIGS7_Equimolar_Correction_Truseq_and_NEB_to_4N_biasReduce_violin.pdf")
save_plot(filename = this.outfile, plot = g.both, base_height = 4, base_width = 4, ncol = 2, nrow = 1)



# Get values for Mann-W test ---
all.intermethod.comparison.groups <- unique(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset$correction)
# UPDATE 1-10-17 -- Through several iterations, it looks like normalization issues aren't actually a big issue. 
# Filter for count>5 in all samples, and do the normalization up-front only on the sample we're applying the correction factors to
# No need for inter-lab normalization as it is.
MANN.WHITNEY.SHIFTS.PLASMA <- sapply(all.intermethod.comparison.groups, USE.NAMES=TRUE, simplify = FALSE, FUN=function(x.group){
  correct.dt <- data.frame(subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, correction==x.group), row.names=1)
  x.ids <- strsplit2(x.group, split=".V.")
  print(x.group)
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.detail%in%x.ids]
  keep.mirs <- (rowSums(dge.test$counts>=5)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  query.ids <- grep(x.ids[1], colnames(dge.test.filt), value=TRUE, ignore.case=TRUE)
  ref.ids <- grep(x.ids[1], invert = TRUE, colnames(dge.test.filt), value=TRUE, ignore.case=TRUE)
  
  dge.query <- dge.test.filt[, query.ids]
  #design.query <- model.matrix(~0+lab.libMethod, dge.query$samples)
  dge.query <- calcNormFactors(dge.query, method="RLE")
  #v.query <- voom(dge.query, normalize.method="quantile")
  
  dge.ref <- dge.test.filt[, ref.ids]
  #design.ref <- model.matrix(~0+lab.libMethod, dge.ref$samples)
  dge.ref <- calcNormFactors(dge.ref, method="RLE")
  #v.ref <- voom(dge.ref, normalize.method="quantile")
  
  merged.miR.ids <- row.names(dge.test.filt)[row.names(dge.test.filt)%in%row.names(correct.dt)]
  #cpm.test.norm.filt.md1 <- v.query$E[merged.miR.ids,]
  cpm.test.norm.filt.md1 <- cpm.DGEList(dge.query, log = TRUE, normalized.lib.sizes = TRUE)[merged.miR.ids,]
  #cpm.test.norm.filt.ref <- v.ref$E[merged.miR.ids,]
  cpm.test.norm.filt.ref <- cpm.DGEList(dge.ref, log = TRUE, normalized.lib.sizes = TRUE)[merged.miR.ids,]
  tt.comp.scl.factor.filt <- correct.dt[merged.miR.ids,]
  
  cpm.test.norm.filt.md.1.scaled <-
    apply(
      cpm.test.norm.filt.md1,
      2,
      FUN = function(x) {
        x.log.scl <- x - tt.comp.scl.factor.filt$logFC
      })
  # combine for re-centering
  cpm.scaled.and.ref <- cbind(cpm.test.norm.filt.ref, cpm.test.norm.filt.md.1.scaled)[, colnames(dge.test.filt)]
  cpm.test.scaled <- voom(2^cpm.scaled.and.ref, design.test, normalize.method="quantile")$E
  cpm.test.scaled.group.means <- data.frame(mean.scaled.query=rowMeans(cpm.test.scaled[, query.ids]), mean.ref=rowMeans(cpm.test.scaled[, ref.ids]), row.names=row.names(cpm.test.scaled))
  cpm.test.scaled.group.means$abs.log.dif <- with(cpm.test.scaled.group.means, abs(mean.scaled.query-mean.ref))
  cpm.test.scaled.group.means$comparison <- "CORRECTED"
  
  # re-do the normalization for the original data to "quantile" so that it's normalized the same way as the scaled data.
  cpm.unscaled.and.ref <- cbind(cpm.test.norm.filt.ref, cpm.test.norm.filt.md1)[, colnames(dge.test.filt)]
  v.test.norm <- voom(2^cpm.unscaled.and.ref, design.test, normalize.method="quantile")
  cpm.test.norm <- v.test.norm$E
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  dge.test.filt <- calcNormFactors(dge.test.filt)
  v.test.norm <- voom(dge.test.filt, design.test, normalize.method="quantile")
  cpm.test.norm <- v.test.norm$E
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  cpm.test.original.group.means <- data.frame(mean.unscaled.query=rowMeans(cpm.test.norm.filt[, query.ids]), mean.ref=rowMeans(cpm.test.norm.filt[, ref.ids]), row.names=row.names(cpm.test.norm.filt))
  cpm.test.original.group.means$abs.log.dif <- with(cpm.test.original.group.means, abs(mean.unscaled.query-mean.ref))
  cpm.test.original.group.means$comparison <- "UNCORRECTED"
  cpm.both <- data.table(keep.rownames = TRUE, rbind(cpm.test.original.group.means[, c("comparison", "abs.log.dif")], cpm.test.scaled.group.means[, c("comparison", "abs.log.dif")]))
  cpm.both[, `:=`(group.id=paste0(x.ids[1], ".V.", x.ids[2]), from.method=x.ids[1], to.method=x.ids[2])]
  w.test.out <- wilcox.test(x = cpm.test.original.group.means$abs.log.dif, y=cpm.test.scaled.group.means$abs.log.dif, alternative = c("greater"),  conf.int = TRUE)
  cpm.both[, `:=`(W=w.test.out$statistic, p.value=w.test.out$p.value, alternative="greater", estimate=w.test.out$estimate)]
  #get matrix to include as well
  colnames(cpm.test.norm.filt.md.1.scaled) <- paste0(colnames(cpm.test.norm.filt.md.1.scaled), ".transf")
  colnames(cpm.unscaled.and.ref) <- paste0(colnames(cpm.unscaled.and.ref), ".orig")
  cpm.full.matr.with.correct <- cbind(cpm.unscaled.and.ref[merged.miR.ids,], cpm.test.norm.filt.md.1.scaled)
  cpm.full.matr.with.correct.norm <- normalizeBetweenArrays(cpm.full.matr.with.correct, method = "quantile")
  return(list(cpm.both=cpm.both, cpm.scaled.matr=cpm.full.matr.with.correct.norm))
}); MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF <- do.call(rbind, lapply(MANN.WHITNEY.SHIFTS.PLASMA, FUN=function(x) x[[1]])); MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF[, comparison:=factor(comparison, levels=c("UNCORRECTED", "CORRECTED"))]


MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS <- unique(subset(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF, select=colnames(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF)[!colnames(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF)%in%c("rn", "comparison", "abs.log.dif")]))
setorder(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS, p.value)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF[, group.id:=factor(group.id, levels=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS$group.id)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS[, group.id:=factor(group.id, levels=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS$group.id)]

ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF, aes(x=abs.log.dif, color=comparison)) + 
  geom_density(size=1.5) + 
  geom_text(data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS, aes(label=paste0("p = ", format(p.value, digits=3)), x=8.5, y=0.6),inherit.aes = FALSE, check_overlap=TRUE) + facet_wrap(~group.id) + theme_bw() + scale_x_continuous(breaks=seq(0, 16, 2), expand=c(0,0), limits=c(-0.1,14)) + scale_y_continuous(breaks=seq(0, 1, 0.1), expand=c(0,0))

# Supp Fig 7 Plasma Denstity Plot Correction Factors  -------------------------------------------------------------------
g.s7.all <- ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF, aes(x=abs.log.dif, color=comparison)) + geom_density(size=1, adjust=1, trim=TRUE) + facet_wrap(~group.id) + theme_bw() + 
  scale_y_continuous(expand=c(0.0, 0.0), breaks=seq(0.0, 1, 0.1), limits=c(0.0, 0.8)) + 
  scale_x_continuous(expand=c(0.01,0), breaks=seq(0,9, 1.0)) + 
  theme(
    axis.text=element_text(color="black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color="black"),
    axis.ticks = element_line(color="black"),
    strip.text = element_text(size=10, color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(y="density",
       x="|Difference VS Reference|(log2)")
g1 <- g.s7.all %+% MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF[group.id =="TruSeq.V.4N_B"]; g1
g2 <- g.s7.all %+% MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF[group.id =="TruSeq.V.NEBNext"]; g2
g3 <- g.s7.all %+% MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF[group.id =="NEBNext.V.4N_B"]; g3
this.outfile <- paste0(outdirs["supplemental_figures"], "/FigS7_DensityPlot_PlasmaPool_TS.V.4NB.pdf"); this.outfile
ggsave(filename = this.outfile, plot = g1, units = "in", width=4, height=3)
this.outfile <- paste0(outdirs["supplemental_figures"], "/FigS7_DensityPlot_PlasmaPool_TS.V.NEB.pdf"); this.outfile
ggsave(filename = this.outfile, plot = g2, units = "in", width=4, height=3)
this.outfile <- paste0(outdirs["supplemental_figures"], "/FigS7_DensityPlot_PlasmaPool_NEB.V.4NB.pdf"); this.outfile
ggsave(filename = this.outfile, plot = g3, units = "in", width=4, height=3)

# Supp Fig 7 Plasma boxpot Plot Correction Factors with mann W -------------------------------------------------------------------
show.groups <- c("TruSeq.V.4N_B", "TruSeq.V.NEBNext", "NEBNext.V.4N_B")
g.s7.all.mannW <- ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF[group.id %in% show.groups], aes(y=abs.log.dif, x=comparison, fill=comparison)) +
  geom_boxplot() + 
  facet_wrap(~group.id) + 
  theme_bw() + 
  geom_text(
    data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS[group.id %in% show.groups], 
    aes(
      label=paste0("p = ", format(p.value, digits=3)),
      x=2, y=8),
    inherit.aes = FALSE, check_overlap=TRUE) +
  scale_y_continuous(breaks=seq(0.0, 8.0, 0.5), limits=c(0.0, 8)) + 
  #scale_x_continuous(expand=c(0.01,0), limits=c(0, 8.0), breaks=seq(0,9, 1.0)) + 
  theme(
    axis.text=element_text(size=10, color="black"),
    axis.title=element_text(size=10),
    plot.title=element_text(size=10),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=10, color="black"),
    #legend.position = "top",
    #legend.direction = "horizontal",
    legend.text=element_text(size=10),
    legend.title=element_text(size=10),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(x=NULL,
       y="|Difference VS Reference|(log2)"); g.s7.all.mannW
this.outfile <- paste0(outdirs["supplemental_figures"], "/FigS7_boxplot_aggregate_PlasmaPool.pdf"); this.outfile
ggsave(this.outfile, g.s7.all.mannW, units = "in", width=8, height=2)
groups.plasm <- unique(apply(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF[, c("to.method", "from.method")], 1, FUN=function(x) ifelse(x[1]==x[2], paste0(x[1], ".orig"), paste0(x[2], ".to.", x[1]))))

truseq.v.4N.matr.plasm <- MANN.WHITNEY.SHIFTS.PLASMA[["TruSeq.V.4N_B"]][[2]]
row.mean.abund.df.plasm <- data.frame(FUNCTION.aveTechReps(truseq.v.4N.matr.plasm, groups.plasm))
ann.colors.correction.facts <- list(from.method=ann_colors$lib.method.simple, to.method=ann_colors$lib.method.simple)


# OLD FIG6A -------------------------------------------------------------------
library(dendsort)
dcols.truseq.v.4N.matr.plasm <- dist(t(truseq.v.4N.matr.plasm))
hc.cols <- hclust(dcols.truseq.v.4N.matr.plasm, method="complete")
# To reorder the plot as much as possible by cols
callback <- function(hc, mat){
  label.split <- strsplit2(hc$labels, split=".", fixed=TRUE)
  wt1 <- ifelse(label.split[, dim(label.split)[2]]=="transf", 1, 0)
  ref.grp <- unique(label.split[wt1>0, 1])
  wt2 <- ifelse(label.split[,1]==ref.grp, 1, 0)
  grp.wt.df <- data.frame(wt1=wt1, wt2=wt2, row.names=hc$labels)
  group.wts <- rowSums(grp.wt.df)
  dend <- reorder(as.dendrogram(hc), wts=group.wts)
  as.hclust(dend)
}
pheatmap(truseq.v.4N.matr.plasm[order(-row.mean.abund.df.plasm$TruSeq.orig),], legend_breaks = seq(0, 16, by=2)  , 
         breaks = seq(0, max(truseq.v.4N.matr.plasm), length.out = 101),
         cluster_rows = FALSE,
         show_rownames = FALSE,
         clustering_distance_cols = dcols.truseq.v.4N.matr.plasm,
         fontsize = 10,
         fontsize_col = 10,
         treeheight_col = 15,
         
         labels_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot$transf.LibMethod.id,
         annotation_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot[, 1:2], 
         annotation_row = row.mean.abund.df.plasm[, c(1, 3, 2)],
         annotation_colors=FUNCTION.filter.ann_colors(ann.colors.correction.facts, TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot),
         clustering_callback=callback, filename="FIG_6C_TruSeq_V_4N_Correction_HEATMAP_PLASMAPOOL_V2.pdf")
pheatmap(truseq.v.4N.matr.plasm[order(-row.mean.abund.df.plasm$TruSeq.orig),], legend_breaks = seq(0, 16, by=2)  , 
         breaks = seq(0, max(truseq.v.4N.matr.plasm), length.out = 101),
         cluster_rows = FALSE,
         show_rownames = FALSE,
         clustering_distance_cols = dcols.truseq.v.4N.matr.plasm,
         fontsize = 10,
         fontsize_col = 10,
         treeheight_col = 15,
         
         labels_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot$transf.LibMethod.id,
         annotation_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot[, 1:2], 
         annotation_row = row.mean.abund.df.plasm[, c(1, 3, 2)],
         annotation_colors=FUNCTION.filter.ann_colors(ann.colors.correction.facts, TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot),
         clustering_callback=callback)


# Mann-whitney tests for equimolar pool
MANN.WHITNEY.SHIFTS.EQ <- sapply(all.intermethod.comparison.groups, USE.NAMES=TRUE, simplify = FALSE, FUN=function(x.group){
  correct.dt <- data.frame(subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, correction==x.group), row.names=1)
  x.ids <- strsplit2(x.group, split=".V.")
  print(x.group)
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[,ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.detail%in%x.ids]
  keep.mirs <- (rowSums(dge.test$counts>=5)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  dge.test.filt <- calcNormFactors(dge.test.filt, method="RLE")
  v.test.norm <- voom(dge.test.filt, design.test, normalize.method="scale")
  cpm.test.norm <- v.test.norm$E
  
  merged.miR.ids <- merge(cpm.test.norm, correct.dt, by=0)$Row.names
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  query.ids <- grep(x.ids[1], colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)
  ref.ids <- grep(x.ids[1], invert = TRUE, colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)
  cpm.test.norm.filt.md1 <- cpm.test.norm.filt[, query.ids]
  cpm.test.norm.filt.ref <- cpm.test.norm.filt[, ref.ids]
  tt.comp.scl.factor.filt <- correct.dt[merged.miR.ids,]
  cpm.test.norm.filt.md.1.scaled <-
    apply(
      cpm.test.norm.filt.md1,
      2,
      FUN = function(x) {
        x.log.scl <- x - tt.comp.scl.factor.filt$logFC
      })
  # combine for re-centering
  cpm.scaled.and.ref <- cbind(cpm.test.norm.filt.ref, cpm.test.norm.filt.md.1.scaled)[, colnames(dge.test)]
  cpm.test.scaled <- voom(2^cpm.scaled.and.ref, design.test, normalize.method="scale")$E
  cpm.test.scaled.group.means <- data.frame(mean.scaled.query=rowMeans(cpm.test.scaled[, query.ids]), mean.ref=rowMeans(cpm.test.scaled[, ref.ids]), row.names=row.names(cpm.test.scaled))
  cpm.test.scaled.group.means$abs.log.dif <- with(cpm.test.scaled.group.means, abs(mean.scaled.query-mean.ref))
  cpm.test.scaled.group.means$comparison <- "CORRECTED"
  cpm.test.original.group.means <- data.frame(mean.unscaled.query=rowMeans(cpm.test.norm.filt[, query.ids]), mean.ref=rowMeans(cpm.test.norm.filt[, ref.ids]), row.names=row.names(cpm.test.norm.filt))
  cpm.test.original.group.means$abs.log.dif <- with(cpm.test.original.group.means, abs(mean.unscaled.query-mean.ref))
  cpm.test.original.group.means$comparison <- "UNCORRECTED"
  cpm.both <- data.table(keep.rownames = TRUE, rbind(cpm.test.original.group.means[, c("comparison", "abs.log.dif")], cpm.test.scaled.group.means[, c("comparison", "abs.log.dif")]))
  cpm.both[, `:=`(group.id=paste0(x.ids[1], ".V.", x.ids[2]), from.method=x.ids[1], to.method=x.ids[2])]
  w.test.out <- wilcox.test(x = cpm.test.original.group.means$abs.log.dif, y=cpm.test.scaled.group.means$abs.log.dif, alternative = c("greater"),  conf.int = TRUE)
  cpm.both[, `:=`(W=w.test.out$statistic, p.value=w.test.out$p.value, alternative="greater", estimate=w.test.out$estimate)]
  
  #get matrix to include as well
  colnames(cpm.test.norm.filt.md.1.scaled) <- paste0(colnames(cpm.test.norm.filt.md.1.scaled), ".transf")
  colnames(cpm.test.norm.filt) <- paste0(colnames(cpm.test.norm.filt), ".orig")
  cpm.full.matr.with.correct <- cbind(cpm.test.norm.filt, cpm.test.norm.filt.md.1.scaled)
  cpm.full.matr.with.correct.norm <- normalizeBetweenArrays(cpm.full.matr.with.correct, method = "scale")
  return(list(cpm.both=cpm.both, cpm.scaled.matr=cpm.full.matr.with.correct.norm))
}); MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF <- do.call(rbind, lapply(MANN.WHITNEY.SHIFTS.EQ, FUN=function(x) x[[1]])); MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF[, comparison:=factor(comparison, levels=c("UNCORRECTED", "CORRECTED"))]

MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF.STATS <- MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF[comparison=="UNCORRECTED", .(W=max(W), p.value=min(p.value), estimate=min(estimate), median.abs.log.dif=quantile(abs.log.dif, 0.9)), by=.(group.id, from.method, to.method, alternative)]
setorder(MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF.STATS, p.value)
MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF[, group.id:=factor(group.id, levels=MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF.STATS$group.id)]
MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF.STATS[, group.id:=factor(group.id, levels=MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF.STATS$group.id)]

g.s7.all.eq <- ggplot(MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF, aes(x=abs.log.dif, color=comparison)) + geom_density(adjust=1, trim=TRUE) + 
  scale_x_continuous(expand=c(0.01,0), breaks=seq(0,20, 1.0)) + 
  theme(
    axis.text=element_text(color="black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color="black"),
    axis.ticks = element_line(color="black"),
    legend.position = "none",
    legend.direction = "horizontal",
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(y="density",
       x="|Difference VS Reference|(log2)")
g.all <- plot_grid(plotlist = list(g1, g2, g3), nrow = 1, labels = paste0("Correction: ", sub(".V.", " to ", show.groups)), label_size = 8)
this.outfile <- paste0(outdirs["supplemental_figures"], "/FigS7_DensityPlot_EquimolarPool_ALL.pdf"); this.outfile
save_plot(filename = this.outfile, ncol = 3, nrow = 1, base_width = 2.5, base_height = 2, plot = g.all)

# Supp Fig 7 Plasma boxpot Plot Correction Factors with mann W -------------------------------------------------------------------
show.groups <- c("TruSeq.V.4N_B", "TruSeq.V.NEBNext", "NEBNext.V.4N_B")
g.s7.all.eq.mannW <- ggplot(MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF[group.id %in% show.groups], aes(y=abs.log.dif, x=comparison, fill=comparison)) +
  geom_boxplot() + 
  facet_wrap(~group.id) + 
  theme_bw() + 
  geom_text(
    data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS[group.id %in% show.groups], 
    aes(
      label=paste0("p = ", format(p.value, digits=3)),
      x=2, y=8),
    inherit.aes = FALSE, check_overlap=TRUE) +
  scale_y_continuous(breaks=seq(0.0, 16, 0.5), limits=c(0.0, 13)) + 
  #scale_x_continuous(expand=c(0.01,0), limits=c(0, 8.0), breaks=seq(0,9, 1.0)) + 
  theme(
    axis.text=element_text(size=10, color="black"),
    axis.title=element_text(size=10),
    plot.title=element_text(size=10),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=10, color="black"),
    #legend.position = "top",
    #legend.direction = "horizontal",
    legend.text=element_text(size=10),
    legend.title=element_text(size=10),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(x=NULL,
       y="|Difference VS Reference|(log2)"); g.s7.all.eq.mannW
this.outfile <- paste0(outdirs["supplemental_figures"], "/FigS7_boxplot_aggregate_EquimolarPool.pdf"); this.outfile
ggsave(this.outfile, g.s7.all.mannW, units = "in", width=6.5, height=3)


# 1-10-17 Now do pairwise by sample--same way as standart method, but dont average over labs
MANN.WHITNEY.SHIFTS.PLASMA.FISHER.BY.LAB <- sapply(all.intermethod.comparison.groups, USE.NAMES=TRUE, simplify = FALSE, FUN=function(x.group){
  correct.dt <- data.frame(subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, correction==x.group), row.names=1)
  x.ids <- strsplit2(x.group, split=".V.")
  print(x.group)
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.detail%in%x.ids]
  keep.mirs <- (rowSums(dge.test$counts>=5)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  query.ids <- grep(x.ids[1], colnames(dge.test.filt), value=TRUE, ignore.case=TRUE)
  ref.ids <- grep(x.ids[1], invert = TRUE, colnames(dge.test.filt), value=TRUE, ignore.case=TRUE)
  
  dge.query <- dge.test.filt[, query.ids]
  #design.query <- model.matrix(~0+lab.libMethod, dge.query$samples)
  dge.query <- calcNormFactors(dge.query, method="RLE")
  #v.query <- voom(dge.query, normalize.method="quantile")
  
  dge.ref <- dge.test.filt[, ref.ids]
  #design.ref <- model.matrix(~0+lab.libMethod, dge.ref$samples)
  dge.ref <- calcNormFactors(dge.ref, method="RLE")
  #v.ref <- voom(dge.ref, normalize.method="quantile")
  
  merged.miR.ids <- row.names(dge.test.filt)[row.names(dge.test.filt)%in%row.names(correct.dt)]
  #cpm.test.norm.filt.md1 <- v.query$E[merged.miR.ids,]
  cpm.test.norm.filt.md1 <- cpm.DGEList(dge.query, log = TRUE, normalized.lib.sizes = TRUE)[merged.miR.ids,]
  #cpm.test.norm.filt.ref <- v.ref$E[merged.miR.ids,]
  cpm.test.norm.filt.ref <- cpm.DGEList(dge.ref, log = TRUE, normalized.lib.sizes = TRUE)[merged.miR.ids,]
  tt.comp.scl.factor.filt <- correct.dt[merged.miR.ids,]
  
  cpm.test.norm.filt.md.1.scaled <-
    apply(
      cpm.test.norm.filt.md1,
      2,
      FUN = function(x) {
        x.log.scl <- x - tt.comp.scl.factor.filt$logFC
      })
  # combine for re-centering
  cpm.scaled.and.ref <- cbind(cpm.test.norm.filt.ref, cpm.test.norm.filt.md.1.scaled)[, colnames(dge.test.filt)]
  cpm.test.scaled <- voom(2^cpm.scaled.and.ref, design.test, normalize.method="quantile")$E
  mean.cpm.test.scaled.ref <- rowMeans(cpm.test.scaled[, ref.ids])
  query.cpm.test.scaled <- cpm.test.scaled[, query.ids]
  abs.dif.scaled.vs.ref <- abs(query.cpm.test.scaled-mean.cpm.test.scaled.ref)
  
  # re-do the normalization for the original data to "quantile" so that it's normalized the same way as the scaled data.
  cpm.unscaled.and.ref <- cbind(cpm.test.norm.filt.ref, cpm.test.norm.filt.md1)[, colnames(dge.test.filt)]
  v.test.norm <- voom(2^cpm.unscaled.and.ref, design.test, normalize.method="quantile")
  cpm.test.norm <- v.test.norm$E
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  query.cpm.test.norm.filt <- cpm.test.norm.filt[, query.ids]
  mean.cpm.test.norm.filt.ref <- rowMeans(cpm.test.norm.filt[, ref.ids])
  abs.dif.unscaled.vs.ref <- abs(query.cpm.test.norm.filt-mean.cpm.test.norm.filt.ref)
  
  # Now get lab-specific differences
  lab.ids <- unique(dge.query$samples$Lab)
  
  lab.wilcox.tests <- sapply(lab.ids, simplify=FALSE, USE.NAMES = TRUE, FUN=function(x){
    print(paste0("\tWILCOX.TEST: ", x))
    lab.query.ids <- query.ids[(strsplit2(query.ids, split=".", fixed=TRUE)[,2])==x]
    
    lab.cpm.test.scaled.group.means <- data.frame(mean.query.cpm.lab = rowMeans(query.cpm.test.scaled[, lab.query.ids]),
                                                  mean.ref.cpm.all.labs = mean.cpm.test.scaled.ref,  
                                                  abs.log.dif = rowMeans(abs.dif.scaled.vs.ref[, lab.query.ids]),
                                                  comparison = "CORRECTED",
                                                  query.Lab = x,
                                                  row.names=row.names(query.cpm.test.scaled))
    
    lab.cpm.test.original.group.means <- data.frame(mean.query.cpm.lab = rowMeans(query.cpm.test.norm.filt[, lab.query.ids]),
                                                    mean.ref.cpm.all.labs = mean.cpm.test.norm.filt.ref,  
                                                    abs.log.dif = rowMeans(abs.dif.unscaled.vs.ref[, lab.query.ids]),
                                                    comparison = "UNCORRECTED",
                                                    query.Lab = x,
                                                    row.names=row.names(query.cpm.test.scaled))
    lab.cpm.both <- data.table(keep.rownames = TRUE, rbind(lab.cpm.test.scaled.group.means, lab.cpm.test.original.group.means))
    lab.cpm.both[, `:=`(group.id=paste0(x.ids[1], ".V.", x.ids[2]), from.method=x.ids[1], to.method=x.ids[2])]
    
    lab.w.test.out <- wilcox.test(x = lab.cpm.test.original.group.means$abs.log.dif, y=lab.cpm.test.scaled.group.means$abs.log.dif, alternative = c("greater"),  conf.int = TRUE)
    lab.cpm.both[, `:=`(W=lab.w.test.out$statistic, p.value=lab.w.test.out$p.value, alternative="greater", estimate=lab.w.test.out$estimate)]
    return(list(full.table.stats=lab.cpm.both, wilcox.test.only=as.list(lab.w.test.out)))
  })
  lab.wilcox.tests.dt <- do.call(rbind, lapply(lab.wilcox.tests, FUN=function(x){x[[1]]}))
  wilcox.tests <- lapply(lab.wilcox.tests, FUN=function(x){x[[2]]})
  
  wilcox.pvals <- unlist(lapply(wilcox.tests, FUN=function(x){x$p.value}))
  df <- 2*length(wilcox.pvals)
  p.combined <- pchisq(-2*sum(log(wilcox.pvals)), df, lower.tail=FALSE)
  
  lab.wilcox.tests.dt[, `:=`(fisher.p.val=p.combined, fisher.df=df)]
  return(lab.wilcox.tests.dt)
});

MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB <- do.call(rbind, MANN.WHITNEY.SHIFTS.PLASMA.FISHER.BY.LAB)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB[, comparison:=factor(comparison, levels=c("UNCORRECTED", "CORRECTED"))]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.STATS <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB[, .(abs.log.dif=max(abs.log.dif)), by=c("group.id", "from.method", "to.method", "query.Lab", "W", "p.value", "alternative", "estimate", "fisher.p.val", "fisher.df", "comparison")]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.STATS[, max.abs.log.dif.group:=max(abs.log.dif), by=group.id]

sapply(show.groups, FUN=function(x){
  gid <- x
  lab1 <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.STATS[group.id == gid & comparison=="CORRECTED", unique(query.Lab)][1]
g <- ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB[group.id == gid],
       aes(y=abs.log.dif,
           x=query.Lab,
           fill=comparison)
) +
  geom_boxplot(position = "dodge") + 
  facet_wrap(~group.id) + 
  geom_text(size=1,
    data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.STATS[group.id == gid & comparison=="CORRECTED"], 
    aes(
      label=paste0("p = ", format(p.value, digits=3)),
      y=max.abs.log.dif.group*1.08
    ),
    inherit.aes = TRUE, check_overlap=TRUE) +
  geom_text(size=1,
    data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.STATS[group.id == gid & comparison=="CORRECTED"], 
    aes(
      label=paste0("w = ", W),
      y=max.abs.log.dif.group*1.05
    ),
    inherit.aes = TRUE, check_overlap=TRUE) +
  geom_text(size=1,
    data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.STATS[group.id == gid & comparison=="CORRECTED" & query.Lab==lab1], 
    aes(
      label=paste0("FDR = ", format(fisher.p.val, digits=3), "; DF = ", fisher.df, "; estimate = ", format(estimate, digits=3)),
      y=max.abs.log.dif.group*1.12
    ),
    inherit.aes = TRUE, check_overlap=TRUE) +
  scale_y_continuous(expand=c(0.01,0), limits = c(0, 12), breaks=seq(0, 12, 1.0)) + 
  theme(text = element_text(size=8),
    axis.text=element_text(color="black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size=1, color="black"),
    axis.ticks = element_line(color="black", size=1.0),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(x=NULL,
       y="|Difference VS Reference|(log2)")

this.outfile <- paste0(outdirs["supplemental_figures"], "/FigS7_BOXPLOT_PLASMAPool_", gid, ".pdf")
ggsave(filename = this.outfile, plot = g, units = "in", width=4, height=3)
})


MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.plot <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB[group.id %in% show.groups]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.plot <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF[group.id %in% show.groups]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.plot[, `:=`(query.Lab="Aggregate", fisher.p.val=NA, fisher.df=NA)]

MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.Add.Agg <- rbind(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.plot, MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.plot, fill=TRUE)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.Add.Agg.STATS <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.Add.Agg[, .(abs.log.dif=max(abs.log.dif)), by=c("group.id", "from.method", "to.method", "query.Lab", "W", "p.value", "alternative", "estimate", "fisher.p.val", "fisher.df", "comparison")]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.Add.Agg.STATS[, max.abs.log.dif.group:=max(abs.log.dif), by=group.id]

sapply(show.groups, FUN=function(x){
  gid <- x
  lab1 <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.Add.Agg.STATS[group.id == gid & comparison=="CORRECTED", unique(query.Lab)][1]
  g <- ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.Add.Agg[group.id == gid],
              aes(y=abs.log.dif,
                  x=query.Lab,
                  fill=comparison)
  ) +
    geom_boxplot(position = "dodge") + 
    facet_wrap(~group.id) + 
    geom_text(size=1,
              data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.Add.Agg.STATS[group.id == gid & comparison=="CORRECTED"], 
              aes(
                label=paste0("p = ", format(p.value, digits=3)),
                y=max.abs.log.dif.group*1.08
              ),
              inherit.aes = TRUE, check_overlap=TRUE) +
    geom_text(size=1,
              data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.Add.Agg.STATS[group.id == gid & comparison=="CORRECTED"], 
              aes(
                label=paste0("w = ", W),
                y=max.abs.log.dif.group*1.05
              ),
              inherit.aes = TRUE, check_overlap=TRUE) +
    geom_text(size=1,
              data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.Add.Agg.STATS[group.id == gid & comparison=="CORRECTED" & query.Lab==lab1], 
              aes(
                label=paste0("FDR = ", format(fisher.p.val, digits=3), "; DF = ", fisher.df, "; estimate = ", format(estimate, digits=3)),
                y=max.abs.log.dif.group*1.12
              ),
              inherit.aes = TRUE, check_overlap=TRUE) +
    scale_y_continuous(expand=c(0.01,0), limits = c(0, 12), breaks=seq(0, 12, 1.0)) + 
    theme(          axis.text=element_text(color="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color="black"),
          axis.ticks = element_line(color="black"),
          legend.position = "top",
          legend.direction = "horizontal",
          legend.key.size=unit(1.2, "lines")
    ) +
    labs(x=NULL,
         y="|Difference VS Reference|(log2)")
  
  this.outfile <- paste0(outdirs["supplemental_figures"], "/FigS8_BOXPLOT_PLASMAPool_withAgg", gid, ".pdf")
  ggsave(filename = this.outfile, plot = g, units = "in", width=4, height=3)
})




# 1/19/2017--Check the rank-ordering too. D. Erle suggested looking at top N miRs before and after downbiasing
downbias.comparison.groups <- grep("4N", grep("DOWNBIAS", unique(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset$correction), value=TRUE), invert=TRUE, value=TRUE)
DOWNBIAS.BY.LAB <- sapply(downbias.comparison.groups, USE.NAMES=TRUE, simplify = FALSE, FUN=function(x.group){
  
  correct.dt <- data.frame(subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, correction==x.group), row.names=1)
  x.ids <- strsplit2(x.group, split=".V.")
  print(x.group)
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple%in%x.ids | ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.detail%in%x.ids | ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$GROUND.TRUTH.SAMPLE%in%x.ids ]
  keep.mirs <- (rowSums(dge.test$counts>=5)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  query.ids <- grep(x.ids[1], colnames(dge.test.filt), value=TRUE, ignore.case=TRUE)
  ref.ids <- grep(x.ids[1], invert = TRUE, colnames(dge.test.filt), value=TRUE, ignore.case=TRUE)
  
  dge.query <- dge.test.filt[, query.ids]
  #design.query <- model.matrix(~0+lab.libMethod, dge.query$samples)
  dge.query <- calcNormFactors(dge.query, method="RLE")
  #v.query <- voom(dge.query, normalize.method="quantile")
  
  dge.ref <- dge.test.filt[, ref.ids]
  #design.ref <- model.matrix(~0+lab.libMethod, dge.ref$samples)
  dge.ref <- calcNormFactors(dge.ref, method="RLE")
  #v.ref <- voom(dge.ref, normalize.method="quantile")
  
  merged.miR.ids <- row.names(dge.test.filt)[row.names(dge.test.filt)%in%row.names(correct.dt)]
  #cpm.test.norm.filt.md1 <- v.query$E[merged.miR.ids,]
  cpm.test.norm.filt.md1 <- cpm.DGEList(dge.query, log = TRUE, normalized.lib.sizes = TRUE)[merged.miR.ids,]
  #cpm.test.norm.filt.ref <- v.ref$E[merged.miR.ids,]
  cpm.test.norm.filt.ref <- cpm.DGEList(dge.ref, log = TRUE, normalized.lib.sizes = TRUE)[merged.miR.ids,]
  tt.comp.scl.factor.filt <- correct.dt[merged.miR.ids,]
  
  cpm.test.norm.filt.md.1.scaled <-
    apply(
      cpm.test.norm.filt.md1,
      2,
      FUN = function(x) {
        x.log.scl <- x - tt.comp.scl.factor.filt$logFC
      })
  # combine for re-centering
  cpm.scaled.and.ref <- cbind(cpm.test.norm.filt.ref, cpm.test.norm.filt.md.1.scaled)[, colnames(dge.test.filt)]
  cpm.test.scaled <- voom(2^cpm.scaled.and.ref, design.test, normalize.method="quantile")$E
  mean.cpm.test.scaled.ref <- rowMeans(cpm.test.scaled[, ref.ids])
  query.cpm.test.scaled <- cpm.test.scaled[, query.ids]
  abs.dif.scaled.vs.ref <- abs(query.cpm.test.scaled-mean.cpm.test.scaled.ref)
  
  # re-do the normalization for the original data to "quantile" so that it's normalized the same way as the scaled data.
  cpm.unscaled.and.ref <- cbind(cpm.test.norm.filt.ref, cpm.test.norm.filt.md1)[, colnames(dge.test.filt)]
  v.test.norm <- voom(2^cpm.unscaled.and.ref, design.test, normalize.method="quantile")
  cpm.test.norm <- v.test.norm$E
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  query.cpm.test.norm.filt <- cpm.test.norm.filt[, query.ids]
  mean.cpm.test.norm.filt.ref <- rowMeans(cpm.test.norm.filt[, ref.ids])
  abs.dif.unscaled.vs.ref <- abs(query.cpm.test.norm.filt-mean.cpm.test.norm.filt.ref)
  
  # Now get lab-specific differences
  lab.ids <- unique(dge.query$samples$Lab)
  
  lab.wilcox.tests <- sapply(lab.ids, simplify=FALSE, USE.NAMES = TRUE, FUN=function(x){
    print(paste0("\tWILCOX.TEST: ", x))
    lab.query.ids <- query.ids[(strsplit2(query.ids, split=".", fixed=TRUE)[,2])==x]
    
    lab.cpm.test.scaled.group.means <- data.frame(mean.query.cpm.lab = rowMeans(query.cpm.test.scaled[, lab.query.ids]),
                                                  mean.ref.cpm.all.labs = mean.cpm.test.scaled.ref,  
                                                  abs.log.dif = rowMeans(abs.dif.scaled.vs.ref[, lab.query.ids]),
                                                  comparison = "CORRECTED",
                                                  query.Lab = x,
                                                  row.names=row.names(query.cpm.test.scaled))
    
    lab.cpm.test.original.group.means <- data.frame(mean.query.cpm.lab = rowMeans(query.cpm.test.norm.filt[, lab.query.ids]),
                                                    mean.ref.cpm.all.labs = mean.cpm.test.norm.filt.ref,  
                                                    abs.log.dif = rowMeans(abs.dif.unscaled.vs.ref[, lab.query.ids]),
                                                    comparison = "UNCORRECTED",
                                                    query.Lab = x,
                                                    row.names=row.names(query.cpm.test.scaled))
    lab.cpm.both <- data.table(keep.rownames = TRUE, rbind(lab.cpm.test.scaled.group.means, lab.cpm.test.original.group.means))
    lab.cpm.both[, `:=`(group.id=paste0(x.ids[1], ".V.", x.ids[2]), from.method=x.ids[1], to.method=x.ids[2])]
    
    lab.w.test.out <- wilcox.test(x = lab.cpm.test.original.group.means$abs.log.dif, y=lab.cpm.test.scaled.group.means$abs.log.dif, alternative = c("greater"),  conf.int = TRUE)
    lab.cpm.both[, `:=`(W=lab.w.test.out$statistic, p.value=lab.w.test.out$p.value, alternative="greater", estimate=lab.w.test.out$estimate)]
    return(list(full.table.stats=lab.cpm.both, wilcox.test.only=as.list(lab.w.test.out)))
  })
  lab.wilcox.tests.dt <- do.call(rbind, lapply(lab.wilcox.tests, FUN=function(x){x[[1]]}))
  wilcox.tests <- lapply(lab.wilcox.tests, FUN=function(x){x[[2]]})
  
  wilcox.pvals <- unlist(lapply(wilcox.tests, FUN=function(x){x$p.value}))
  df <- 2*length(wilcox.pvals)
  p.combined <- pchisq(-2*sum(log(wilcox.pvals)), df, lower.tail=FALSE)
  
  lab.wilcox.tests.dt[, `:=`(fisher.p.val=p.combined, fisher.df=df)]
  return(lab.wilcox.tests.dt)
});

MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS <- do.call(rbind, DOWNBIAS.BY.LAB)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.UNCORRECTED <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[comparison=="UNCORRECTED", .(mean.query.cpm.lab.UNCORRECTED=mean(mean.query.cpm.lab), mean.ref.cpm.all.labs.UNCORRECTED=mean(mean.ref.cpm.all.labs)), by=.(rn, query.Lab, group.id, from.method)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.UNCORRECTED[, rn:=sub("1$", "", rn)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.CORRECTED <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[comparison=="CORRECTED", .(mean.query.cpm.lab.CORRECTED=mean(mean.query.cpm.lab), mean.ref.cpm.all.labs.CORRECTED=mean(mean.ref.cpm.all.labs)), by=.(rn, query.Lab, group.id, from.method)]
setkeyv(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.UNCORRECTED, c("rn", "query.Lab", "group.id", "from.method"))
setkeyv(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.CORRECTED, c("rn", "query.Lab", "group.id", "from.method"))
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY <- merge(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.UNCORRECTED, MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.CORRECTED)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt <- melt.data.table(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY, id.vars =  c("rn", "query.Lab", "group.id", "from.method"), variable.name = "CPM.MEASURE", value.name = "CPM")
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt[, rank:=rank(-CPM, ties.method = "min", na.last = TRUE), by=.(query.Lab, group.id, CPM.MEASURE)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt[, min.rank:=min(rank), by=.(rn, from.method)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt[, miR.id.short:=tstrsplit(rn, split=":", fixed=TRUE)[1]]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt[, CPM.MEASURE.2:=ifelse(CPM.MEASURE=="mean.query.cpm.lab.UNCORRECTED", paste0(from.method, ".UNCORRECTED"), ifelse(CPM.MEASURE=="mean.query.cpm.lab.CORRECTED", paste0(from.method, ".CORRECTED"), "REFERENCE")), by=1:nrow(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.mean <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt[, .(rank=min(rank)), by=.(rn, query.Lab, group.id, from.method, miR.id.short, CPM.MEASURE.2)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.mean[, min.rank:=min(rank), by=.(rn, from.method)]


MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.truseq <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.mean[min.rank<=10 & from.method=="TruSeq"]
setorderv(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.truseq, c("CPM.MEASURE.2", "rank"))
mir.id.order <- unique(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.truseq$miR.id.short)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.truseq[, miR.id.short:=factor(miR.id.short, mir.id.order)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.truseq[, CPM.MEASURE.2:=factor(CPM.MEASURE.2, c("REFERENCE", "TruSeq.UNCORRECTED", "TruSeq.CORRECTED"))]
ids.ref.lt.10 <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.truseq[rank<=10 & CPM.MEASURE.2=="REFERENCE"]$miR.id.short

ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.truseq[miR.id.short%in%ids.ref.lt.10], aes(x=miR.id.short, pos=CPM.MEASURE.2, y=rank, fill=CPM.MEASURE.2)) + geom_boxplot() + theme_bw() + scale_y_continuous(breaks = seq(0,200,20)) + theme(
  axis.text=element_text(size=10, color="black"),
  axis.text.x = element_text(hjust=1, angle = 50),
  axis.title=element_text(size=10),
  plot.title=element_text(size=10),
  panel.grid.major = element_line(color=NA),
  panel.border = element_rect(size=1.5, color="black"),
  axis.ticks = element_line(color="black", size=1.5),
  strip.text = element_text(size=10, color="black"),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.text=element_text(size=10),
  legend.title=element_text(size=10),
  legend.key.size=unit(1.2, "lines")
) 
ggsave("UNDERBIASED_MIRS_TRUSEQ_DOWNBIAS.pdf")

ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.truseq[!miR.id.short%in%ids.ref.lt.10], aes(x=miR.id.short, pos=CPM.MEASURE.2, y=rank, fill=CPM.MEASURE.2)) + geom_boxplot() + theme_bw() + scale_y_continuous(breaks = seq(0,200,20)) + theme(
  axis.text=element_text(size=10, color="black"),
  axis.text.x = element_text(hjust=1, angle = 50),
  axis.title=element_text(size=10),
  plot.title=element_text(size=10),
  panel.grid.major = element_line(color=NA),
  panel.border = element_rect(size=1.5, color="black"),
  axis.ticks = element_line(color="black", size=1.5),
  strip.text = element_text(size=10, color="black"),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.text=element_text(size=10),
  legend.title=element_text(size=10),
  legend.key.size=unit(1.2, "lines")
) 
ggsave("OVERBIASED_MIRS_TRUSEQ_DOWNBIAS.pdf")


# add mann-whitney results to the original scaling factors table
setkey(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, correction)
write.table(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, "SUPPLEMENTARY_TABLE_CORRECTION_FACTORS.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

#
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS <- do.call(rbind, DOWNBIAS.BY.LAB)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[, rn:=sub("1$", "", rn)]
setorder(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS, rn)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[, comparison:=factor(comparison, levels=c("UNCORRECTED", "CORRECTED"))]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[, `:=`(query.rank.in.sample=rank(-mean.query.cpm.lab, ties.method = "ave"),
                                                                            ref.rank.in.sample=rank(-mean.ref.cpm.all.labs, ties.method = "ave")),
                                                                     by=.(query.Lab, comparison, group.id, from.method, to.method)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[, ref.rank.in.sample.uncorrected:=min(ref.rank.in.sample), by=.(rn, group.id)]

MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.ref.top10 <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[ref.rank.in.sample.uncorrected<=10]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.ref.top10[, miR.id.short:=tstrsplit(rn, split=":", fixed=TRUE)[1]]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.ref.top10[, miR.id.short.rank:=paste0(miR.id.short, "(", ref.rank.in.sample.uncorrected, ")")]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.ref.top10.truseq <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.ref.top10[from.method=="TruSeq"]
setorder(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.ref.top10.truseq, ref.rank.in.sample.uncorrected)
mir.id.sort.truseq <- unique(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.ref.top10.truseq$miR.id.short)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.ref.top10.truseq[, miR.id.short:=factor(miR.id.short, levels=mir.id.sort.truseq)]

ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.ref.top10.truseq, aes(x=miR.id.short, pos=comparison, y=query.rank.in.sample, color=comparison)) + geom_jitter(height=0, width=0.1) + theme_bw() + theme(
  axis.text=element_text(size=10, color="black"),
  axis.text.x = element_text(hjust=1, angle = 50),
  axis.title=element_text(size=10),
  plot.title=element_text(size=10),
  panel.grid.major = element_line(color=NA),
  panel.border = element_rect(size=1.5, color="black"),
  axis.ticks = element_line(color="black", size=1.5),
  strip.text = element_text(size=10, color="black"),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.text=element_text(size=10),
  legend.title=element_text(size=10),
  legend.key.size=unit(1.2, "lines")
) 

MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[, mean.query.cpm:=mean(mean.query.cpm.lab), by=.(group.id, comparison, rn)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[, query.rank.in.sample:=rank(-mean.query.cpm, ties.method = "ave"), by=.(group.id, comparison, query.Lab)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[from.method=="TruSeq"]

MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq.top10.ids <- unique(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq[query.rank.in.sample<=10 & comparison=="UNCORRECTED"]$rn)

MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq.query.top10 <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq[rn%in%MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq.top10.ids]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq.query.top10[, miR.id.short:=tstrsplit(rn, split=":", fixed=TRUE)[1]]
setorder(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq.query.top10, query.rank.in.sample)
mir.id.sort.truseq <- unique(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq.query.top10$miR.id.short)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq.query.top10[, miR.id.short:=factor(miR.id.short, levels=mir.id.sort.truseq)]

ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq.query.top10, aes(x=miR.id.short, pos=comparison, y=query.rank.in.sample, color=comparison)) + geom_jitter(height=0, width=0.1) + theme_bw() + theme(
  axis.text=element_text(size=10, color="black"),
  axis.text.x = element_text(hjust=1, angle = 50),
  axis.title=element_text(size=10),
  plot.title=element_text(size=10),
  panel.grid.major = element_line(color=NA),
  panel.border = element_rect(size=1.5, color="black"),
  axis.ticks = element_line(color="black", size=1.5),
  strip.text = element_text(size=10, color="black"),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.text=element_text(size=10),
  legend.title=element_text(size=10),
  legend.key.size=unit(1.2, "lines")
) 







######## NOW ALSO SHOW CORRECTIONS FOR EQUIMOLAR POOL WITH EQUIMOLAR SEQS
# Now this is where we need to diverge to include the other methods. 
# The thought is that scaling factors are going to convert between two different methods, so all of our filtering, normalization, centering and testing should 
# be done without the potential influence of the other library prep types.

lib.comparison.groups <- list(TruSeq.V.4N=c("TruSeq", "4N"), 
                              TruSeq.V.NEB=c("TruSeq", "NEBNext"),
                              TruSeq.V.CleanTag=c("TruSeq", "CleanTag"),
                              NEB.V.CleanTag=c("NEBNext", "CleanTag"),
                              NEB.V.4N=c("NEBNext", "4N"),
                              CleanTag.V.4N=c("CleanTag", "4N")
)


ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.EQUIMOLAR <- do.call(rbind, sapply(names(lib.comparison.groups), USE.NAMES=TRUE, simplify = FALSE, FUN=function(x.id.list){
  x.ids <- unlist(lib.comparison.groups[x.id.list], recursive = TRUE, use.names = FALSE)
  dge.comp <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[, ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.simple%in%x.ids]
  design.comp <- model.matrix(~0+lab.libMethod, dge.comp$samples)
  colnames(design.comp) <- unique(dge.comp$samples$lab.libMethod)
  lib.methods.comp <- strsplit2(strsplit2(colnames(design.comp), split=".", fixed=TRUE)[,1], split="_")[,1]
  comp.contrasts.fodder <- sapply(unique(lib.methods.comp), USE.NAMES=FALSE, FUN=function(lib.m){
    mat <- design.comp
    y.sum <- rowSums(apply(mat, 1, FUN=function(x){
      x*ifelse(lib.methods.comp==lib.m, 1, 0)
    }))
    y.sum/sum(y.sum)
  })
  colnames(comp.contrasts.fodder) <- unique(lib.methods.comp)
  out.contrast <- as.numeric(comp.contrasts.fodder[, x.ids[1]]-comp.contrasts.fodder[, x.ids[2]])
  # Now filter miRs
  dge.comp.filt <- dge.comp[rowSums(dge.comp$counts>0)==dim(dge.comp)[2],]
  v.comp <- voom(dge.comp.filt, design.comp, normalize.method = "scale")
  fit.comp <- lmFit(v.comp, design.comp)
  fit.comp <- eBayes(fit.comp, robust = TRUE, trend = TRUE)
  fit2.comp <- contrasts.fit(fit.comp, out.contrast)
  fit2.comp <- eBayes(fit2.comp, trend = TRUE, robust = TRUE)
  tt.comp <- data.frame(topTable(fit2.comp, n = Inf, confint=0.95, sort.by = "none"))
  tt.comp.dt <- data.table(data.frame(miR.id.full=row.names(tt.comp), tt.comp))
  tt.comp.dt[, comparison:=x.id.list]
  tt.comp.scl.factor <- tt.comp[, c("logFC", "CI.L", "CI.R")]
  
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[, ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.simple%in%x.ids]
  keep.mirs <- (rowSums(dge.test$counts>=1)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  v.test.norm <- voom(dge.test.filt, design.test, normalize.method="scale")
  cpm.test.norm <- v.test.norm$E
  
  merged.miR.ids <- merge(cpm.test.norm, tt.comp.scl.factor, by=0)$Row.names
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  cpm.test.norm.filt.md1 <- cpm.test.norm.filt[, grep(x.ids[1], colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)]
  tt.comp.scl.factor.filt <- tt.comp.scl.factor[merged.miR.ids,]
  
  # set a floor min. Logic here is that if we had a miR that had a count of 1 and the scaling factor 
  # indicates a drop of -5.0 fold, we do not want to make the new log2 count log2(1) - 5.0, since our scale will be way off
  # for plotting. So, if the value is below the min, set it to the min value. In other words. Limit the amount 
  # of downward scaling to the minimum value of the sample.
  #cpm.test.norm.filt.scaled <- apply(cpm.test.norm.filt, 2, FUN=function(x){
  #  x.floor <- min(x[x>0])
  #  x.log.scl <- x-tt.comp.scl.factor.filt$logFC
  #  return(ifelse(x.log.scl<x.floor, x.floor, x.log.scl))
  #})
  
  # this does it without a floor
  # Decided to use median scaling to get the scaling factors
  # Then after scaling, use quantile normalization to match up the distributions
  # Thinking that this might be a more appropriate way to handle the fact that
  # some of the scaling at the low end will lead to values less than the minimum. 
  cpm.test.norm.filt.md.1.scaled <- apply(cpm.test.norm.filt.md1, 2, FUN=function(x){
    #x.floor <- min(x[x>0])
    x.log.scl <- x-tt.comp.scl.factor.filt$logFC
    #return(ifelse(x.log.scl<x.floor, x.floor, x.log.scl))
    return(x.log.scl)
  })
  colnames(cpm.test.norm.filt) <- paste0(colnames(cpm.test.norm.filt), ".orig")
  colnames(cpm.test.norm.filt.md.1.scaled) <- paste0(colnames(cpm.test.norm.filt.md.1.scaled), ".transf")
  cpm.test.orig.and.scaled <- cbind(cpm.test.norm.filt, cpm.test.norm.filt.md.1.scaled)
  cpm.test.orig.and.scaled.rescaled <- normalizeBetweenArrays(cpm.test.orig.and.scaled, method = "quantile")
  ids.reference <- grep(x.ids[1], colnames(cpm.test.orig.and.scaled.rescaled), value=TRUE, ignore.case=TRUE, invert = TRUE)
  ids.test.orig <- grep(x.ids[1], colnames(cpm.test.norm.filt), value=TRUE, ignore.case=TRUE)
  ids.test.scaled <- colnames(cpm.test.norm.filt.md.1.scaled)
  # Return this to get a matrix of transformed and original values
  # return(cpm.test.orig.and.scaled.rescaled)
  # Now get an average of the values we did not transform
  cpm.rescaled.test.orig <- cpm.test.orig.and.scaled.rescaled[, ids.test.orig]
  cpm.rescaled.test.scaled <- cpm.test.orig.and.scaled.rescaled[, ids.test.scaled]
  cpm.rescaled.ref.orig <- cpm.test.orig.and.scaled.rescaled[, ids.reference]
  cpm.mean.reference <- rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.reference])
  
  # Check scale norm too. Didnt look quite as good, so removed.
  cpm.test.orig.and.scaled.rescaled.med <- normalizeBetweenArrays(cpm.test.orig.and.scaled, method = "scale")
  cpm.rescaled.test.orig.med <- cpm.test.orig.and.scaled.rescaled.med[, ids.test.orig]
  cpm.rescaled.test.scaled.med <- cpm.test.orig.and.scaled.rescaled.med[, ids.test.scaled]
  cpm.rescaled.ref.orig.med <- cpm.test.orig.and.scaled.rescaled.med[, ids.reference]
  cpm.mean.reference.med <- rowMeans(cpm.test.orig.and.scaled.rescaled.med[, ids.reference])
  
  
  cpm.test.orig.and.scaled.rescaled.mean <- data.table(data.frame(miRNA.full.ID=row.names(cpm.test.orig.and.scaled.rescaled), 
                                                                  scaling.factor.logFC=tt.comp.scl.factor.filt$logFC,
                                                                  cpm.mean.reference.quantile = cpm.mean.reference,      
                                                                  cpm.mean.test.orig.quantile = rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.test.orig]),
                                                                  cpm.mean.test.scaled.quantile = rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.test.scaled]),
                                                                  abs.dif.orig.quantile = rowMeans(apply(cpm.rescaled.test.orig, 2, FUN=function(x){ abs(x-cpm.mean.reference)})),
                                                                  abs.dif.transf.quantile = rowMeans(apply(cpm.rescaled.test.scaled, 2, FUN=function(x){ abs(x-cpm.mean.reference)}))
  ))
  cpm.test.orig.and.scaled.rescaled.mean[, comparison:=x.id.list]
  return(cpm.test.orig.and.scaled.rescaled.mean)
})); #ggplot(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.LIST[[3]]) + geom_density(aes(color="VS.ORIG", x=abs.dif.orig)) + geom_density(aes(color="VS.TRANSF", x=abs.dif.transf))
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR <- melt(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.EQUIMOLAR, measure.vars=c("abs.dif.orig.quantile", "abs.dif.transf.quantile"), variable.name = "ORIG.or.TRANSF.and.SCL", value.name = "mean.abs.diff")
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR[, SCALE.METHOD:=tstrsplit(ORIG.or.TRANSF.and.SCL, split=".", fixed=TRUE)[4]]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR[, ORIG.or.TRANSF:=sub(paste0(".", SCALE.METHOD), "", ORIG.or.TRANSF.and.SCL)]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR[, `:=`(SCALE.METHOD=ifelse(SCALE.METHOD=="quantile", "QUANTILE", "MEDIAN"), ORIG.or.TRANSF=ifelse(ORIG.or.TRANSF=="abs.dif.orig", "ORIGINAL", "TRANSFORMED"))]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR.quantile <- subset(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR, SCALE.METHOD=="QUANTILE")
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR.quantile[, pool:="SynthEQ"]
# Filter the equimolar results to include only the miRs we looked at in the plasma pool
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR.quantile.filt <- subset(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR.quantile, miRNA.full.ID%in%ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile$miRNA.full.ID)
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile[, pool:="PlasmaPool"]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile.EQandPlasm <- rbind(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile, ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR.quantile.filt)
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile.EQandPlasm[, format.alpha:=ifelse(ORIG.or.TRANSF=="ORIGINAL", 1, 0.5)]

# First just plot the equimolar pool. Include all miRs; not just the ones in the plasma pool
ggplot(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR.quantile, aes(x=mean.abs.diff, color=ORIG.or.TRANSF)) + geom_density(size=1.5) + facet_wrap(~comparison) + theme_bw() + 
  scale_y_continuous(expand=c(0.01, 0.0), breaks=seq(0.0, 10.0, 0.5), limits=c(0,3.5)) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0.0, 9.0), breaks=seq(0,9, 1.0)) + 
  theme(
    axis.text=element_text(size=14, face="bold", color="black"),
    axis.title=element_text(face="bold", size=16),
    plot.title=element_text(face="bold", size=18),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=16, face="bold", color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=14, face="bold"),
    legend.title=element_text(size=15, face="bold"),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(title="SynthEQ",
       y="probability density",
       x="|Difference VS Reference|(log2)")



# Plot the equimolar and plasma overlaid
ggplot(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile.EQandPlasm, aes(x=mean.abs.diff, color=pool)) + geom_density(size=1.5, aes(linetype=ORIG.or.TRANSF)) + facet_wrap(~comparison) + theme_bw() + 
  scale_linetype_manual(values=c("dotted", "solid")) +
  scale_y_continuous(expand=c(0.0, 0.0), breaks=seq(0.0, 10.0, 0.5), limits=c(0,3.5)) + 
  scale_x_continuous(expand=c(0.01,0), limits=c(0.0, 13.5), breaks=seq(0,15, 1.0)) + 
  theme(
    axis.text=element_text(size=14, face="bold", color="black"),
    axis.title=element_text(face="bold", size=16),
    plot.title=element_text(face="bold", size=18),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=16, face="bold", color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=14, face="bold"),
    legend.title=element_text(size=15, face="bold"),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(y="probability density",
       x="|Difference VS Reference|(log2)")







#######################################

### ABOVE: Used equimolar pool scaling factors to correct the plasma pool. Now get plasma pool scaling factors.
# Use the filtered set where we required >=1 count in all samples and filtered low count libs, but before we required it to have a corresponding match in the eq pool
# This is what was run:
# ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR.FILT.ZEROS <- ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR[rowSums(ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR$counts>=1)==dim(ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR)[2], ]

# design.plasma <- model.matrix(~0 + lab.libMethod, ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR.FILT.ZEROS$samples)
# colnames(design.plasma) <- unique(ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR.FILT.ZEROS$samples$lab.libMethod)

# NEW--now do the weighted voom like we did with the equimolar pool
v.plasma <- voom(ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR.FILT.ZEROS, design.plasma, normalize.method="scale", plot=TRUE)
lib.methods.design.plasma <- strsplit2(colnames(design.plasma), split=".", fixed=TRUE)[,1]

my.contrasts.fodder.plasma <- sapply(unique(lib.methods.design.plasma), USE.NAMES=FALSE, FUN=function(lib.m){
  mat <- design.plasma
  y.sum <- rowSums(apply(mat, 1, FUN=function(x){
    x*ifelse(lib.methods.design.plasma==lib.m, 1, 0)
  }))
  y.sum/sum(y.sum)
})

all.4n.sum <- as.numeric(rowSums(apply(design.plasma, 1, FUN=function(x){
  x*ifelse(strsplit2(lib.methods.design.plasma, split="_")[,1]=="4N", 1, 0)
})))
my.contrasts.fodder.plasma <- cbind(my.contrasts.fodder.plasma, all.4n.sum/sum(all.4n.sum))
colnames(my.contrasts.fodder) <- c(unique(lib.methods.design), "4N_ABCD")
my.contrasts.plasma <- cbind(truseq.to.4NA=my.contrasts.fodder[, "TruSeq"]-my.contrasts.fodder[,"4N_A"],
                             truseq.to.4NB=my.contrasts.fodder[, "TruSeq"]-my.contrasts.fodder[,"4N_B"],
                             truseq.to.4NC=my.contrasts.fodder[, "TruSeq"]-my.contrasts.fodder[,"4N_C"],
                             truseq.to.4ND=my.contrasts.fodder[, "TruSeq"]-my.contrasts.fodder[,"4N_D"],
                             truseq.to.avg.4N=my.contrasts.fodder[,"TruSeq"]-my.contrasts.fodder[,"4N_ABCD"])

fit.plasm <- lmFit(v.plasma, design.plasma)
fit.plasm <- eBayes(fit.plasm, robust = TRUE, trend = TRUE)
fit2.plasm <- contrasts.fit(fit.plasm, my.contrasts.plasma)
fit2.plasm <- eBayes(fit2.plasm, trend = TRUE, robust = TRUE)
tt.truseq.vs.4N.ABCD.plasm <- topTable(fit2.plasm, 5, n = Inf, confint=0.95, sort.by = "none")
tt.truseq.vs.4N.ALL.plasm <- topTable(fit2.plasm, 1:4, n = Inf)

v.noqual.plasm <- voom(ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR.FILT.ZEROS, design.plasma, normalize.method="scale", plot=TRUE)
lcpm.plasm <- v.noqual.plasm$E
lcpm.truseq.plasm <- lcpm.plasm[, grep("truseq", colnames(lcpm.plasm), ignore.case=TRUE)]
lcpm.4N.plasm <- lcpm.plasm[, grep("4N", colnames(lcpm.plasm), ignore.case=TRUE)]
lcpm.4N.merge.tt.4NALL.plasm <- merge(lcpm.4N.plasm, tt.truseq.vs.4N.ABCD.plasm[, c("logFC", "CI.L", "CI.R")], by=0)
mean.lcpm.truseq.plasm <- rowMeans(lcpm.truseq.plasm)
mean.lcpm.4N.plasm <- rowMeans(lcpm.4N.plasm)
lcpm.4N.transform.plasm <- apply(lcpm.4N.plasm, 2, FUN=function(x){
  x+lcpm.4N.merge.tt.4NALL.plasm$logFC
})
mean.lcpm.4N.transform.plasm <- rowMeans(lcpm.4N.transform.plasm)
merge.mean.plasm <- cbind(mean.lcpm.truseq=mean.lcpm.truseq.plasm, mean.lcpm.4N=mean.lcpm.4N.plasm, mean.lcpm.4N.transform=mean.lcpm.4N.transform.plasm )


colnames(lcpm.4N.transform.plasm) <- paste0(colnames(lcpm.4N.plasm), ".transf")
colnames(lcpm.plasm) <- paste0(colnames(lcpm.plasm), ".orig")
lcpm.merge.transf.plasm <- cbind(lcpm.plasm, lcpm.4N.transform.plasm)
name.split.plasm <- strsplit2(colnames(lcpm.merge.transf.plasm), split=".", fixed=TRUE)
sample.info.plasm <- data.frame(row.names=colnames(lcpm.merge.transf.plasm), lib.method.simple=paste(strsplit2(name.split.plasm[,1], split="_")[,1], name.split.plasm[,5], sep="."))
pheatmap(lcpm.merge.transf.plasm, show_rownames = FALSE, annotation_col = sample.info.plasm)

lcpm.truseq.merge.tt.4NALL.plasm <- merge(lcpm.truseq.plasm, tt.truseq.vs.4N.ABCD.plasm[, c("logFC", "CI.L", "CI.R")], by=0)
lcpm.truseq.transform.plasm <- apply(lcpm.truseq.plasm, 2, FUN=function(x){
  x-lcpm.truseq.merge.tt.4NALL.plasm$logFC
})
colnames(lcpm.truseq.transform.plasm) <- paste0(colnames(lcpm.truseq.transform.plasm), ".transf")
lcpm.merge.transf.plasm <- cbind(lcpm.plasm, lcpm.truseq.transform.plasm)
name.split.plasm <- strsplit2(colnames(lcpm.merge.transf.plasm), split=".", fixed=TRUE)
sample.info.plasm <- data.frame(row.names=colnames(lcpm.merge.transf.plasm), lib.method.simple=paste(strsplit2(name.split.plasm[,1], split="_")[,1], name.split.plasm[,5], sep="."))
pheatmap(lcpm.merge.transf.plasm, show_rownames = FALSE, annotation_col = sample.info.plasm)










# END Cleanup 9-11-17 ----





ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE <- ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE[, ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE$samples$pool=="SynthEQ"]
# filter to get counts found in all samples
ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.FILT <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE[rowSums(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE$counts>0)==dim(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE)[2],]
ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.FILT <- calcNormFactors(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.FILT)
design.eq <- model.matrix(~0+lab.libMethod, ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.FILT$samples)
colnames(design.eq) <- (unique(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.FILT$samples$lab.libMethod))
lib.methods.design <- strsplit2(colnames(design.eq), split=".", fixed=TRUE)[,1]

my.contrasts.fodder <- sapply(unique(lib.methods.design), USE.NAMES=FALSE, FUN=function(lib.m){
  mat <- design.eq
  y.sum <- rowSums(apply(mat, 1, FUN=function(x){
    x*ifelse(lib.methods.design==lib.m, 1, 0)
  }))
  y.sum/sum(y.sum)
})

all.4n.sum <- as.numeric(rowSums(apply(design.eq, 1, FUN=function(x){
  x*ifelse(strsplit2(lib.methods.design, split="_")[,1]=="4N", 1, 0)
})))
my.contrasts.fodder <- cbind(my.contrasts.fodder, all.4n.sum/sum(all.4n.sum))
colnames(my.contrasts.fodder) <- c(unique(lib.methods.design), "4N_ABCD")
my.contrasts <- cbind(truseq.to.4NA=my.contrasts.fodder[, "TruSeq"]-my.contrasts.fodder[,"4N_A"],
                      truseq.to.4NB=my.contrasts.fodder[, "TruSeq"]-my.contrasts.fodder[,"4N_B"],
                      truseq.to.4NC=my.contrasts.fodder[, "TruSeq"]-my.contrasts.fodder[,"4N_C"],
                      truseq.to.4ND=my.contrasts.fodder[, "TruSeq"]-my.contrasts.fodder[,"4N_D"],
                      truseq.to.avg.4N=my.contrasts.fodder[,"TruSeq"]-my.contrasts.fodder[,"4N_ABCD"])


v <- voom(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.FILT, design.eq, normalize.method = "scale",  plot=TRUE)
fit <- lmFit(v, design.eq)
fit <- eBayes(fit, robust = TRUE, trend = TRUE)
fit2 <- contrasts.fit(fit, my.contrasts)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)
tt.truseq.vs.4N.ABCD <- topTable(fit2, 5, n = Inf, confint=0.95, sort.by = "none")
tt.truseq.vs.4N.ALL <- topTable(fit2, 1:4, n = Inf)

v.noqual <- voom(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.FILT, design.eq, normalize.method="scale")
lcpm <- v.noqual$E
lcpm.truseq <- lcpm[, grep("truseq", colnames(lcpm), ignore.case=TRUE)]
lcpm.4N <- lcpm[, grep("4N", colnames(lcpm), ignore.case=TRUE)]
lcpm.4N.merge.tt.4NALL <- merge(lcpm.4N, tt.truseq.vs.4N.ABCD[, c("logFC", "CI.L", "CI.R")], by=0)
mean.lcpm.truseq <- rowMeans(lcpm.truseq)
mean.lcpm.4N <- rowMeans(lcpm.4N)
lcpm.4N.transform <- apply(lcpm.4N, 2, FUN=function(x){
  x+lcpm.4N.merge.tt.4NALL$logFC
})
mean.lcpm.4N.transform <- rowMeans(lcpm.4N.transform)
merge.mean <- cbind(mean.lcpm.truseq=mean.lcpm.truseq, mean.lcpm.4N=mean.lcpm.4N, mean.lcpm.4N.transform=mean.lcpm.4N.transform )


colnames(lcpm.4N.transform) <- paste0(colnames(lcpm.4N.transform), ".transf")
colnames(lcpm) <- paste0(colnames(lcpm), ".orig")
lcpm.merge.transf <- cbind(lcpm, lcpm.4N.transform)
name.split <- strsplit2(colnames(lcpm.merge.transf), split=".", fixed=TRUE)
sample.info <- data.frame(row.names=colnames(lcpm.merge.transf), lib.method.simple=paste(strsplit2(name.split[,1], split="_")[,1], name.split[,5], sep="."))
pheatmap(lcpm.merge.transf, show_rownames = FALSE, annotation_col = sample.info)

lcpm.truseq.merge.tt.4NALL <- merge(lcpm.truseq, tt.truseq.vs.4N.ABCD[, c("logFC", "CI.L", "CI.R")], by=0)
lcpm.truseq.transform <- apply(lcpm.truseq, 2, FUN=function(x){
  x-lcpm.truseq.merge.tt.4NALL$logFC
})
colnames(lcpm.truseq.transform) <- paste0(colnames(lcpm.truseq.transform), ".transf")
lcpm.merge.transf <- cbind(lcpm, lcpm.truseq.transform)
name.split <- strsplit2(colnames(lcpm.merge.transf), split=".", fixed=TRUE)
sample.info <- data.frame(row.names=colnames(lcpm.merge.transf), lib.method.simple=paste(strsplit2(name.split[,1], split="_")[,1], name.split[,5], sep="."))
pheatmap(lcpm.merge.transf, show_rownames = FALSE, annotation_col = sample.info)

# now get plasma
# Drop those with total miR counts <100K and keep only in-house 4N and truseq, as we did with equimolar pools
ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.detail%in%ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.FILT$samples$lib.method.detail & ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.size>100000]
# Can only look at those with counts across the board
ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR.FILT.ZEROS <- ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR[rowSums(ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR$counts>=1)==dim(ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR)[2], ]

# do same transformation as above
design.plasma <- model.matrix(~0 + lab.libMethod, ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR.FILT.ZEROS$samples)
colnames(design.plasma) <- unique(ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR.FILT.ZEROS$samples$lab.libMethod)
v.plasma <- voom(ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR.FILT.ZEROS, design.plasma, normalize.method="scale")
lcpm.plasma.filt <- v.plasma$E[row.names(v.plasma$E)%in%row.names(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.FILT),]

lcpm.truseq.plasma <- lcpm.plasma.filt[, grep("truseq", colnames(lcpm.plasma.filt), ignore.case=TRUE)]

lcpm.truseq.plasma.merge.tt.4NALL <- merge(lcpm.truseq.plasma, tt.truseq.vs.4N.ABCD[, c("logFC", "CI.L", "CI.R")], by=0)

lcpm.truseq.plasma.transform <- apply(lcpm.truseq.plasma, 2, FUN=function(x){
  x.floor <- min(x[x>0])
  x.log.scl <- x-lcpm.truseq.plasma.merge.tt.4NALL$logFC
  return(ifelse(x.log.scl<x.floor, x.floor, x.log.scl))
})
colnames(lcpm.truseq.plasma.transform) <- paste0(colnames(lcpm.truseq.plasma.transform), ".transf")
colnames(lcpm.plasma.filt) <- paste0(colnames(lcpm.plasma.filt), ".orig")
lcpm.plasma.merge.transf <- data.frame(merge(lcpm.plasma.filt, lcpm.truseq.plasma.transform, by=0), row.names=1)
name.split.plasma <- strsplit2(colnames(lcpm.plasma.merge.transf), split=".", fixed=TRUE)
sample.info.plasm <- data.frame(row.names=colnames(lcpm.plasma.merge.transf), lib.method.simple=paste(strsplit2(name.split.plasma[,1], split="_")[,1], name.split.plasma[,5], sep="."))
lcpm.plasma.merge.transf.scale <- normalizeBetweenArrays(lcpm.plasma.merge.transf, method = "scale")
lcpm.plasma.merge.transf.quantile <- normalizeBetweenArrays(lcpm.plasma.merge.transf, method = "quantile")
pheatmap(lcpm.plasma.merge.transf.scale, show_rownames = FALSE, annotation_col = sample.info.plasm, main = "PLASMA POOL CORRECTIONS: MEDIAN CENTRED", labels_col=sub("^X", "", sub("PlasmaPool.", "", row.names(sample.info.plasm))))
pheatmap(lcpm.plasma.merge.transf.quantile, show_rownames = FALSE, annotation_col = sample.info.plasm, main = "PLASMA POOL CORRECTIONS: QUANTILE NORM", labels_col=sub("^X", "", sub("PlasmaPool.", "", row.names(sample.info.plasm))))

group.id.split <- strsplit2(colnames(lcpm.plasma.merge.transf.quantile), split=".", fixed=TRUE)
group.id <- paste(strsplit2(group.id.split[,1], split="_")[,1], group.id.split[,dim(group.id.split)[2]], sep="_")
n.samples.per.group <- as.numeric(table(group.id))
lcpm.plasma.merge.transf.quantile.sum <- sumTechReps(lcpm.plasma.merge.transf.quantile, ID = group.id)
lcpm.plasma.merge.transf.quantile.mean <- sapply(1:dim(lcpm.plasma.merge.transf.quantile.sum)[2], FUN=function(x){ y<-lcpm.plasma.merge.transf.quantile.sum[,x]/n.samples.per.group[x]; return(y) })
colnames(lcpm.plasma.merge.transf.quantile.mean ) <- colnames(lcpm.plasma.merge.transf.quantile.sum)
lcpm.plasma.merge.transf.quantile.mean.dt <- data.table(data.frame(miR.id=row.names(lcpm.plasma.merge.transf.quantile.mean), lcpm.plasma.merge.transf.quantile.mean))
lcpm.plasma.merge.transf.quantile.mean.dt[, `:=`(orig.4N.vs.orig.TruSeq = ((X4N_orig-TruSeq_orig)),
                                                 orig.4N.vs.transf.Truseq = ((X4N_orig-TruSeq_transf)) ,
                                                 orig.TruSeq.vs.trans.Truseq = ((TruSeq_transf-TruSeq_orig)) )]
lcpm.plasma.merge.transf.quantile.mean.dt.melt <- melt(lcpm.plasma.merge.transf.quantile.mean.dt, measure.vars = c("orig.4N.vs.orig.TruSeq", "orig.4N.vs.transf.Truseq", "orig.TruSeq.vs.trans.Truseq"), variable.name = "comparison", value.name = "abs.difference")
ggplot(lcpm.plasma.merge.transf.quantile.mean.dt.melt, aes(y=abs.difference, x=comparison, fill=comparison)) + geom_boxplot()


lcpm.plasma.merge.transf.scale.sum <- sumTechReps(lcpm.plasma.merge.transf.scale, ID = group.id)
lcpm.plasma.merge.transf.scale.mean <- sapply(1:dim(lcpm.plasma.merge.transf.scale.sum)[2], FUN=function(x){ y<-lcpm.plasma.merge.transf.scale.sum[,x]/n.samples.per.group[x]; return(y) })
colnames(lcpm.plasma.merge.transf.scale.mean ) <- colnames(lcpm.plasma.merge.transf.scale.sum)
lcpm.plasma.merge.transf.scale.mean.dt <- data.table(data.frame(miR.id=row.names(lcpm.plasma.merge.transf.scale.mean), lcpm.plasma.merge.transf.scale.mean))
lcpm.plasma.merge.transf.scale.mean.dt[, `:=`(orig.4N.vs.orig.TruSeq = ((X4N_orig-TruSeq_orig)),
                                              orig.4N.vs.transf.Truseq = ((X4N_orig-TruSeq_transf)) ,
                                              orig.TruSeq.vs.trans.Truseq = ((TruSeq_transf-TruSeq_orig)) )]
lcpm.plasma.merge.transf.scale.mean.dt.melt <- melt(lcpm.plasma.merge.transf.scale.mean.dt, measure.vars = c("orig.4N.vs.orig.TruSeq", "orig.4N.vs.transf.Truseq", "orig.TruSeq.vs.trans.Truseq"), variable.name = "comparison", value.name = "abs.difference")
ggplot(lcpm.plasma.merge.transf.scale.mean.dt.melt, aes(y=abs.difference, x=comparison, fill=comparison)) + geom_boxplot()


lcpm.plasma.merge.transf.noNorm <- lcpm.plasma.merge.transf
lcpm.plasma.merge.transf.noNorm.sum <- sumTechReps(lcpm.plasma.merge.transf.noNorm, ID = group.id)
lcpm.plasma.merge.transf.noNorm.mean <- sapply(1:dim(lcpm.plasma.merge.transf.noNorm.sum)[2], FUN=function(x){ y<-lcpm.plasma.merge.transf.noNorm.sum[,x]/n.samples.per.group[x]; return(y) })
colnames(lcpm.plasma.merge.transf.noNorm.mean ) <- colnames(lcpm.plasma.merge.transf.noNorm.sum)
lcpm.plasma.merge.transf.noNorm.mean.dt <- data.table(data.frame(miR.id=row.names(lcpm.plasma.merge.transf.noNorm.mean), lcpm.plasma.merge.transf.noNorm.mean))
lcpm.plasma.merge.transf.noNorm.mean.dt[, `:=`(orig.4N.vs.orig.TruSeq = abs(X4N_orig-TruSeq_orig)/((X4N_orig+TruSeq_orig)/2),
                                               orig.4N.vs.transf.Truseq = abs(X4N_orig-TruSeq_transf)/((X4N_orig+TruSeq_transf)/2) ,
                                               orig.TruSeq.vs.trans.Truseq = abs(TruSeq_transf-TruSeq_orig)/((TruSeq_transf+TruSeq_orig)/2))]
lcpm.plasma.merge.transf.noNorm.mean.dt.melt <- melt(lcpm.plasma.merge.transf.noNorm.mean.dt, measure.vars = c("orig.4N.vs.orig.TruSeq", "orig.4N.vs.transf.Truseq", "orig.TruSeq.vs.trans.Truseq"), variable.name = "comparison", value.name = "abs.percent.difference")
ggplot(lcpm.plasma.merge.transf.noNorm.mean.dt.melt, aes(y=abs.percent.difference, x=comparison, fill=comparison)) + geom_boxplot() + theme_bw()
ggplot(lcpm.plasma.merge.transf.noNorm.mean.dt.melt, aes(x=abs.percent.difference, color=comparison)) + geom_density(size=1.5) + theme_bw() + scale_y_continuous(expand=c(0.0, 0.0), breaks=seq(0.0, 10, 0.5), limits=c(0.0, 6.0)) + scale_x_continuous(expand=c(0,0), limits=c(0.0,1.05), breaks=seq(0,1,.2), labels = scales::percent) + theme(
  axis.text=element_text(size=14, face="bold", color="black"),
  axis.title=element_text(face="bold", size=16),
  plot.title=element_text(face="bold", size=18),
  panel.grid.major = element_line(color=NA),
  panel.border = element_rect(size=1.5, color="black"),
  axis.ticks = element_line(color="black", size=1.5),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.text=element_text(size=14, face="bold"),
  legend.title=element_text(size=15, face="bold"),
  legend.key.size=unit(1.2, "lines")
) +
  labs(y="probability density",
       x="|% Difference|")
ggsave("Comparison_ScalingFactors_DensityPlot.png")

ggplot(lcpm.plasma.merge.transf.noNorm.mean.dt.melt[comparison!="orig.TruSeq.vs.trans.Truseq"], aes(x=abs.percent.difference, color=comparison)) + geom_density(size=1.5) + theme_bw() + scale_y_continuous(expand=c(0.0, 0.0), breaks=seq(0.0, 6, 0.5), limits=c(0.0, 6.0)) + scale_x_continuous(expand=c(0,0), limits=c(0.0,1.05), breaks=seq(0,1,.2), labels = scales::percent) + theme(
  axis.text=element_text(size=14, face="bold", color="black"),
  axis.title=element_text(face="bold", size=16),
  plot.title=element_text(face="bold", size=18),
  panel.grid.major = element_line(color=NA),
  panel.border = element_rect(size=1.5, color="black"),
  axis.ticks = element_line(color="black", size=1.5),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.text=element_text(size=14, face="bold"),
  legend.title=element_text(size=15, face="bold"),
  legend.key.size=unit(1.2, "lines")
) +
  labs(y="probability density",
       x="|% Difference|")
ggsave("Comparison_ScalingFactors_DensityPlot.png")
ggsave("Comparison_ScalingFactors_DensityPlot_drop_truseq_compare.png")

#####################################################
### Now we want to see what the scaling factors look like for the other methods if you start with the equimolar pool. 
### Filtered out the NEB and CleanTag samples early, so need to re-make here... 

ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.CAST.ALL <- data.frame(dcast.data.table(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE,  miR.ID~lab.libMethod.pool.replicate, value.var="multimapAdjustedReadCount", fill=0, fun.aggregate = "max"), row.names=1)

sample.info.all <- data.frame(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE[, .(n.miRs=.N), by=.(lab.libMethod.pool.replicate, lab.libMethod.pool, lab.libMethod, Lab, lib.method.detail, lib.method.simple, lib.method.old2, pool, replicate)], row.names=1)
row.names(sample.info.all) <- make.names(row.names(sample.info.all))

sample.info.sort.all <- sample.info.all[colnames(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.CAST.ALL),]
ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL <- DGEList(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.CAST.ALL)
ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL$samples <- cbind(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL$samples, sample.info.sort.all)
# Now repeat what we did above.
drop.4n.libs <- c("4N_NEXTflex", "4N_Xu")
# Get equimolar pool samples and remove 4N NextFlex and 4N-Xu samples, as well as any library with <100k miRNA reads.
ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR <- ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL[, ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL$samples$pool=="SynthEQ" & (!ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL$samples$lib.method.detail%in%drop.4n.libs) & ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL$samples$lib.size>=100000]
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA <- ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL[, ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL$samples$pool=="PlasmaPool" & (!ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL$samples$lib.method.detail%in%drop.4n.libs) & ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL$samples$lib.size>=100000]
# Now this is where we need to diverge to include the other methods. 
# The thought is that scaling factors are going to convert between two different methods, so all of our filtering, normalization, centering and testing should 
# be done without the potential influence of the other library prep types.






lib.comparison.groups <- list(TruSeq.V.4N=c("TruSeq", "4N"), 
                              TruSeq.V.NEB=c("TruSeq", "NEBNext"),
                              TruSeq.V.CleanTag=c("TruSeq", "CleanTag"),
                              NEB.V.CleanTag=c("NEBNext", "CleanTag"),
                              NEB.V.4N=c("NEBNext", "4N"),
                              CleanTag.V.4N=c("CleanTag", "4N")
)
# Do the first half just to get the scaling factors into a data table format
ALL.COMPARISON.SCALING.FACTOR.DT <- do.call(rbind, 
                                            sapply(names(lib.comparison.groups), USE.NAMES=TRUE, simplify = FALSE, FUN=function(x.id.list){
                                              x.ids <- unlist(lib.comparison.groups[x.id.list], recursive = TRUE, use.names = FALSE)
                                              dge.comp <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[, ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.simple%in%x.ids]
                                              design.comp <- model.matrix(~0+lab.libMethod, dge.comp$samples)
                                              colnames(design.comp) <- unique(dge.comp$samples$lab.libMethod)
                                              lib.methods.comp <- strsplit2(strsplit2(colnames(design.comp), split=".", fixed=TRUE)[,1], split="_")[,1]
                                              comp.contrasts.fodder <- sapply(unique(lib.methods.comp), USE.NAMES=FALSE, FUN=function(lib.m){
                                                mat <- design.comp
                                                y.sum <- rowSums(apply(mat, 1, FUN=function(x){
                                                  x*ifelse(lib.methods.comp==lib.m, 1, 0)
                                                }))
                                                y.sum/sum(y.sum)
                                              })
                                              colnames(comp.contrasts.fodder) <- unique(lib.methods.comp)
                                              out.contrast <- as.numeric(comp.contrasts.fodder[, x.ids[1]]-comp.contrasts.fodder[, x.ids[2]])
                                              # Now filter miRs
                                              dge.comp.filt <- dge.comp[rowSums(dge.comp$counts>0)==dim(dge.comp)[2],]
                                              # v.comp <- voomWithQualityWeights(dge.comp.filt, design.comp, normalize.method = "scale")
                                              # Drop the quality weighting, just for simplicity
                                              v.comp <- voom(dge.comp.filt, design.comp, normalize.method = "scale")
                                              fit.comp <- lmFit(v.comp, design.comp)
                                              fit.comp <- eBayes(fit.comp, robust = TRUE, trend = TRUE)
                                              fit2.comp <- contrasts.fit(fit.comp, out.contrast)
                                              fit2.comp <- eBayes(fit2.comp, trend = TRUE, robust = TRUE)
                                              tt.comp <- data.frame(topTable(fit2.comp, n = Inf, confint=0.95, sort.by = "none"))
                                              tt.comp.dt <- data.table(data.frame(miR.id.full=row.names(tt.comp), tt.comp))
                                              tt.comp.dt[, comparison:=x.id.list]
                                              return(tt.comp.dt)
                                            })
)

# Repeat the first half, but continue to merge with the plasma data
library(ModelMetrics)
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT <- do.call(rbind,  sapply(names(lib.comparison.groups), USE.NAMES=TRUE, simplify = FALSE, FUN=function(x.id.list){
  x.ids <- unlist(lib.comparison.groups[x.id.list], recursive = TRUE, use.names = FALSE)
  dge.comp <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[, ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.simple%in%x.ids]
  design.comp <- model.matrix(~0+lab.libMethod, dge.comp$samples)
  colnames(design.comp) <- unique(dge.comp$samples$lab.libMethod)
  lib.methods.comp <- strsplit2(strsplit2(colnames(design.comp), split=".", fixed=TRUE)[,1], split="_")[,1]
  comp.contrasts.fodder <- sapply(unique(lib.methods.comp), USE.NAMES=FALSE, FUN=function(lib.m){
    mat <- design.comp
    y.sum <- rowSums(apply(mat, 1, FUN=function(x){
      x*ifelse(lib.methods.comp==lib.m, 1, 0)
    }))
    y.sum/sum(y.sum)
  })
  colnames(comp.contrasts.fodder) <- unique(lib.methods.comp)
  out.contrast <- as.numeric(comp.contrasts.fodder[, x.ids[1]]-comp.contrasts.fodder[, x.ids[2]])
  # Now filter miRs
  dge.comp.filt <- dge.comp[rowSums(dge.comp$counts>0)==dim(dge.comp)[2],]
  # v.comp <- voomWithQualityWeights(dge.comp.filt, design.comp, normalize.method = "scale")
  # Drop the quality weighting, just for simplicity
  v.comp <- voom(dge.comp.filt, design.comp, normalize.method = "scale")
  fit.comp <- lmFit(v.comp, design.comp)
  fit.comp <- eBayes(fit.comp, robust = TRUE, trend = TRUE)
  fit2.comp <- contrasts.fit(fit.comp, out.contrast)
  fit2.comp <- eBayes(fit2.comp, trend = TRUE, robust = TRUE)
  tt.comp <- data.frame(topTable(fit2.comp, n = Inf, confint=0.95, sort.by = "none"))
  tt.comp.dt <- data.table(data.frame(miR.id.full=row.names(tt.comp), tt.comp))
  tt.comp.dt[, comparison:=x.id.list]
  tt.comp.scl.factor <- tt.comp[, c("logFC", "CI.L", "CI.R")]
  
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple%in%x.ids]
  keep.mirs <- (rowSums(dge.test$counts>=1)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  v.test.norm <- voom(dge.test.filt, design.test, normalize.method="scale")
  cpm.test.norm <- v.test.norm$E
  
  merged.miR.ids <- merge(cpm.test.norm, tt.comp.scl.factor, by=0)$Row.names
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  cpm.test.norm.filt.md1 <- cpm.test.norm.filt[, grep(x.ids[1], colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)]
  tt.comp.scl.factor.filt <- tt.comp.scl.factor[merged.miR.ids,]
  
  # set a floor min. Logic here is that if we had a miR that had a count of 1 and the scaling factor 
  # indicates a drop of -5.0 fold, we do not want to make the new log2 count log2(1) - 5.0, since our scale will be way off
  # for plotting. So, if the value is below the min, set it to the min value. In other words. Limit the amount 
  # of downward scaling to the minimum value of the sample.
  #cpm.test.norm.filt.scaled <- apply(cpm.test.norm.filt, 2, FUN=function(x){
  #  x.floor <- min(x[x>0])
  #  x.log.scl <- x-tt.comp.scl.factor.filt$logFC
  #  return(ifelse(x.log.scl<x.floor, x.floor, x.log.scl))
  #})
  
  # this does it without a floor
  # Decided to use median scaling to get the scaling factors
  # Then after scaling, use quantile normalization to match up the distributions
  # Thinking that this might be a more appropriate way to handle the fact that
  # some of the scaling at the low end will lead to values less than the minimum. 
  cpm.test.norm.filt.md.1.scaled <- apply(cpm.test.norm.filt.md1, 2, FUN=function(x){
    #x.floor <- min(x[x>0])
    x.log.scl <- x-tt.comp.scl.factor.filt$logFC
    #return(ifelse(x.log.scl<x.floor, x.floor, x.log.scl))
    return(x.log.scl)
  })
  colnames(cpm.test.norm.filt) <- paste0(colnames(cpm.test.norm.filt), ".orig")
  colnames(cpm.test.norm.filt.md.1.scaled) <- paste0(colnames(cpm.test.norm.filt.md.1.scaled), ".transf")
  cpm.test.orig.and.scaled <- cbind(cpm.test.norm.filt, cpm.test.norm.filt.md.1.scaled)
  cpm.test.orig.and.scaled.rescaled <- normalizeBetweenArrays(cpm.test.orig.and.scaled, method = "quantile")
  ids.reference <- grep(x.ids[1], colnames(cpm.test.orig.and.scaled.rescaled), value=TRUE, ignore.case=TRUE, invert = TRUE)
  ids.test.orig <- grep(x.ids[1], colnames(cpm.test.norm.filt), value=TRUE, ignore.case=TRUE)
  ids.test.scaled <- colnames(cpm.test.norm.filt.md.1.scaled)
  # Return this to get a matrix of transformed and original values
  # return(cpm.test.orig.and.scaled.rescaled)
  # Now get an average of the values we did not transform
  cpm.rescaled.test.orig <- cpm.test.orig.and.scaled.rescaled[, ids.test.orig]
  cpm.rescaled.test.scaled <- cpm.test.orig.and.scaled.rescaled[, ids.test.scaled]
  cpm.rescaled.ref.orig <- cpm.test.orig.and.scaled.rescaled[, ids.reference]
  cpm.mean.reference <- rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.reference])
  
  # Check scale norm too. Didnt look quite as good, so removed.
  cpm.test.orig.and.scaled.rescaled.med <- normalizeBetweenArrays(cpm.test.orig.and.scaled, method = "scale")
  cpm.rescaled.test.orig.med <- cpm.test.orig.and.scaled.rescaled.med[, ids.test.orig]
  cpm.rescaled.test.scaled.med <- cpm.test.orig.and.scaled.rescaled.med[, ids.test.scaled]
  cpm.rescaled.ref.orig.med <- cpm.test.orig.and.scaled.rescaled.med[, ids.reference]
  cpm.mean.reference.med <- rowMeans(cpm.test.orig.and.scaled.rescaled.med[, ids.reference])
  
  
  cpm.test.orig.and.scaled.rescaled.mean <- data.table(data.frame(miRNA.full.ID=row.names(cpm.test.orig.and.scaled.rescaled), 
                                                                  scaling.factor.logFC=tt.comp.scl.factor.filt$logFC,
                                                                  cpm.mean.reference.quantile = cpm.mean.reference,      
                                                                  cpm.mean.test.orig.quantile = rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.test.orig]),
                                                                  cpm.mean.test.scaled.quantile = rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.test.scaled]),
                                                                  abs.dif.orig.quantile = rowMeans(apply(cpm.rescaled.test.orig, 2, FUN=function(x){ abs(x-cpm.mean.reference)})),
                                                                  abs.dif.transf.quantile = rowMeans(apply(cpm.rescaled.test.scaled, 2, FUN=function(x){ abs(x-cpm.mean.reference)}))
  ))
  cpm.test.orig.and.scaled.rescaled.mean[, comparison:=x.id.list]
  return(cpm.test.orig.and.scaled.rescaled.mean)
  #return(cpm.test.orig.and.scaled.rescaled)
})); #ggplot(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.LIST[[3]]) + geom_density(aes(color="VS.ORIG", x=abs.dif.orig)) + geom_density(aes(color="VS.TRANSF", x=abs.dif.transf))
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif <- melt(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT, measure.vars=c("abs.dif.orig.quantile", "abs.dif.transf.quantile"), variable.name = "ORIG.or.TRANSF.and.SCL", value.name = "mean.abs.diff")
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif[, SCALE.METHOD:=tstrsplit(ORIG.or.TRANSF.and.SCL, split=".", fixed=TRUE)[4]]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif[, ORIG.or.TRANSF:=sub(paste0(".", SCALE.METHOD), "", ORIG.or.TRANSF.and.SCL)]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif[, `:=`(SCALE.METHOD=ifelse(SCALE.METHOD=="quantile", "QUANTILE", "MEDIAN"), ORIG.or.TRANSF=ifelse(ORIG.or.TRANSF=="abs.dif.orig", "ORIGINAL", "TRANSFORMED"))]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile <- subset(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif, SCALE.METHOD=="QUANTILE")
ggplot(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile, aes(x=mean.abs.diff, color=ORIG.or.TRANSF)) + geom_density(size=1.5) + facet_wrap(~comparison) + theme_bw() + 
  scale_y_continuous(expand=c(0.0, 0.0), breaks=seq(0.0, 1.0, 0.1), limits=c(0.0, 0.9)) + 
  scale_x_continuous(expand=c(0.01,0), limits=c(0.0, 9.0), breaks=seq(0,9, 1.0)) + 
  theme(
    axis.text=element_text(size=14, face="bold", color="black"),
    axis.title=element_text(face="bold", size=16),
    plot.title=element_text(face="bold", size=18),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=16, face="bold", color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=14, face="bold"),
    legend.title=element_text(size=15, face="bold"),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(y="probability density",
       x="|Difference VS Reference|(log2)")

# Get the matrices for plotting

ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.LIST <- sapply(names(lib.comparison.groups), USE.NAMES=TRUE, simplify = FALSE, FUN=function(x.id.list){
  x.ids <- unlist(lib.comparison.groups[x.id.list], recursive = TRUE, use.names = FALSE)
  dge.comp <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[, ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.simple%in%x.ids]
  design.comp <- model.matrix(~0+lab.libMethod, dge.comp$samples)
  colnames(design.comp) <- unique(dge.comp$samples$lab.libMethod)
  lib.methods.comp <- strsplit2(strsplit2(colnames(design.comp), split=".", fixed=TRUE)[,1], split="_")[,1]
  comp.contrasts.fodder <- sapply(unique(lib.methods.comp), USE.NAMES=FALSE, FUN=function(lib.m){
    mat <- design.comp
    y.sum <- rowSums(apply(mat, 1, FUN=function(x){
      x*ifelse(lib.methods.comp==lib.m, 1, 0)
    }))
    y.sum/sum(y.sum)
  })
  colnames(comp.contrasts.fodder) <- unique(lib.methods.comp)
  out.contrast <- as.numeric(comp.contrasts.fodder[, x.ids[1]]-comp.contrasts.fodder[, x.ids[2]])
  # Now filter miRs
  dge.comp.filt <- dge.comp[rowSums(dge.comp$counts>0)==dim(dge.comp)[2],]
  v.comp <- voom(dge.comp.filt, design.comp, normalize.method = "scale")
  fit.comp <- lmFit(v.comp, design.comp)
  fit.comp <- eBayes(fit.comp, robust = TRUE, trend = TRUE)
  fit2.comp <- contrasts.fit(fit.comp, out.contrast)
  fit2.comp <- eBayes(fit2.comp, trend = TRUE, robust = TRUE)
  tt.comp <- data.frame(topTable(fit2.comp, n = Inf, confint=0.95, sort.by = "none"))
  tt.comp.dt <- data.table(data.frame(miR.id.full=row.names(tt.comp), tt.comp))
  tt.comp.dt[, comparison:=x.id.list]
  tt.comp.scl.factor <- tt.comp[, c("logFC", "CI.L", "CI.R")]
  
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple%in%x.ids]
  
  keep.mirs <- (rowSums(dge.test$counts>=1)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  v.test.norm <- voom(dge.test.filt, design.test, normalize.method="scale")
  cpm.test.norm <- v.test.norm$E
  
  merged.miR.ids <- merge(cpm.test.norm, tt.comp.scl.factor, by=0)$Row.names
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  cpm.test.norm.filt.md1 <- cpm.test.norm.filt[, grep(x.ids[1], colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)]
  tt.comp.scl.factor.filt <- tt.comp.scl.factor[merged.miR.ids,]
  
  # set a floor min. Logic here is that if we had a miR that had a count of 1 and the scaling factor 
  # indicates a drop of -5.0 fold, we do not want to make the new log2 count log2(1) - 5.0, since our scale will be way off
  # for plotting. So, if the value is below the min, set it to the min value. In other words. Limit the amount 
  # of downward scaling to the minimum value of the sample.
  #cpm.test.norm.filt.scaled <- apply(cpm.test.norm.filt, 2, FUN=function(x){
  #  x.floor <- min(x[x>0])
  #  x.log.scl <- x-tt.comp.scl.factor.filt$logFC
  #  return(ifelse(x.log.scl<x.floor, x.floor, x.log.scl))
  #})
  
  # this does it without a floor
  # Decided to use median scaling to get the scaling factors
  # Then after scaling, use quantile normalization to match up the distributions
  # Thinking that this might be a more appropriate way to handle the fact that
  # some of the scaling at the low end will lead to values less than the minimum. 
  cpm.test.norm.filt.md.1.scaled <- apply(cpm.test.norm.filt.md1, 2, FUN=function(x){
    #x.floor <- min(x[x>0])
    x.log.scl <- x-tt.comp.scl.factor.filt$logFC
    #return(ifelse(x.log.scl<x.floor, x.floor, x.log.scl))
    return(x.log.scl)
  })
  ap <- paste(x.ids, collapse=".to.")
  colnames(cpm.test.norm.filt.md.1.scaled) <- paste(colnames(cpm.test.norm.filt.md.1.scaled), ap, sep="_corrected_")
  
  # return the scaled matrices (not rescaled)
  return(cpm.test.norm.filt.md.1.scaled)
})
# Now check to make sure miRs are in all the files
miRs.in.all <- table(unlist(lapply(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.LIST, row.names)))==length(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.LIST)
# get miR ids found in all of the scaled count matrices
miRs.in.all.ids <- names(miRs.in.all)[miRs.in.all==TRUE]
# filter and merge combine the matrices based on those miR IDS.
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR <- do.call(cbind, lapply(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.LIST, FUN=function(x){ x[miRs.in.all.ids, ]}))
# re-merge with unscaled values
# first make sure we have voom-transformed values
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.FILT.SCALED.MIRS <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[miRs.in.all.ids,]
design.plasm <- model.matrix(~0+lab.libMethod, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.FILT.SCALED.MIRS$samples)
v.ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.FILT.SCALED.MIRS.norm <- voom(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.FILT.SCALED.MIRS, design.plasm)
cpm.ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.FILT.SCALED.MIRS.norm <- v.ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.FILT.SCALED.MIRS.norm$E

ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG <- cbind(cpm.ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.FILT.SCALED.MIRS.norm, ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR)
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE <- normalizeBetweenArrays(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG, method="quantile")
# Now get the annotation information
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO <- data.frame(strsplit2(colnames(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE), split="_corrected_"), stringsAsFactors = FALSE)
colnames(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO) <- c("lab.libMethod.replicate", "correctionType")
row.names(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO) <- colnames(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE)

ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO$correctionType <- ifelse(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO$correctionType=="", "original", as.character(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO$correctionType))
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO$lab.libMethod.replicate.transf <- row.names(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO)
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG <- data.table(merge(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.FILT.SCALED.MIRS$samples, by.x=1, by.y=0, all.x = TRUE))
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG[, lib.method.simple.orig:=lib.method.simple]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG[, lib.method.simple.transf:=tstrsplit(correctionType, split=".to.", fixed=TRUE)[2]]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG[, lib.method.simple.transf:=ifelse(is.na(lib.method.simple.transf), lib.method.simple, lib.method.simple.transf)]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG[, is.corrected:=ifelse(correctionType=="original", "ORIGINAL", "TRANSFORMED")]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG[, `:=`(lib.method.simple.transf=factor(lib.method.simple.transf, levels=lib.simple.levels), lib.method.simple.orig=factor(lib.method.simple.orig, levels=lib.simple.levels))]

ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG[, transf.LibMethod.id:=sub("X4N", "4N", sub(".PlasmaPool", "", make.names(paste(lab.libMethod.replicate, substr(is.corrected, 1, 4), sep="."), unique = TRUE)))]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df <- data.frame(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG, row.names=ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG$lab.libMethod.replicate.transf)
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df <- ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df[colnames(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE),]                                                                                     



ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.for.plot <- ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df[, c("lib.method.simple.transf", "lib.method.simple.orig", "transf.LibMethod.id")]

ann.colors.plasm = list(
  lib.method.simple.orig = hue_pal()(4),
  lib.method.simple.transf = hue_pal()(4),
  is.corrected = grey_pal()(2)
)

names(ann.colors.plasm$lib.method.simple.orig) <- lib.simple.levels
names(ann.colors.plasm$lib.method.simple.transf) <- lib.simple.levels
names(ann.colors.plasm$is.corrected) <- c("TRANSFORMED", "ORIGINAL")
pheatmap(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE,
         show_rownames = FALSE,
         fontsize = 14,
         fontsize_col = 12,
         labels_col = ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.for.plot$transf.LibMethod.id,
         annotation_col = ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.for.plot[, 1:2], 
         annotation_colors=ann.colors.plasm 
)


pheatmap(
  cor(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE, use="pair", method = "spearman")^2,
  show_rownames = TRUE, 
  fontsize = 8,
  fontsize_row=8,
  annotation_col = ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.for.plot[, 1:2],
  annotation_row = ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.for.plot[, 1:2],
  annotation_colors=ann.colors.plasm, 
  labels_col = ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.for.plot$transf.LibMethod.id,
  labels_row =ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.for.plot$transf.LibMethod.id)

# Now sort by row abundance and don't cluster
comparison.groups <- apply(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.for.plot[, 1:2], 1, paste, collapse=".to.")
unique.pairs.groups <- unique(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.for.plot[, 1:2])
row.names(unique.pairs.groups) <- apply(unique.pairs.groups[, 1:2], 1, paste, collapse=".to.")
row.annot.scaling.factor.miRs <- FUNCTION.medianTechReps(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE, ID = comparison.groups)
row.annot.scaling.factor.miRs.sort <- row.annot.scaling.factor.miRs[order(-rowMeans(row.annot.scaling.factor.miRs)),]
colnames(unique.pairs.groups) <- c("Original", "Corrected")
pheatmap(row.annot.scaling.factor.miRs,
         border_color = "black",
         show_rownames = FALSE,
         fontsize = 14,
         fontsize_col = 12,
         annotation_col = unique.pairs.groups, 
         annotation_colors=ann.colors.plasm 
)


# subset the truseq 4n comparison
truseq.4n.comparison.cols <- grep("NEBNext|CleanTag", row.names(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.for.plot), invert=TRUE, value=TRUE, ignore.case=TRUE)
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.truSeq.4N <- ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE[, truseq.4n.comparison.cols] 
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.for.plot.truseq.4N <- ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.for.plot[truseq.4n.comparison.cols, ]

orig.4n.ids <- row.names(subset(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.for.plot.truseq.4N, lib.method.simple.transf=="4N" & lib.method.simple.orig=="4N"))
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.truSeq.4N.sort <- ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.truSeq.4N[order(-rowMeans(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.truSeq.4N[,orig.4n.ids])), ]
comparison.groups.4N.Truseq <- apply(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.for.plot.truseq.4N[, 1:2], 1, FUN=function(x){
  ifelse(x[1]==x[2], paste0(x[1], ".orig"), paste(x[1], x[2], sep=".to."))})

row.annot.miR.4N.truseq <- FUNCTION.medianTechReps(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.truSeq.4N.sort, ID = comparison.groups.4N.Truseq)
colnames(row.annot.miR.4N.truseq) <- make.names(colnames(row.annot.miR.4N.truseq))
row.annot.miR.4N.truseq.2 <- ifelse(row.annot.miR.4N.truseq>0, 1, 0)
row.names(row.annot.miR.4N.truseq.2) <- make.names(row.names(row.annot.miR.4N.truseq.2))
row.annot.miR.4N.truseq.2 <- ifelse(row.annot.miR.4N.truseq>0, 1, 0)
pheatmap(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.truSeq.4N.sort,
         cluster_rows = FALSE,
         show_rownames = FALSE,
         fontsize = 14,
         fontsize_col = 12,
         labels_col = sub(".[12]$", "SF", ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.for.plot.truseq.4N$transf.LibMethod.id),
         annotation_col = ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.for.plot.truseq.4N[, 1:2], 
         annotation_colors=ann.colors.plasm 
)





######## NOW ALSO SHOW CORRECTIONS FOR EQUIMOLAR POOL WITH EQUIMOLAR SEQS
# Now this is where we need to diverge to include the other methods. 
# The thought is that scaling factors are going to convert between two different methods, so all of our filtering, normalization, centering and testing should 
# be done without the potential influence of the other library prep types.

lib.comparison.groups <- list(TruSeq.V.4N=c("TruSeq", "4N"), 
                              TruSeq.V.NEB=c("TruSeq", "NEBNext"),
                              TruSeq.V.CleanTag=c("TruSeq", "CleanTag"),
                              NEB.V.CleanTag=c("NEBNext", "CleanTag"),
                              NEB.V.4N=c("NEBNext", "4N"),
                              CleanTag.V.4N=c("CleanTag", "4N")
)


ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.EQUIMOLAR <- do.call(rbind, sapply(names(lib.comparison.groups), USE.NAMES=TRUE, simplify = FALSE, FUN=function(x.id.list){
  x.ids <- unlist(lib.comparison.groups[x.id.list], recursive = TRUE, use.names = FALSE)
  dge.comp <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[, ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.simple%in%x.ids]
  design.comp <- model.matrix(~0+lab.libMethod, dge.comp$samples)
  colnames(design.comp) <- unique(dge.comp$samples$lab.libMethod)
  lib.methods.comp <- strsplit2(strsplit2(colnames(design.comp), split=".", fixed=TRUE)[,1], split="_")[,1]
  comp.contrasts.fodder <- sapply(unique(lib.methods.comp), USE.NAMES=FALSE, FUN=function(lib.m){
    mat <- design.comp
    y.sum <- rowSums(apply(mat, 1, FUN=function(x){
      x*ifelse(lib.methods.comp==lib.m, 1, 0)
    }))
    y.sum/sum(y.sum)
  })
  colnames(comp.contrasts.fodder) <- unique(lib.methods.comp)
  out.contrast <- as.numeric(comp.contrasts.fodder[, x.ids[1]]-comp.contrasts.fodder[, x.ids[2]])
  # Now filter miRs
  dge.comp.filt <- dge.comp[rowSums(dge.comp$counts>0)==dim(dge.comp)[2],]
  v.comp <- voom(dge.comp.filt, design.comp, normalize.method = "scale")
  fit.comp <- lmFit(v.comp, design.comp)
  fit.comp <- eBayes(fit.comp, robust = TRUE, trend = TRUE)
  fit2.comp <- contrasts.fit(fit.comp, out.contrast)
  fit2.comp <- eBayes(fit2.comp, trend = TRUE, robust = TRUE)
  tt.comp <- data.frame(topTable(fit2.comp, n = Inf, confint=0.95, sort.by = "none"))
  tt.comp.dt <- data.table(data.frame(miR.id.full=row.names(tt.comp), tt.comp))
  tt.comp.dt[, comparison:=x.id.list]
  tt.comp.scl.factor <- tt.comp[, c("logFC", "CI.L", "CI.R")]
  
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[, ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.simple%in%x.ids]
  keep.mirs <- (rowSums(dge.test$counts>=1)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  v.test.norm <- voom(dge.test.filt, design.test, normalize.method="scale")
  cpm.test.norm <- v.test.norm$E
  
  merged.miR.ids <- merge(cpm.test.norm, tt.comp.scl.factor, by=0)$Row.names
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  cpm.test.norm.filt.md1 <- cpm.test.norm.filt[, grep(x.ids[1], colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)]
  tt.comp.scl.factor.filt <- tt.comp.scl.factor[merged.miR.ids,]
  
  # set a floor min. Logic here is that if we had a miR that had a count of 1 and the scaling factor 
  # indicates a drop of -5.0 fold, we do not want to make the new log2 count log2(1) - 5.0, since our scale will be way off
  # for plotting. So, if the value is below the min, set it to the min value. In other words. Limit the amount 
  # of downward scaling to the minimum value of the sample.
  #cpm.test.norm.filt.scaled <- apply(cpm.test.norm.filt, 2, FUN=function(x){
  #  x.floor <- min(x[x>0])
  #  x.log.scl <- x-tt.comp.scl.factor.filt$logFC
  #  return(ifelse(x.log.scl<x.floor, x.floor, x.log.scl))
  #})
  
  # this does it without a floor
  # Decided to use median scaling to get the scaling factors
  # Then after scaling, use quantile normalization to match up the distributions
  # Thinking that this might be a more appropriate way to handle the fact that
  # some of the scaling at the low end will lead to values less than the minimum. 
  cpm.test.norm.filt.md.1.scaled <- apply(cpm.test.norm.filt.md1, 2, FUN=function(x){
    #x.floor <- min(x[x>0])
    x.log.scl <- x-tt.comp.scl.factor.filt$logFC
    #return(ifelse(x.log.scl<x.floor, x.floor, x.log.scl))
    return(x.log.scl)
  })
  colnames(cpm.test.norm.filt) <- paste0(colnames(cpm.test.norm.filt), ".orig")
  colnames(cpm.test.norm.filt.md.1.scaled) <- paste0(colnames(cpm.test.norm.filt.md.1.scaled), ".transf")
  cpm.test.orig.and.scaled <- cbind(cpm.test.norm.filt, cpm.test.norm.filt.md.1.scaled)
  cpm.test.orig.and.scaled.rescaled <- normalizeBetweenArrays(cpm.test.orig.and.scaled, method = "quantile")
  ids.reference <- grep(x.ids[1], colnames(cpm.test.orig.and.scaled.rescaled), value=TRUE, ignore.case=TRUE, invert = TRUE)
  ids.test.orig <- grep(x.ids[1], colnames(cpm.test.norm.filt), value=TRUE, ignore.case=TRUE)
  ids.test.scaled <- colnames(cpm.test.norm.filt.md.1.scaled)
  # Return this to get a matrix of transformed and original values
  # return(cpm.test.orig.and.scaled.rescaled)
  # Now get an average of the values we did not transform
  cpm.rescaled.test.orig <- cpm.test.orig.and.scaled.rescaled[, ids.test.orig]
  cpm.rescaled.test.scaled <- cpm.test.orig.and.scaled.rescaled[, ids.test.scaled]
  cpm.rescaled.ref.orig <- cpm.test.orig.and.scaled.rescaled[, ids.reference]
  cpm.mean.reference <- rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.reference])
  
  # Check scale norm too. Didnt look quite as good, so removed.
  cpm.test.orig.and.scaled.rescaled.med <- normalizeBetweenArrays(cpm.test.orig.and.scaled, method = "scale")
  cpm.rescaled.test.orig.med <- cpm.test.orig.and.scaled.rescaled.med[, ids.test.orig]
  cpm.rescaled.test.scaled.med <- cpm.test.orig.and.scaled.rescaled.med[, ids.test.scaled]
  cpm.rescaled.ref.orig.med <- cpm.test.orig.and.scaled.rescaled.med[, ids.reference]
  cpm.mean.reference.med <- rowMeans(cpm.test.orig.and.scaled.rescaled.med[, ids.reference])
  
  
  cpm.test.orig.and.scaled.rescaled.mean <- data.table(data.frame(miRNA.full.ID=row.names(cpm.test.orig.and.scaled.rescaled), 
                                                                  scaling.factor.logFC=tt.comp.scl.factor.filt$logFC,
                                                                  cpm.mean.reference.quantile = cpm.mean.reference,      
                                                                  cpm.mean.test.orig.quantile = rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.test.orig]),
                                                                  cpm.mean.test.scaled.quantile = rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.test.scaled]),
                                                                  abs.dif.orig.quantile = rowMeans(apply(cpm.rescaled.test.orig, 2, FUN=function(x){ abs(x-cpm.mean.reference)})),
                                                                  abs.dif.transf.quantile = rowMeans(apply(cpm.rescaled.test.scaled, 2, FUN=function(x){ abs(x-cpm.mean.reference)}))
  ))
  cpm.test.orig.and.scaled.rescaled.mean[, comparison:=x.id.list]
  return(cpm.test.orig.and.scaled.rescaled.mean)
})); #ggplot(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.LIST[[3]]) + geom_density(aes(color="VS.ORIG", x=abs.dif.orig)) + geom_density(aes(color="VS.TRANSF", x=abs.dif.transf))
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR <- melt(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.EQUIMOLAR, measure.vars=c("abs.dif.orig.quantile", "abs.dif.transf.quantile"), variable.name = "ORIG.or.TRANSF.and.SCL", value.name = "mean.abs.diff")
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR[, SCALE.METHOD:=tstrsplit(ORIG.or.TRANSF.and.SCL, split=".", fixed=TRUE)[4]]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR[, ORIG.or.TRANSF:=sub(paste0(".", SCALE.METHOD), "", ORIG.or.TRANSF.and.SCL)]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR[, `:=`(SCALE.METHOD=ifelse(SCALE.METHOD=="quantile", "QUANTILE", "MEDIAN"), ORIG.or.TRANSF=ifelse(ORIG.or.TRANSF=="abs.dif.orig", "ORIGINAL", "TRANSFORMED"))]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR.quantile <- subset(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR, SCALE.METHOD=="QUANTILE")
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR.quantile[, pool:="SynthEQ"]
# Filter the equimolar results to include only the miRs we looked at in the plasma pool
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR.quantile.filt <- subset(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR.quantile, miRNA.full.ID%in%ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile$miRNA.full.ID)
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile[, pool:="PlasmaPool"]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile.EQandPlasm <- rbind(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile, ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR.quantile.filt)
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile.EQandPlasm[, format.alpha:=ifelse(ORIG.or.TRANSF=="ORIGINAL", 1, 0.5)]

# First just plot the equimolar pool. Include all miRs; not just the ones in the plasma pool
ggplot(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR.quantile, aes(x=mean.abs.diff, color=ORIG.or.TRANSF)) + geom_density(size=1.5) + facet_wrap(~comparison) + theme_bw() + 
  scale_y_continuous(expand=c(0.01, 0.0), breaks=seq(0.0, 10.0, 0.5), limits=c(0,3.5)) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0.0, 9.0), breaks=seq(0,9, 1.0)) + 
  theme(
    axis.text=element_text(size=14, face="bold", color="black"),
    axis.title=element_text(face="bold", size=16),
    plot.title=element_text(face="bold", size=18),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=16, face="bold", color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=14, face="bold"),
    legend.title=element_text(size=15, face="bold"),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(title="SynthEQ",
       y="probability density",
       x="|Difference VS Reference|(log2)")



# Plot the equimolar and plasma overlaid
ggplot(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile.EQandPlasm, aes(x=mean.abs.diff, color=pool)) + geom_density(size=1.5, aes(linetype=ORIG.or.TRANSF)) + facet_wrap(~comparison) + theme_bw() + 
  scale_linetype_manual(values=c("dotted", "solid")) +
  scale_y_continuous(expand=c(0.0, 0.0), breaks=seq(0.0, 10.0, 0.5), limits=c(0,3.5)) + 
  scale_x_continuous(expand=c(0.01,0), limits=c(0.0, 13.5), breaks=seq(0,15, 1.0)) + 
  theme(
    axis.text=element_text(size=14, face="bold", color="black"),
    axis.title=element_text(face="bold", size=16),
    plot.title=element_text(face="bold", size=18),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=16, face="bold", color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=14, face="bold"),
    legend.title=element_text(size=15, face="bold"),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(y="probability density",
       x="|Difference VS Reference|(log2)")

#########################################
########################################
# 12-10-16 SCALING FACTORS AFTER DISCUSSION WITH GALAS
# GET SEPARETE ONES FOR EACH 4N GROUP, BUT ONLY PLOT 4N_ALL VS TRUSEQ
#####################################################
### Now we want to see what the scaling factors look like for the other methods if you start with the equimolar pool. 
### Filtered out the NEB and CleanTag samples early, so need to re-make here... 

ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.CAST.ALL <- data.frame(dcast.data.table(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE, ReferenceID~lab.libMethod.pool.replicate, value.var="multimapAdjustedReadCount", fill=0, fun.aggregate = "max"), row.names=1)

sample.info.all <- data.frame(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE[, .(n.miRs=.N), by=.(lab.libMethod.pool.replicate, lab.libMethod.pool, lab.libMethod, Lab, lib.method.detail, lib.method.simple, lib.method.old2, pool, replicate)], row.names=1)
row.names(sample.info.all) <- make.names(row.names(sample.info.all))

sample.info.sort.all <- sample.info.all[colnames(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.CAST.ALL),]
ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL <- DGEList(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.CAST.ALL)
ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL$samples <- cbind(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL$samples, sample.info.sort.all)
# Now repeat what we did above.
drop.4n.libs <- c("4N_NEXTflex", "4N_Xu")
# Get equimolar pool samples and remove 4N NextFlex and 4N-Xu samples, as well as any library with <100k miRNA reads.
ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR <- ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL[, ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL$samples$pool=="SynthEQ" & (!ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL$samples$lib.method.detail%in%drop.4n.libs) & ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL$samples$lib.size>=100000]
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA <- ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL[, ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL$samples$pool=="PlasmaPool" & (!ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL$samples$lib.method.detail%in%drop.4n.libs) & ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL$samples$lib.size>=100000]
# Now this is where we need to diverge to include the other methods. 
# The thought is that scaling factors are going to convert between two different methods, so all of our filtering, normalization, centering and testing should 
# be done without the potential influence of the other library prep types.

# NEW 12-11-16 -- ADD GROUND TRUTH CALCS AND USE THE GALAS 4N LIBS (4NB LAB2) FOR THIS PURPOSE
# add a column to the eq samples to indicate these
# We'll call this downbiasing 
#ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$GROUND.TRUTH.SAMPLE <- ifelse(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$Lab=="Lab2" & ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.detail=="4N_B", "DOWNBIAS", "NOT_DOWNBIAS_SAMPLE")
# ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$GROUND.TRUTH.SAMPLE <- ifelse(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$Lab=="Lab2" & ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.detail=="4N_B", "DOWNBIAS", "NOT_DOWNBIAS_SAMPLE")

# Based on reviewers, we added more 4N_B libraries, so use all of them here for now
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$GROUND.TRUTH.SAMPLE <- ifelse(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.detail=="4N_B", "DOWNBIAS", "NOT_DOWNBIAS_SAMPLE")
ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$GROUND.TRUTH.SAMPLE <- ifelse(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.detail=="4N_B", "DOWNBIAS", "NOT_DOWNBIAS_SAMPLE")


lib.comparison.groups <- list(TruSeq.V.4N=c("TruSeq", "4N"), 
                              TruSeq.V.NEBNext=c("TruSeq", "NEBNext"),
                              TruSeq.V.CleanTag=c("TruSeq", "CleanTag"),
                              NEBNext.V.CleanTag=c("NEBNext", "CleanTag"),
                              NEBNext.V.4N=c("NEBNext", "4N"),
                              CleanTag.V.4N=c("CleanTag", "4N"),
                              X4N.V.DOWNBIAS=c("4N", "DOWNBIAS")
)
downbias.labLibMethod <- unique(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples[ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$GROUND.TRUTH.SAMPLE=="DOWNBIAS", ]$lab.libMethod)

# Do the first half just to get the scaling factors into a data table format
# UPDATE 1/10 -- explicitly filter and normalize ahead of time, in case it wasn't dont before.
ALL.COMPARISON.SCALING.FACTOR.DT <- do.call(rbind, 
                                            sapply(names(lib.comparison.groups), USE.NAMES=TRUE, simplify = FALSE, FUN=function(x.id.list){
                                              print(x.id.list)
                                              x.ids <- unlist(lib.comparison.groups[x.id.list], recursive = TRUE, use.names = FALSE)
                                              dge.comp <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[, ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.simple%in%x.ids | ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$GROUND.TRUTH.SAMPLE%in%x.ids]
                                              design.comp <- model.matrix(~0+lab.libMethod, dge.comp$samples)
                                              colnames(design.comp) <- ifelse(unique(dge.comp$samples$lab.libMethod)%in%downbias.labLibMethod, "DOWNBIAS", unique(dge.comp$samples$lab.libMethod))
                                              
                                              lib.methods.comp <- strsplit2(strsplit2(colnames(design.comp), split=".", fixed=TRUE)[,1], split="_")[,1]
                                              comp.contrasts.fodder <- sapply(unique(lib.methods.comp), USE.NAMES=FALSE, FUN=function(lib.m){
                                                mat <- design.comp
                                                y.sum <- rowSums(apply(mat, 1, FUN=function(x){
                                                  x*ifelse(lib.methods.comp==lib.m, 1, 0)
                                                }))
                                                y.sum/sum(y.sum)
                                              })
                                              colnames(comp.contrasts.fodder) <- unique(lib.methods.comp)
                                              out.contrast <- as.numeric(comp.contrasts.fodder[, x.ids[1]]-comp.contrasts.fodder[, x.ids[2]])
                                              # Now filter miRs
                                              dge.comp.filt <- dge.comp[rowSums(dge.comp$counts>0)==dim(dge.comp)[2],]
                                              dge.comp.filt <- calcNormFactors(dge.comp.filt, method = "RLE")
                                              # v.comp <- voomWithQualityWeights(dge.comp.filt, design.comp, normalize.method = "scale")
                                              # Drop the quality weighting, just for simplicity
                                              v.comp <- voom(dge.comp.filt, design.comp, normalize.method = "scale")
                                              fit.comp <- lmFit(v.comp, design.comp)
                                              fit.comp <- eBayes(fit.comp, robust = TRUE, trend = TRUE)
                                              fit2.comp <- contrasts.fit(fit.comp, out.contrast)
                                              fit2.comp <- eBayes(fit2.comp, trend = TRUE, robust = TRUE)
                                              tt.comp <- data.frame(topTable(fit2.comp, n = Inf, confint=0.95, sort.by = "none"))
                                              tt.comp.dt <- data.table(data.frame(miR.id.full=row.names(tt.comp), tt.comp))
                                              tt.comp.dt[, comparison:=x.id.list]
                                              return(tt.comp.dt)
                                            })
)

lib.detail.for.comp <- levels(factor(unique(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.detail), levels=lib.detail.levels, exclude = "4N_Xu"))
lib.detail.matr <- matrix(data = 1, nrow = length(lib.detail.for.comp), ncol= length(lib.detail.for.comp))
dimnames(lib.detail.matr) <- list(lib.detail.for.comp, lib.detail.for.comp)
comparisons.detail.tmp <- data.frame(melt(lib.detail.matr)[melt(upper.tri(lib.detail.matr, diag = FALSE))[,3],1:2], stringsAsFactors = FALSE)
ground.truth.comps <- data.frame(Var1=lib.detail.for.comp, Var2="DOWNBIAS")
comparisons.detail <- rbind(comparisons.detail.tmp, ground.truth.comps)
row.names(comparisons.detail) <- NULL
ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS <- do.call(rbind, 
                                                    sapply(1:dim(comparisons.detail)[1], USE.NAMES=FALSE, simplify = FALSE, FUN=function(x.row){
                                                      
                                                      x.ids <- c(as.character(comparisons.detail[x.row,]$Var1), as.character(comparisons.detail[x.row,]$Var2))
                                                      #print(x.ids)
                                                      x.id.list <- paste0(x.ids[1], ".V.", x.ids[2])
                                                      dge.comp <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[, ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.detail%in%x.ids | ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$GROUND.TRUTH.SAMPLE%in%x.ids]
                                                      
                                                      design.comp <- model.matrix(~0+lab.libMethod, dge.comp$samples)
                                                      colnames(design.comp) <- ifelse(any(x.ids=="DOWNBIAS") & unique(dge.comp$samples$lab.libMethod)%in%downbias.labLibMethod, "DOWNBIAS", unique(dge.comp$samples$lab.libMethod))
                                                      if(length(unique(colnames(design.comp)))<2){
                                                        warning(paste0("skipping ", x.ids[1], " ", x.ids[2]))
                                                        return(NULL)
                                                      }
                                                      lib.methods.comp <- strsplit2(colnames(design.comp), split=".", fixed=TRUE)[,1]
                                                      
                                                      comp.contrasts.fodder <- sapply(unique(lib.methods.comp), USE.NAMES=FALSE, FUN=function(lib.m){
                                                        mat <- design.comp
                                                        y.sum <- rowSums(apply(mat, 1, FUN=function(x){
                                                          x*ifelse(lib.methods.comp==lib.m, 1, 0)
                                                        }))
                                                        y.sum/sum(y.sum)
                                                      })
                                                      colnames(comp.contrasts.fodder) <- unique(lib.methods.comp)
                                                      out.contrast <- as.numeric(comp.contrasts.fodder[, x.ids[1]]-comp.contrasts.fodder[, x.ids[2]])
                                                      # Now filter miRs
                                                      dge.comp.filt <- dge.comp[rowSums(dge.comp$counts>0)==dim(dge.comp)[2],]
                                                      # v.comp <- voomWithQualityWeights(dge.comp.filt, design.comp, normalize.method = "scale")
                                                      # Drop the quality weighting, just for simplicity
                                                      dge.comp.filt <- calcNormFactors(dge.comp.filt, method = "RLE")
                                                      v.comp <- voom(dge.comp.filt, design.comp, normalize.method = "scale")
                                                      fit.comp <- lmFit(v.comp, design.comp)
                                                      fit.comp <- eBayes(fit.comp, robust = TRUE, trend = TRUE)
                                                      fit2.comp <- contrasts.fit(fit.comp, out.contrast)
                                                      fit2.comp <- eBayes(fit2.comp, trend = TRUE, robust = TRUE)
                                                      tt.comp <- data.frame(topTable(fit2.comp, n = Inf, confint=0.95, sort.by = "none"))
                                                      tt.comp.dt <- data.table(data.frame(miR.id.full=row.names(tt.comp), tt.comp))
                                                      tt.comp.dt[, comparison:=x.id.list]
                                                      return(tt.comp.dt)
                                                    })
)

# save the data from the original correction factor calculation input so that we can try alternatives
save(file = "20170810_CrossU01_ScalingFactor_Calculation_Input_Data.RData",
     list = Hmisc::Cs(comparisons.detail,
                      ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR,
                      downbias.labLibMethod,
                      ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA,
                      CROSS.U01.PLASM.METADATA,
                      CROSS.U01.SYNTH.METADATA,
                      CROSS.U01.SYNTH.METADATA.WITH.MIR.SENSE
     ))



ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG <- rbind(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS, ALL.COMPARISON.SCALING.FACTOR.DT[grepl("4N", comparison)])
ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG[, c("miR.name", "miRBase.accession", "mature.family"):=tstrsplit(miR.id.full, split=":")[c(1,2,5)]]
ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG[, c("from.method", "to.method"):=tstrsplit(comparison, split=".V.")]
ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset <- subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG, select=c("miR.id.full", "miR.name", "miRBase.accession", "mature.family", "comparison", "from.method", "to.method", "logFC", "CI.L", "CI.R"))
setnames(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, c("CI.L", "CI.R", "comparison"), c("CI.95L", "CI.95R", "correction"))
ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset[, from.method:=sub("^X4N", "4N", from.method)]
ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset[, to.method:=sub("^X4N", "4N", to.method)]
ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset[, correction:=sub("^X4N", "4N", correction)]
write.table(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset[to.method!="DOWNBIAS",], "20170110_INTER-METHOD_Correction_Factors_equimolar_pool_voom_limma.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset[to.method=="DOWNBIAS",], "20170110_DOWNBIAS_Correction_Factors_equimolar_pool_voom_limma.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


# Now remove clean tag & 4N sub-methods from plots to show
plot.methods <- c("TruSeq", "NEBNext", "4N")
ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset.FOR.PLOTTING <- subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, from.method%in%plot.methods & (to.method%in%plot.methods | to.method=="DOWNBIAS"))
lib.comparison.groups.for.plot <- unique(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset.FOR.PLOTTING$correction)

# Apply to plasma data
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.QUANTILENORM <- do.call(rbind,  sapply(lib.comparison.groups.for.plot, USE.NAMES=FALSE, simplify = FALSE, FUN=function(x.group){
  
  correct.dt <- data.frame(subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset.FOR.PLOTTING, correction==x.group), row.names=1)
  x.ids <- strsplit2(x.group, split=".V.")
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple%in%x.ids | ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$GROUND.TRUTH.SAMPLE%in%x.ids]
  keep.mirs <- (rowSums(dge.test$counts>=5)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  v.test.norm <- voom(dge.test.filt, design.test, normalize.method="scale")
  cpm.test.norm <- v.test.norm$E
  
  merged.miR.ids <- merge(cpm.test.norm, correct.dt, by=0)$Row.names
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  cpm.test.norm.filt.md1 <- cpm.test.norm.filt[, grep(x.ids[1], colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)]
  tt.comp.scl.factor.filt <- correct.dt[merged.miR.ids,]
  cpm.test.norm.filt.md.1.scaled <-
    apply(
      cpm.test.norm.filt.md1,
      2,
      FUN = function(x) {
        x.log.scl <- x - tt.comp.scl.factor.filt$logFC
      })
  colnames(cpm.test.norm.filt) <-
    paste0(colnames(cpm.test.norm.filt), ".orig")
  colnames(cpm.test.norm.filt.md.1.scaled) <-
    paste0(colnames(cpm.test.norm.filt.md.1.scaled), ".transf")
  cpm.test.orig.and.scaled <-
    cbind(cpm.test.norm.filt, cpm.test.norm.filt.md.1.scaled)
  col.names.combined <- colnames(cpm.test.orig.and.scaled)
  group.des.df <- ifelse(grepl("transf$", col.names.combined), paste0(x.ids[1], ".to.", x.ids[2]), ifelse(grepl(x.ids[1], col.names.combined), paste0(x.ids[1], ".orig"), paste0(x.ids[2], ".orig")))
  group.design <- model.matrix(~0+group.des.df)
  colnames(group.design) <- sub("group.des.df", "", colnames(group.design))
  cpm.test.orig.and.scaled.rescaled <-
    voom(2^cpm.test.orig.and.scaled, group.design, normalize.method=ifelse(x.ids[2]=="DOWNBIAS", "scale", "quantile"))$E
  ids.reference <- if(x.ids[2]=="DOWNBIAS"){
    paste0(row.names(dge.test.filt$samples[dge.test.filt$samples$GROUND.TRUTH.SAMPLE=="DOWNBIAS",]), ".orig")
  } else{
    grep(
      x.ids[1],
      colnames(cpm.test.orig.and.scaled.rescaled),
      value = TRUE,
      ignore.case = TRUE,
      invert = TRUE
    )
  }
  
  ids.test.orig <-
    grep(x.ids[1],
         colnames(cpm.test.norm.filt),
         value = TRUE,
         ignore.case = TRUE)
  ids.test.scaled <- colnames(cpm.test.norm.filt.md.1.scaled)
  # Return this to get a matrix of transformed and original values
  # return(cpm.test.orig.and.scaled.rescaled)
  # Now get an average of the values only did intra-lab scaling
  cpm.rescaled.test.orig <-
    cpm.test.orig.and.scaled.rescaled[, ids.test.orig]
  cpm.rescaled.test.scaled <-
    cpm.test.orig.and.scaled.rescaled[, ids.test.scaled]
  cpm.rescaled.ref.orig <-
    cpm.test.orig.and.scaled.rescaled[, ids.reference]
  cpm.mean.reference <-
    rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.reference])
  
  # Check scale norm too. Didnt look quite as good, so removed.
  cpm.test.orig.and.scaled.rescaled.med <-
    voom(2^cpm.test.orig.and.scaled, design = group.design, normalize.method = "scale")$E
  cpm.rescaled.test.orig.med <-
    cpm.test.orig.and.scaled.rescaled.med[, ids.test.orig]
  cpm.rescaled.test.scaled.med <-
    cpm.test.orig.and.scaled.rescaled.med[, ids.test.scaled]
  cpm.rescaled.ref.orig.med <-
    cpm.test.orig.and.scaled.rescaled.med[, ids.reference] 
  cpm.mean.reference.med <-
    rowMeans(cpm.test.orig.and.scaled.rescaled.med[, ids.reference])
  
  
  cpm.test.orig.and.scaled.rescaled.mean <-
    data.table(
      data.frame(
        miRNA.full.ID = row.names(cpm.test.orig.and.scaled.rescaled),
        scaling.factor.logFC =
          tt.comp.scl.factor.filt$logFC,
        cpm.mean.reference = cpm.mean.reference,
        cpm.mean.test.orig = rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.test.orig]),
        cpm.mean.test.scaled = rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.test.scaled]),
        abs.dif.orig = rowMeans(apply(
          cpm.rescaled.test.orig,
          2,
          FUN = function(x) {
            abs(x - cpm.mean.reference)
          }
        )),
        abs.dif.transf = rowMeans(apply(
          cpm.rescaled.test.scaled,
          2,
          FUN = function(x) {
            abs(x - cpm.mean.reference)
          }
        ))
      )
    )
  cpm.test.orig.and.scaled.rescaled.mean[, comparison := x.group]
  return(cpm.test.orig.and.scaled.rescaled.mean)
  #return(cpm.test.orig.and.scaled.rescaled)
}))
#ggplot(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.LIST[[3]]) + geom_density(aes(color="VS.ORIG", x=abs.dif.orig)) + geom_density(aes(color="VS.TRANSF", x=abs.dif.transf))


ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif <- melt(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.QUANTILENORM, measure.vars=c("abs.dif.orig", "abs.dif.transf"), variable.name = "ORIG.or.TRANSF.and.SCL", value.name = "mean.abs.diff")
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif[, SCALE.METHOD:="VOOM.NO.INTERLAB"]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif[, ORIG.or.TRANSF:=sub(paste0(".", SCALE.METHOD), "", ORIG.or.TRANSF.and.SCL), by=1:nrow(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif)]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif[, `:=`(SCALE.METHOD=toupper(SCALE.METHOD), ORIG.or.TRANSF=ifelse(ORIG.or.TRANSF=="abs.dif.orig", "ORIGINAL", "TRANSFORMED"))]
#ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile <- subset(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif, SCALE.METHOD=="QUANTILE")
ggplot(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif, aes(x=mean.abs.diff, color=ORIG.or.TRANSF)) + geom_density(size=1.5) + facet_wrap(~comparison) + theme_bw() + 
  scale_y_continuous(expand=c(0.0, 0.0), breaks=seq(0.0, 1.5, 0.1), limits=c(0.0, 1.4)) + 
  scale_x_continuous(expand=c(0.01,0), limits=c(-0.1, 8.0), breaks=seq(0,9, 1.0)) + 
  theme(
    axis.text=element_text(size=14, face="bold", color="black"),
    axis.title=element_text(face="bold", size=16),
    plot.title=element_text(face="bold", size=18),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=16, face="bold", color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=14, face="bold"),
    legend.title=element_text(size=15, face="bold"),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(y="probability density",
       x="|Difference VS Reference|(log2)")





# Get the matrices for plotting
ALL.COMPARISON.SCALING.FACTOR.CORRECTED.MATR.LIST.UNRESCALED <- sapply(lib.comparison.groups.for.plot, USE.NAMES=FALSE, simplify = FALSE, FUN=function(x.group){
  
  correct.dt <- data.frame(subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset.FOR.PLOTTING, correction==x.group), row.names=1)
  x.ids <- strsplit2(x.group, split=".V.")
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple%in%x.ids | ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$GROUND.TRUTH.SAMPLE%in%x.ids]
  keep.mirs <- (rowSums(dge.test$counts>=1)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  v.test.norm <- voom(dge.test.filt, design.test, normalize.method="scale")
  cpm.test.norm <- v.test.norm$E
  
  merged.miR.ids <- merge(cpm.test.norm, correct.dt, by=0)$Row.names
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  cpm.test.norm.filt.md1 <- cpm.test.norm.filt[, grep(x.ids[1], colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)]
  tt.comp.scl.factor.filt <- correct.dt[merged.miR.ids,]
  cpm.test.norm.filt.md.1.scaled <-
    apply(
      cpm.test.norm.filt.md1,
      2,
      FUN = function(x) {
        x.log.scl <- x - tt.comp.scl.factor.filt$logFC
      })
  colnames(cpm.test.norm.filt) <-
    paste0(colnames(cpm.test.norm.filt), ".orig")
  colnames(cpm.test.norm.filt.md.1.scaled) <-
    paste0(colnames(cpm.test.norm.filt.md.1.scaled), ".transf")
  cpm.test.orig.and.scaled <-
    cbind(cpm.test.norm.filt, cpm.test.norm.filt.md.1.scaled)
})
names(ALL.COMPARISON.SCALING.FACTOR.CORRECTED.MATR.LIST.UNRESCALED) <- lib.comparison.groups.for.plot
### 12-10-2016 Only include truseq to 4N for initial plots
TruSeq.V.4N.Scaled.And.Orig <- ALL.COMPARISON.SCALING.FACTOR.CORRECTED.MATR.LIST.UNRESCALED[["TruSeq.V.4N"]]

# Now get the annotation information
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO <- data.frame(strsplit2(colnames(TruSeq.V.4N.Scaled.And.Orig), split=".", fixed=TRUE), stringsAsFactors = FALSE, row.names=colnames(TruSeq.V.4N.Scaled.And.Orig))
colnames(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO) <- c("lib.method.detail", "Lab", "pool", "replicate", "orig.or.transf")
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO$lab.libMethod.pool.replicate <- apply(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO[, c("lib.method.detail", "Lab", "pool", "replicate")], 1, paste, collapse=".")
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO$lab.libMethod.pool.replicate.transf <- row.names(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO)

TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG <- data.table(merge(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.FILT.SCALED.MIRS$samples, by.x=which(colnames(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO)=="lab.libMethod.pool.replicate"), by.y=0, all.x = TRUE))
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG[, is.corrected:=ifelse(orig.or.transf=="orig", "ORIGINAL", "CORRECTED")]
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG[, transf.LibMethod.id:=sub("X4N", "4N", sub(".PlasmaPool", "", make.names(paste(lab.libMethod.pool.replicate, substr(is.corrected, 1, 4), sep="."), unique = TRUE)))]
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG[, `:=`(from.method=lib.method.simple, to.method=ifelse(is.corrected=="ORIGINAL", as.character(lib.method.simple), "4N"))]

TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df <- data.frame(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG, row.names=which(colnames(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG)=="lab.libMethod.pool.replicate.transf"))
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df <- TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df[colnames(TruSeq.V.4N.Scaled.And.Orig),]                                                                                     
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot <- TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df[, c("to.method", "from.method", "transf.LibMethod.id")]
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot$from.method <- factor(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot$from.method, lib.simple.levels)
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot$to.method <- factor(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot$to.method, lib.simple.levels)

groups <- apply(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot[,1:2], 1, FUN=function(x) ifelse(x[1]==x[2], paste0(x[1], ".orig"), paste0(x[1], ".to.", x[2])))
design.scl <- model.matrix(~0+groups)
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE <- voom(2^TruSeq.V.4N.Scaled.And.Orig, design.scl, normalize.method="quantile", plot=TRUE)$E
library(scales)
ann.colors.plasm = list(
  from.method = hue_pal()(4),
  to.method = hue_pal()(4),
  is.corrected = grey_pal()(2)
)
names(ann.colors.plasm$from.method) <- lib.simple.levels
names(ann.colors.plasm$to.method) <- lib.simple.levels
names(ann.colors.plasm$is.corrected) <- c("TRANSFORMED", "ORIGINAL")
row.mean.abund.df <- data.frame(FUNCTION.aveTechReps(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE, groups))
# Full heatmap with row  clustering
pheatmap(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE,
         show_rownames = FALSE,
         fontsize = 14,
         fontsize_col = 12,
         labels_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot$transf.LibMethod.id,
         annotation_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot[, 1:2], 
         annotation_row = row.mean.abund.df,
         annotation_colors=ann.colors.plasm 
)
# Full heatmap with truseq sorting


TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.sort.truseq.abund <- TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE[order(-row.mean.abund.df$TruSeq.orig),]
pheatmap(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.sort.truseq.abund,
         cluster_rows = FALSE,
         show_rownames = FALSE,
         fontsize = 10,
         fontsize_col = 10,
         treeheight_col = 15,
         labels_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot$transf.LibMethod.id,
         annotation_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot[, 1:2], 
         annotation_row = row.mean.abund.df,
         annotation_colors=ann.colors.plasm, 
         filename = "FIG_6B_TruSeq_4N_Correction_Heatmap_PLASMAPOOL.pdf"
)

# Muneesh wanted to see it if we filter mIRs with < 5 counts
these.sample.ids <- unique(sub(".orig|.transf", "", row.names(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot)))
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.this.subset <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, these.sample.ids]
mir.id.filt.gt.5 <- row.names(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.this.subset[rowSums(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.this.subset$counts>=5)==dim(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.this.subset)[2],])

TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.filt.5count <- TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE[row.names(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE)%in%mir.id.filt.gt.5,]
row.mean.abund.df.filt.5count <- row.mean.abund.df[row.names(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.filt.5count),]

TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.sort.truseq.abund.5count <- TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.filt.5count[order(-row.mean.abund.df.filt.5count$TruSeq.orig),]
pheatmap(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.sort.truseq.abund.5count,
         cluster_rows = FALSE,
         show_rownames = FALSE,
         fontsize = 14,
         fontsize_col = 12,
         labels_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot$transf.LibMethod.id,
         annotation_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot[, 1:2], 
         annotation_row = row.mean.abund.df.filt.5count,
         annotation_colors=ann.colors.plasm 
)


#####################################
### look at the downsampling biased corrections
#Density plot of diff to ref 
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif[, comparison:=sub(".V.DOWNBIAS", ".DOWNBIASED", comparison)]
ggplot(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif[comparison=="TruSeq.DOWNBIASED"], aes(x=mean.abs.diff, color=ORIG.or.TRANSF)) + geom_density(size=1.5) + facet_wrap(~comparison) + theme_bw() + 
  scale_y_continuous(expand=c(0.0, 0.0), breaks=seq(0.0, 1.0, 0.1), limits=c(0.0, 1.0)) + 
  scale_x_continuous(expand=c(0.01,0), limits=c(-0.1, 8.0), breaks=seq(0,9, 1.0)) + 
  theme(
    axis.text=element_text(size=14, face="bold", color="black"),
    axis.title=element_text(face="bold", size=16),
    plot.title=element_text(face="bold", size=18),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=16, face="bold", color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=14, face="bold"),
    legend.title=element_text(size=15, face="bold"),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(y="probability density",
       x="|Difference VS Reference|(log2)")


## And for the matrices
TruSeq.DOWNBIAS.Scaled.And.Orig <- ALL.COMPARISON.SCALING.FACTOR.CORRECTED.MATR.LIST.UNRESCALED[["TruSeq.V.DOWNBIAS"]]

TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO <- data.frame(strsplit2(colnames(TruSeq.DOWNBIAS.Scaled.And.Orig), split=".", fixed=TRUE), stringsAsFactors = FALSE, row.names=colnames(TruSeq.DOWNBIAS.Scaled.And.Orig))
colnames(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO) <- c("lib.method.detail", "Lab", "pool", "replicate", "orig.or.transf")
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO$lab.libMethod.pool.replicate <- apply(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO[, c("lib.method.detail", "Lab", "pool", "replicate")], 1, paste, collapse=".")
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO$lab.libMethod.pool.replicate.transf <- row.names(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO)

TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG <- data.table(merge(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.FILT.SCALED.MIRS$samples, by.x=which(colnames(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO)=="lab.libMethod.pool.replicate"), by.y=0, all.x = TRUE))
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG[, is.corrected:=ifelse(orig.or.transf=="orig", "ORIGINAL", "CORRECTED")]
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG[, transf.LibMethod.id:=sub("X4N", "4N", sub(".PlasmaPool", "", make.names(paste(lab.libMethod.pool.replicate, substr(is.corrected, 1, 4), sep="."), unique = TRUE)))]
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG[, lib.method.detail:=sub("^X4N", "4N", lib.method.detail.x)]
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG[, c("lib.method.detail.x", "lib.method.detail.y"):=NULL]
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG[, `:=`(from.method=ifelse(lab.libMethod%in%downbias.labLibMethod, "Downbias.ref", lib.method.detail), to.method=ifelse(is.corrected=="CORRECTED" | lab.libMethod%in%downbias.labLibMethod, "Downbias.ref", lib.method.detail)), by=1:nrow(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG)]

TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df <- data.frame(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG, row.names=which(colnames(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG)=="lab.libMethod.pool.replicate.transf"))
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df <- TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df[colnames(TruSeq.DOWNBIAS.Scaled.And.Orig),]                                                                                     
lib.simple.and.detail.levels.with.downbias <- c(lib.simple.levels, "Downbias.ref", lib.detail.levels[!lib.detail.levels%in%lib.simple.levels])
lib.levels.use <- lib.simple.and.detail.levels.with.downbias[lib.simple.and.detail.levels.with.downbias%in%unique(c(as.character(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df$from.method), as.character(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df$to.method)))]
lib.levels.use.pos <- match(lib.levels.use, lib.simple.and.detail.levels.with.downbias)

TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df$from.method <- factor(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df$from.method, lib.levels.use)
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df$to.method <- factor(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df$to.method, lib.levels.use)
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df.for.plot <- TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df[, c("to.method", "from.method", "transf.LibMethod.id")]
groups <- apply(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df[,c("from.method", "to.method")], 1, FUN=function(x){
  g.tmp <- ifelse(x[1]==x[2], paste0(x[1], ".orig"), paste0(x[1], ".to.", x[2]))
  sub(".to.Downbias.ref", ".DOWNBIASED", g.tmp)
} )
design.scl <- model.matrix(~0+groups)
TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE <- voom(2^TruSeq.DOWNBIAS.Scaled.And.Orig, design.scl, normalize.method="quantile", plot=TRUE)$E

# Custom palatte.. Does Truseq, CleanTag, NEB and 4N as normal. Then adds a greenish color for DOWNBIAS, and then takes the hue_pal_4 value and shades them for the different 4N methods
custom.lib.method.pallate <- c(hue_pal()(4), "#3FA1FF", hcl(285, c = 100, l = 65, alpha = 0.9*seq(4,1)/4))

ann.colors.plasm = list(
  from.method = custom.lib.method.pallate[lib.levels.use.pos],
  to.method = custom.lib.method.pallate[lib.levels.use.pos],
  is.corrected = grey_pal()(2)
)
names(ann.colors.plasm$from.method) <- lib.simple.and.detail.levels.with.downbias[lib.levels.use.pos]
names(ann.colors.plasm$to.method) <- lib.simple.and.detail.levels.with.downbias[lib.levels.use.pos]
names(ann.colors.plasm$is.corrected) <- c("TRANSFORMED", "ORIGINAL")
row.mean.abund.df <- data.frame(FUNCTION.aveTechReps(TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE, groups))
# Full heatmap with row  clustering
pheatmap(TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE,
         show_rownames = FALSE,
         fontsize = 14,
         fontsize_col = 12,
         labels_col = TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df.for.plot$transf.LibMethod.id,
         annotation_col = TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df.for.plot[, 1:2], 
         annotation_row = row.mean.abund.df,
         annotation_colors=ann.colors.plasm 
)

# Sort heatmap
row.mean.abund.df.sort <- row.mean.abund.df[order(-row.mean.abund.df$TruSeq.orig),]
# Full heatmap with row  clustering
pheatmap(TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE[row.names(row.mean.abund.df.sort),],
         show_rownames = FALSE,
         cluster_rows = FALSE,
         fontsize = 14,
         fontsize_col = 12,
         labels_col = TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df.for.plot$transf.LibMethod.id,
         annotation_col = TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df.for.plot[, 1:2], 
         annotation_row = row.mean.abund.df.sort,
         annotation_colors=ann.colors.plasm 
)



############################################
#### REPEAT FOR EQUIMOLAR POOL 12-11-16


# Now remove clean tag & 4N sub-methods from plots to show
plot.methods <- c("TruSeq", "NEBNext", "4N")
ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset.FOR.PLOTTING <- subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, from.method%in%plot.methods & (to.method%in%plot.methods | to.method=="DOWNBIAS"))
lib.comparison.groups.for.plot <- unique(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset.FOR.PLOTTING$correction)

# Apply to equimolar data

ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.QUANTILENORM.EQ <- do.call(rbind,  sapply(lib.comparison.groups.for.plot, USE.NAMES=FALSE, simplify = FALSE, FUN=function(x.group){
  correct.dt <- data.frame(subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset.FOR.PLOTTING, correction==x.group), row.names=1)
  x.ids <- strsplit2(x.group, split=".V.")
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[, (ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.simple%in%x.ids | ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$GROUND.TRUTH.SAMPLE%in%x.ids)]
  keep.mirs <- (rowSums(dge.test$counts>=1)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  v.test.norm <- voom(dge.test.filt, design.test, normalize.method="scale")
  cpm.test.norm <- v.test.norm$E
  
  merged.miR.ids <- merge(cpm.test.norm, correct.dt, by=0)$Row.names
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  cpm.test.norm.filt.md1 <- cpm.test.norm.filt[, grep(x.ids[1], colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)]
  tt.comp.scl.factor.filt <- correct.dt[merged.miR.ids,]
  cpm.test.norm.filt.md.1.scaled <-
    apply(
      cpm.test.norm.filt.md1,
      2,
      FUN = function(x) {
        x.log.scl <- x - tt.comp.scl.factor.filt$logFC
      })
  colnames(cpm.test.norm.filt) <-
    paste0(colnames(cpm.test.norm.filt), ".orig")
  colnames(cpm.test.norm.filt.md.1.scaled) <-
    paste0(colnames(cpm.test.norm.filt.md.1.scaled), ".transf")
  cpm.test.orig.and.scaled <-
    cbind(cpm.test.norm.filt, cpm.test.norm.filt.md.1.scaled)
  col.names.combined <- colnames(cpm.test.orig.and.scaled)
  group.des.df <- ifelse(grepl("transf$", col.names.combined), paste0(x.ids[1], ".to.", x.ids[2]), ifelse(grepl(x.ids[1], col.names.combined), paste0(x.ids[1], ".orig"), paste0(x.ids[2], ".orig")))
  group.design <- model.matrix(~0+group.des.df)
  colnames(group.design) <- sub("group.des.df", "", colnames(group.design))
  cpm.test.orig.and.scaled.rescaled <-
    voom(2^cpm.test.orig.and.scaled, group.design, normalize.method=ifelse(x.ids[2]=="DOWNBIAS", "scale", "quantile"))$E
  ids.reference <- if(x.ids[2]=="DOWNBIAS"){
    grep(
      downbias.labLibMethod,
      colnames(cpm.test.orig.and.scaled.rescaled),
      value=TRUE,
      ignore.case=TRUE
    )
  } else{
    grep(
      x.ids[1],
      colnames(cpm.test.orig.and.scaled.rescaled),
      value = TRUE,
      ignore.case = TRUE,
      invert = TRUE
    )
  }
  
  ids.test.orig <-
    grep(x.ids[1],
         colnames(cpm.test.norm.filt),
         value = TRUE,
         ignore.case = TRUE)
  ids.test.scaled <- colnames(cpm.test.norm.filt.md.1.scaled)
  # Return this to get a matrix of transformed and original values
  # return(cpm.test.orig.and.scaled.rescaled)
  # Now get an average of the values only did intra-lab scaling
  cpm.rescaled.test.orig <-
    cpm.test.orig.and.scaled.rescaled[, ids.test.orig]
  cpm.rescaled.test.scaled <-
    cpm.test.orig.and.scaled.rescaled[, ids.test.scaled]
  cpm.rescaled.ref.orig <-
    cpm.test.orig.and.scaled.rescaled[, ids.reference]
  cpm.mean.reference <-
    rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.reference])
  
  # Check scale norm too. Didnt look quite as good, so removed.
  cpm.test.orig.and.scaled.rescaled.med <-
    voom(2^cpm.test.orig.and.scaled, design = group.design, normalize.method = "scale")$E
  cpm.rescaled.test.orig.med <-
    cpm.test.orig.and.scaled.rescaled.med[, ids.test.orig]
  cpm.rescaled.test.scaled.med <-
    cpm.test.orig.and.scaled.rescaled.med[, ids.test.scaled]
  cpm.rescaled.ref.orig.med <-
    cpm.test.orig.and.scaled.rescaled.med[, ids.reference] 
  cpm.mean.reference.med <-
    rowMeans(cpm.test.orig.and.scaled.rescaled.med[, ids.reference])
  
  
  cpm.test.orig.and.scaled.rescaled.mean <-
    data.table(
      data.frame(
        miRNA.full.ID = row.names(cpm.test.orig.and.scaled.rescaled),
        scaling.factor.logFC =
          tt.comp.scl.factor.filt$logFC,
        cpm.mean.reference = cpm.mean.reference,
        cpm.mean.test.orig = rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.test.orig]),
        cpm.mean.test.scaled = rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.test.scaled]),
        abs.dif.orig = rowMeans(apply(
          cpm.rescaled.test.orig,
          2,
          FUN = function(x) {
            abs(x - cpm.mean.reference)
          }
        )),
        abs.dif.transf = rowMeans(apply(
          cpm.rescaled.test.scaled,
          2,
          FUN = function(x) {
            abs(x - cpm.mean.reference)
          }
        ))
      )
    )
  cpm.test.orig.and.scaled.rescaled.mean[, comparison := x.group]
  return(cpm.test.orig.and.scaled.rescaled.mean)
  #return(cpm.test.orig.and.scaled.rescaled)
}))

ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQ <- melt(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.QUANTILENORM.EQ, measure.vars=c("abs.dif.orig", "abs.dif.transf"), variable.name = "ORIG.or.TRANSF.and.SCL", value.name = "mean.abs.diff")
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQ[, SCALE.METHOD:="QUANTILE"]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQ[, ORIG.or.TRANSF:=sub(paste0(".", SCALE.METHOD), "", ORIG.or.TRANSF.and.SCL), by=1:nrow(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQ)]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQ[, `:=`(SCALE.METHOD=toupper(SCALE.METHOD), ORIG.or.TRANSF=ifelse(ORIG.or.TRANSF=="abs.dif.orig", "ORIGINAL", "TRANSFORMED"))]
#ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile <- subset(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif, SCALE.METHOD=="QUANTILE")
ggplot(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQ[comparison=="TruSeq.V.4N"], aes(x=mean.abs.diff, color=ORIG.or.TRANSF)) + geom_density(size=1.5) + facet_wrap(~comparison) + theme_bw() + 
  scale_y_continuous(expand=c(0.0, 0.0), breaks=seq(0.0, 2.0, 0.2), limits=c(0.0, 1.6)) + 
  scale_x_continuous(expand=c(0.01,0), limits=c(-0.1, 8.0), breaks=seq(0,9, 1.0)) + 
  theme(
    axis.text=element_text(size=10, color="black"),
    axis.title=element_text(size=10),
    plot.title=element_text(size=10),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=10, color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=10),
    legend.title=element_text(size=10),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(y="probability density",
       x="|Difference VS Reference|(log2)")



# PLOT THE EQUIMOLAR DOWNBIASED TOO
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQ[, comparison:=sub(".V.DOWNBIAS", ".DOWNBIASED", comparison)]
ggplot(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQ[comparison=="TruSeq.DOWNBIASED"], aes(x=mean.abs.diff, color=ORIG.or.TRANSF)) + geom_density(size=1.5) + facet_wrap(~comparison) + theme_bw() + 
  scale_y_continuous(expand=c(0.0, 0.0), breaks=seq(0.0, 3.0, 0.2), limits=c(0,3)) + 
  scale_x_continuous(expand=c(0.01,0), limits=c(-0.1, 9.5), breaks=seq(0,9, 1.0)) + 
  theme(
    axis.text=element_text(size=14, face="bold", color="black"),
    axis.title=element_text(face="bold", size=16),
    plot.title=element_text(face="bold", size=18),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=16, face="bold", color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=14, face="bold"),
    legend.title=element_text(size=15, face="bold"),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(y="probability density",
       x="|Difference VS Reference|(log2)")




# Get the matrices for plotting EQUIMOLAR
ALL.COMPARISON.SCALING.FACTOR.CORRECTED.MATR.LIST.UNRESCALED.EQ <- sapply(lib.comparison.groups.for.plot, USE.NAMES=FALSE, simplify = FALSE, FUN=function(x.group){
  
  correct.dt <- data.frame(subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset.FOR.PLOTTING, correction==x.group), row.names=1)
  x.ids <- strsplit2(x.group, split=".V.")
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[, ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.simple%in%x.ids | ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$GROUND.TRUTH.SAMPLE%in%x.ids]
  keep.mirs <- (rowSums(dge.test$counts>=1)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  v.test.norm <- voom(dge.test.filt, design.test, normalize.method="scale")
  cpm.test.norm <- v.test.norm$E
  
  merged.miR.ids <- merge(cpm.test.norm, correct.dt, by=0)$Row.names
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  cpm.test.norm.filt.md1 <- cpm.test.norm.filt[, grep(x.ids[1], colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)]
  tt.comp.scl.factor.filt <- correct.dt[merged.miR.ids,]
  cpm.test.norm.filt.md.1.scaled <-
    apply(
      cpm.test.norm.filt.md1,
      2,
      FUN = function(x) {
        x.log.scl <- x - tt.comp.scl.factor.filt$logFC
      })
  colnames(cpm.test.norm.filt) <-
    paste0(colnames(cpm.test.norm.filt), ".orig")
  colnames(cpm.test.norm.filt.md.1.scaled) <-
    paste0(colnames(cpm.test.norm.filt.md.1.scaled), ".transf")
  cpm.test.orig.and.scaled <-
    cbind(cpm.test.norm.filt, cpm.test.norm.filt.md.1.scaled)
})
names(ALL.COMPARISON.SCALING.FACTOR.CORRECTED.MATR.LIST.UNRESCALED.EQ) <- lib.comparison.groups.for.plot
### 12-10-2016 Only include truseq to 4N for initial plots
TruSeq.V.4N.Scaled.And.Orig.eq <- ALL.COMPARISON.SCALING.FACTOR.CORRECTED.MATR.LIST.UNRESCALED.EQ[["TruSeq.V.4N"]]

# Now get the annotation information
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ <- data.frame(strsplit2(colnames(TruSeq.V.4N.Scaled.And.Orig.eq), split=".", fixed=TRUE), stringsAsFactors = FALSE, row.names=colnames(TruSeq.V.4N.Scaled.And.Orig.eq))
colnames(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ) <- c("lib.method.detail", "Lab", "pool", "replicate", "orig.or.transf")
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ$lab.libMethod.pool.replicate <- apply(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ[, c("lib.method.detail", "Lab", "pool", "replicate")], 1, paste, collapse=".")
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ$lab.libMethod.pool.replicate.transf <- row.names(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ)

TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG <- data.table(merge(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ, ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples, by.x=which(colnames(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ)=="lab.libMethod.pool.replicate"), by.y=0, all.x = TRUE))
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG[, is.corrected:=ifelse(orig.or.transf=="orig", "ORIGINAL", "CORRECTED")]
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG[, transf.LibMethod.id:=sub("X4N", "4N", sub(".SynthEQ", "", make.names(paste(lab.libMethod.pool.replicate, substr(is.corrected, 1, 4), sep="."), unique = TRUE)))]
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG[, `:=`(from.method=lib.method.simple, to.method=ifelse(is.corrected=="ORIGINAL", lib.method.simple, "4N"))]

TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df <- data.frame(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG, row.names=which(colnames(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG)=="lab.libMethod.pool.replicate.transf"))
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df <- TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df[colnames(TruSeq.V.4N.Scaled.And.Orig.eq),]                                                                                     
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df.for.plot <- TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df[, c("to.method", "from.method", "transf.LibMethod.id")]
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df.for.plot$from.method <- factor(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df.for.plot$from.method, lib.simple.levels)
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df.for.plot$to.method <- factor(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df.for.plot$to.method, lib.simple.levels)
ann.colors.plasm = list(
  from.method = hue_pal()(4),
  to.method = hue_pal()(4),
  is.corrected = grey_pal()(2)
)
names(ann.colors.plasm$from.method) <- lib.simple.levels
names(ann.colors.plasm$to.method) <- lib.simple.levels
names(ann.colors.plasm$is.corrected) <- c("TRANSFORMED", "ORIGINAL")

#now rescale
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df <- TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df[colnames(TruSeq.V.4N.Scaled.And.Orig.eq),]
groups.eq <- apply(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df[, c("to.method", "from.method")], 1, FUN=function(x) ifelse(x[1]==x[2], paste0(x[1], ".orig"), paste0(x[1], ".to.", x[2])))
design.scl <- model.matrix(~0+groups.eq)
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ <- voom(2^TruSeq.V.4N.Scaled.And.Orig.eq, design.scl, plot=TRUE)$E
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df.for.plot <- TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df.for.plot[colnames(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ),]


row.mean.abund.df.eq <- data.frame(FUNCTION.aveTechReps(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ, groups.eq))
# Full heatmap with row  clustering
pheatmap(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ,
         show_rownames = FALSE,
         fontsize = 14,
         fontsize_col = 12,
         labels_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df.for.plot$transf.LibMethod.id,
         annotation_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df.for.plot[, 1:2], 
         annotation_row = row.mean.abund.df.eq,
         annotation_colors=ann.colors.plasm 
)
# Full heatmap with truseq sorting


TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ.sort.truseq.abund <- TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ[order(-row.mean.abund.df.eq$TruSeq.orig),]
pheatmap(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ.sort.truseq.abund,
         cluster_rows = FALSE,
         show_rownames = FALSE,
         fontsize = 10,
         fontsize_col = 10,
         treeheight_col = 15,
         labels_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df.for.plot$transf.LibMethod.id,
         annotation_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df.for.plot[, 1:2], 
         annotation_row = row.mean.abund.df.eq,
         annotation_colors=ann.colors.plasm#,
         #filename = "FIG_6A_TruSeq_4N_Correction_Heatmap_EquimolarPool.pdf"
)


# 1-19-17  Melt to see if we can see the reduction in bias
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ.melt <- data.table(melt(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ.sort.truseq.abund))
setnames(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ.melt, c("miR.id.full", "sample.ID", "CPM"))
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ.melt[, c("lib.method.detail", "Lab", "Pool", "replicate", "correction"):=tstrsplit(sample.ID, split=".", fixed=TRUE)]
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ.melt[, lib.method.detail:=sub("^X", "", lib.method.detail)]
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ.melt[, lib.method.simple:=tstrsplit(lib.method.detail, split="_")[1]]
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ.melt[, `:=`(lib.method.detail=factor(lib.method.detail, lib.detail.levels), lib.method.simple=factor(lib.method.simple, lib.simple.levels))]
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ.melt.summary <- TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ.melt[, .(ave.cpm=mean(CPM)), by=.(miR.id.full, lib.method.detail, lib.method.simple, Lab, Pool, correction)]
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ.melt.summary[, lab.libMethod:=paste(lib.method.detail, Lab, sep="."), by=1:nrow(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ.melt.summary)]
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ.melt.summary[, lab.libMethod.correction:=paste(lib.method.detail, Lab, correction, sep="."), by=1:nrow(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ.melt.summary)]
lab.lib.levels <- unique(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ.melt.summary$lab.libMethod.correction)
lab.lib.levels.ordered <- c(grep("orig", grep("TruSeq" , lab.lib.levels, value = TRUE), value=TRUE),
                            grep("transf", grep("TruSeq" , lab.lib.levels, value = TRUE), value=TRUE),
                            grep("4N" , lab.lib.levels, value = TRUE))
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ.melt.summary[, lab.libMethod.correction:=factor(lab.libMethod.correction, levels=lab.lib.levels.ordered)]
TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ.melt.summary[, correction:=factor(correction, levels=c("orig", "transf"))]
ggplot(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.EQ.melt.summary, aes(x=lab.libMethod.correction, alpha=correction, fill=lib.method.simple, y=ave.cpm)) + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + scale_alpha_discrete(range=c(1,0.5)) + theme_bw() + theme(
  axis.text=element_text(size=10, color="black"),
  axis.text.x=element_text(hjust=1, angle=50),
  axis.title=element_text(size=10),
  plot.title=element_text(size=10),
  panel.grid.major = element_line(color=NA),
  panel.border = element_rect(size=1.5, color="black"),
  axis.ticks = element_line(color="black", size=1.5),
  strip.text = element_text(size=16, face="bold", color="black"),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.text=element_text(size=14, face="bold"),
  legend.title=element_text(size=15, face="bold"),
  legend.key.size=unit(1.2, "lines")
) + labs(y="Mean CPM (log2)", x=NULL)


# Muneesh wanted to see it if we filter mIRs with < 5 counts
these.sample.ids <- unique(sub(".orig|.transf", "", row.names(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot)))
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.this.subset <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, these.sample.ids]
mir.id.filt.gt.5 <- row.names(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.this.subset[rowSums(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.this.subset$counts>=5)==dim(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.this.subset)[2],])

TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.filt.5count <- TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE[row.names(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE)%in%mir.id.filt.gt.5,]
row.mean.abund.df.filt.5count <- row.mean.abund.df[row.names(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.filt.5count),]

TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.sort.truseq.abund.5count <- TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.filt.5count[order(-row.mean.abund.df.filt.5count$TruSeq.orig),]
pheatmap(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.sort.truseq.abund.5count,
         cluster_rows = FALSE,
         show_rownames = FALSE,
         fontsize = 14,
         fontsize_col = 12,
         labels_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot$transf.LibMethod.id,
         annotation_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot[, 1:2], 
         annotation_row = row.mean.abund.df.filt.5count,
         annotation_colors=ann.colors.plasm 
)


### PLOT THE DOWNBIASED EQUIMOLAR HEATMAPS TOO
TruSeq.DOWNBIAS.Scaled.And.Orig.eq <- ALL.COMPARISON.SCALING.FACTOR.CORRECTED.MATR.LIST.UNRESCALED.EQ[["TruSeq.V.DOWNBIAS"]]

TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO <- data.frame(strsplit2(colnames(TruSeq.DOWNBIAS.Scaled.And.Orig.eq), split=".", fixed=TRUE), stringsAsFactors = FALSE, row.names=colnames(TruSeq.DOWNBIAS.Scaled.And.Orig.eq))
colnames(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO) <- c("lib.method.detail", "Lab", "pool", "replicate", "orig.or.transf")
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO$lab.libMethod.pool.replicate <- apply(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO[, c("lib.method.detail", "Lab", "pool", "replicate")], 1, paste, collapse=".")
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO$lab.libMethod.pool.replicate.transf <- row.names(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO)

TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG <- data.table(merge(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO, ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples, by.x=which(colnames(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO)=="lab.libMethod.pool.replicate"), by.y=0, all.x = TRUE))
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG[, is.corrected:=ifelse(orig.or.transf=="orig", "ORIGINAL", "CORRECTED")]
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG[, transf.LibMethod.id:=sub("X4N", "4N", sub(".SynthEQ", "", make.names(paste(lab.libMethod.pool.replicate, substr(is.corrected, 1, 4), sep="."), unique = TRUE)))]
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG[, lib.method.detail:=sub("^X4N", "4N", lib.method.detail.x)]
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG[, c("lib.method.detail.x", "lib.method.detail.y"):=NULL]
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG[, `:=`(from.method=ifelse(lab.libMethod%in%downbias.labLibMethod, "Downbias.ref", lib.method.detail), to.method=ifelse(is.corrected=="CORRECTED" | lab.libMethod%in%downbias.labLibMethod, "Downbias.ref", lib.method.detail)), by=1:nrow(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG)]

TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df <- data.frame(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG, row.names=which(colnames(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG)=="lab.libMethod.pool.replicate.transf"))
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df <- TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df[colnames(TruSeq.DOWNBIAS.Scaled.And.Orig.eq),]                                                                                     
lib.simple.and.detail.levels.with.downbias <- c(lib.simple.levels, "Downbias.ref", lib.detail.levels[!lib.detail.levels%in%lib.simple.levels])
lib.levels.use <- lib.simple.and.detail.levels.with.downbias[lib.simple.and.detail.levels.with.downbias%in%unique(c(as.character(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df$from.method), as.character(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df$to.method)))]
lib.levels.use.pos <- match(lib.levels.use, lib.simple.and.detail.levels.with.downbias)

TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df$from.method <- factor(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df$from.method, lib.levels.use)
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df$to.method <- factor(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df$to.method, lib.levels.use)
TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df.for.plot <- TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df[, c("to.method", "from.method", "transf.LibMethod.id")]
groups <- apply(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df[,c("from.method", "to.method")], 1, FUN=function(x){
  g.tmp <- ifelse(x[1]==x[2], paste0(x[1], ".orig"), paste0(x[1], ".to.", x[2]))
  sub(".to.Downbias.ref", ".DOWNBIASED", g.tmp)
} )
design.scl <- model.matrix(~0+groups)
TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE <- voom(2^TruSeq.DOWNBIAS.Scaled.And.Orig.eq, design.scl, plot=TRUE)$E

# Custom palatte.. Does Truseq, CleanTag, NEB and 4N as normal. Then adds a greenish color for DOWNBIAS, and then takes the hue_pal_4 value and shades them for the different 4N methods
custom.lib.method.pallate <- c(hue_pal()(4), "#3FA1FF", hcl(285, c = 100, l = 65, alpha = 0.9*seq(4,1)/4))

ann.colors.plasm = list(
  from.method = custom.lib.method.pallate[lib.levels.use.pos],
  to.method = custom.lib.method.pallate[lib.levels.use.pos],
  is.corrected = grey_pal()(2)
)
names(ann.colors.plasm$from.method) <- lib.simple.and.detail.levels.with.downbias[lib.levels.use.pos]
names(ann.colors.plasm$to.method) <- lib.simple.and.detail.levels.with.downbias[lib.levels.use.pos]
names(ann.colors.plasm$is.corrected) <- c("TRANSFORMED", "ORIGINAL")
row.mean.abund.df <- data.frame(FUNCTION.aveTechReps(TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE, groups))
# Full heatmap with row  clustering
pheatmap(TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE,
         show_rownames = FALSE,
         fontsize = 14,
         fontsize_col = 12,
         labels_col = TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df.for.plot$transf.LibMethod.id,
         annotation_col = TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df.for.plot[, 1:2], 
         annotation_row = row.mean.abund.df,
         annotation_colors=ann.colors.plasm 
)

# Sort heatmap
row.mean.abund.df.sort <- row.mean.abund.df[order(-row.mean.abund.df$TruSeq.orig),]
# Full heatmap with row  clustering
pheatmap(TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE[row.names(row.mean.abund.df.sort),],
         show_rownames = FALSE,
         cluster_rows = FALSE,
         fontsize = 14,
         fontsize_col = 12,
         labels_col = TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df.for.plot$transf.LibMethod.id,
         annotation_col = TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG.df.for.plot[, 1:2], 
         annotation_row = row.mean.abund.df.sort,
         annotation_colors=ann.colors.plasm 
)


TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.MELT <- data.table(melt(TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE))
setnames(TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.MELT, c("miR.id", "lab.libMethod.pool.replicate.transf", "cpm"))
setkey(TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.MELT, lab.libMethod.pool.replicate.transf)
setkey(TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG, lab.libMethod.pool.replicate.transf)
TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.MELT.WITH.SAMPLE.INFO <- merge(TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.MELT, TruSeq.DOWNBIAS.Scaled.And.Orig.SAMPLE.INFO.WITH.ORIG)
TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.MELT.WITH.SAMPLE.INFO[, `:=`(from.method=factor(from.method, lib.levels.use), to.method=factor(to.method, lib.levels.use))]
TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.MELT.WITH.SAMPLE.INFO.summary.cpm <- TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.MELT.WITH.SAMPLE.INFO[, .(mean.cpm=mean(cpm)), by=.(miR.id, Lab.x, pool.x, orig.or.transf, lib.size, lab.libMethod.pool, lab.libMethod, lib.method.simple, GROUND.TRUTH.SAMPLE, is.corrected)]
TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.MELT.WITH.SAMPLE.INFO.summary.cpm[, lib.method.simple:=factor(lib.method.simple, levels=lib.simple.levels)]
ann.colors.viol = list(
  lib.method.simple = hue_pal()(4)
)
names(ann.colors.viol$lib.method.simple) <- levels(TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.MELT.WITH.SAMPLE.INFO.summary.cpm$lib.method.simple)
TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.MELT.WITH.SAMPLE.INFO.summary.cpm.filt <- TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.MELT.WITH.SAMPLE.INFO.summary.cpm[lab.libMethod%in%sample.levels]
TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.MELT.WITH.SAMPLE.INFO.summary.cpm.filt[, lab.libMethod:=factor(lab.libMethod, sample.levels)]
TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.MELT.WITH.SAMPLE.INFO.summary.cpm.filt[, is.corrected:=ifelse(as.character(lab.libMethod)%in%downbias.labLibMethod, "ORIGINAL", as.character(is.corrected))]
TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.MELT.WITH.SAMPLE.INFO.summary.cpm.filt[, is.corrected:=factor(is.corrected, levels=c("ORIGINAL", "CORRECTED"))]
ggplot(
  TruSeq.DOWNBIAS.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.MELT.WITH.SAMPLE.INFO.summary.cpm.filt, 
  aes(x=lab.libMethod, y=2^mean.cpm, fill=lib.method.simple, alpha=is.corrected)) + geom_violin(
    draw_quantiles = c(0.25,0.5,0.75),
    size=1.0
  ) + scale_fill_manual(values = ann.colors.viol$lib.method.simple) +
  scale_alpha_discrete(range=c(1,0.5)) +
  scale_y_log10(
    limits=c(10^-2,10^5),
    breaks = scales::trans_breaks("log10", n = 7,  function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) + 
  facet_wrap(~is.corrected, strip.position = "bottom", nrow = 1, scale="free_x") + #annotation_logticks(sides = "l", size = 1, color="black") + 
  theme_bw() +
  theme(
    axis.text=element_text(size=14, face="bold", color="black"),
    axis.text.x=element_text(angle=50, hjust=1),
    axis.line=element_line(color="black", size=1.5),
    axis.title=element_text(face="bold", size=16),
    plot.title=element_text(face="bold", size=18),
    panel.grid.major = element_line(color=NA),
    panel.border = element_blank(),
    panel.spacing = unit(0, "lines"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text=element_text(size=14, face="bold", color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=14, face="bold"),
    legend.title=element_text(size=15, face="bold"),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(
    title="EQUIMOLAR POOL:\nDOWNBIAS CORRECTION FACTORS",
    x="Lab",
    y="Average Counts-Per-Million",
    fill="Library Prep Method"
  ) 

ggplot(ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.DETECTED.ALL.BY.ITER[downsampled.to<downsample.to.vec[length(downsample.to.vec)],],
       aes(
         x=factor(round(log10(downsampled.to), digits=1)),
         y=n.miRs.all.detected.by.iter, fill=lib.method.simple)) + 
  geom_boxplot() + 
  theme_bw() + 
  labs(x="log10 Sequencing Depth (total miRNA-mapping reads)") +
  facet_wrap(~lib.method.simple, strip.position = "bottom", nrow = 1) +
  theme(
    axis.text=element_text(size=14, face="bold", color="black"),
    axis.line=element_line(color="black", size=1.5),
    axis.title=element_text(face="bold", size=16),
    plot.title=element_text(face="bold", size=18),
    panel.grid.major = element_line(color=NA),
    panel.border = element_blank(),
    panel.spacing = unit(0, "lines"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text=element_text(size=14, face="bold", color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=14, face="bold"),
    legend.title=element_text(size=15, face="bold"),
    legend.key.size=unit(1.2, "lines")) 



############################
###
# Get values for Mann-W test
all.intermethod.comparison.groups <- grep("DOWNBIAS", unique(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset$correction), invert=TRUE, value=TRUE)
# UPDATE 1-10-17 -- Through several iterations, it looks like normalization issues aren't actually a big issue. 
# Filter for count>5 in all samples, and do the normalization up-front only on the sample we're applying the correction factors to
# No need for inter-lab normalization as it is.
MANN.WHITNEY.SHIFTS.PLASMA <- sapply(all.intermethod.comparison.groups, USE.NAMES=TRUE, simplify = FALSE, FUN=function(x.group){
  correct.dt <- data.frame(subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, correction==x.group), row.names=1)
  x.ids <- strsplit2(x.group, split=".V.")
  print(x.group)
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple%in%x.ids | ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.detail%in%x.ids]
  keep.mirs <- (rowSums(dge.test$counts>=5)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  query.ids <- grep(x.ids[1], colnames(dge.test.filt), value=TRUE, ignore.case=TRUE)
  ref.ids <- grep(x.ids[1], invert = TRUE, colnames(dge.test.filt), value=TRUE, ignore.case=TRUE)
  
  dge.query <- dge.test.filt[, query.ids]
  #design.query <- model.matrix(~0+lab.libMethod, dge.query$samples)
  dge.query <- calcNormFactors(dge.query, method="RLE")
  #v.query <- voom(dge.query, normalize.method="quantile")
  
  dge.ref <- dge.test.filt[, ref.ids]
  #design.ref <- model.matrix(~0+lab.libMethod, dge.ref$samples)
  dge.ref <- calcNormFactors(dge.ref, method="RLE")
  #v.ref <- voom(dge.ref, normalize.method="quantile")
  
  merged.miR.ids <- row.names(dge.test.filt)[row.names(dge.test.filt)%in%row.names(correct.dt)]
  #cpm.test.norm.filt.md1 <- v.query$E[merged.miR.ids,]
  cpm.test.norm.filt.md1 <- cpm.DGEList(dge.query, log = TRUE, normalized.lib.sizes = TRUE)[merged.miR.ids,]
  #cpm.test.norm.filt.ref <- v.ref$E[merged.miR.ids,]
  cpm.test.norm.filt.ref <- cpm.DGEList(dge.ref, log = TRUE, normalized.lib.sizes = TRUE)[merged.miR.ids,]
  tt.comp.scl.factor.filt <- correct.dt[merged.miR.ids,]
  
  cpm.test.norm.filt.md.1.scaled <-
    apply(
      cpm.test.norm.filt.md1,
      2,
      FUN = function(x) {
        x.log.scl <- x - tt.comp.scl.factor.filt$logFC
      })
  # combine for re-centering
  cpm.scaled.and.ref <- cbind(cpm.test.norm.filt.ref, cpm.test.norm.filt.md.1.scaled)[, colnames(dge.test.filt)]
  cpm.test.scaled <- voom(2^cpm.scaled.and.ref, design.test, normalize.method="quantile")$E
  cpm.test.scaled.group.means <- data.frame(mean.scaled.query=rowMeans(cpm.test.scaled[, query.ids]), mean.ref=rowMeans(cpm.test.scaled[, ref.ids]), row.names=row.names(cpm.test.scaled))
  cpm.test.scaled.group.means$abs.log.dif <- with(cpm.test.scaled.group.means, abs(mean.scaled.query-mean.ref))
  cpm.test.scaled.group.means$comparison <- "CORRECTED"
  
  # re-do the normalization for the original data to "quantile" so that it's normalized the same way as the scaled data.
  cpm.unscaled.and.ref <- cbind(cpm.test.norm.filt.ref, cpm.test.norm.filt.md1)[, colnames(dge.test.filt)]
  v.test.norm <- voom(2^cpm.unscaled.and.ref, design.test, normalize.method="quantile")
  cpm.test.norm <- v.test.norm$E
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  dge.test.filt <- calcNormFactors(dge.test.filt)
  v.test.norm <- voom(dge.test.filt, design.test, normalize.method="quantile")
  cpm.test.norm <- v.test.norm$E
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  cpm.test.original.group.means <- data.frame(mean.unscaled.query=rowMeans(cpm.test.norm.filt[, query.ids]), mean.ref=rowMeans(cpm.test.norm.filt[, ref.ids]), row.names=row.names(cpm.test.norm.filt))
  cpm.test.original.group.means$abs.log.dif <- with(cpm.test.original.group.means, abs(mean.unscaled.query-mean.ref))
  cpm.test.original.group.means$comparison <- "UNCORRECTED"
  cpm.both <- data.table(keep.rownames = TRUE, rbind(cpm.test.original.group.means[, c("comparison", "abs.log.dif")], cpm.test.scaled.group.means[, c("comparison", "abs.log.dif")]))
  cpm.both[, `:=`(group.id=paste0(x.ids[1], ".V.", x.ids[2]), from.method=x.ids[1], to.method=x.ids[2])]
  w.test.out <- wilcox.test(x = cpm.test.original.group.means$abs.log.dif, y=cpm.test.scaled.group.means$abs.log.dif, alternative = c("greater"),  conf.int = TRUE)
  cpm.both[, `:=`(W=w.test.out$statistic, p.value=w.test.out$p.value, alternative="greater", estimate=w.test.out$estimate)]
  #get matrix to include as well
  colnames(cpm.test.norm.filt.md.1.scaled) <- paste0(colnames(cpm.test.norm.filt.md.1.scaled), ".transf")
  colnames(cpm.unscaled.and.ref) <- paste0(colnames(cpm.unscaled.and.ref), ".orig")
  cpm.full.matr.with.correct <- cbind(cpm.unscaled.and.ref[merged.miR.ids,], cpm.test.norm.filt.md.1.scaled)
  cpm.full.matr.with.correct.norm <- normalizeBetweenArrays(cpm.full.matr.with.correct, method = "quantile")
  return(list(cpm.both=cpm.both, cpm.scaled.matr=cpm.full.matr.with.correct.norm))
}); MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF <- do.call(rbind, lapply(MANN.WHITNEY.SHIFTS.PLASMA, FUN=function(x) x[[1]])); MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF[, comparison:=factor(comparison, levels=c("UNCORRECTED", "CORRECTED"))]


MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS <- unique(subset(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF, select=colnames(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF)[!colnames(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF)%in%c("rn", "comparison", "abs.log.dif")]))
setorder(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS, p.value)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF[, group.id:=factor(group.id, levels=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS$group.id)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS[, group.id:=factor(group.id, levels=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS$group.id)]

ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF, aes(x=abs.log.dif, color=comparison)) + 
  geom_density(size=1.5) + 
  geom_text(data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS, aes(label=paste0("p = ", format(p.value, digits=3)), x=8.5, y=0.6),inherit.aes = FALSE, check_overlap=TRUE) + facet_wrap(~group.id) + theme_bw() + scale_x_continuous(breaks=seq(0, 16, 2), expand=c(0,0), limits=c(-0.1,14)) + scale_y_continuous(breaks=seq(0, 1, 0.1), expand=c(0,0))

# FIG6B -------------------------------------------------------------------


ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF[group.id =="TruSeq.V.4N"], aes(x=abs.log.dif, color=comparison)) + geom_density(size=1.5, adjust=1, trim=TRUE) + facet_wrap(~group.id) + theme_bw() + 
  scale_y_continuous(expand=c(0.0, 0.0), breaks=seq(0.0, 1, 0.1), limits=c(0.0, 0.8)) + 
  scale_x_continuous(expand=c(0.01,0), breaks=seq(0,9, 1.0)) + 
  theme(
    axis.text=element_text(size=10, color="black"),
    axis.title=element_text(size=10),
    plot.title=element_text(size=10),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=10, color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=10),
    legend.title=element_text(size=10),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(y="density",
       x="|Difference VS Reference|(log2)")
ggsave("FIG_6D_TruSeq_V_4N_Correction_DensityPlot_PLASMAPOOL.svg", units = "in", width=4, height=3)

# FIG6C -------------------------------------------------------------------


ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF[group.id =="TruSeq.V.4N"], aes(y=abs.log.dif, x=comparison, fill=comparison)) +
  geom_boxplot() + 
  facet_wrap(~group.id) + 
  theme_bw() + 
  geom_text(
    data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS[group.id =="TruSeq.V.4N"], 
    aes(
      label=paste0("p = ", format(p.value, digits=3)),
      x=2, y=8),
    inherit.aes = FALSE, check_overlap=TRUE) +
  scale_y_continuous(breaks=seq(0.0, 8.0, 0.5), limits=c(0.0, 8)) + 
  #scale_x_continuous(expand=c(0.01,0), limits=c(0, 8.0), breaks=seq(0,9, 1.0)) + 
  theme(
    axis.text=element_text(size=10, color="black"),
    axis.title=element_text(size=10),
    plot.title=element_text(size=10),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=10, color="black"),
    #legend.position = "top",
    #legend.direction = "horizontal",
    legend.text=element_text(size=10),
    legend.title=element_text(size=10),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(x=NULL,
       y="|Difference VS Reference|(log2)")
ggsave("FIG_6D_TruSeq_V_4N_Correction_BOXPLOT_PLASMAPOOL.svg", units = "in", width=4, height=2)
groups.plasm <- apply(TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot[, c("to.method", "from.method")], 1, FUN=function(x) ifelse(x[1]==x[2], paste0(x[1], ".orig"), paste0(x[2], ".to.", x[1])))

truseq.v.4N.matr.plasm <- MANN.WHITNEY.SHIFTS.PLASMA[["TruSeq.V.4N"]][[2]]
row.mean.abund.df.plasm <- data.frame(FUNCTION.aveTechReps(truseq.v.4N.matr.plasm, groups.plasm))
ann.colors.correction.facts <- list(from.method=ann_colors$lib.method.simple, to.method=ann_colors$lib.method.simple)


# FIG6A -------------------------------------------------------------------

library(dendsort)
dcols.truseq.v.4N.matr.plasm <- dist(t(truseq.v.4N.matr.plasm))
hc.cols <- hclust(dcols.truseq.v.4N.matr.plasm, method="complete")

# To reorder the plot as much as possible by cols
callback <- function(hc, mat){
  label.split <- strsplit2(hc$labels, split=".", fixed=TRUE)
  wt1 <- ifelse(label.split[, dim(label.split)[2]]=="transf", 1, 0)
  ref.grp <- unique(label.split[wt1>0, 1])
  wt2 <- ifelse(label.split[,1]==ref.grp, 1, 0)
  grp.wt.df <- data.frame(wt1=wt1, wt2=wt2, row.names=hc$labels)
  group.wts <- rowSums(grp.wt.df)
  dend <- reorder(as.dendrogram(hc), wts=group.wts)
  as.hclust(dend)
}
pheatmap(truseq.v.4N.matr.plasm[order(-row.mean.abund.df.plasm$TruSeq.orig),], legend_breaks = seq(0, 16, by=2)  , 
         breaks = seq(0, max(truseq.v.4N.matr.plasm), length.out = 101),
         cluster_rows = FALSE,
         show_rownames = FALSE,
         clustering_distance_cols = dcols.truseq.v.4N.matr.plasm,
         fontsize = 10,
         fontsize_col = 10,
         treeheight_col = 15,
         
         labels_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot$transf.LibMethod.id,
         annotation_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot[, 1:2], 
         annotation_row = row.mean.abund.df.plasm[, c(1, 3, 2)],
         annotation_colors=FUNCTION.filter.ann_colors(ann.colors.correction.facts, TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot),
         clustering_callback=callback, filename="FIG_6C_TruSeq_V_4N_Correction_HEATMAP_PLASMAPOOL_V2.pdf")
pheatmap(truseq.v.4N.matr.plasm[order(-row.mean.abund.df.plasm$TruSeq.orig),], legend_breaks = seq(0, 16, by=2)  , 
         breaks = seq(0, max(truseq.v.4N.matr.plasm), length.out = 101),
         cluster_rows = FALSE,
         show_rownames = FALSE,
         clustering_distance_cols = dcols.truseq.v.4N.matr.plasm,
         fontsize = 10,
         fontsize_col = 10,
         treeheight_col = 15,
         
         labels_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot$transf.LibMethod.id,
         annotation_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot[, 1:2], 
         annotation_row = row.mean.abund.df.plasm[, c(1, 3, 2)],
         annotation_colors=FUNCTION.filter.ann_colors(ann.colors.correction.facts, TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.WITH.ORIG.df.for.plot),
         clustering_callback=callback)


# Mann-whitney tests for equimolar pool
MANN.WHITNEY.SHIFTS.EQ <- sapply(all.intermethod.comparison.groups, USE.NAMES=TRUE, simplify = FALSE, FUN=function(x.group){
  correct.dt <- data.frame(subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, correction==x.group), row.names=1)
  x.ids <- strsplit2(x.group, split=".V.")
  print(x.group)
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[, ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.simple%in%x.ids | ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.detail%in%x.ids]
  keep.mirs <- (rowSums(dge.test$counts>=5)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  dge.test.filt <- calcNormFactors(dge.test.filt, method="RLE")
  v.test.norm <- voom(dge.test.filt, design.test, normalize.method="scale")
  cpm.test.norm <- v.test.norm$E
  
  merged.miR.ids <- merge(cpm.test.norm, correct.dt, by=0)$Row.names
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  query.ids <- grep(x.ids[1], colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)
  ref.ids <- grep(x.ids[1], invert = TRUE, colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)
  cpm.test.norm.filt.md1 <- cpm.test.norm.filt[, query.ids]
  cpm.test.norm.filt.ref <- cpm.test.norm.filt[, ref.ids]
  tt.comp.scl.factor.filt <- correct.dt[merged.miR.ids,]
  cpm.test.norm.filt.md.1.scaled <-
    apply(
      cpm.test.norm.filt.md1,
      2,
      FUN = function(x) {
        x.log.scl <- x - tt.comp.scl.factor.filt$logFC
      })
  # combine for re-centering
  cpm.scaled.and.ref <- cbind(cpm.test.norm.filt.ref, cpm.test.norm.filt.md.1.scaled)[, colnames(dge.test)]
  cpm.test.scaled <- voom(2^cpm.scaled.and.ref, design.test, normalize.method="scale")$E
  cpm.test.scaled.group.means <- data.frame(mean.scaled.query=rowMeans(cpm.test.scaled[, query.ids]), mean.ref=rowMeans(cpm.test.scaled[, ref.ids]), row.names=row.names(cpm.test.scaled))
  cpm.test.scaled.group.means$abs.log.dif <- with(cpm.test.scaled.group.means, abs(mean.scaled.query-mean.ref))
  cpm.test.scaled.group.means$comparison <- "CORRECTED"
  cpm.test.original.group.means <- data.frame(mean.unscaled.query=rowMeans(cpm.test.norm.filt[, query.ids]), mean.ref=rowMeans(cpm.test.norm.filt[, ref.ids]), row.names=row.names(cpm.test.norm.filt))
  cpm.test.original.group.means$abs.log.dif <- with(cpm.test.original.group.means, abs(mean.unscaled.query-mean.ref))
  cpm.test.original.group.means$comparison <- "UNCORRECTED"
  cpm.both <- data.table(keep.rownames = TRUE, rbind(cpm.test.original.group.means[, c("comparison", "abs.log.dif")], cpm.test.scaled.group.means[, c("comparison", "abs.log.dif")]))
  cpm.both[, `:=`(group.id=paste0(x.ids[1], ".V.", x.ids[2]), from.method=x.ids[1], to.method=x.ids[2])]
  w.test.out <- wilcox.test(x = cpm.test.original.group.means$abs.log.dif, y=cpm.test.scaled.group.means$abs.log.dif, alternative = c("greater"),  conf.int = TRUE)
  cpm.both[, `:=`(W=w.test.out$statistic, p.value=w.test.out$p.value, alternative="greater", estimate=w.test.out$estimate)]
  
  #get matrix to include as well
  colnames(cpm.test.norm.filt.md.1.scaled) <- paste0(colnames(cpm.test.norm.filt.md.1.scaled), ".transf")
  colnames(cpm.test.norm.filt) <- paste0(colnames(cpm.test.norm.filt), ".orig")
  cpm.full.matr.with.correct <- cbind(cpm.test.norm.filt, cpm.test.norm.filt.md.1.scaled)
  cpm.full.matr.with.correct.norm <- normalizeBetweenArrays(cpm.full.matr.with.correct, method = "scale")
  return(list(cpm.both=cpm.both, cpm.scaled.matr=cpm.full.matr.with.correct.norm))
}); MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF <- do.call(rbind, lapply(MANN.WHITNEY.SHIFTS.EQ, FUN=function(x) x[[1]])); MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF[, comparison:=factor(comparison, levels=c("UNCORRECTED", "CORRECTED"))]

MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF.STATS <- MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF[comparison=="UNCORRECTED", .(W=max(W), p.value=min(p.value), estimate=min(estimate), median.abs.log.dif=quantile(abs.log.dif, 0.9)), by=.(group.id, from.method, to.method, alternative)]
setorder(MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF.STATS, p.value)
MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF[, group.id:=factor(group.id, levels=MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF.STATS$group.id)]
MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF.STATS[, group.id:=factor(group.id, levels=MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF.STATS$group.id)]

ggplot(MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF, aes(x=abs.log.dif, color=comparison)) + 
  geom_density(size=1.5) + 
  geom_text(data=MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF.STATS, aes(label=paste0("p = ", format(p.value, digits=3)), x=8.5, y=0.6),inherit.aes = FALSE, check_overlap=TRUE) + facet_wrap(~group.id) + theme_bw() + scale_x_continuous(breaks=seq(0, 16, 2), expand=c(0,0), limits=c(-0.1,14)) + scale_y_continuous(breaks=seq(0, 1, 0.1), expand=c(0,0))

ggplot(MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF[group.id=="TruSeq.V.4N"], aes(x=abs.log.dif, color=comparison)) + 
  geom_density(size=1.5) + 
  geom_text(data=MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF.STATS, aes(label=paste0("p = ", format(p.value, digits=3)), x=6, y=2),inherit.aes = FALSE, check_overlap=TRUE) #+ facet_wrap(~group.id) + theme_bw() #+ scale_x_continuous(breaks=seq(0, 16, 2), expand=c(0,0), limits=c(-0.1,14)) + scale_y_continuous(breaks=seq(0, 1, 0.1), expand=c(0,0))
ggplot(MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF[group.id =="TruSeq.V.4N"], aes(x=abs.log.dif, color=comparison)) + geom_density(size=1.5, trim=TRUE) + facet_wrap(~group.id) + theme_bw() + 
  scale_y_continuous(expand=c(0.0, 0.0), breaks=seq(0.0, 10.0, 0.25)) + 
  scale_x_continuous(expand=c(0.0,0.01), breaks=seq(0,9, 1.0)) + 
  theme(
    axis.text=element_text(size=10, color="black"),
    axis.title=element_text(size=10),
    plot.title=element_text(size=10),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=10, color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=10),
    legend.title=element_text(size=10),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(y="probability density",
       x="|Difference VS Reference|(log2)")
ggsave("OLD_FIG_6B_TruSeq_V_4N_Correction_DensityPlot_EQPOOL_SUPPLEMENT.svg", units = "in", width=4, height=3)

ggplot(MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF[group.id =="TruSeq.V.4N"], aes(y=abs.log.dif, x=comparison, fill=comparison)) +
  geom_boxplot() + 
  facet_wrap(~group.id) + 
  theme_bw() + 
  geom_text(
    data=MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF.STATS[group.id =="TruSeq.V.4N"], 
    aes(
      label=paste0("p = ", format(p.value, digits=3)),
      x=2, y=8),
    inherit.aes = FALSE, check_overlap=TRUE) +
  scale_y_continuous(breaks=seq(0.0, 8.0, 0.5), limits=c(0.0, 8)) + 
  #scale_x_continuous(expand=c(0.01,0), limits=c(0, 8.0), breaks=seq(0,9, 1.0)) + 
  theme(
    axis.text=element_text(size=10, color="black"),
    axis.title=element_text(size=10),
    plot.title=element_text(size=10),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=10, color="black"),
    #legend.position = "top",
    #legend.direction = "horizontal",
    legend.text=element_text(size=10),
    legend.title=element_text(size=10),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(x=NULL,
       y="|Difference VS Reference|(log2)")
ggsave("OLD_FIG_6B_TruSeq_V_4N_Correction_BOXPLOT_EQPOOL_SUPPLEMENT.svg", units = "in", width=4, height=2)


MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF.TRUSEQ.V.4N <- subset(MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF, group.id=="TruSeq.V.4N")



truseq.v.4N.matr <- MANN.WHITNEY.SHIFTS.EQ[["TruSeq.V.4N"]][[2]]

row.mean.abund.df.eq <- data.frame(FUNCTION.aveTechReps(truseq.v.4N.matr, groups.eq))

pheatmap(truseq.v.4N.matr[order(-row.mean.abund.df.eq$TruSeq.orig),],
         cluster_rows = FALSE,
         show_rownames = FALSE,
         fontsize = 10,
         fontsize_col = 10,
         treeheight_col = 15,
         labels_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df.for.plot$transf.LibMethod.id,
         annotation_col = TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df.for.plot[, 1:2], 
         annotation_row = row.mean.abund.df.eq,
         annotation_colors= FUNCTION.filter.ann_colors(ann.colors.correction.facts, TruSeq.V.4N.Scaled.And.Orig.SCALING.FACTOR.SCALED.MATR.MERGE.ORIG.RESCALE.SAMPLE.INFO.EQ.WITH.ORIG.df.for.plot), filename="OLD_FIG_6A_TruSeq_V_4N_Correction_HEATMAP_EQPOOL_SUPPLEMENT.pdf")



######
# Same Mann-Whitney tests, except do the tests pairwise for all replicates and then summarize the p-values. This summarizes pvalues for all miRNAs in the set, doing wilcox tests row-wise.
library(matrixStats)
MANN.WHITNEY.SHIFTS.PLASMA.FISHER <- sapply(all.intermethod.comparison.groups, USE.NAMES=TRUE, simplify = FALSE, FUN=function(x.group){
  #x.group <- "TruSeq.V.4N"
  correct.dt <- data.frame(subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, correction==x.group), row.names=1)
  x.ids <- strsplit2(x.group, split=".V.")
  print(x.group)
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple%in%x.ids | ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.detail%in%x.ids]
  keep.mirs <- (rowSums(dge.test$counts>=5)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  dge.test.filt <- calcNormFactors(dge.test.filt, method="RLE")
  cpm.test.norm <- cpm.DGEList(dge.test.filt, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 0)
  
  merged.miR.ids <- merge(cpm.test.norm, correct.dt, by=0)$Row.names
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  query.ids <- grep(x.ids[1], colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)
  ref.ids <- grep(x.ids[1], invert = TRUE, colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)
  cpm.test.norm.filt.md1 <- cpm.test.norm.filt[, query.ids]
  cpm.test.norm.filt.ref <- cpm.test.norm.filt[, ref.ids]
  tt.comp.scl.factor.filt <- correct.dt[merged.miR.ids,]
  cpm.test.norm.filt.md.1.scaled <-
    apply(
      cpm.test.norm.filt.md1,
      2,
      FUN = function(x) {
        x.log.scl <- x - tt.comp.scl.factor.filt$logFC
      })
  # combine for re-centering
  cpm.scaled.and.ref <- cbind(cpm.test.norm.filt.ref, cpm.test.norm.filt.md.1.scaled)[, colnames(cpm.test.norm.filt)]
  cpm.test.scaled <- voom(2^cpm.scaled.and.ref, design.test, normalize.method="scale")$E
  
  
  cpm.test.scaled.group.means <- abs(cpm.test.scaled[, query.ids] - rowMeans(cpm.test.scaled[, ref.ids]))
  cpm.test.original.group.means <- abs(cpm.test.norm.filt[, query.ids] - rowMeans(cpm.test.norm.filt[, ref.ids]))
  wilcox.tests <- sapply(row.names(cpm.test.scaled.group.means), simplify = FALSE, USE.NAMES=TRUE, FUN=function(i){
    wilcox.test(x=cpm.test.original.group.means[i, ], y=cpm.test.scaled.group.means[i, ], alternative="greater")
  })
  wilcox.pvals <- unlist(lapply(wilcox.tests, FUN=function(x){x$p.value}))
  df <- 2*length(wilcox.pvals)
  p.combined <- pchisq(-2*sum(log(wilcox.pvals)), df, lower.tail=FALSE)
  
  cpm.scaled.df <- data.frame(comparison="CORRECTED", mean.abs.log.dif=rowMeans(cpm.test.scaled.group.means), sd.abs.log.dif=rowSds(cpm.test.scaled.group.means), pval=wilcox.pvals[row.names(cpm.test.scaled.group.means)])
  cpm.original.df <- data.frame(comparison="UNCORRECTED", mean.abs.log.dif=rowMeans(cpm.test.original.group.means), sd.abs.log.dif=rowSds(cpm.test.original.group.means), pval=wilcox.pvals[row.names(cpm.test.original.group.means)])
  
  cpm.both <- data.table(keep.rownames = TRUE, rbind(cpm.original.df, cpm.scaled.df))
  cpm.both[, `:=`(group.id=paste0(x.ids[1], ".V.", x.ids[2]), from.method=x.ids[1], to.method=x.ids[2])]
  cpm.both[, `:=`(fisher.p.val=p.combined, alternative="greater")]
  return(cpm.both)
}); MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER <- do.call(rbind, MANN.WHITNEY.SHIFTS.PLASMA.FISHER); MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER[, comparison:=factor(comparison, levels=c("UNCORRECTED", "CORRECTED"))]

MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS.PLASMA.FISHER.STATS <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER[comparison=="UNCORRECTED", .(p.value=min(fisher.p.val)), by=.(group.id, from.method, to.method, alternative)]
setorder(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS.PLASMA.FISHER.STATS, p.value)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER[, group.id:=factor(group.id, levels=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS.PLASMA.FISHER.STATS$group.id)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS.PLASMA.FISHER.STATS[, group.id:=factor(group.id, levels=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS.PLASMA.FISHER.STATS$group.id)]

ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER, aes(x=mean.abs.log.dif, color=comparison)) + 
  geom_density(size=1.5) + 
  geom_text(data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS.PLASMA.FISHER.STATS, aes(label=paste0("p = ", format(p.value, digits=3)), x=8.5, y=0.6),inherit.aes = FALSE, check_overlap=TRUE) + facet_wrap(~group.id) + theme_bw() + scale_x_continuous(breaks=seq(0, 16, 2), expand=c(0,0), limits=c(-0.1,14)) + scale_y_continuous(breaks=seq(0, 1, 0.1), expand=c(0,0))
ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER[comparison=="CORRECTED"], aes(y=-1*log10(pval), x=group.id, fill=from.method)) + geom_violin()



# 1-10-17 Now do pairwise by sample--same way as standart method, but dont average over labs
library(limma)
MANN.WHITNEY.SHIFTS.PLASMA.FISHER.BY.LAB <- sapply(all.intermethod.comparison.groups, USE.NAMES=TRUE, simplify = FALSE, FUN=function(x.group){
  correct.dt <- data.frame(subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, correction==x.group), row.names=1)
  x.ids <- strsplit2(x.group, split=".V.")
  print(x.group)
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple%in%x.ids | ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.detail%in%x.ids]
  keep.mirs <- (rowSums(dge.test$counts>=5)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  query.ids <- grep(x.ids[1], colnames(dge.test.filt), value=TRUE, ignore.case=TRUE)
  ref.ids <- grep(x.ids[1], invert = TRUE, colnames(dge.test.filt), value=TRUE, ignore.case=TRUE)
  
  dge.query <- dge.test.filt[, query.ids]
  #design.query <- model.matrix(~0+lab.libMethod, dge.query$samples)
  dge.query <- calcNormFactors(dge.query, method="RLE")
  #v.query <- voom(dge.query, normalize.method="quantile")
  
  dge.ref <- dge.test.filt[, ref.ids]
  #design.ref <- model.matrix(~0+lab.libMethod, dge.ref$samples)
  dge.ref <- calcNormFactors(dge.ref, method="RLE")
  #v.ref <- voom(dge.ref, normalize.method="quantile")
  
  merged.miR.ids <- row.names(dge.test.filt)[row.names(dge.test.filt)%in%row.names(correct.dt)]
  #cpm.test.norm.filt.md1 <- v.query$E[merged.miR.ids,]
  cpm.test.norm.filt.md1 <- cpm.DGEList(dge.query, log = TRUE, normalized.lib.sizes = TRUE)[merged.miR.ids,]
  #cpm.test.norm.filt.ref <- v.ref$E[merged.miR.ids,]
  cpm.test.norm.filt.ref <- cpm.DGEList(dge.ref, log = TRUE, normalized.lib.sizes = TRUE)[merged.miR.ids,]
  tt.comp.scl.factor.filt <- correct.dt[merged.miR.ids,]
  
  cpm.test.norm.filt.md.1.scaled <-
    apply(
      cpm.test.norm.filt.md1,
      2,
      FUN = function(x) {
        x.log.scl <- x - tt.comp.scl.factor.filt$logFC
      })
  # combine for re-centering
  cpm.scaled.and.ref <- cbind(cpm.test.norm.filt.ref, cpm.test.norm.filt.md.1.scaled)[, colnames(dge.test.filt)]
  cpm.test.scaled <- voom(2^cpm.scaled.and.ref, design.test, normalize.method="quantile")$E
  mean.cpm.test.scaled.ref <- rowMeans(cpm.test.scaled[, ref.ids])
  query.cpm.test.scaled <- cpm.test.scaled[, query.ids]
  abs.dif.scaled.vs.ref <- abs(query.cpm.test.scaled-mean.cpm.test.scaled.ref)
  
  # re-do the normalization for the original data to "quantile" so that it's normalized the same way as the scaled data.
  cpm.unscaled.and.ref <- cbind(cpm.test.norm.filt.ref, cpm.test.norm.filt.md1)[, colnames(dge.test.filt)]
  v.test.norm <- voom(2^cpm.unscaled.and.ref, design.test, normalize.method="quantile")
  cpm.test.norm <- v.test.norm$E
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  query.cpm.test.norm.filt <- cpm.test.norm.filt[, query.ids]
  mean.cpm.test.norm.filt.ref <- rowMeans(cpm.test.norm.filt[, ref.ids])
  abs.dif.unscaled.vs.ref <- abs(query.cpm.test.norm.filt-mean.cpm.test.norm.filt.ref)
  
  # Now get lab-specific differences
  lab.ids <- unique(dge.query$samples$Lab)
  
  lab.wilcox.tests <- sapply(lab.ids, simplify=FALSE, USE.NAMES = TRUE, FUN=function(x){
    print(paste0("\tWILCOX.TEST: ", x))
    lab.query.ids <- query.ids[(strsplit2(query.ids, split=".", fixed=TRUE)[,2])==x]
    
    lab.cpm.test.scaled.group.means <- data.frame(mean.query.cpm.lab = rowMeans(query.cpm.test.scaled[, lab.query.ids]),
                                                  mean.ref.cpm.all.labs = mean.cpm.test.scaled.ref,  
                                                  abs.log.dif = rowMeans(abs.dif.scaled.vs.ref[, lab.query.ids]),
                                                  comparison = "CORRECTED",
                                                  query.Lab = x,
                                                  row.names=row.names(query.cpm.test.scaled))
    
    lab.cpm.test.original.group.means <- data.frame(mean.query.cpm.lab = rowMeans(query.cpm.test.norm.filt[, lab.query.ids]),
                                                    mean.ref.cpm.all.labs = mean.cpm.test.norm.filt.ref,  
                                                    abs.log.dif = rowMeans(abs.dif.unscaled.vs.ref[, lab.query.ids]),
                                                    comparison = "UNCORRECTED",
                                                    query.Lab = x,
                                                    row.names=row.names(query.cpm.test.scaled))
    lab.cpm.both <- data.table(keep.rownames = TRUE, rbind(lab.cpm.test.scaled.group.means, lab.cpm.test.original.group.means))
    lab.cpm.both[, `:=`(group.id=paste0(x.ids[1], ".V.", x.ids[2]), from.method=x.ids[1], to.method=x.ids[2])]
    
    lab.w.test.out <- wilcox.test(x = lab.cpm.test.original.group.means$abs.log.dif, y=lab.cpm.test.scaled.group.means$abs.log.dif, alternative = c("greater"),  conf.int = TRUE)
    lab.cpm.both[, `:=`(W=lab.w.test.out$statistic, p.value=lab.w.test.out$p.value, alternative="greater", estimate=lab.w.test.out$estimate)]
    return(list(full.table.stats=lab.cpm.both, wilcox.test.only=as.list(lab.w.test.out)))
  })
  lab.wilcox.tests.dt <- do.call(rbind, lapply(lab.wilcox.tests, FUN=function(x){x[[1]]}))
  wilcox.tests <- lapply(lab.wilcox.tests, FUN=function(x){x[[2]]})
  
  wilcox.pvals <- unlist(lapply(wilcox.tests, FUN=function(x){x$p.value}))
  df <- 2*length(wilcox.pvals)
  p.combined <- pchisq(-2*sum(log(wilcox.pvals)), df, lower.tail=FALSE)
  
  lab.wilcox.tests.dt[, `:=`(fisher.p.val=p.combined, fisher.df=df)]
  return(lab.wilcox.tests.dt)
});

MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB <- do.call(rbind, MANN.WHITNEY.SHIFTS.PLASMA.FISHER.BY.LAB)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB[, comparison:=factor(comparison, levels=c("UNCORRECTED", "CORRECTED"))]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.STATS <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB[, .(abs.log.dif=max(abs.log.dif)), by=c("group.id", "from.method", "to.method", "query.Lab", "W", "p.value", "alternative", "estimate", "fisher.p.val", "fisher.df", "comparison")]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.STATS[, max.abs.log.dif.group:=max(abs.log.dif), by=group.id]

ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB[group.id =="TruSeq.V.4N"],
       aes(y=abs.log.dif,
           x=query.Lab,
           fill=comparison)
) +
  geom_boxplot(position = "dodge") + 
  facet_wrap(~group.id) + 
  theme_bw() + 
  geom_text(
    data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.STATS[group.id =="TruSeq.V.4N" & comparison=="CORRECTED"], 
    aes(
      label=paste0("p = ", format(p.value, digits=3)),
      y=max.abs.log.dif.group*1.08
    ),
    inherit.aes = TRUE, check_overlap=TRUE) +
  geom_text(
    data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.STATS[group.id =="TruSeq.V.4N" & comparison=="CORRECTED"], 
    aes(
      label=paste0("w = ", W),
      y=max.abs.log.dif.group*1.05
    ),
    inherit.aes = TRUE, check_overlap=TRUE) +
  geom_text(
    data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.STATS[group.id =="TruSeq.V.4N" & comparison=="CORRECTED" & query.Lab=="Lab3"], 
    aes(
      label=paste0("FDR = ", format(fisher.p.val, digits=3), "; DF = ", fisher.df, "; estimate = ", format(estimate, digits=3)),
      y=max.abs.log.dif.group*1.12
    ),
    inherit.aes = TRUE, check_overlap=TRUE) +
  scale_y_continuous(expand=c(0.01,0), limits = c(0,10), breaks=seq(0, 12, 1.0)) + 
  theme(
    axis.text=element_text(size=10, color="black"),
    axis.title=element_text(size=10),
    plot.title=element_text(size=10),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=10, color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=10),
    legend.title=element_text(size=10),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(x=NULL,
       y="|Difference VS Reference|(log2)")
ggsave("FIG_6D_TruSeq_V_4N_Correction_BOXPLOT_ALL_LABS_PLASMAPOOL.svg", units = "in", width=5, height=2.5)

gid <- "TruSeq.V.4N_B"
ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB[group.id == gid],
       aes(y=abs.log.dif,
           x=query.Lab,
           fill=comparison)
) +
  geom_boxplot(position = "dodge") + 
  facet_wrap(~group.id) + 
  theme_bw() + 
  geom_text(
    data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.STATS[group.id == gid & comparison=="CORRECTED"], 
    aes(
      label=paste0("p = ", format(p.value, digits=3)),
      y=max.abs.log.dif.group*1.08
    ),
    inherit.aes = TRUE, check_overlap=TRUE) +
  geom_text(
    data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.STATS[group.id == gid & comparison=="CORRECTED"], 
    aes(
      label=paste0("w = ", W),
      y=max.abs.log.dif.group*1.05
    ),
    inherit.aes = TRUE, check_overlap=TRUE) +
  geom_text(
    data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.STATS[group.id == gid & comparison=="CORRECTED" & query.Lab=="Lab3"], 
    aes(
      label=paste0("FDR = ", format(fisher.p.val, digits=3), "; DF = ", fisher.df, "; estimate = ", format(estimate, digits=3)),
      y=max.abs.log.dif.group*1.12
    ),
    inherit.aes = TRUE, check_overlap=TRUE) +
  scale_y_continuous(expand=c(0.01,0), limits = c(0, 12), breaks=seq(0, 12, 1.0)) + 
  theme(
    axis.text=element_text(size=10, color="black"),
    axis.title=element_text(size=10),
    plot.title=element_text(size=10),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=10, color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=10),
    legend.title=element_text(size=10),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(x=NULL,
       y="|Difference VS Reference|(log2)")

MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.TruSeq.V.4N <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB[group.id =="TruSeq.V.4N"]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.TruSeq.V.4N <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF[group.id =="TruSeq.V.4N"]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.TruSeq.V.4N[, `:=`(query.Lab="Aggregate", fisher.p.val=NA, fisher.df=NA)]

MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.TruSeq.V.4N.Add.Agg <- rbind(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.TruSeq.V.4N, MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.TruSeq.V.4N, fill=TRUE)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.TruSeq.V.4N.Add.Agg.STATS <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.TruSeq.V.4N.Add.Agg[, .(abs.log.dif=max(abs.log.dif)), by=c("group.id", "from.method", "to.method", "query.Lab", "W", "p.value", "alternative", "estimate", "fisher.p.val", "fisher.df", "comparison")]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.TruSeq.V.4N.Add.Agg.STATS[, max.abs.log.dif.group:=max(abs.log.dif), by=group.id]

ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.TruSeq.V.4N.Add.Agg,
       aes(y=abs.log.dif,
           x=query.Lab,
           fill=comparison)
) +
  geom_boxplot(position = "dodge") + 
  facet_wrap(~group.id) + 
  theme_bw() + 
  geom_text(
    data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.TruSeq.V.4N.Add.Agg.STATS[group.id =="TruSeq.V.4N" & comparison=="CORRECTED"], 
    aes(
      label=format(p.value, digits=3),
      y=max.abs.log.dif.group*1.08
    ),
    inherit.aes = TRUE, check_overlap=TRUE) +
  geom_text(
    data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.TruSeq.V.4N.Add.Agg.STATS[group.id =="TruSeq.V.4N" & comparison=="CORRECTED"], 
    aes(
      label=W,
      y=max.abs.log.dif.group*1.05
    ),
    inherit.aes = TRUE, check_overlap=TRUE) +
  geom_text(
    data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.TruSeq.V.4N.Add.Agg.STATS[group.id =="TruSeq.V.4N" & comparison=="CORRECTED" & query.Lab=="Lab3"], 
    aes(
      label=paste0("FDR = ", format(fisher.p.val, digits=3), "; DF = ", fisher.df, "; estimate = ", format(estimate, digits=3)),
      y=max.abs.log.dif.group*1.12
    ),
    inherit.aes = TRUE, check_overlap=TRUE) +
  scale_y_continuous(expand=c(0.01,0), limits = c(0,10), breaks=seq(0, 12, 1.0)) + 
  theme(
    axis.text=element_text(size=10, color="black"),
    axis.title=element_text(size=10),
    plot.title=element_text(size=10),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=10, color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=10),
    legend.title=element_text(size=10),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(x=NULL,
       y="|Difference VS Reference|(log2)")
ggsave("FIG_6D_TruSeq_V_4N_Correction_BOXPLOT_ALL_LABS_PLASMAPOOL_With_agg_data.svg", units = "in", width=8, height=3)


# 2/9/17 Do mann-whitney by miR
#MANN.WHITNEY.SHIFTS.PLASMA.BY.MIR <- sapply(all.intermethod.comparison.groups, USE.NAMES=TRUE, simplify = FALSE, FUN=function(x.group){
x.group <- "TruSeq.V.4N"
correct.dt <- data.frame(subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, correction==x.group), row.names=1)
x.ids <- strsplit2(x.group, split=".V.")
print(x.group)
#Finally, merge back with plasma data.
dge.test <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple%in%x.ids | ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.detail%in%x.ids]
keep.mirs <- (rowSums(dge.test$counts>=5)==dim(dge.test)[2]) 
design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
dge.test.filt <- dge.test[keep.mirs,]
dge.test.filt <- calcNormFactors(dge.test.filt, method="RLE")
cpm.test.norm <- cpm.DGEList(dge.test.filt, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 0)

merged.miR.ids <- merge(cpm.test.norm, correct.dt, by=0)$Row.names
cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
query.ids <- grep(x.ids[1], colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)
ref.ids <- grep(x.ids[1], invert = TRUE, colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)
cpm.test.norm.filt.md1 <- cpm.test.norm.filt[, query.ids]
cpm.test.norm.filt.ref <- cpm.test.norm.filt[, ref.ids]
tt.comp.scl.factor.filt <- correct.dt[merged.miR.ids,]
cpm.test.norm.filt.md.1.scaled <-
  apply(
    cpm.test.norm.filt.md1,
    2,
    FUN = function(x) {
      x.log.scl <- x - tt.comp.scl.factor.filt$logFC
    })
# combine for re-centering
cpm.scaled.and.ref <- cbind(cpm.test.norm.filt.ref, cpm.test.norm.filt.md.1.scaled)[, colnames(cpm.test.norm.filt)]
cpm.test.scaled <- voom(2^cpm.scaled.and.ref, design.test, normalize.method="scale")$E

# do differential expression
design.test <- model.matrix(~0+lib.method.simple, data = dge.test.filt$samples)
colnames(design.test) <- c("x4N", "TruSeq")
contrasts <- makeContrasts(TruSeq-x4N, levels=design.test)
voom.test.scaled <- voom(cpm.scaled.and.ref, design.test, normalize.method="scale")
fit.scaled <- lmFit(voom.test.scaled, design.test)
fit.scaled.2 <- contrasts.fit(fit.scaled, contrasts)
fit.scaled.2 <- treat(fit.scaled.2, lfc=log2(1.1), robust=TRUE)
tt.scaled <- topTreat(fit.scaled.2, confint=TRUE, number=Inf)
voom.test.orig <- voom(dge.test.filt[merged.miR.ids,], design.test, normalize.method="scale")
fit.orig <- lmFit(voom.test.orig, design.test)
fit.orig.2 <- contrasts.fit(fit.orig, contrasts)
fit.orig.2 <- treat(fit.orig.2, lfc=log2(1.1), robust=TRUE)
tt.orig <- topTreat(fit.orig.2, confint=TRUE, number=Inf)
tt.merged <- data.table(merge(tt.scaled, tt.orig, by=0), keep.rownames = TRUE)
setnames(tt.merged, colnames(tt.merged), sub(".y$", ".orig", sub(".x$", ".scaled", colnames(tt.merged))))
tt.merged[, xtabs(~(adj.P.Val.orig<=0.01)+(adj.P.Val.scaled<=0.01))]
ggplot(tt.merged, aes(x=-1*log10(adj.P.Val.orig), y=-1*log10(adj.P.Val.scaled))) + geom_hline(yintercept = -1*log10(0.01)) + geom_vline(xintercept = -1*log10(0.01)) + geom_point() + theme_bw() + scale_y_continuous(breaks=seq(0,40,2)) + scale_x_continuous(breaks=seq(0,40,2))
ggsave("limma_before_after_adj_1.1fold.pdf")
fit.scaled <- lmFit(voom.test.scaled, design.test)
fit.scaled.2 <- contrasts.fit(fit.scaled, contrasts)
fit.scaled.2 <- treat(fit.scaled.2, lfc=log2(1.2), robust=TRUE)
tt.scaled <- topTreat(fit.scaled.2, confint=TRUE, number=Inf)
voom.test.orig <- voom(dge.test.filt[merged.miR.ids,], design.test, normalize.method="scale")
fit.orig <- lmFit(voom.test.orig, design.test)
fit.orig.2 <- contrasts.fit(fit.orig, contrasts)
fit.orig.2 <- treat(fit.orig.2, lfc=log2(1.2), robust=TRUE)
tt.orig <- topTreat(fit.orig.2, confint=TRUE, number=Inf)
tt.merged <- data.table(merge(tt.scaled, tt.orig, by=0), keep.rownames = TRUE)
setnames(tt.merged, colnames(tt.merged), sub(".y$", ".orig", sub(".x$", ".scaled", colnames(tt.merged))))
tt.merged[, xtabs(~(adj.P.Val.orig<=0.01)+(adj.P.Val.scaled<=0.01))]
ggplot(tt.merged, aes(x=-1*log10(adj.P.Val.orig), y=-1*log10(adj.P.Val.scaled))) + geom_hline(yintercept = -1*log10(0.01)) + geom_vline(xintercept = -1*log10(0.01)) + geom_point() + theme_bw() + scale_y_continuous(breaks=seq(0,40,2)) + scale_x_continuous(breaks=seq(0,40,2))
ggsave("limma_before_after_adj_1.2fold.pdf")
fit.scaled.2 <- treat(fit.scaled.2, lfc=log2(1.0), robust=TRUE)
tt.scaled <- topTreat(fit.scaled.2, confint=TRUE, number=Inf)
fit.orig.2 <- treat(fit.orig.2, lfc=log2(1.0), robust=TRUE)
tt.orig <- topTreat(fit.orig.2, confint=TRUE, number=Inf)
tt.merged <- data.table(merge(tt.scaled, tt.orig, by=0), keep.rownames = TRUE)
setnames(tt.merged, colnames(tt.merged), sub(".y$", ".orig", sub(".x$", ".scaled", colnames(tt.merged))))
ggplot(tt.merged, aes(x=-1*log10(adj.P.Val.orig), y=-1*log10(adj.P.Val.scaled))) + geom_hline(yintercept = -1*log10(0.01)) + geom_vline(xintercept = -1*log10(0.01)) + geom_point() + theme_bw() + scale_y_continuous(breaks=seq(0,40,2)) + scale_x_continuous(breaks=seq(0,40,2))
tt.merged[, xtabs(~(adj.P.Val.orig<=0.01)+(adj.P.Val.scaled<=0.01))]
ggsave("limma_before_after_adj_1.0fold.pdf")


plot(density(((2^rowMeans(cpm.test.scaled[, query.ids])) - (2^rowMeans(cpm.test.scaled[, ref.ids])))/((2^rowMeans(cpm.test.scaled[, ref.ids])))))

cpm.test.scaled.group.means <- abs(cpm.test.scaled[, query.ids] - rowMeans(cpm.test.scaled[, ref.ids]))
cpm.test.original.group.means <- abs(cpm.test.norm.filt[, query.ids] - rowMeans(cpm.test.norm.filt[, ref.ids]))
pair.ttest.orig.vs.scl <- sapply(row.names(cpm.test.scaled.group.means), simplify = FALSE, USE.NAMES=TRUE, FUN=function(i){
  t.test(y = cpm.test.original.group.means[i, ], x = cpm.test.scaled.group.means[i, ], paired = TRUE, var.equal = FALSE, conf.level = 0.95, alternative = "greater")
})
pair.wilcox.orig.vs.scl <- sapply(row.names(cpm.test.scaled.group.means), simplify = FALSE, USE.NAMES=TRUE, FUN=function(i){
  wilcox.test(y = cpm.test.original.group.means[i, ], x = cpm.test.scaled.group.means[i, ], paired = TRUE, conf.level = 0.95, alternative = "greater")
})

wilcox.tests.original <- sapply(row.names(cpm.test.norm.filt), simplify = FALSE, USE.NAMES=TRUE, FUN=function(i){
  wilcox.test(x=cpm.test.norm.filt[i, ref.ids], y=cpm.test.norm.filt[i, query.ids], conf.int = TRUE)
})

wilcox.tests.scaled <- sapply(row.names(cpm.test.scaled), simplify = FALSE, USE.NAMES=TRUE, FUN=function(i){
  wilcox.test(x=cpm.test.scaled[i, ref.ids], y=cpm.test.scaled[i, query.ids], conf.int = TRUE)
})
row.mean.orig <- rowMeans(cpm.test.norm.filt)
row.mean.scaled <- rowMeans(cpm.test.scaled)
lab.wilcox.tests.dt.orig <- data.table(do.call(rbind, lapply(wilcox.tests.original, 
                                                             FUN=function(x){
                                                               data.frame(statistic=x$statistic, p.value=x$p.value, alternative=x$alternative, conf.int.l=x$conf.int[1], conf.int.h=x$conf.int[2], estimate=x$estimate)
                                                             })), keep.rownames = TRUE)
lab.wilcox.tests.dt.scaled <- data.table(do.call(rbind, lapply(wilcox.tests.scaled, 
                                                               FUN=function(x){
                                                                 data.frame(statistic=x$statistic, p.value=x$p.value, alternative=x$alternative, conf.int.l=x$conf.int[1], conf.int.h=x$conf.int[2], estimate=x$estimate)
                                                               })), keep.rownames = TRUE)
setkeyv(lab.wilcox.tests.dt.scaled, c("rn", "alternative"))
setkeyv(lab.wilcox.tests.dt.orig, c("rn", "alternative"))
lab.wilcox.tests.dt.MERGED <- merge(lab.wilcox.tests.dt.orig, lab.wilcox.tests.dt.scaled)
setnames(lab.wilcox.tests.dt.MERGED, grep(".x$", colnames(lab.wilcox.tests.dt.MERGED), value = TRUE), sub(".x$", ".ORIGvsREF", grep(".x$", colnames(lab.wilcox.tests.dt.MERGED), value = TRUE)))
setnames(lab.wilcox.tests.dt.MERGED, grep(".y$", colnames(lab.wilcox.tests.dt.MERGED), value = TRUE), sub(".y$", ".SCALEDvsREF", grep(".y$", colnames(lab.wilcox.tests.dt.MERGED), value = TRUE)))

sig.level=0.001
lab.wilcox.tests.dt.MERGED[, `:=`(ORIG.SIGNIF=ifelse(p.value.ORIGvsREF<=sig.level, 1, 0), SCALED.SIGNIF=ifelse(p.value.SCALEDvsREF<=sig.level, 1, 0))]
lab.wilcox.tests.dt.MERGED[, WHICH.SIGNIF:=paste0(ORIG.SIGNIF, SCALED.SIGNIF)]
lab.wilcox.tests.dt.MERGED[, `:=`(negLog.p.value.ORIGvsREF=-1*log10(p.value.ORIGvsREF), negLog.p.value.SCALEDvsREF=-1*log10(p.value.SCALEDvsREF))]
lab.wilcox.tests.dt.MERGED[, `:=`(adj.p.value.ORIGvsREF=p.adjust(p.value.ORIGvsREF, method = "fdr"), adj.p.value.SCALEDvsREF=p.adjust(p.value.SCALEDvsREF, method = "fdr"))]
lab.wilcox.tests.dt.MERGED[, `:=`(negLog.p.adj.ORIGvsREF=-1*log10(adj.p.value.ORIGvsREF), negLog.p.adj.SCALEDvsREF=-1*log10(adj.p.value.SCALEDvsREF))]
ggplot(lab.wilcox.tests.dt.MERGED, aes(x=negLog.p.value.ORIGvsREF, y=negLog.p.value.SCALEDvsREF)) + geom_jitter() 
ggplot(lab.wilcox.tests.dt.MERGED, aes(x=negLog.p.adj.ORIGvsREF, y=negLog.p.adj.SCALEDvsREF)) + geom_jitter() 
lab.wilcox.tests.dt.MERGED[, `:=`(ORIG.ADJ.SIGNIF=ifelse(adj.p.value.ORIGvsREF<=0.01, 1, 0), SCALED.ADJ.SIGNIF=ifelse(adj.p.value.SCALEDvsREF<=0.01, 1, 0))]
lab.wilcox.tests.dt.MERGED[, `:=`(abs.dist.zero.ORIG=abs(estimate.ORIGvsREF), abs.dist.zero.SCALED=abs(estimate.SCALEDvsREF))]
lab.wilcox.tests.dt.MERGED[, dif.abs.dist.orig.vs.scaled:=abs.dist.zero.ORIG-abs.dist.zero.SCALED]
lab.wilcox.tests.dt.MERGED[, mean.cpm.orig:=row.mean.orig[rn]]
lab.wilcox.tests.dt.MERGED[, mean.cpm.scaled:=row.mean.scaled[rn]]
setkey(lab.wilcox.tests.dt.MERGED, rn)


cpm.scaled.df <- data.frame(comparison="CORRECTED", mean.abs.log.dif=rowMeans(cpm.test.scaled.group.means), sd.abs.log.dif=rowSds(cpm.test.scaled.group.means), pval=wilcox.pvals[row.names(cpm.test.scaled.group.means)])
cpm.original.df <- data.frame(comparison="UNCORRECTED", mean.abs.log.dif=rowMeans(cpm.test.original.group.means), sd.abs.log.dif=rowSds(cpm.test.original.group.means), pval=wilcox.pvals[row.names(cpm.test.original.group.means)])

cpm.both <- data.table(keep.rownames = TRUE, rbind(cpm.original.df, cpm.scaled.df))
cpm.both[, `:=`(group.id=paste0(x.ids[1], ".V.", x.ids[2]), from.method=x.ids[1], to.method=x.ids[2])]
cpm.both[, `:=`(fisher.p.val=p.combined, alternative="greater")]
return(cpm.both)
}); MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER <- do.call(rbind, MANN.WHITNEY.SHIFTS.PLASMA.FISHER); MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER[, comparison:=factor(comparison, levels=c("UNCORRECTED", "CORRECTED"))]

MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS.PLASMA.FISHER.STATS <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER[comparison=="UNCORRECTED", .(p.value=min(fisher.p.val)), by=.(group.id, from.method, to.method, alternative)]
setorder(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS.PLASMA.FISHER.STATS, p.value)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER[, group.id:=factor(group.id, levels=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS.PLASMA.FISHER.STATS$group.id)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS.PLASMA.FISHER.STATS[, group.id:=factor(group.id, levels=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS.PLASMA.FISHER.STATS$group.id)]

ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER, aes(x=mean.abs.log.dif, color=comparison)) + 
  geom_density(size=1.5) + 
  geom_text(data=MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.STATS.PLASMA.FISHER.STATS, aes(label=paste0("p = ", format(p.value, digits=3)), x=8.5, y=0.6),inherit.aes = FALSE, check_overlap=TRUE) + facet_wrap(~group.id) + theme_bw() + scale_x_continuous(breaks=seq(0, 16, 2), expand=c(0,0), limits=c(-0.1,14)) + scale_y_continuous(breaks=seq(0, 1, 0.1), expand=c(0,0))
ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER[comparison=="CORRECTED"], aes(y=-1*log10(pval), x=group.id, fill=from.method)) + geom_violin()

##


# 1/19/2017--Check the rank-ordering too. D. Erle suggested looking at top N miRs before and after downbiasing
downbias.comparison.groups <- grep("4N", grep("DOWNBIAS", unique(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset$correction), value=TRUE), invert=TRUE, value=TRUE)
DOWNBIAS.BY.LAB <- sapply(downbias.comparison.groups, USE.NAMES=TRUE, simplify = FALSE, FUN=function(x.group){
  
  correct.dt <- data.frame(subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, correction==x.group), row.names=1)
  x.ids <- strsplit2(x.group, split=".V.")
  print(x.group)
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple%in%x.ids | ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.detail%in%x.ids | ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$GROUND.TRUTH.SAMPLE%in%x.ids ]
  keep.mirs <- (rowSums(dge.test$counts>=5)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  query.ids <- grep(x.ids[1], colnames(dge.test.filt), value=TRUE, ignore.case=TRUE)
  ref.ids <- grep(x.ids[1], invert = TRUE, colnames(dge.test.filt), value=TRUE, ignore.case=TRUE)
  
  dge.query <- dge.test.filt[, query.ids]
  #design.query <- model.matrix(~0+lab.libMethod, dge.query$samples)
  dge.query <- calcNormFactors(dge.query, method="RLE")
  #v.query <- voom(dge.query, normalize.method="quantile")
  
  dge.ref <- dge.test.filt[, ref.ids]
  #design.ref <- model.matrix(~0+lab.libMethod, dge.ref$samples)
  dge.ref <- calcNormFactors(dge.ref, method="RLE")
  #v.ref <- voom(dge.ref, normalize.method="quantile")
  
  merged.miR.ids <- row.names(dge.test.filt)[row.names(dge.test.filt)%in%row.names(correct.dt)]
  #cpm.test.norm.filt.md1 <- v.query$E[merged.miR.ids,]
  cpm.test.norm.filt.md1 <- cpm.DGEList(dge.query, log = TRUE, normalized.lib.sizes = TRUE)[merged.miR.ids,]
  #cpm.test.norm.filt.ref <- v.ref$E[merged.miR.ids,]
  cpm.test.norm.filt.ref <- cpm.DGEList(dge.ref, log = TRUE, normalized.lib.sizes = TRUE)[merged.miR.ids,]
  tt.comp.scl.factor.filt <- correct.dt[merged.miR.ids,]
  
  cpm.test.norm.filt.md.1.scaled <-
    apply(
      cpm.test.norm.filt.md1,
      2,
      FUN = function(x) {
        x.log.scl <- x - tt.comp.scl.factor.filt$logFC
      })
  # combine for re-centering
  cpm.scaled.and.ref <- cbind(cpm.test.norm.filt.ref, cpm.test.norm.filt.md.1.scaled)[, colnames(dge.test.filt)]
  cpm.test.scaled <- voom(2^cpm.scaled.and.ref, design.test, normalize.method="quantile")$E
  mean.cpm.test.scaled.ref <- rowMeans(cpm.test.scaled[, ref.ids])
  query.cpm.test.scaled <- cpm.test.scaled[, query.ids]
  abs.dif.scaled.vs.ref <- abs(query.cpm.test.scaled-mean.cpm.test.scaled.ref)
  
  # re-do the normalization for the original data to "quantile" so that it's normalized the same way as the scaled data.
  cpm.unscaled.and.ref <- cbind(cpm.test.norm.filt.ref, cpm.test.norm.filt.md1)[, colnames(dge.test.filt)]
  v.test.norm <- voom(2^cpm.unscaled.and.ref, design.test, normalize.method="quantile")
  cpm.test.norm <- v.test.norm$E
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  query.cpm.test.norm.filt <- cpm.test.norm.filt[, query.ids]
  mean.cpm.test.norm.filt.ref <- rowMeans(cpm.test.norm.filt[, ref.ids])
  abs.dif.unscaled.vs.ref <- abs(query.cpm.test.norm.filt-mean.cpm.test.norm.filt.ref)
  
  # Now get lab-specific differences
  lab.ids <- unique(dge.query$samples$Lab)
  
  lab.wilcox.tests <- sapply(lab.ids, simplify=FALSE, USE.NAMES = TRUE, FUN=function(x){
    print(paste0("\tWILCOX.TEST: ", x))
    lab.query.ids <- query.ids[(strsplit2(query.ids, split=".", fixed=TRUE)[,2])==x]
    
    lab.cpm.test.scaled.group.means <- data.frame(mean.query.cpm.lab = rowMeans(query.cpm.test.scaled[, lab.query.ids]),
                                                  mean.ref.cpm.all.labs = mean.cpm.test.scaled.ref,  
                                                  abs.log.dif = rowMeans(abs.dif.scaled.vs.ref[, lab.query.ids]),
                                                  comparison = "CORRECTED",
                                                  query.Lab = x,
                                                  row.names=row.names(query.cpm.test.scaled))
    
    lab.cpm.test.original.group.means <- data.frame(mean.query.cpm.lab = rowMeans(query.cpm.test.norm.filt[, lab.query.ids]),
                                                    mean.ref.cpm.all.labs = mean.cpm.test.norm.filt.ref,  
                                                    abs.log.dif = rowMeans(abs.dif.unscaled.vs.ref[, lab.query.ids]),
                                                    comparison = "UNCORRECTED",
                                                    query.Lab = x,
                                                    row.names=row.names(query.cpm.test.scaled))
    lab.cpm.both <- data.table(keep.rownames = TRUE, rbind(lab.cpm.test.scaled.group.means, lab.cpm.test.original.group.means))
    lab.cpm.both[, `:=`(group.id=paste0(x.ids[1], ".V.", x.ids[2]), from.method=x.ids[1], to.method=x.ids[2])]
    
    lab.w.test.out <- wilcox.test(x = lab.cpm.test.original.group.means$abs.log.dif, y=lab.cpm.test.scaled.group.means$abs.log.dif, alternative = c("greater"),  conf.int = TRUE)
    lab.cpm.both[, `:=`(W=lab.w.test.out$statistic, p.value=lab.w.test.out$p.value, alternative="greater", estimate=lab.w.test.out$estimate)]
    return(list(full.table.stats=lab.cpm.both, wilcox.test.only=as.list(lab.w.test.out)))
  })
  lab.wilcox.tests.dt <- do.call(rbind, lapply(lab.wilcox.tests, FUN=function(x){x[[1]]}))
  wilcox.tests <- lapply(lab.wilcox.tests, FUN=function(x){x[[2]]})
  
  wilcox.pvals <- unlist(lapply(wilcox.tests, FUN=function(x){x$p.value}))
  df <- 2*length(wilcox.pvals)
  p.combined <- pchisq(-2*sum(log(wilcox.pvals)), df, lower.tail=FALSE)
  
  lab.wilcox.tests.dt[, `:=`(fisher.p.val=p.combined, fisher.df=df)]
  return(lab.wilcox.tests.dt)
});

MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS <- do.call(rbind, DOWNBIAS.BY.LAB)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.UNCORRECTED <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[comparison=="UNCORRECTED", .(mean.query.cpm.lab.UNCORRECTED=mean(mean.query.cpm.lab), mean.ref.cpm.all.labs.UNCORRECTED=mean(mean.ref.cpm.all.labs)), by=.(rn, query.Lab, group.id, from.method)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.UNCORRECTED[, rn:=sub("1$", "", rn)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.CORRECTED <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[comparison=="CORRECTED", .(mean.query.cpm.lab.CORRECTED=mean(mean.query.cpm.lab), mean.ref.cpm.all.labs.CORRECTED=mean(mean.ref.cpm.all.labs)), by=.(rn, query.Lab, group.id, from.method)]
setkeyv(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.UNCORRECTED, c("rn", "query.Lab", "group.id", "from.method"))
setkeyv(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.CORRECTED, c("rn", "query.Lab", "group.id", "from.method"))
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY <- merge(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.UNCORRECTED, MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.CORRECTED)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt <- melt.data.table(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY, id.vars =  c("rn", "query.Lab", "group.id", "from.method"), variable.name = "CPM.MEASURE", value.name = "CPM")
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt[, rank:=rank(-CPM, ties.method = "min", na.last = TRUE), by=.(query.Lab, group.id, CPM.MEASURE)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt[, min.rank:=min(rank), by=.(rn, from.method)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt[, miR.id.short:=tstrsplit(rn, split=":", fixed=TRUE)[1]]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt[, CPM.MEASURE.2:=ifelse(CPM.MEASURE=="mean.query.cpm.lab.UNCORRECTED", paste0(from.method, ".UNCORRECTED"), ifelse(CPM.MEASURE=="mean.query.cpm.lab.CORRECTED", paste0(from.method, ".CORRECTED"), "REFERENCE")), by=1:nrow(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.mean <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt[, .(rank=min(rank)), by=.(rn, query.Lab, group.id, from.method, miR.id.short, CPM.MEASURE.2)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.mean[, min.rank:=min(rank), by=.(rn, from.method)]


MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.truseq <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.mean[min.rank<=10 & from.method=="TruSeq"]
setorderv(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.truseq, c("CPM.MEASURE.2", "rank"))
mir.id.order <- unique(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.truseq$miR.id.short)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.truseq[, miR.id.short:=factor(miR.id.short, mir.id.order)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.truseq[, CPM.MEASURE.2:=factor(CPM.MEASURE.2, c("REFERENCE", "TruSeq.UNCORRECTED", "TruSeq.CORRECTED"))]
ids.ref.lt.10 <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.truseq[rank<=10 & CPM.MEASURE.2=="REFERENCE"]$miR.id.short

ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.truseq[miR.id.short%in%ids.ref.lt.10], aes(x=miR.id.short, pos=CPM.MEASURE.2, y=rank, fill=CPM.MEASURE.2)) + geom_boxplot() + theme_bw() + scale_y_continuous(breaks = seq(0,200,20)) + theme(
  axis.text=element_text(size=10, color="black"),
  axis.text.x = element_text(hjust=1, angle = 50),
  axis.title=element_text(size=10),
  plot.title=element_text(size=10),
  panel.grid.major = element_line(color=NA),
  panel.border = element_rect(size=1.5, color="black"),
  axis.ticks = element_line(color="black", size=1.5),
  strip.text = element_text(size=10, color="black"),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.text=element_text(size=10),
  legend.title=element_text(size=10),
  legend.key.size=unit(1.2, "lines")
) 
ggsave("UNDERBIASED_MIRS_TRUSEQ_DOWNBIAS.pdf")

ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.SUMMARY.melt.truseq[!miR.id.short%in%ids.ref.lt.10], aes(x=miR.id.short, pos=CPM.MEASURE.2, y=rank, fill=CPM.MEASURE.2)) + geom_boxplot() + theme_bw() + scale_y_continuous(breaks = seq(0,200,20)) + theme(
  axis.text=element_text(size=10, color="black"),
  axis.text.x = element_text(hjust=1, angle = 50),
  axis.title=element_text(size=10),
  plot.title=element_text(size=10),
  panel.grid.major = element_line(color=NA),
  panel.border = element_rect(size=1.5, color="black"),
  axis.ticks = element_line(color="black", size=1.5),
  strip.text = element_text(size=10, color="black"),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.text=element_text(size=10),
  legend.title=element_text(size=10),
  legend.key.size=unit(1.2, "lines")
) 
ggsave("OVERBIASED_MIRS_TRUSEQ_DOWNBIAS.pdf")


# add mann-whitney results to the original scaling factors table
setkey(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, correction)
write.table(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, "SUPPLEMENTARY_TABLE_CORRECTION_FACTORS.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

#
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS <- do.call(rbind, DOWNBIAS.BY.LAB)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[, rn:=sub("1$", "", rn)]
setorder(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS, rn)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[, comparison:=factor(comparison, levels=c("UNCORRECTED", "CORRECTED"))]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[, `:=`(query.rank.in.sample=rank(-mean.query.cpm.lab, ties.method = "ave"),
                                                                            ref.rank.in.sample=rank(-mean.ref.cpm.all.labs, ties.method = "ave")),
                                                                     by=.(query.Lab, comparison, group.id, from.method, to.method)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[, ref.rank.in.sample.uncorrected:=min(ref.rank.in.sample), by=.(rn, group.id)]

MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.ref.top10 <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[ref.rank.in.sample.uncorrected<=10]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.ref.top10[, miR.id.short:=tstrsplit(rn, split=":", fixed=TRUE)[1]]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.ref.top10[, miR.id.short.rank:=paste0(miR.id.short, "(", ref.rank.in.sample.uncorrected, ")")]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.ref.top10.truseq <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.ref.top10[from.method=="TruSeq"]
setorder(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.ref.top10.truseq, ref.rank.in.sample.uncorrected)
mir.id.sort.truseq <- unique(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.ref.top10.truseq$miR.id.short)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.ref.top10.truseq[, miR.id.short:=factor(miR.id.short, levels=mir.id.sort.truseq)]

ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.ref.top10.truseq, aes(x=miR.id.short, pos=comparison, y=query.rank.in.sample, color=comparison)) + geom_jitter(height=0, width=0.1) + theme_bw() + theme(
  axis.text=element_text(size=10, color="black"),
  axis.text.x = element_text(hjust=1, angle = 50),
  axis.title=element_text(size=10),
  plot.title=element_text(size=10),
  panel.grid.major = element_line(color=NA),
  panel.border = element_rect(size=1.5, color="black"),
  axis.ticks = element_line(color="black", size=1.5),
  strip.text = element_text(size=10, color="black"),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.text=element_text(size=10),
  legend.title=element_text(size=10),
  legend.key.size=unit(1.2, "lines")
) 

MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[, mean.query.cpm:=mean(mean.query.cpm.lab), by=.(group.id, comparison, rn)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[, query.rank.in.sample:=rank(-mean.query.cpm, ties.method = "ave"), by=.(group.id, comparison, query.Lab)]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS[from.method=="TruSeq"]

MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq.top10.ids <- unique(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq[query.rank.in.sample<=10 & comparison=="UNCORRECTED"]$rn)

MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq.query.top10 <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq[rn%in%MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq.top10.ids]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq.query.top10[, miR.id.short:=tstrsplit(rn, split=":", fixed=TRUE)[1]]
setorder(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq.query.top10, query.rank.in.sample)
mir.id.sort.truseq <- unique(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq.query.top10$miR.id.short)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq.query.top10[, miR.id.short:=factor(miR.id.short, levels=mir.id.sort.truseq)]

ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.DOWNBIAS.truseq.query.top10, aes(x=miR.id.short, pos=comparison, y=query.rank.in.sample, color=comparison)) + geom_jitter(height=0, width=0.1) + theme_bw() + theme(
  axis.text=element_text(size=10, color="black"),
  axis.text.x = element_text(hjust=1, angle = 50),
  axis.title=element_text(size=10),
  plot.title=element_text(size=10),
  panel.grid.major = element_line(color=NA),
  panel.border = element_rect(size=1.5, color="black"),
  axis.ticks = element_line(color="black", size=1.5),
  strip.text = element_text(size=10, color="black"),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.text=element_text(size=10),
  legend.title=element_text(size=10),
  legend.key.size=unit(1.2, "lines")
) 







######## NOW ALSO SHOW CORRECTIONS FOR EQUIMOLAR POOL WITH EQUIMOLAR SEQS
# Now this is where we need to diverge to include the other methods. 
# The thought is that scaling factors are going to convert between two different methods, so all of our filtering, normalization, centering and testing should 
# be done without the potential influence of the other library prep types.

lib.comparison.groups <- list(TruSeq.V.4N=c("TruSeq", "4N"), 
                              TruSeq.V.NEB=c("TruSeq", "NEBNext"),
                              TruSeq.V.CleanTag=c("TruSeq", "CleanTag"),
                              NEB.V.CleanTag=c("NEBNext", "CleanTag"),
                              NEB.V.4N=c("NEBNext", "4N"),
                              CleanTag.V.4N=c("CleanTag", "4N")
)


ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.EQUIMOLAR <- do.call(rbind, sapply(names(lib.comparison.groups), USE.NAMES=TRUE, simplify = FALSE, FUN=function(x.id.list){
  x.ids <- unlist(lib.comparison.groups[x.id.list], recursive = TRUE, use.names = FALSE)
  dge.comp <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[, ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.simple%in%x.ids]
  design.comp <- model.matrix(~0+lab.libMethod, dge.comp$samples)
  colnames(design.comp) <- unique(dge.comp$samples$lab.libMethod)
  lib.methods.comp <- strsplit2(strsplit2(colnames(design.comp), split=".", fixed=TRUE)[,1], split="_")[,1]
  comp.contrasts.fodder <- sapply(unique(lib.methods.comp), USE.NAMES=FALSE, FUN=function(lib.m){
    mat <- design.comp
    y.sum <- rowSums(apply(mat, 1, FUN=function(x){
      x*ifelse(lib.methods.comp==lib.m, 1, 0)
    }))
    y.sum/sum(y.sum)
  })
  colnames(comp.contrasts.fodder) <- unique(lib.methods.comp)
  out.contrast <- as.numeric(comp.contrasts.fodder[, x.ids[1]]-comp.contrasts.fodder[, x.ids[2]])
  # Now filter miRs
  dge.comp.filt <- dge.comp[rowSums(dge.comp$counts>0)==dim(dge.comp)[2],]
  v.comp <- voom(dge.comp.filt, design.comp, normalize.method = "scale")
  fit.comp <- lmFit(v.comp, design.comp)
  fit.comp <- eBayes(fit.comp, robust = TRUE, trend = TRUE)
  fit2.comp <- contrasts.fit(fit.comp, out.contrast)
  fit2.comp <- eBayes(fit2.comp, trend = TRUE, robust = TRUE)
  tt.comp <- data.frame(topTable(fit2.comp, n = Inf, confint=0.95, sort.by = "none"))
  tt.comp.dt <- data.table(data.frame(miR.id.full=row.names(tt.comp), tt.comp))
  tt.comp.dt[, comparison:=x.id.list]
  tt.comp.scl.factor <- tt.comp[, c("logFC", "CI.L", "CI.R")]
  
  #Finally, merge back with plasma data.
  dge.test <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[, ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.simple%in%x.ids]
  keep.mirs <- (rowSums(dge.test$counts>=1)==dim(dge.test)[2]) 
  design.test <- model.matrix(~0+lab.libMethod, dge.test$samples)
  dge.test.filt <- dge.test[keep.mirs,]
  v.test.norm <- voom(dge.test.filt, design.test, normalize.method="scale")
  cpm.test.norm <- v.test.norm$E
  
  merged.miR.ids <- merge(cpm.test.norm, tt.comp.scl.factor, by=0)$Row.names
  cpm.test.norm.filt <- cpm.test.norm[merged.miR.ids, ]
  cpm.test.norm.filt.md1 <- cpm.test.norm.filt[, grep(x.ids[1], colnames(cpm.test.norm), value=TRUE, ignore.case=TRUE)]
  tt.comp.scl.factor.filt <- tt.comp.scl.factor[merged.miR.ids,]
  
  # set a floor min. Logic here is that if we had a miR that had a count of 1 and the scaling factor 
  # indicates a drop of -5.0 fold, we do not want to make the new log2 count log2(1) - 5.0, since our scale will be way off
  # for plotting. So, if the value is below the min, set it to the min value. In other words. Limit the amount 
  # of downward scaling to the minimum value of the sample.
  #cpm.test.norm.filt.scaled <- apply(cpm.test.norm.filt, 2, FUN=function(x){
  #  x.floor <- min(x[x>0])
  #  x.log.scl <- x-tt.comp.scl.factor.filt$logFC
  #  return(ifelse(x.log.scl<x.floor, x.floor, x.log.scl))
  #})
  
  # this does it without a floor
  # Decided to use median scaling to get the scaling factors
  # Then after scaling, use quantile normalization to match up the distributions
  # Thinking that this might be a more appropriate way to handle the fact that
  # some of the scaling at the low end will lead to values less than the minimum. 
  cpm.test.norm.filt.md.1.scaled <- apply(cpm.test.norm.filt.md1, 2, FUN=function(x){
    #x.floor <- min(x[x>0])
    x.log.scl <- x-tt.comp.scl.factor.filt$logFC
    #return(ifelse(x.log.scl<x.floor, x.floor, x.log.scl))
    return(x.log.scl)
  })
  colnames(cpm.test.norm.filt) <- paste0(colnames(cpm.test.norm.filt), ".orig")
  colnames(cpm.test.norm.filt.md.1.scaled) <- paste0(colnames(cpm.test.norm.filt.md.1.scaled), ".transf")
  cpm.test.orig.and.scaled <- cbind(cpm.test.norm.filt, cpm.test.norm.filt.md.1.scaled)
  cpm.test.orig.and.scaled.rescaled <- normalizeBetweenArrays(cpm.test.orig.and.scaled, method = "quantile")
  ids.reference <- grep(x.ids[1], colnames(cpm.test.orig.and.scaled.rescaled), value=TRUE, ignore.case=TRUE, invert = TRUE)
  ids.test.orig <- grep(x.ids[1], colnames(cpm.test.norm.filt), value=TRUE, ignore.case=TRUE)
  ids.test.scaled <- colnames(cpm.test.norm.filt.md.1.scaled)
  # Return this to get a matrix of transformed and original values
  # return(cpm.test.orig.and.scaled.rescaled)
  # Now get an average of the values we did not transform
  cpm.rescaled.test.orig <- cpm.test.orig.and.scaled.rescaled[, ids.test.orig]
  cpm.rescaled.test.scaled <- cpm.test.orig.and.scaled.rescaled[, ids.test.scaled]
  cpm.rescaled.ref.orig <- cpm.test.orig.and.scaled.rescaled[, ids.reference]
  cpm.mean.reference <- rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.reference])
  
  # Check scale norm too. Didnt look quite as good, so removed.
  cpm.test.orig.and.scaled.rescaled.med <- normalizeBetweenArrays(cpm.test.orig.and.scaled, method = "scale")
  cpm.rescaled.test.orig.med <- cpm.test.orig.and.scaled.rescaled.med[, ids.test.orig]
  cpm.rescaled.test.scaled.med <- cpm.test.orig.and.scaled.rescaled.med[, ids.test.scaled]
  cpm.rescaled.ref.orig.med <- cpm.test.orig.and.scaled.rescaled.med[, ids.reference]
  cpm.mean.reference.med <- rowMeans(cpm.test.orig.and.scaled.rescaled.med[, ids.reference])
  
  
  cpm.test.orig.and.scaled.rescaled.mean <- data.table(data.frame(miRNA.full.ID=row.names(cpm.test.orig.and.scaled.rescaled), 
                                                                  scaling.factor.logFC=tt.comp.scl.factor.filt$logFC,
                                                                  cpm.mean.reference.quantile = cpm.mean.reference,      
                                                                  cpm.mean.test.orig.quantile = rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.test.orig]),
                                                                  cpm.mean.test.scaled.quantile = rowMeans(cpm.test.orig.and.scaled.rescaled[, ids.test.scaled]),
                                                                  abs.dif.orig.quantile = rowMeans(apply(cpm.rescaled.test.orig, 2, FUN=function(x){ abs(x-cpm.mean.reference)})),
                                                                  abs.dif.transf.quantile = rowMeans(apply(cpm.rescaled.test.scaled, 2, FUN=function(x){ abs(x-cpm.mean.reference)}))
  ))
  cpm.test.orig.and.scaled.rescaled.mean[, comparison:=x.id.list]
  return(cpm.test.orig.and.scaled.rescaled.mean)
})); #ggplot(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.LIST[[3]]) + geom_density(aes(color="VS.ORIG", x=abs.dif.orig)) + geom_density(aes(color="VS.TRANSF", x=abs.dif.transf))
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR <- melt(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.EQUIMOLAR, measure.vars=c("abs.dif.orig.quantile", "abs.dif.transf.quantile"), variable.name = "ORIG.or.TRANSF.and.SCL", value.name = "mean.abs.diff")
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR[, SCALE.METHOD:=tstrsplit(ORIG.or.TRANSF.and.SCL, split=".", fixed=TRUE)[4]]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR[, ORIG.or.TRANSF:=sub(paste0(".", SCALE.METHOD), "", ORIG.or.TRANSF.and.SCL)]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR[, `:=`(SCALE.METHOD=ifelse(SCALE.METHOD=="quantile", "QUANTILE", "MEDIAN"), ORIG.or.TRANSF=ifelse(ORIG.or.TRANSF=="abs.dif.orig", "ORIGINAL", "TRANSFORMED"))]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR.quantile <- subset(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR, SCALE.METHOD=="QUANTILE")
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR.quantile[, pool:="SynthEQ"]
# Filter the equimolar results to include only the miRs we looked at in the plasma pool
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR.quantile.filt <- subset(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR.quantile, miRNA.full.ID%in%ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile$miRNA.full.ID)
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile[, pool:="PlasmaPool"]
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile.EQandPlasm <- rbind(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile, ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR.quantile.filt)
ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile.EQandPlasm[, format.alpha:=ifelse(ORIG.or.TRANSF=="ORIGINAL", 1, 0.5)]

# First just plot the equimolar pool. Include all miRs; not just the ones in the plasma pool
ggplot(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.EQUIMOLAR.quantile, aes(x=mean.abs.diff, color=ORIG.or.TRANSF)) + geom_density(size=1.5) + facet_wrap(~comparison) + theme_bw() + 
  scale_y_continuous(expand=c(0.01, 0.0), breaks=seq(0.0, 10.0, 0.5), limits=c(0,3.5)) +
  scale_x_continuous(expand=c(0.01,0), limits=c(0.0, 9.0), breaks=seq(0,9, 1.0)) + 
  theme(
    axis.text=element_text(size=14, face="bold", color="black"),
    axis.title=element_text(face="bold", size=16),
    plot.title=element_text(face="bold", size=18),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=16, face="bold", color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=14, face="bold"),
    legend.title=element_text(size=15, face="bold"),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(title="SynthEQ",
       y="probability density",
       x="|Difference VS Reference|(log2)")



# Plot the equimolar and plasma overlaid
ggplot(ALL.COMPARISON.SCALING.FACTOR.SCALED.MATR.DT.melt.abs.dif.quantile.EQandPlasm, aes(x=mean.abs.diff, color=pool)) + geom_density(size=1.5, aes(linetype=ORIG.or.TRANSF)) + facet_wrap(~comparison) + theme_bw() + 
  scale_linetype_manual(values=c("dotted", "solid")) +
  scale_y_continuous(expand=c(0.0, 0.0), breaks=seq(0.0, 10.0, 0.5), limits=c(0,3.5)) + 
  scale_x_continuous(expand=c(0.01,0), limits=c(0.0, 13.5), breaks=seq(0,15, 1.0)) + 
  theme(
    axis.text=element_text(size=14, face="bold", color="black"),
    axis.title=element_text(face="bold", size=16),
    plot.title=element_text(face="bold", size=18),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    strip.text = element_text(size=16, face="bold", color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=14, face="bold"),
    legend.title=element_text(size=15, face="bold"),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(y="probability density",
       x="|Difference VS Reference|(log2)")







#######################################

### ABOVE: Used equimolar pool scaling factors to correct the plasma pool. Now get plasma pool scaling factors.
# Use the filtered set where we required >=1 count in all samples and filtered low count libs, but before we required it to have a corresponding match in the eq pool
# This is what was run:
# ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR.FILT.ZEROS <- ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR[rowSums(ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR$counts>=1)==dim(ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR)[2], ]

# design.plasma <- model.matrix(~0 + lab.libMethod, ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR.FILT.ZEROS$samples)
# colnames(design.plasma) <- unique(ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR.FILT.ZEROS$samples$lab.libMethod)

# NEW--now do the weighted voom like we did with the equimolar pool
v.plasma <- voom(ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR.FILT.ZEROS, design.plasma, normalize.method="scale", plot=TRUE)
lib.methods.design.plasma <- strsplit2(colnames(design.plasma), split=".", fixed=TRUE)[,1]

my.contrasts.fodder.plasma <- sapply(unique(lib.methods.design.plasma), USE.NAMES=FALSE, FUN=function(lib.m){
  mat <- design.plasma
  y.sum <- rowSums(apply(mat, 1, FUN=function(x){
    x*ifelse(lib.methods.design.plasma==lib.m, 1, 0)
  }))
  y.sum/sum(y.sum)
})

all.4n.sum <- as.numeric(rowSums(apply(design.plasma, 1, FUN=function(x){
  x*ifelse(strsplit2(lib.methods.design.plasma, split="_")[,1]=="4N", 1, 0)
})))
my.contrasts.fodder.plasma <- cbind(my.contrasts.fodder.plasma, all.4n.sum/sum(all.4n.sum))
colnames(my.contrasts.fodder) <- c(unique(lib.methods.design), "4N_ABCD")
my.contrasts.plasma <- cbind(truseq.to.4NA=my.contrasts.fodder[, "TruSeq"]-my.contrasts.fodder[,"4N_A"],
                             truseq.to.4NB=my.contrasts.fodder[, "TruSeq"]-my.contrasts.fodder[,"4N_B"],
                             truseq.to.4NC=my.contrasts.fodder[, "TruSeq"]-my.contrasts.fodder[,"4N_C"],
                             truseq.to.4ND=my.contrasts.fodder[, "TruSeq"]-my.contrasts.fodder[,"4N_D"],
                             truseq.to.avg.4N=my.contrasts.fodder[,"TruSeq"]-my.contrasts.fodder[,"4N_ABCD"])

fit.plasm <- lmFit(v.plasma, design.plasma)
fit.plasm <- eBayes(fit.plasm, robust = TRUE, trend = TRUE)
fit2.plasm <- contrasts.fit(fit.plasm, my.contrasts.plasma)
fit2.plasm <- eBayes(fit2.plasm, trend = TRUE, robust = TRUE)
tt.truseq.vs.4N.ABCD.plasm <- topTable(fit2.plasm, 5, n = Inf, confint=0.95, sort.by = "none")
tt.truseq.vs.4N.ALL.plasm <- topTable(fit2.plasm, 1:4, n = Inf)

v.noqual.plasm <- voom(ENDOG.MIR.COUNTS.MATURE.FULL.DGE.FILT.PLASMA.ONLY.FOR.SCL.FACTOR.FILT.ZEROS, design.plasma, normalize.method="scale", plot=TRUE)
lcpm.plasm <- v.noqual.plasm$E
lcpm.truseq.plasm <- lcpm.plasm[, grep("truseq", colnames(lcpm.plasm), ignore.case=TRUE)]
lcpm.4N.plasm <- lcpm.plasm[, grep("4N", colnames(lcpm.plasm), ignore.case=TRUE)]
lcpm.4N.merge.tt.4NALL.plasm <- merge(lcpm.4N.plasm, tt.truseq.vs.4N.ABCD.plasm[, c("logFC", "CI.L", "CI.R")], by=0)
mean.lcpm.truseq.plasm <- rowMeans(lcpm.truseq.plasm)
mean.lcpm.4N.plasm <- rowMeans(lcpm.4N.plasm)
lcpm.4N.transform.plasm <- apply(lcpm.4N.plasm, 2, FUN=function(x){
  x+lcpm.4N.merge.tt.4NALL.plasm$logFC
})
mean.lcpm.4N.transform.plasm <- rowMeans(lcpm.4N.transform.plasm)
merge.mean.plasm <- cbind(mean.lcpm.truseq=mean.lcpm.truseq.plasm, mean.lcpm.4N=mean.lcpm.4N.plasm, mean.lcpm.4N.transform=mean.lcpm.4N.transform.plasm )


colnames(lcpm.4N.transform.plasm) <- paste0(colnames(lcpm.4N.plasm), ".transf")
colnames(lcpm.plasm) <- paste0(colnames(lcpm.plasm), ".orig")
lcpm.merge.transf.plasm <- cbind(lcpm.plasm, lcpm.4N.transform.plasm)
name.split.plasm <- strsplit2(colnames(lcpm.merge.transf.plasm), split=".", fixed=TRUE)
sample.info.plasm <- data.frame(row.names=colnames(lcpm.merge.transf.plasm), lib.method.simple=paste(strsplit2(name.split.plasm[,1], split="_")[,1], name.split.plasm[,5], sep="."))
pheatmap(lcpm.merge.transf.plasm, show_rownames = FALSE, annotation_col = sample.info.plasm)

lcpm.truseq.merge.tt.4NALL.plasm <- merge(lcpm.truseq.plasm, tt.truseq.vs.4N.ABCD.plasm[, c("logFC", "CI.L", "CI.R")], by=0)
lcpm.truseq.transform.plasm <- apply(lcpm.truseq.plasm, 2, FUN=function(x){
  x-lcpm.truseq.merge.tt.4NALL.plasm$logFC
})
colnames(lcpm.truseq.transform.plasm) <- paste0(colnames(lcpm.truseq.transform.plasm), ".transf")
lcpm.merge.transf.plasm <- cbind(lcpm.plasm, lcpm.truseq.transform.plasm)
name.split.plasm <- strsplit2(colnames(lcpm.merge.transf.plasm), split=".", fixed=TRUE)
sample.info.plasm <- data.frame(row.names=colnames(lcpm.merge.transf.plasm), lib.method.simple=paste(strsplit2(name.split.plasm[,1], split="_")[,1], name.split.plasm[,5], sep="."))
pheatmap(lcpm.merge.transf.plasm, show_rownames = FALSE, annotation_col = sample.info.plasm)


########################################################################################################
##### GROUND TRUTH SCALING FACTORS 12-7-2016
### Instead of scaling from one library prep method to another, use the equimolar pools to correct to the "expected" value. 



####################################
#### # miRs detected using downsampling 12-7-2016
# Try to get the # of miRs detected at different read depths. Use downsampling to different levels.
library(metaseqR)
# start with ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL
downsample.to <- round(10^seq(from = 4, to = 6.5, by = 0.5), digits=0)
n.miRs <- dim(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL)[1]
ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA <- ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL[, ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL$samples$pool=="PlasmaPool"]
ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS <- ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA[rowSums(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA$counts)>0,]
# Idea: Downsample to the different levels in downsample.to variable and count miRs detected.
# Method: 
# 1. Remove libraries with total counts < the indicated threshold
# 2. Append a dummy sample where counts sum to the "downsample.to" value. 
# 3. Run the downsampling and drop the dummy column
# 4. Get the number of samples and # of miRs detected after downsampling
downsample.miRs.above.cutoff.matr <- sapply(downsample.to, FUN=function(x){
  # get the count matrix
  dge.COUNTS <- ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS$counts
  # generate dummy counts
  # the c(row.names(dge.COUNTS)) addition makes sure all miRs are present in the dummy counts
  dummy.counts <- as.numeric(table(c(row.names(dge.COUNTS), sample(row.names(dge.COUNTS), size = x-dim(dge.COUNTS)[1], replace = TRUE))))
  samples.below.cutoff <- ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS$samples$lib.size < x
  # replace libraries with counts lt the cutoff with the dummy counts
  dge.COUNTS[, samples.below.cutoff] <- dummy.counts
  # downsample those after rounding
  dge.COUNTS.downsampled <- round(downsample.counts(dge.COUNTS), digits=0)
  n.miRs.detected <- apply(dge.COUNTS.downsampled, 2, FUN=function(x){ sum(ifelse(x>0, 1, 0)) })
  n.miRs.detected.na.filt <- ifelse(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS$samples$lib.size<x, NA, n.miRs.detected)
})

colnames(downsample.miRs.above.cutoff.matr) <- paste0("downampleLevel.", 1:dim(downsample.miRs.above.cutoff.matr)[2])
downsample.miRs.above.cutoff.matr.add.sample.info <- data.table(data.frame(lab.libMethod.pool.replicate=row.names(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS$samples), cbind(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS$samples, downsample.miRs.above.cutoff.matr)))
downsample.miRs.above.cutoff.matr.add.sample.info.melt <- melt.data.table(downsample.miRs.above.cutoff.matr.add.sample.info, id.vars=c("lab.libMethod.pool.replicate", colnames(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS$samples)), variable.name = "downsampleLevel", value.name = "n.miRs.detected")
downsample.miRs.above.cutoff.matr.add.sample.info.melt[, downsampled.to:=downsample.to[as.numeric(tstrsplit(downsampleLevel, split=".", fixed=TRUE)[2])], by=1:nrow(downsample.miRs.above.cutoff.matr.add.sample.info.melt)]

downsample.miRs.above.cutoff.matr.add.sample.info.melt.filt.missing <- subset(downsample.miRs.above.cutoff.matr.add.sample.info.melt, !is.na(n.miRs.detected))

downsample.miRs.above.cutoff.matr.add.sample.info.melt.filt.missing[, `:=`(lab.libMethod=factor(lab.libMethod, levels=sample.levels), lib.method.simple=factor(lib.method.simple, levels=lib.simple.levels))]
ggplot(downsample.miRs.above.cutoff.matr.add.sample.info.melt.filt.missing, aes(x=factor(round(log10(downsampled.to), digits=1)), pos=lib.method.simple, y=n.miRs.detected, fill=lib.method.simple)) + geom_boxplot(position="dodge") + theme_bw() + labs(x="log10 Sequencing Depth (total miRNA-mapping reads)") + theme(
  axis.text=element_text(size=14, face="bold", color="black"),
  axis.title=element_text(face="bold", size=16),
  plot.title=element_text(face="bold", size=18),
  panel.grid.major = element_line(color=NA),
  panel.border = element_rect(size=1.5, color="black"),
  axis.ticks = element_line(color="black", size=1.5),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.text=element_text(size=14, face="bold"),
  legend.title=element_text(size=15, face="bold"),
  legend.key.size=unit(1.2, "lines")
)

# now test consistency of detection
downsample.miRs.above.cutoff.consistency <- sapply(downsample.to, simplify=FALSE, USE.NAMES = TRUE, FUN=function(x){
  # get the count matrix
  dge.COUNTS <- ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS$counts
  # generate dummy counts
  # the c(row.names(dge.COUNTS)) addition makes sure all miRs are present in the dummy counts
  dummy.counts <- as.numeric(table(c(row.names(dge.COUNTS), sample(row.names(dge.COUNTS), size = x-dim(dge.COUNTS)[1], replace = TRUE))))
  samples.below.cutoff <- ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS$samples$lib.size < x
  # replace libraries with counts lt the cutoff with the dummy counts
  dge.COUNTS[, samples.below.cutoff] <- dummy.counts
  # downsample those after rounding
  dge.COUNTS.downsampled <- round(downsample.counts(dge.COUNTS), digits=0)
  dge.COUNTS.downsampled[, samples.below.cutoff] <- -1
  dge.COUNTS.downsampled.detected <- ifelse(dge.COUNTS.downsampled>0, 1, 0)
  dge.COUNTS.downsampled.total <- ifelse(dge.COUNTS.downsampled>=0, 1, 0)
  dge.COUNTS.downsampled.detected.tot <- sumTechReps(dge.COUNTS.downsampled.detected, ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS$samples$lib.method.simple)
  dge.COUNTS.downsampled.total.tot <- sumTechReps(dge.COUNTS.downsampled.total, ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS$samples$lib.method.simple)
  dge.COUNTS.downsampled.detected.perc.tot <- dge.COUNTS.downsampled.detected.tot/dge.COUNTS.downsampled.total.tot
  n.miRs.detected <- apply(dge.COUNTS.downsampled.detected.perc.tot, 2, FUN=function(x) sum(ifelse(x==1, 1, 0)))
  n.samples.remaining <- apply(dge.COUNTS.downsampled.total.tot, 2, max, na.rm=TRUE)
  list(n.miRs.detected=n.miRs.detected, n.samples.remaining=n.samples.remaining)
})
names(downsample.miRs.above.cutoff.consistency) <- paste0("downampleLevel.", 1:dim(downsample.miRs.above.cutoff.matr)[2])
downsample.miRs.above.cutoff.consistency.matr <- do.call(rbind, lapply(downsample.miRs.above.cutoff.consistency, unlist))

downsample.miRs.above.cutoff.consistency.dt <- data.table(data.frame(downsampleLevel=row.names(downsample.miRs.above.cutoff.consistency.matr), downsample.miRs.above.cutoff.consistency.matr, stringsAsFactors = FALSE))
value.cols.n.miRs <- grep("miRs", colnames(downsample.miRs.above.cutoff.consistency.dt), value = TRUE)
value.cols.n.samples <- grep("samples", colnames(downsample.miRs.above.cutoff.consistency.dt), value = TRUE)
downsample.miRs.above.cutoff.consistency.dt.melt.samples <- melt.data.table(downsample.miRs.above.cutoff.consistency.dt, id.vars="downsampleLevel", measure = value.cols.n.samples, variable.name = "measure.lib.method.simple", value.name = c("n.samples.remaining"))
downsample.miRs.above.cutoff.consistency.dt.melt.miRs <- melt.data.table(downsample.miRs.above.cutoff.consistency.dt, id.vars="downsampleLevel", measure = value.cols.n.miRs, variable.name = "measure.lib.method.simple", value.name = c("n.miRs.detected.in.all"))
downsample.miRs.above.cutoff.consistency.dt.melt.samples[, lib.method.simple:=tstrsplit(measure.lib.method.simple, split=".remaining.")[2]]
downsample.miRs.above.cutoff.consistency.dt.melt.miRs[, lib.method.simple:=tstrsplit(measure.lib.method.simple, split=".detected.")[2]]
downsample.miRs.above.cutoff.consistency.dt.melt.samples[, measure.lib.method.simple:=NULL]
downsample.miRs.above.cutoff.consistency.dt.melt.miRs[, measure.lib.method.simple:=NULL]
setkeyv(downsample.miRs.above.cutoff.consistency.dt.melt.samples, c("downsampleLevel", "lib.method.simple"))
setkeyv(downsample.miRs.above.cutoff.consistency.dt.melt.miRs, c("downsampleLevel", "lib.method.simple"))

downsample.miRs.above.cutoff.consistency.dt.melt <- merge(downsample.miRs.above.cutoff.consistency.dt.melt.miRs, downsample.miRs.above.cutoff.consistency.dt.melt.samples)
downsample.miRs.above.cutoff.consistency.dt.melt[is.na(n.miRs.detected.in.all), n.miRs.detected.in.all:=0]
downsample.miRs.above.cutoff.consistency.dt.melt[, downsampled.to:=downsample.to[as.numeric(tstrsplit(downsampleLevel, split=".", fixed=TRUE)[2])], by=1:nrow(downsample.miRs.above.cutoff.consistency.dt.melt)]

downsample.miRs.above.cutoff.consistency.dt.melt[, `:=`(lib.method.simple=factor(lib.method.simple, levels=lib.simple.levels))]
ggplot(downsample.miRs.above.cutoff.consistency.dt.melt,
       aes(x=factor(round(log10(downsampled.to), digits=1)), pos=lib.method.simple, y=n.miRs.detected.in.all, fill=lib.method.simple)) +
  geom_bar(position="dodge", stat="identity", width=0.75) + 
  theme_bw() + 
  scale_y_continuous(breaks=seq(0,1000, 100), limits=c(0,900), expand = c(0, 0)) +
  geom_text(
    aes(
      label=n.samples.remaining,
      y=n.miRs.detected.in.all+25,
      hjust="center",
      size=2,
      fontface="bold"), 
    position = position_dodge(width=0.75)) + 
  labs(x="log10 Sequencing Depth (total miRNA-mapping reads)") + theme(
    axis.text=element_text(size=14, face="bold", color="black"),
    axis.title=element_text(face="bold", size=16),
    plot.title=element_text(face="bold", size=18),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=14, face="bold"),
    legend.title=element_text(size=15, face="bold"),
    legend.key.size=unit(1.2, "lines")
  )

### DOWNSAMPLING BRUTE FORCE.. THE METASEQR IS NOT RANDOM

ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS
test.count <- as.numeric(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS$counts[,1])
test.count.round <- round(test.count, digits=0)
nonzero.index <- which(test.count.round>0)


nonzero.index.rep <- unlist(sapply(nonzero.index, FUN=function(x){ rep.int(x, times=test.count.round[x]) } ))
nonzero.index.subs <- table(sample(nonzero.index.rep, size = 100000, replace = FALSE))
nonzero.index.subs2 <- table(c(1:length(test.count), sample(nonzero.index.rep, size = 100000, replace = FALSE)))-1

ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts <- round(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS$counts, digits=0)
library(plyr)
library(foreach)
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl) 
downsample.to.vec <- round(10^seq(from = 4, to = 6.5, by = 0.5), digits=0)


ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts.EXPANDED <- do.call(rbind, sapply(1:dim(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts)[2], USE.NAMES=FALSE, simplify=FALSE,
                                                                                                                        FUN=function(x.col){
                                                                                                                          sample.id <- colnames(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts)[x.col]
                                                                                                                          x <- ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts[,x.col]
                                                                                                                          idx.nonzero <- which(x>0)
                                                                                                                          idx.nonzero.rep <- unlist(sapply(idx.nonzero, 
                                                                                                                                                           FUN=function(idxs){
                                                                                                                                                             rep.int(idxs, times=x[idxs])
                                                                                                                                                           }
                                                                                                                          )
                                                                                                                          )
                                                                                                                          idx.nonzero.rep.dt <- data.table(data.frame(lab.libMethod.replicate=sample.id, mir.idx=idx.nonzero.rep))
                                                                                                                        }))

ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts.EXPANDED[, lib.size:=.N, by=lab.libMethod.replicate]


ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts.EXPANDED.downsampled <- do.call(rbind, sapply(downsample.to.vec, USE.NAMES=FALSE, simplify=FALSE, FUN=function(downsample.to) ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts.EXPANDED[lib.size>=downsample.to, .(mir.idx.new=c(mir.idx, sample(mir.idx, downsample.to, replace = FALSE))), by=lab.libMethod.replicate][, .(downsampled.count=.N-1, downsample.to=downsample.to), by=.(lab.libMethod.replicate, mir.idx.new)]))


downsample.iter.df <- data.frame(expand.grid(1:100, downsample.to.vec))
colnames(downsample.iter.df) <- c("iter", "downsample.to")
max.mir.idx <- ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts.EXPANDED[, max(mir.idx)]
do.downsampling <- function(){
  ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts.EXPANDED.downsampled <<- do.call(rbind, sapply(1:dim(downsample.iter.df)[1], USE.NAMES=FALSE, simplify=FALSE, FUN=function(x){
    iter=downsample.iter.df[x, 1]
    downsample.to=downsample.iter.df[x, 2]
    print(paste0(x, " of ", dim(downsample.iter.df)[1]))
    ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts.EXPANDED[lib.size>=downsample.to, .(mir.idx=sample(mir.idx, downsample.to, replace = FALSE)), by=lab.libMethod.replicate][, .(downsampled.count=.N, iter=iter, downsample.to=downsample.to), by=.(lab.libMethod.replicate, mir.idx)]
  }))
  
  write.table(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts.EXPANDED.downsampled, "ALL_ITERATIONS_DOWNSAMPLING_LEVELS_PLASMA.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  
}


# This will take ~ half hour to run
# do.downsampling
ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts.EXPANDED.downsampled <- fread("ALL_ITERATIONS_DOWNSAMPLING_LEVELS_PLASMA.txt")
setnames(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts.EXPANDED.downsampled, "downsample.to", "downsampled.to")
all.plasm.miRNA.dt <- data.table(data.frame(mir.idx=1:dim(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS)[1], mir.seq.id=row.names(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS)))

ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts.EXPANDED.downsampled
sample.info.plasma.for.downsampling.dt <- data.table(data.frame(lab.libMethod.replicate=row.names(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS$samples), ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS$samples))

setkey(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts.EXPANDED.downsampled, lab.libMethod.replicate)
setkey(sample.info.plasma.for.downsampling.dt, lab.libMethod.replicate)
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO <- merge(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts.EXPANDED.downsampled, sample.info.plasma.for.downsampling.dt)
setkey(ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO, mir.idx)
setkey(all.plasm.miRNA.dt, mir.idx)
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO[all.plasm.miRNA.dt, mir.seq.id:=mir.seq.id]
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO[, total.libs.for.method.simple:=length(unique(lab.libMethod.replicate)), by=.(lib.method.simple, downsampled.to)]
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY <- ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO[downsampled.count>0, .(n.libs.with.mir.detected=.N), by=.(mir.seq.id, lib.method.simple, iter, downsampled.to, total.libs.for.method.simple)]


ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.DETECTED.ALL.BY.ITER <- ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY[n.libs.with.mir.detected>=total.libs.for.method.simple, .(n.miRs.all.detected.by.iter=.N), by=.(lib.method.simple, iter, downsampled.to, total.libs.for.method.simple)]

downsample.to.vec <- downsample.to.vec
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.DETECTED.ALL.BY.ITER[, lib.method.simple:=factor(lib.method.simple, lib.simple.levels)]

ggplot(ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.DETECTED.ALL.BY.ITER[downsampled.to<downsample.to.vec[length(downsample.to.vec)],],
       aes(
         x=factor(round(log10(downsampled.to), digits=1)),
         y=n.miRs.all.detected.by.iter, fill=lib.method.simple)) + 
  geom_boxplot() + 
  theme_bw() + 
  labs(x="log10 Sequencing Depth (total miRNA-mapping reads)") +
  facet_wrap(~lib.method.simple, strip.position = "bottom", nrow = 1) +
  theme(
    axis.text=element_text(size=14, face="bold", color="black"),
    axis.line=element_line(color="black", size=1.5),
    axis.title=element_text(face="bold", size=16),
    plot.title=element_text(face="bold", size=18),
    panel.grid.major = element_line(color=NA),
    panel.border = element_blank(),
    panel.spacing = unit(0, "lines"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text=element_text(size=14, face="bold", color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=14, face="bold"),
    legend.title=element_text(size=15, face="bold"),
    legend.key.size=unit(1.2, "lines")) 

n.libs.by.method.downsample.level <- ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.DETECTED.ALL.BY.ITER[downsampled.to<downsample.to.vec[length(downsample.to.vec)], .(max.n.detect.pos=max(n.miRs.all.detected.by.iter)), .(lib.method.simple, downsampled.to, total.libs.for.method.simple)]
n.libs.by.method.downsample.level[, downsampled.to.group:=factor(round(log10(downsampled.to), digits=1))]
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.DETECTED.ALL.BY.ITER[, downsampled.to.group:=factor(round(log10(downsampled.to), digits=1))]
ggplot(ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.DETECTED.ALL.BY.ITER[downsampled.to<downsample.to.vec[length(downsample.to.vec)],],
       aes(
         x=lib.method.simple,
         y=n.miRs.all.detected.by.iter, fill=lib.method.simple)) + 
  geom_boxplot() + 
  geom_text(data=n.libs.by.method.downsample.level, aes(label=paste0("n=", total.libs.for.method.simple), y=max.n.detect.pos+15)) +
  theme_bw() +
  facet_wrap(~downsampled.to.group, strip.position = "top", nrow = 1) +
  labs(x="library prep method",
       y="# miRNAs detected in all samples",
       strip="log10 Sequencing Depth (total miRNA-mapping reads)") +
  scale_y_continuous(breaks = seq(0,500,50)) +
  theme(panel.grid.minor = element_blank(),
        axis.text=element_text(size=14, face="bold", color="black"),
        axis.text.x = element_text(angle=50, hjust=1, face="bold", color="black"),
        axis.line=element_line(color="black", size=1.5),
        axis.title=element_text(face="bold", size=16),
        plot.title=element_text(face="bold", size=18),
        panel.grid.major = element_line(color=NA),
        #panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        #strip.background = element_blank(),
        strip.text=element_text(size=14, face="bold", color="black"),
        axis.ticks = element_line(color="black", size=1.5),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=15, face="bold"),
        legend.key.size=unit(1.2, "lines")) 



ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE <- ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO[downsampled.count>0, .(n.mirs.detected.by.sample.iter=length(unique(mir.idx))), by=.(lab.libMethod.replicate, lab.libMethod, lib.method.simple, iter, downsampled.to)]
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE[, downsampled.to.factor:=factor(round(log10(downsampled.to), digits=1))]
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE[, lab.libMethod:=factor(lab.libMethod, levels=sample.levels)]
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE[, lib.method.simple:=factor(lib.method.simple, levels=lib.simple.levels)]
ggplot(ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE, aes(x=lab.libMethod, fill=lib.method.simple, y=n.mirs.detected.by.sample.iter)) + geom_boxplot() + theme_bw() + labs(x="log10 Sequencing Depth (total miRNA-mapping reads)") + theme(
  axis.text=element_text(size=14, face="bold", color="black"),
  axis.title=element_text(face="bold", size=16),
  plot.title=element_text(face="bold", size=18),
  panel.grid.major = element_line(color=NA),
  panel.border = element_rect(size=1.5, color="black"),
  axis.ticks = element_line(color="black", size=1.5),
  legend.position = "top",
  legend.direction = "horizontal",
  legend.text=element_text(size=14, face="bold"),
  legend.title=element_text(size=15, face="bold"),
  legend.key.size=unit(1.2, "lines")
) + facet_wrap(~downsampled.to.factor)

ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE[, downsampled.to.group:=factor(round(log10(downsampled.to), digits=1))]
n.libs.by.method.downsample.level2 <- ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE[downsampled.to<downsample.to.vec[length(downsample.to.vec)], .(max.n.detect.pos=max(n.mirs.detected.by.sample.iter)), .(lib.method.simple, downsampled.to)]
n.libs.by.method.downsample.level2[, downsampled.to.group:=factor(round(log10(downsampled.to), digits=1))]

ggplot(ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE[downsampled.to<downsample.to.vec[length(downsample.to.vec)-1],],
       aes(
         x=lib.method.simple,
         y=n.mirs.detected.by.sample.iter, fill=lib.method.simple)) + 
  geom_boxplot() + 
  theme_bw() +
  facet_wrap(~downsampled.to.group, strip.position = "top", nrow = 1) +
  labs(x="library prep method",
       y="# miRNAs detected in all samples",
       strip="log10 Sequencing Depth (total miRNA-mapping reads)") +
  scale_y_continuous(breaks = seq(0,1000,100)) +
  theme(panel.grid.minor = element_blank(),
        axis.text=element_text(size=14, face="bold", color="black"),
        axis.text.x = element_text(angle=50, hjust=1, face="bold", color="black"),
        axis.line=element_line(color="black", size=1.5),
        axis.title=element_text(face="bold", size=16),
        plot.title=element_text(face="bold", size=18),
        panel.grid.major = element_line(color=NA),
        #panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        #strip.background = element_blank(),
        strip.text=element_text(size=14, face="bold", color="black"),
        axis.ticks = element_line(color="black", size=1.5),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=15, face="bold"),
        legend.key.size=unit(1.2, "lines")) 

total.iterations <- 100
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER <- ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO[downsampled.count>0, .(n.iterations.detected=.N), by=.(lab.libMethod.replicate, mir.idx, downsampled.to, lib.size, lab.libMethod.pool, lab.libMethod, Lab, lib.method.detail, lib.method.simple, pool, mir.seq.id, total.libs.for.method.simple)]
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER[, perc.iter.detect:=n.iterations.detected/total.iterations]
perc.iterations.cutoff=0.9

ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.SUMMARY.BY.SAMPLE <- ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER[perc.iter.detect>perc.iterations.cutoff, .(n.mirs.detected.by.sample.cutoff.perc=length(unique(mir.idx))), by=.(lab.libMethod.replicate, lab.libMethod, lib.method.simple, downsampled.to)]
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.SUMMARY.BY.SAMPLE[, downsampled.to.group:=factor(round(log10(downsampled.to), digits=1))]
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.SUMMARY.BY.SAMPLE[, lab.libMethod:=factor(lab.libMethod, levels=sample.levels)]
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.SUMMARY.BY.SAMPLE[, lib.method.simple:=factor(lib.method.simple, levels=lib.simple.levels)]
ggplot(ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.SUMMARY.BY.SAMPLE[downsampled.to<downsample.to.vec[length(downsample.to.vec)-1],],
       aes(
         x=lib.method.simple,
         y=n.mirs.detected.by.sample.cutoff.perc, fill=lib.method.simple)) + 
  geom_boxplot() + 
  theme_bw() +
  facet_wrap(~downsampled.to.group, strip.position = "top", nrow = 1) +
  labs(x="library prep method",
       y="# miRNAs detected in all samples",
       strip="log10 Sequencing Depth (total miRNA-mapping reads)") +
  scale_y_continuous(breaks = seq(0,1000,100)) +
  theme(panel.grid.minor = element_blank(),
        axis.text=element_text(size=14, face="bold", color="black"),
        axis.text.x = element_text(angle=50, hjust=1, face="bold", color="black"),
        axis.line=element_line(color="black", size=1.5),
        axis.title=element_text(face="bold", size=16),
        plot.title=element_text(face="bold", size=18),
        panel.grid.major = element_line(color=NA),
        #panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        #strip.background = element_blank(),
        strip.text=element_text(size=14, face="bold", color="black"),
        axis.ticks = element_line(color="black", size=1.5),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.text=element_text(size=14, face="bold"),
        legend.title=element_text(size=15, face="bold"),
        legend.key.size=unit(1.2, "lines")) 

#########################################
# START: NEW 12/16/2016 WE WANT TO SEE THE HEATMAPS FOR THOSE OTHER THAN 4N. START HERE SINCE THIS WILL BE CONSISTENT WITH OTHER DOWNSTREAM STEPS.
# COULD TRY USING THE MIRS FOUND TO BE "DETECTED" IN THE ANALYSIS, ABOVE
ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts.EXPANDED.RECOLLAPSED <- ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts.EXPANDED[, .(count.orig=.N), by=.(lab.libMethod.replicate, mir.idx, lib.size)]
setkeyv(ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts.EXPANDED.RECOLLAPSED, c("lab.libMethod.replicate", "mir.idx"))
setkeyv(ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER, c("lab.libMethod.replicate", "mir.idx"))


ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.AND.ORIG.COUNTS <- subset(merge(ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER, ENDOG.MIR.COUNTS.PLASMA.AND.EQUIMOLAR.MATURE.FULL.DGE.ALL.PLASMA.NOZEROS.round.counts.EXPANDED.RECOLLAPSED), lab.libMethod.replicate%in%row.names(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples))
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.AND.ORIG.COUNTS[, downsampled.to.group:=factor(round(log10(downsampled.to), digits=1))]
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.AND.ORIG.COUNTS[, cpm.orig:=(count.orig*(10^6))/lib.size.x]
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.AND.ORIG.COUNTS[, repl:=tstrsplit(lab.libMethod.replicate, split=".", fixed=TRUE)[4]]
ggplot(ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.AND.ORIG.COUNTS, aes(x=count.orig, y=perc.iter.detect, color=downsampled.to.group)) + geom_point(aes(alpha=0.25)) + geom_hline(yintercept = 0.5) + facet_wrap(~lab.libMethod.replicate) + scale_x_log10() + scale_color_brewer("Blues", direction = -1)

ggplot(ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.AND.ORIG.COUNTS[perc.iter.detect<0.9,], aes(color=lib.method.simple, x=cpm.orig)) + stat_ecdf() + facet_wrap(~downsampled.to.group, scale="free_x")

ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.AND.ORIG.COUNTS[, total.libs.for.method.simple:=length(unique(lab.libMethod.replicate)), by=.(lib.method.simple, downsampled.to.group)]
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.AND.ORIG.COUNTS.SUMMARIZE.OVER.LIB.METHOD <- ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.AND.ORIG.COUNTS[, .(tot.detected.lib.method.simple=sum(n.iterations.detected)), by=.(mir.idx, mir.seq.id, downsampled.to, lib.method.simple, total.libs.for.method.simple)]
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.AND.ORIG.COUNTS.SUMMARIZE.OVER.LIB.METHOD[, total.iterations.method.simple:=total.libs.for.method.simple*100]
ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.AND.ORIG.COUNTS.SUMMARIZE.OVER.LIB.METHOD[, perc.iter.detect.over.samples:= tot.detected.lib.method.simple/total.iterations.method.simple]

ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.AND.ORIG.COUNTS.SUMMARIZE.OVER.LIB.METHOD[, n.lib.methods.over.90.perc.detected:=sum(ifelse(perc.iter.detect.over.samples>0.9, 1, 0)), by=.(mir.seq.id, downsampled.to)]
common.mirs <- unique(as.character(ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.AND.ORIG.COUNTS.SUMMARIZE.OVER.LIB.METHOD[downsampled.to==100000 & n.lib.methods.over.90.perc.detected==4]$mir.seq.id))
common.mirs.to.two.protocols <- unique(as.character(ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.AND.ORIG.COUNTS.SUMMARIZE.OVER.LIB.METHOD[downsampled.to==100000 & n.lib.methods.over.90.perc.detected>=2]$mir.seq.id))
common.mirs.to.three.protocols <- unique(as.character(ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.AND.ORIG.COUNTS.SUMMARIZE.OVER.LIB.METHOD[downsampled.to==100000 & n.lib.methods.over.90.perc.detected>=3]$mir.seq.id))
common.mirs.to.1.protocols <- unique(as.character(ALL.DOWNSAMPLED.PLAMA.WITH.SAMPLE.INFO.SUMMARY.BY.SAMPLE.SUMMARIZE.ITER.AND.ORIG.COUNTS.SUMMARIZE.OVER.LIB.METHOD[downsampled.to==100000 & perc.iter.detect.over.samples>=0.95]$mir.seq.id))
protocol.specific <- common.mirs.to.1.protocols

ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA <- calcNormFactors(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA, method="RLE")
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm <- cpm(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA, normalized.lib.sizes = TRUE, prior.count=1, log=TRUE)


plasma.col.info <- data.frame(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples[, c("lib.size", "Lab", "lib.method.detail", "lib.method.simple")], row.names=row.names(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples))
plasma.col.info$lib.method.simple <- factor(plasma.col.info$lib.method.simple, levels=lib.simple.levels)
plasma.col.info$lib.method.detail <- factor(plasma.col.info$lib.method.detail, levels=lib.detail.levels)

ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.common.miRs <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm[common.mirs,]
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.common.miRs.2.methods <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm[common.mirs.to.two.protocols,]
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.common.miRs.3.methods <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm[common.mirs.to.three.protocols,]
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.common.miRs.1.methods <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm[common.mirs.to.1.protocols,]

drop.common.mirs <- common.mirs.to.1.protocols[!common.mirs.to.1.protocols%in%common.mirs]
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.common.miRs.not.common <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm[drop.common.mirs,]

pheatmap(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.common.miRs, main = "Plasma Pool:\nRLE-Normalized, Log2 CPM", show_rownames = FALSE, show_colnames=FALSE, cluster_rows = TRUE, annotation_col = plasma.col.info, fontsize = 12, annotation_colors = FUNCTION.filter.ann_colors(ann_colors, plasma.col.info), width=14, height=10)#, filename = "PlasmaPool_heatmap_logCPM_CPM100_cutoff.png")
pheatmap(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.common.miRs.2.methods, main = "Plasma Pool:\nRLE-Normalized, Log2 CPM", show_rownames = FALSE, show_colnames=FALSE, cluster_rows = TRUE, annotation_col = plasma.col.info, fontsize = 12, annotation_colors = FUNCTION.filter.ann_colors(ann_colors, plasma.col.info), width=14, height=10)#, filename = "PlasmaPool_heatmap_logCPM_CPM100_cutoff.png")
pheatmap(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.common.miRs.not.common, main = "Plasma Pool:\nRLE-Normalized, Log2 CPM", show_rownames = FALSE, show_colnames=FALSE, cluster_rows = TRUE, annotation_col = plasma.col.info, fontsize = 12, annotation_colors = FUNCTION.filter.ann_colors(ann_colors, plasma.col.info), width=14, height=10)#, filename = "PlasmaPool_heatmap_logCPM_CPM100_cutoff.png")
pheatmap(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.common.miRs.1.methods, main = "Plasma Pool:\nRLE-Normalized, Log2 CPM", show_rownames = FALSE, show_colnames=FALSE, cluster_rows = TRUE, annotation_col = plasma.col.info, fontsize = 12, annotation_colors = FUNCTION.filter.ann_colors(ann_colors, plasma.col.info), width=14, height=10)#, filename = "PlasmaPool_heatmap_logCPM_CPM100_cutoff.png")

ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.gt.cutoff <- ifelse((10^ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm)>=100, 1, 0)
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.gt.cutoff.tot <- sumTechReps(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.gt.cutoff,  ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple)
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.gt.cutoff.tot.max <- apply(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.gt.cutoff.tot, 2, max)
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.gt.cutoff.tot.max.matr <- matrix(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.gt.cutoff.tot.max, nrow = dim(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.gt.cutoff.tot)[1], ncol = 4, byrow = TRUE)
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.gt.cutoff.tot.rowMax <- rowMax(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.gt.cutoff.tot/ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.gt.cutoff.tot.max.matr)
mirs.all.libs.1.method <- row.names(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.gt.cutoff[ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.gt.cutoff.tot.rowMax==1,])


pheatmap(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.common.miRs, main = "Plasma Pool:\nRLE-Normalized, Log2 CPM", show_rownames = FALSE, show_colnames=FALSE, cluster_rows = TRUE, annotation_col = plasma.col.info, fontsize = 12, annotation_colors = FUNCTION.filter.ann_colors(ann_colors, plasma.col.info), width=14, height=10)#, filename = "PlasmaPool_heatmap_logCPM_CPM100_cutoff.png")

pheatmap(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm[mirs.all.libs.1.method,], main = "Plasma Pool:\nRLE-Normalized, Log2 CPM", show_rownames = FALSE, show_colnames=TRUE, cluster_rows = TRUE, annotation_col = plasma.col.info[, c(1,4)], fontsize = 12, annotation_colors = FUNCTION.filter.ann_colors(ann_colors, plasma.col.info[, c(1,4)]), width=14, height=10, labels_col = sub(".PlasmaPool", "", sub("^X", "", colnames(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm))))#, filename = "PlasmaPool_heatmap_logCPM_CPM100_cutoff.png")


pheatmap(cor(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.filt, method = "spearman", use = "pair")^2,  annotation_col = plasma.col.info, annotation_colors = FUNCTION.filter.ann_colors(ann_colors, plasma.col.info), annotation_row=plasma.col.info, show_rownames=FALSE, show_colnames=FALSE, fontsize = 14, filename="PlasmaPool_CORRELATION_SPEARMAN_heatmap_logCPM_CPM100_cutoff.jpeg")

# for manuscrpt
plasma.col.info.man <- plasma.col.info
plasma.col.info.man$log2.lib.size <- log2(plasma.col.info.man$lib.size)
plasma.col.info.man <- plasma.col.info.man[, !colnames(plasma.col.info.man)%in%c("Lab", "lib.size")]

# 1-8-17 I think it makes more sense to use the CPM>100 cutoff based on the manusript progression. Comment this out until we figure out for sure what to use
#pheatmap(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.common.miRs,
#         show_rownames = FALSE,
#         show_colnames=TRUE,
#         cluster_rows = TRUE, 
#         labels_col = sub("PlasmaPool.", "", sub("^X", "", colnames(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.DROP.cpm.common.miRs))),
#         annotation_col = plasma.col.info.man[, c(3,1,2)],
##         fontsize = 10,
#         treeheight_row = 15,
#         treeheight_col = 15,
#         annotation_colors = FUNCTION.filter.ann_colors(ann_colors, plasma.col.info.man),
#         width=14, 
#         height=10, filename = "FIG_5A_PlasmaPool_heatmap_logCPM_CPM100_cutoff.pdf")








# END: NEW 12/16/2016 WE WANT TO SEE THE HEATMAPS FOR THOSE OTHER THAN 4N. START HERE SINCE THIS WILL BE CONSISTENT WITH OTHER DOWNSTREAM STEPS
#################################################


##################
#### NEW %CV INTRA-LAB FOR MANUSCRIPT!! 11-29-16
### DIDNT INCLUDE CLEANTAG/NEBNEXT IN ORIGINAL CALCULATION SO DO IT AGAIN HERE AND GET SUMMARY INFO
# START WITH THIS: 
#ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA
#Intralab CV
lib.method.group <- unique(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lab.libMethod)
cpm.cutoff=100


## 1-8-17 Think we should use a CPM>100 cutoff. The main concern is that it makes a messy plot when we use our standard cutoff of cpm>100 in all libraries of >= 1 protocol. Therefore, 
# let's try a little different one. On the downsampling cutoff, we chose to plot only commonly-detected miRs for the main text and supplement a heatmap with those detected accross at least one protocol.
# So lets to a cutoff of >= 90% libraries detecting the miR (at >=100cpm), and then plot the others (>100CPM in at least one method) as a supplement
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.norm <- calcNormFactors(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA, method="RLE")
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm <- cpm(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.norm, normalized.lib.sizes = TRUE, prior.count = 1, log = FALSE)
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.logcpm <- cpm(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.norm, normalized.lib.sizes = TRUE, prior.count = 1, log = TRUE)

mirs.passing.90perc.all <- (rowSums(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm>=100)/dim(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm)[2])>=0.9
LOG.CPM.mirs.passing.90perc.all <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.logcpm[mirs.passing.90perc.all,]
plasma.col.info <- data.frame(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples[, c("lib.size", "Lab", "lib.method.detail", "lib.method.simple")], row.names=row.names(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples))
plasma.col.info$lib.method.simple <- factor(plasma.col.info$lib.method.simple, levels=lib.simple.levels)
plasma.col.info$lib.method.detail <- factor(plasma.col.info$lib.method.detail, levels=lib.detail.levels)
plasma.col.info$log.lib.size <- log2(plasma.col.info$lib.size)
plasma.col.info.man <- plasma.col.info[, !colnames(plasma.col.info)%in%c("lib.size", "Lab")]
pheatmap(LOG.CPM.mirs.passing.90perc.all,
         show_rownames = FALSE,
         show_colnames=TRUE,
         cluster_rows = TRUE, 
         labels_col = sub("PlasmaPool.", "", sub("^X", "", colnames(LOG.CPM.mirs.passing.90perc.all))),
         annotation_col = plasma.col.info.man[, c(3,1,2)],
         fontsize = 10,
         treeheight_row = 15,
         treeheight_col = 15,
         annotation_colors = FUNCTION.filter.ann_colors(ann_colors, plasma.col.info.man),
         width=14, 
         height=10)

# Now get the protocol-specific ones
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.gt100 <- ifelse(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm>=100, 1, 0)
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.all <- ifelse(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm>=-Inf, 1, 0)
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple <- factor(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple, lib.simple.levels)
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.n.libs.gt100 <- sumTechReps(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.gt100, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple)
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.n.libs <- sumTechReps(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.all, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple)

ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100 <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.n.libs.gt100/ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.n.libs
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.filt <- rowMaxs(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100)==1

LOG.CPM.mirs.passing.100cpm.by.method <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.logcpm[ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.filt,]
LOG.CPM.mirs.passing.100cpm.by.method.drop.common.miRs <- LOG.CPM.mirs.passing.100cpm.by.method[!row.names(LOG.CPM.mirs.passing.100cpm.by.method)%in%row.names(LOG.CPM.mirs.passing.90perc.all),]
# FIG5A -------------------------------------------------------------------

# Use this instead
# To reorder the plot as much as possible by highest-->lowest expression
LOG.CPM.mirs.passing.100cpm.by.method.row.max <- rowMaxs(LOG.CPM.mirs.passing.100cpm.by.method)
rmin <- min(LOG.CPM.mirs.passing.100cpm.by.method.row.max)
rmax <- max(LOG.CPM.mirs.passing.100cpm.by.method.row.max)
LOG.CPM.mirs.passing.100cpm.by.method.row.max.scl <- 1-(LOG.CPM.mirs.passing.100cpm.by.method.row.max-rmin)/(rmax-rmin)

dist.mirs <- dist(LOG.CPM.mirs.passing.100cpm.by.method)
hc.orig <- hclust(dist.mirs, method = "complete")

callback.mirs <- function(tree_row, mat){
  mir.wts <- LOG.CPM.mirs.passing.100cpm.by.method.row.max.scl[tree_row$labels]
  dend <- reorder(as.dendrogram(tree_row), wts=mir.wts)
  as.hclust(dend)
}

pheatmap(LOG.CPM.mirs.passing.100cpm.by.method, 
         show_rownames = FALSE,
         breaks=seq(max(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.logcpm), min(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.logcpm), length.out=101),
         show_colnames=TRUE,
         cluster_rows = TRUE, 
         annotation_colors=FUNCTION.filter.ann_colors(ann_colors, plasma.col.info.man[, c(3,1,2)]),
         labels_col = sub("PlasmaPool.", "", sub("^X", "", colnames(LOG.CPM.mirs.passing.100cpm.by.method))),
         annotation_col = plasma.col.info.man[, c(3,1,2)],
         fontsize = 10,
         treeheight_row = 15,
         treeheight_col = 15,
         width=14, 
         height=10, 
         clustering_callback = callback.mirs)#,  filename = "FIG_5A_PlasmaPool_heatmap_logCPM_CPM100_cutoff.pdf"

# Muneesh want's to see if the bias has anything to do with the protocol-specific factors we see.
ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR <- calcNormFactors(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR[rowSums(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$counts>0)==dim(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR)[2],], method="RLE")
ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm <- cpm(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR, log=TRUE, normalized.lib.sizes = TRUE, prior.count = 1)
ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.simple <- factor(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.simple, levels=lib.simple.levels)
ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm.mean <- FUNCTION.aveTechReps(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm, ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples$lib.method.simple)


ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm.mean.df <- data.frame(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm.mean)
pheatmap(LOG.CPM.mirs.passing.100cpm.by.method[merge.mean.eq.with.perc.passing$Row.names,], 
         annotation_colors = FUNCTION.filter.ann_colors(ann_colors, plasma.col.info.man),
         breaks=seq(min(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm), max(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm), length.out = 101),
         show_rownames = FALSE,
         show_colnames=TRUE,
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         labels_col = sub("PlasmaPool.", "", sub("^X", "", colnames(LOG.CPM.mirs.passing.100cpm.by.method))),
         annotation_col = plasma.col.info.man[, c(3,1,2)],
         annotation_row = data.frame(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100[merge.mean.eq.with.perc.passing$Row.names,]),
         fontsize = 10,
         treeheight_row = 15,
         treeheight_col = 15,
         width=14, 
         height=10)
ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm.mean.dt <- data.table(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm.mean.df, keep.rownames = TRUE)
ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm.mean.dt.melt <- melt(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm.mean.dt, id.vars = "rn", variable.name = "lib.method.simple", value.name = "mean.logCPM.eqPool")

ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.dt <- data.table(data.frame(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100), keep.rownames=TRUE)
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.dt.melt <- melt(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.dt, id.vars = "rn", variable.name = "lib.method.simple", value.name = "perc.plasma.libs.gt.100cpm")
setkeyv(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm.mean.dt.melt, c("rn", "lib.method.simple"))
setkeyv(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.dt.melt, c("rn", "lib.method.simple"))

merge.melt.mean.eq.plasma.gt.100 <- merge(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm.mean.dt.melt, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.dt.melt)
merge.melt.mean.eq.plasma.gt.100[, perc.is.gt.90:=ifelse(perc.plasma.libs.gt.100cpm>=0.9, "GT.90.perc", "NOT.GT.90.perc")]
merge.melt.mean.eq.plasma.gt.100[, n.libs.all.passing:=sum(ifelse(perc.plasma.libs.gt.100cpm>=0.9, 1, 0)), by=rn]
merge.melt.mean.eq.plasma.gt.100[, lib.method.simple:=factor(sub("^X", "", lib.method.simple), lib.simple.levels)]
setkey(merge.melt.mean.eq.plasma.gt.100, lib.method.simple)
merge.melt.mean.eq.plasma.gt.100[, which.passing:=paste(ifelse(perc.plasma.libs.gt.100cpm>=0.9, 1, 0), collapse="."), by=rn]
merge.melt.mean.eq.plasma.gt.100[, rank.mean.logCPM:=rank(mean.logCPM.eqPool), by=lib.method.simple]

setkey(merge.melt.mean.eq.plasma.gt.100, which.passing)
merge.melt.mean.eq.plasma.gt.100.filt <- data.frame(subset(merge.melt.mean.eq.plasma.gt.100, n.libs.all.passing>0 & n.libs.all.passing <4)[, .(n.libs.all.passing=max(n.libs.all.passing)), by=.(rn, which.passing)], row.names=1)
merge.melt.mean.eq.plasma.gt.100.filt.in.eq <- merge.melt.mean.eq.plasma.gt.100.filt[row.names(merge.melt.mean.eq.plasma.gt.100.filt)%in%row.names(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm),]
ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm.filt <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm[row.names(merge.melt.mean.eq.plasma.gt.100.filt.in.eq),]
col.info.eq.for.plot <- data.frame(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR$samples[, c("lib.size", "lib.method.detail", "lib.method.simple")])
col.info.eq.for.plot$lib.size <- log2(col.info.eq.for.plot$lib.size)
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.row.annot <- sapply(colnames(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100), FUN=function(x){
  x.val <- ifelse(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100[, x]==1, x, "NOT.PASSING")
  x.val <- factor(x.val, levels=c(x, "NOT.PASSING"))
})
colnames(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.row.annot) <- make.names(colnames(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.row.annot))
ann_colors_2 <- FUNCTION.filter.ann_colors(ann_colors, col.info.eq.for.plot)
ann_colors_2$TruSeq <- c(ann_colors$lib.method.simple, NOT.PASSING="grey")
ann_colors_2$X4N <- c(ann_colors$lib.method.simple, NOT.PASSING="grey")
ann_colors_2$CleanTag <- c(ann_colors$lib.method.simple, NOT.PASSING="grey")
ann_colors_2$NEBNext <- c(ann_colors$lib.method.simple, NOT.PASSING="grey")
ann_colors_2b <- c(FUNCTION.filter.ann_colors(ann_colors_2, col.info.eq.for.plot), FUNCTION.filter.ann_colors(ann_colors_2, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.row.annot.filt))


ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.row.annot.filt <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.row.annot[row.names(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm.filt),]
add.vals.for.sort <- merge(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm.mean, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.row.annot.filt, by=0)
add.vals.for.sort$val.not.passing.over.passing.perc.mean <- sapply(1:dim(add.vals.for.sort)[1], FUN=function(x){
  x.cpm <- as.numeric(add.vals.for.sort[x, 2:5])
  x.pass <- add.vals.for.sort[x, 6:9]
  mean.all <- mean(x.cpm)
  out.val <- (mean(x.cpm[x.pass!="NOT.PASSING"])-mean(x.cpm[x.pass=="NOT.PASSING"]))/mean.all
  return(out.val)
})
add.vals.for.sort.dt <- data.table(add.vals.for.sort)
setnames(add.vals.for.sort.dt, "Row.names", "rn")
setkey(add.vals.for.sort.dt, rn)
setorderv(add.vals.for.sort.dt, c("X4N",  "CleanTag.y", "NEBNext.y", "TruSeq.y", "val.not.passing.over.passing.perc.mean"), c(1,1,1,1,-1))
add.vals.for.sort.dt[, n.passing:=sum(ifelse(.SD=="NOT.PASSING", 0, 1)), by=1:nrow(add.vals.for.sort.dt), .SDcols=c("X4N",  "CleanTag.y", "NEBNext.y", "TruSeq.y")]
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.row.annot.filt2 <- data.frame(add.vals.for.sort, row.names=1)[, c(9:5)]
colnames(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.row.annot.filt2) <- sub(".y$", "", colnames(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.row.annot.filt2))
add.vals.for.sort.dt.SORT.AND.FILT <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.cpm.perc.libs.gt100.row.annot.filt2[add.vals.for.sort.dt$rn,]
add.vals.for.sort.dt.allBut1.pass <- add.vals.for.sort.dt.SORT.AND.FILT[add.vals.for.sort.dt[n.passing==3]$rn,]
add.vals.for.sort.dt.allBut1.pass$which.not.passing <- sapply(1:nrow(add.vals.for.sort.dt.allBut1.pass), FUN=function(x){x.val <- add.vals.for.sort.dt.allBut1.pass[x,]; sub("^X", "", colnames(add.vals.for.sort.dt.allBut1.pass)[x.val=="NOT.PASSING"])})

ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm.filt <- ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm[row.names(add.vals.for.sort.dt.allBut1.pass),]
ann_colors_2$which.not.passing <- ann_colors_2$lib.method.simple
ann_colors_2.inv <- sapply(1:length(ann_colors_2), FUN=function(x){ if(!names(ann_colors_2)[x]%in%c(lib.simple.levels, make.names(lib.simple.levels))){
  return(ann_colors_2[[x]])
} else{
  old.vals <- ann_colors_2[[x]]
  new.vals <- rep("grey", length(old.vals))
  names(new.vals) <- names(old.vals)
  new.vals["NOT.PASSING"] <- old.vals[sub("^X", "", names(ann_colors_2)[x])]
  return(new.vals)
}
})
pheatmap(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm.filt, drop_levels = TRUE,
         annotation_colors = ann_colors_2,
         breaks=seq(min(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm), max(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm), length.out = 101),
         show_rownames = FALSE,
         show_colnames=TRUE,
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         labels_col = sub("SynthEQ.", "", sub("^X", "", colnames(ENDOG.MIR.COUNTS.EQUIMOLAR.MATURE.FULL.DGE.ALL.EQUIMOLAR.cpm))),
         annotation_col = col.info.eq.for.plot,
         annotation_row = add.vals.for.sort.dt.allBut1.pass[, c(1,6)],
         fontsize = 10,
         treeheight_row = 15,
         treeheight_col = 15,
         width=14, 
         height=10, filename = "testing_why_gaps.pdf")



merge.melt.mean.eq.plasma.gt.100.mean.by.eq.grp <- dcast.data.table(merge.melt.mean.eq.plasma.gt.100[n.libs.all.passing>0 & n.libs.all.passing<4], rn+n.libs.all.passing+which.passing~perc.is.gt.90, fun.aggregate = mean, fill = NA, value.var = "mean.logCPM.eqPool")
merge.melt.mean.eq.plasma.gt.100.mean.by.eq.grp[, pass.minus.notPass:=(GT.90.perc-NOT.GT.90.perc), by=1:nrow(merge.melt.mean.eq.plasma.gt.100.mean.by.eq.grp)]

ggplot(merge.melt.mean.eq.plasma.gt.100.mean.by.eq.grp, aes(fill=which.passing, x=which.passing, y=pass.minus.notPass)) + geom_boxplot()



PLASMA.TRUSEQ.INTRALAB.CV.ALL <- do.call(rbind, sapply(lib.method.group, USE.NAMES = FALSE, simplify = FALSE, FUN=function(grp){
  dge.grp <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lab.libMethod==grp]
  dge.grp.filt.zeros <- dge.grp[rowSums(dge.grp$counts>0)==dim(dge.grp)[2],]
  dge.grp.filt.zeros <- calcNormFactors(dge.grp.filt.zeros, method="RLE")
  sample.ids <- colnames(dge.grp.filt.zeros)
  cpm.subs <- cpm(dge.grp.filt.zeros, normalized.lib.sizes=FALSE, log=FALSE)
  cv <- apply(cpm.subs, 1, FUN=function(x){
    sdx <- sd(x)
    meanx <- mean(x)
    cv.ln <- sdx/meanx
    return(cv.ln*100)
  })
  qcd <- apply(cpm.subs, 1, FUN=function(x){
    q1 <- quantile(x, 0.25)
    q3 <- quantile(x, 0.75)
    iqr <- (q3-q1)/2
    mh <- (q1+q3)/2
    return(iqr/mh)
  })
  disp.df <- data.table(data.frame(mature.seqID=row.names(cpm.subs), intralab.cv=cv, intralab.qcd=qcd, lab.libMethod=grp))
  return(disp.df)
}))

PLASMA.TRUSEQ.INTRALAB.CV.ALL.CPM100.FILT <- do.call(rbind, sapply(lib.method.group, USE.NAMES = FALSE, simplify = FALSE, FUN=function(grp){
  dge.grp <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lab.libMethod==grp]
  dge.grp.filt.zeros <- dge.grp[rowSums(dge.grp$counts>0)==dim(dge.grp)[2],]
  dge.grp.filt.zeros <- calcNormFactors(dge.grp.filt.zeros, method="RLE")
  sample.ids <- colnames(dge.grp.filt.zeros)
  cpm.subs.full <- cpm(dge.grp.filt.zeros, normalized.lib.sizes=FALSE, log=FALSE)
  cpm.subs <- cpm.subs.full[rowSums(cpm.subs.full>=cpm.cutoff)==dim(cpm.subs.full)[2],]
  cv <- apply(cpm.subs, 1, FUN=function(x){
    sdx <- sd(x)
    meanx <- mean(x)
    cv.ln <- sdx/meanx
    return(cv.ln*100)
  })
  qcd <- apply(cpm.subs, 1, FUN=function(x){
    q1 <- quantile(x, 0.25)
    q3 <- quantile(x, 0.75)
    iqr <- (q3-q1)/2
    mh <- (q1+q3)/2
    return(iqr/mh)
  })
  disp.df <- data.table(data.frame(mature.seqID=row.names(cpm.subs), intralab.cv=cv, intralab.qcd=qcd, lab.libMethod=grp))
  return(disp.df)
}))

PLASMA.TRUSEQ.INTRALAB.CV.ALL[, c("lib.method.detail", "Lab"):=tstrsplit(lab.libMethod, split=".", fixed=TRUE)]
PLASMA.TRUSEQ.INTRALAB.CV.ALL[, lib.method.simple:=tstrsplit(lib.method.detail, split="_")[1]]
PLASMA.TRUSEQ.INTRALAB.CV.ALL[, `:=`(lib.method.detail=factor(lib.method.detail, levels=lib.detail.levels), lib.method.simple=factor(lib.method.simple, levels=lib.simple.levels), lab.libMethod=factor(lab.libMethod, levels=sample.levels))]
ggplot(PLASMA.TRUSEQ.INTRALAB.CV.ALL,
       aes(x = lab.libMethod,
           y = intralab.cv,
           fill=lib.method.simple))  +
  geom_violin(
    draw_quantiles=c(0.25, 0.5, 0.75)
    #outlier.color = NA,
    #size = 0.8,
    #width = 0.6
  ) + scale_y_continuous(limits = c(0,210), breaks = seq(0, 220, 20)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", colour = "black", angle=50, hjust=1),
    axis.text.y = element_text(size = 14, face = "bold", colour = "black"),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major = element_line(color = NA),
    panel.border = element_rect(size = 1.5, color = "black"),
    axis.ticks = element_line(color = "black", size = 1),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=16, face="bold"),
    legend.title=element_text(size=18, face="bold"),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(title = "PlasmaPool: %CV Within-labs",
       x = "Library Prep Method",
       y = "%CV") 
ggsave("PLASMA_POOL_CV_INTRALAB_noNorm_FiltZeroCountOnly_violin.png")


#Intralab qcd
ggplot(PLASMA.TRUSEQ.INTRALAB.CV.ALL,
       aes(x = lab.libMethod,
           y = intralab.qcd,
           fill=lib.method.simple))  +
  geom_violin(
    draw_quantiles=c(0.25, 0.5, 0.75)
    #outlier.color = NA,
    #size = 0.8,
    #width = 0.6
  ) + scale_y_continuous(limits=c(0,1), breaks = seq(0, 1, 0.1)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14, face = "bold", colour = "black", angle=50, hjust=1),
    axis.text.y = element_text(size = 14, face = "bold", colour = "black"),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.major = element_line(color = NA),
    panel.border = element_rect(size = 1.5, color = "black"),
    axis.ticks = element_line(color = "black", size = 1),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=16, face="bold"),
    legend.title=element_text(size=18, face="bold"),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(title = "PlasmaPool: QCD Within-labs",
       x = "Library Prep Method",
       y = "QCD") 
ggsave("PlasmaPool_QCD_INTRALAB_noNorm_FiltZeroCountOnly_violin.png")

# FIG5B--%CV --------------------------------------------------------------

# Now for CPM>100
PLASMA.TRUSEQ.INTRALAB.CV.ALL.CPM100.FILT[, c("lib.method.detail", "Lab"):=tstrsplit(lab.libMethod, split=".", fixed=TRUE)]
PLASMA.TRUSEQ.INTRALAB.CV.ALL.CPM100.FILT[, lib.method.simple:=tstrsplit(lib.method.detail, split="_")[1]]
PLASMA.TRUSEQ.INTRALAB.CV.ALL.CPM100.FILT[, `:=`(lib.method.detail=factor(lib.method.detail, levels=lib.detail.levels), lib.method.simple=factor(lib.method.simple, levels=lib.simple.levels), lab.libMethod=factor(lab.libMethod, levels=sample.levels))]
ggplot(PLASMA.TRUSEQ.INTRALAB.CV.ALL.CPM100.FILT,
       aes(x = lab.libMethod,
           y = intralab.cv,
           fill=lib.method.simple))  +
  geom_violin(
    draw_quantiles=c(0.25, 0.5, 0.75)
    #outlier.color = NA,
    #size = 0.8,
    #width = 0.6
  ) + scale_y_continuous(limits = c(0,210), breaks = seq(0, 220, 20)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10, colour = "black", angle=50, hjust=1),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 10),
    panel.grid.major = element_line(color = NA),
    panel.border = element_rect(size = 0.5, color = "black"),
    axis.ticks = element_line(color = "black", size = 1),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=10),
    legend.title=element_text(size=10)
    #legend.key.size=unit(1.0, "lines")
  ) +
  labs(#title = "PlasmaPool: %CV Within-labs",
    #x = "Library Prep Method",
    x = NULL,
    y = "%CV") 
ggsave("FIG_5B_PLASMA_POOL_CV_INTRALAB_noNorm_FiltCPM_gt100_violin.svg", units = "in", width = 5.5, height = 4)


# FIG5B--QCD --------------------------------------------------------------


#Intralab qcd
ggplot(PLASMA.TRUSEQ.INTRALAB.CV.ALL.CPM100.FILT,
       aes(x = lab.libMethod,
           y = intralab.qcd,
           fill=lib.method.simple))  +
  geom_violin(
    draw_quantiles=c(0.25, 0.5, 0.75)
    #outlier.color = NA,
    #size = 0.8,
    #width = 0.6
  )  + scale_y_continuous(limits=c(0,1), breaks = seq(0, 1, 0.1)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10, colour = "black", angle=50, hjust=1),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title = element_text(size = 10),
    plot.title = element_text(size = 10),
    panel.grid.major = element_line(color = NA),
    panel.border = element_rect(size = 0.5, color = "black"),
    axis.ticks = element_line(color = "black", size = 1),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=10),
    legend.title=element_text(size=10)
    #legend.key.size=unit(1.0, "lines")
  ) +
  labs(#title = "PlasmaPool: %CV Within-labs",
    #x = "Library Prep Method",
    x = NULL,
    y = "QCD") 
ggsave("FIG_5B_PLASMA_POOL_QCD_INTRALAB_noNorm_FiltCPM_gt100_violin.svg", units = "in", width = 5.5, height = 4)


# Make summary info for %CV/QCD in plasma pools
summary.quantiles <- c(0.02, 0.25, 0.50, 0.75, 0.98)

PLASMA.TRUSEQ.INTRALAB.CV.ALL.CPM100.FILT[, (sub("^0.", "intralab.cv", summary.quantiles)):=lapply(summary.quantiles, USE.NAMES=FALSE, SIMPLIFY=FALSE, quantile, x=intralab.cv), by=.(lib.method.simple)]
PLASMA.TRUSEQ.INTRALAB.CV.ALL.CPM100.FILT[, (sub("^0.", "intralab.qcd", summary.quantiles)):=lapply(summary.quantiles, USE.NAMES=FALSE, SIMPLIFY=FALSE, quantile, x=intralab.qcd), by=.(lib.method.simple)]
PLASMA.TRUSEQ.INTRALAB.CV.ALL.CPM100.FILT.SUMMARY <- PLASMA.TRUSEQ.INTRALAB.CV.ALL.CPM100.FILT[, lapply(.SD, max), by=.(lib.method.simple), .SDcols=c(sub("^0.", "intralab.cv", summary.quantiles), sub("^0.", "intralab.qcd", summary.quantiles))]
#print(PLASMA.TRUSEQ.INTRALAB.CV.ALL.CPM100.FILT.SUMMARY, row.names = FALSE)

PLASMA.TRUSEQ.INTRALAB.CV.ALL[, (sub("^0.", "intralab.cv", summary.quantiles)):=lapply(summary.quantiles, USE.NAMES=FALSE, SIMPLIFY=FALSE, quantile, x=intralab.cv), by=.(lib.method.simple)]
PLASMA.TRUSEQ.INTRALAB.CV.ALL[, (sub("^0.", "intralab.qcd", summary.quantiles)):=lapply(summary.quantiles, USE.NAMES=FALSE, SIMPLIFY=FALSE, quantile, x=intralab.qcd), by=.(lib.method.simple)]
PLASMA.TRUSEQ.INTRALAB.CV.ALL.SUMMARY <- PLASMA.TRUSEQ.INTRALAB.CV.ALL[, lapply(.SD, max), by=.(lib.method.simple), .SDcols=c(sub("^0.", "intralab.cv", summary.quantiles), sub("^0.", "intralab.qcd", summary.quantiles))]

PLASMA.TRUSEQ.INTRALAB.CV.ALL.SUMMARY[, FILTER:="ALL.NONZERO"]
PLASMA.TRUSEQ.INTRALAB.CV.ALL.CPM100.FILT.SUMMARY[, FILTER:="ALL.CPM.GT.100"]

print(rbind(PLASMA.TRUSEQ.INTRALAB.CV.ALL.SUMMARY, PLASMA.TRUSEQ.INTRALAB.CV.ALL.CPM100.FILT.SUMMARY), row.names = FALSE)

# Now TruSeq Inter-Lab for summary info
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.TRUSEQ <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lib.method.simple=="TruSeq"]
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.TRUSEQ.NOZEROS <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.TRUSEQ[rowSums(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.TRUSEQ$counts>0)==dim(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.TRUSEQ)[2],]
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.TRUSEQ.NOZEROS.CPM.noNorm <- cpm(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.TRUSEQ.NOZEROS, normalized.lib.sizes=FALSE, log=FALSE)
MIR.IDS.TruSeq.All.CPM.gt100 <- row.names(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.TRUSEQ.NOZEROS.CPM.noNorm[rowSums(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.TRUSEQ.NOZEROS.CPM.noNorm>=100)==dim(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.TRUSEQ.NOZEROS.CPM.noNorm)[2],])

ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.TRUSEQ.NOZEROS.CPM.noNorm.aveReps <- FUNCTION.aveTechReps(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.TRUSEQ.NOZEROS.CPM.noNorm, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.TRUSEQ$samples$Lab)



ENDOG.MIR.TRUSEQ.FILT.CV.NOZEROCOUNTS <- apply(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.TRUSEQ.NOZEROS.CPM.noNorm.aveReps, 1, FUN=function(x){
  sdx <- sd(x)
  meanx <- mean(x)
  cv.ln <- 100*sdx/meanx
  return(cv.ln)
})
ENDOG.MIR.TRUSEQ.FILT.QCD.NOZEROCOUNTS <- apply(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.TRUSEQ.NOZEROS.CPM.noNorm.aveReps, 1, FUN=function(x){
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr <- (q3-q1)/2
  mh <- (q3+q1)/2
  qcd <- iqr/mh
  return(qcd)
})
ENDOG.MIR.TRUSEQ.FILT.CV.CPM.gt100 <- apply(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.TRUSEQ.NOZEROS.CPM.noNorm.aveReps[MIR.IDS.TruSeq.All.CPM.gt100,], 1, FUN=function(x){
  sdx <- sd(x)
  meanx <- mean(x)
  cv.ln <- 100*sdx/meanx
  return(cv.ln)
})
ENDOG.MIR.TRUSEQ.FILT.QCD.CPM.gt100 <- apply(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.TRUSEQ.NOZEROS.CPM.noNorm.aveReps[MIR.IDS.TruSeq.All.CPM.gt100,], 1, FUN=function(x){
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr <- (q3-q1)/2
  mh <- (q3+q1)/2
  qcd <- iqr/mh
  return(qcd)
})

do.call(rbind, list(ENDOG.MIR.TRUSEQ.FILT.CV.NOZEROCOUNTS=quantile(ENDOG.MIR.TRUSEQ.FILT.CV.NOZEROCOUNTS, summary.quantiles),
                    ENDOG.MIR.TRUSEQ.FILT.QCD.NOZEROCOUNTS=quantile(ENDOG.MIR.TRUSEQ.FILT.QCD.NOZEROCOUNTS, summary.quantiles),
                    ENDOG.MIR.TRUSEQ.FILT.CV.CPM.gt100=quantile(ENDOG.MIR.TRUSEQ.FILT.CV.CPM.gt100, summary.quantiles),
                    ENDOG.MIR.TRUSEQ.FILT.QCD.CPM.gt100=quantile(ENDOG.MIR.TRUSEQ.FILT.QCD.CPM.gt100, summary.quantiles)))


###################################




PLASMA.TRUSEQ.INTRALAB.CV.ALL.CPM100.FILT <- do.call(rbind, sapply(lib.method.group, USE.NAMES = FALSE, simplify = FALSE, FUN=function(grp){
  dge.grp <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[, ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lab.libMethod==grp]
  dge.grp.filt.zeros <- dge.grp[rowSums(dge.grp$counts>0)==dim(dge.grp)[2],]
  dge.grp.filt.zeros <- calcNormFactors(dge.grp.filt.zeros, method="RLE")
  sample.ids <- colnames(dge.grp.filt.zeros)
  cpm.subs.full <- cpm(dge.grp.filt.zeros, normalized.lib.sizes=FALSE, log=FALSE)
  cpm.subs <- cpm.subs.full[rowSums(cpm.subs.full>=cpm.cutoff)==dim(cpm.subs.full)[2],]
  cv <- apply(cpm.subs, 1, FUN=function(x){
    sdx <- sd(x)
    meanx <- mean(x)
    cv.ln <- sdx/meanx
    return(cv.ln*100)
  })
  qcd <- apply(cpm.subs, 1, FUN=function(x){
    q1 <- quantile(x, 0.25)
    q3 <- quantile(x, 0.75)
    iqr <- (q3-q1)/2
    mh <- (q1+q3)/2
    return(iqr/mh)
  })
  disp.df <- data.table(data.frame(mature.seqID=row.names(cpm.subs), intralab.cv=cv, intralab.qcd=qcd, lab.libMethod=grp))
  return(disp.df)
}))

############### SUMMARY CORRELATION COEFFICIENT TABLES FOR PLASMA POOL ##############################
# Convert the cpm.rle matrix to a melted data.frame for data.table calculations of correlations.
lib.method.group <- unique(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$samples$lab.libMethod)
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.FILT.ZEROS <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA[rowSums(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA$counts>0)==dim(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA)[2],]
ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.FILT.ZEROS <- calcNormFactors(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.FILT.ZEROS, method="RLE")
cpm.rle.all.plasm <- cpm(ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.FILT.ZEROS, normalized.lib.sizes = TRUE)
cor.matrix.all.plasma.square <- cor(cpm.rle.all.plasm, use="pair", method="spearman")
cor.matrix.all.plasma <- data.table(cbind(melt(cor.matrix.all.plasma.square), melt(upper.tri(cor.matrix.all.plasma.square, diag=TRUE))))
setnames(cor.matrix.all.plasma, c("lab.libMethod.replicate.A", "lab.libMethod.replicate.B", "rho", "row", "col", "upper.tri"))
cor.matrix.all.plasma[, `:=`(id.A=ifelse(upper.tri==TRUE, as.character(lab.libMethod.replicate.A), as.character(lab.libMethod.replicate.B)), id.B=ifelse(upper.tri==TRUE, as.character(lab.libMethod.replicate.B), as.character(lab.libMethod.replicate.A)))]
cor.matrix.all.plasma[, `:=`(lab.libMethod.replicate.A=id.A, lab.libMethod.replicate.B=id.B)]
cor.matrix.all.plasma[, `:=`(id.A=NULL, id.B=NULL)]
short.sample.info.plasma.A <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.FILT.ZEROS$samples[, c("Lab", "lib.method.detail", "lib.method.simple", "replicate")]
short.sample.info.plasma.B <- ENDOG.MIR.COUNTS.PLASMA.MATURE.FULL.DGE.ALL.PLASMA.FILT.ZEROS$samples[, c("Lab", "lib.method.detail", "lib.method.simple", "replicate")]

colnames(short.sample.info.plasma.A) <- paste0(colnames(short.sample.info.plasma.A), ".A")
colnames(short.sample.info.plasma.B) <- paste0(colnames(short.sample.info.plasma.B), ".B")

cor.matrix.all.plasma.tmp <- merge.data.frame(cor.matrix.all.plasma, short.sample.info.plasma.A, by.x=1, by.y=0)
cor.matrix.all.plasma.sample.info <- data.table(merge.data.frame(cor.matrix.all.plasma.tmp, short.sample.info.plasma.B, by.x=2, by.y=0))
lib.detail.levels <- rev(unique(cor.matrix.all.plasma.sample.info$lib.method.detail.A))
lib.simple.levels <- rev(unique(cor.matrix.all.plasma.sample.info$lib.method.simple.A))
reorder.paste.fun <- function(query.vec, lookup.vec){
  matches <- match(query.vec, lookup.vec)
  comp <- paste(query.vec[order(matches)], collapse=".VS.")
  return(comp)
}

cor.matrix.all.plasma.sample.info[, comparison.detail:=apply(.SD, 1, reorder.paste.fun, lookup.vec=c(lib.detail.levels, "4N_NEXT")), .SDcols=c("lib.method.detail.A", "lib.method.detail.B")]
cor.matrix.all.plasma.sample.info[, comparison.simple:=apply(.SD, 1, reorder.paste.fun, lookup.vec=lib.simple.levels), .SDcols=c("lib.method.simple.A", "lib.method.simple.B")]


setkeyv(cor.matrix.all.plasma.sample.info, c("lab.libMethod.replicate.B", "lab.libMethod.replicate.A"))
cor.matrix.all.plasma.sample.info.unique <- unique(cor.matrix.all.plasma.sample.info)


# form meanrho
FUNCTION.meanRhoFromCoefs <- function(r, ns, nr){
  delr <- length(r) - length(r[(r < 1) & (r > -1)])
  r <- r[(r < 1) & (r > -1)]
  rz <- 1/2 * log((1 + r)/(1 - r))
  mrz <- mean(rz)
  coeff <- (exp(2 * mrz) - 1)/(exp(2 * mrz) + 1)
  SE <- sqrt(1/(ns - 3))
  u <- coeff/SE
  mean.r <- mean(r)
  r02 <- quantile(r, 0.02)
  r98 <- quantile(r, 0.98)
  p.value <- 2 * (1 - pnorm(abs(u)))
  rval <- list(irr.name = "Rho", nlibs= nr, value = coeff, stat.name = "z", statistic = u, p.value = p.value, simple.mean.r=mean.r, r02=r02, r98=r98)
}
n.miRs <- dim(cpm.rle.all.plasm)[1]
cor.matrix.all.plasma.sample.info.cor.summary.details <- cor.matrix.all.plasma.sample.info.unique[, { nr=length(unique(c(lab.libMethod.replicate.A, lab.libMethod.replicate.B))); FUNCTION.meanRhoFromCoefs(rho, ns=n.miRs, nr=nr) }, by=comparison.detail]
cor.matrix.all.plasma.sample.info.cor.summary.simple <- cor.matrix.all.plasma.sample.info.unique[, { nr=length(unique(c(lab.libMethod.replicate.A, lab.libMethod.replicate.B))); FUNCTION.meanRhoFromCoefs(rho, ns=n.miRs, nr=nr) }, by=comparison.simple]
setnames(cor.matrix.all.plasma.sample.info.cor.summary.details, "comparison.detail", "comparison")
setnames(cor.matrix.all.plasma.sample.info.cor.summary.simple, "comparison.simple", "comparison")
cor.matrix.all.plasma.sample.info.cor.summary.all <- rbind(cor.matrix.all.plasma.sample.info.cor.summary.simple, cor.matrix.all.plasma.sample.info.cor.summary.details)
cor.matrix.all.plasma.sample.info.cor.summary.all[, c("MethodA", "MethodB"):=tstrsplit(comparison, split=".VS.")]
cor.matrix.all.plasma.sample.info.cor.summary.all
write.table(cor.matrix.all.plasma.sample.info.cor.summary.all, "Spearman_correlation_plasma_Pool_summaries.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)



##################




pheatmap(log10(mature.counts.dt.cast.cpm.df+1))

mir.sample.info <- mature.counts.dt.zerocounts[, .(total.mir.adjusted.counts=max(total.mir))]
mature.counts.dt.zerocounts[, cpm.rank.in.sample:=rank(-cpm.multimap.adjust, ties.method = "min"), by=.(SAMPLE.ID)]
mature.counts.dt.zerocounts[, SAMPLE.ID:=gsub("-", ".", SAMPLE.ID)]
mir.sample.info <- mature.counts.dt.zerocounts[, .(sum.multimapAdjustedReadCount.top10.miRs=sum(ifelse(cpm.rank.in.sample<=10, multimapAdjustedReadCount, 0))), by=.(SAMPLE.ID, total.mir.adjusted.counts, Lab, lib.method, replicate)]
mir.sample.info[, percent.reads.top10.miRs:=sum.multimapAdjustedReadCount.top10.miRs/total.mir.adjusted.counts]
mir.sample.info[, log.total.adjusted.counts:=log10(total.mir.adjusted.counts)]
mir.sample.info[, lib.method.simple:=tstrsplit(lib.method, split="-")[1]]
col.annot <- data.frame(subset(mir.sample.info, select=c("SAMPLE.ID", "log.total.adjusted.counts", "Lab", "lib.method.simple")), row.names=1)

pheatmap(log10(mature.counts.dt.cast.cpm.df+1), show_rownames = FALSE, annotation_col = col.annot, fontsize = 14, filename = "MANUSCRIPT_Full_Plasma_miRNA_heatmap.tiff")
mir.sample.info[, SAMPLE:=paste(Lab, lib.method, sep=":")]
mir.sample.info[, SAMPLE:=factor(SAMPLE, levels=levels)]
ggplot(mir.sample.info, aes(x=SAMPLE, y=percent.reads.top10.miRs, fill=replicate)) + 
  geom_bar(stat="identity", position="dodge", color="black", width=0.6) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1), expand=c(0, 0), labels = scales::percent) +
  theme(
    axis.text=element_text(size=14, face="bold"),
    axis.text.x=element_text(angle=50, hjust=1),
    axis.title=element_text(face="bold", size=16),
    plot.title=element_text(face="bold", size=18),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=14, face="bold"),
    legend.title=element_text(size=15, face="bold"),
    legend.key.size=unit(1.2, "lines")
  ) +
  labs(y="% Total Counts",
       title="Percent of total miRNA-mapping reads\nmapping to top 10 miRs")
ggsave("MANUSCRIPT_Percent_plasma_reads_from_top10_expressed_mirs.jpeg", width=12, height=9)

mature.counts.dt.zerocounts.by.sample <- mature.counts.dt.zerocounts[, .(mir.tot.adj.read.count=sum(multimapAdjustedReadCount), mean.cpm.multimap.adj=mean(cpm.multimap.adjust)), by=.(mir.ID, miR.Accession, mature.fam.id, Lab, lib.method)]
mature.counts.dt.zerocounts.by.sample[, SAMPLE:=paste(Lab, lib.method, sep=":")]
unique.sample.ids <- unique(mature.counts.dt.zerocounts.by.sample$SAMPLE)
levels <- c(grep("truseq", unique.sample.ids, ignore.case=T, value=TRUE), grep("NebNext", unique.sample.ids, ignore.case=T, value=TRUE), grep("4N", unique.sample.ids, ignore.case=T, value=TRUE))
mature.counts.dt.zerocounts.by.sample[, SAMPLE:=factor(SAMPLE, levels=levels)]

ggplot(mature.counts.dt.zerocounts.by.sample, aes(x=SAMPLE, fill=lib.method, y=mean.cpm.multimap.adj)) + geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) + scale_y_log10()

# get some info on detected mirs
min.adjusted.count=1
unique.mirs.detected <- mature.counts.dt.zerocounts[multimapAdjustedReadCount>=min.adjusted.count, .(n.mirs.detected=.N, total.mir.adj.counts=sum(multimapAdjustedReadCount)), by=.(SAMPLE.ID)]

mirs.n.samples.detected <- mature.counts.dt.zerocounts[multimapAdjustedReadCount>=min.adjusted.count, .(n.samples.detected=.N, median.adjusted.cpm=median(cpm.multimap.adjust)), by=.(mir.ID, miR.Accession, mature.fam.id)]
mirs.n.samples.detected.by.lib <- mature.counts.dt.zerocounts[multimapAdjustedReadCount>=min.adjusted.count, .(n.samples.detected=.N, median.adjusted.cpm=median(cpm.multimap.adjust)), by=.(mir.ID, miR.Accession, mature.fam.id, Lab , lib.method)]
mirs.n.samples.detected.by.lib[, SAMPLE:=paste(Lab, lib.method, sep=":")]
mirs.n.samples.detected.by.lib[, SAMPLE:=factor(SAMPLE, levels=levels)]
ggplot(mirs.n.samples.detected.by.lib, aes(x=SAMPLE, fill=factor(n.samples.detected, levels=c(4,3,2,1)))) + geom_bar() + theme_bw()


mir.max.nreplicates <- mirs.n.samples.detected.by.lib[, .(n.samples.detected.bin=.N), by=.(mir.ID, miR.Accession, mature.fam.id, n.samples.detected)]
mir.id.filter <- unique(subset(mir.max.nreplicates, n.samples.detected==4 & n.samples.detected.bin>1)$mature.fam.id)


pheatmap(log10(mature.counts.dt.cast.cpm.df[mir.id.filter,]+1), show_rownames = FALSE, annotation_col = col.annot, fontsize = 14, filename = "MANUSCRIPT_FILTERED_Plasma_miRNA_heatmap.tiff")
pheatmap(cor(log10(mature.counts.dt.cast.cpm.df[mir.id.filter,]+1), method="spearman", use="pair")^2, annotation_col = col.annot, annotation_row=col.annot,  filename = "MANUSCRIPT_FILTERED_Plasma_miRNA_heatmap_Spearman_Correlation.tiff")
library(edgeR)
plotMDS(log10(mature.counts.dt.cast.cpm.df[mir.id.filter,]+1), col=)

############


# Reads in all synthetic libraries based on the metadata input table, and saves an RData file. Will be transfered to Windows OS for manipulation and plotting
library(data.table)

setwd("/data/ANALYSIS/CROSS-U01/MANUSCRIPT_DATA_CROSS-U01_SEPT2016/MANUSCRIPT_SYNTHETIC_POOL_FILES")
sequence.annot <- fread("/data/ANALYSIS/SYNTHETIC_sRNA_POOLS/20160328_CrossU01_Synthetic_RNA_comparison/DGALA_DATA/SPIKE-IN_SYNTHETIC_POOL_FULL_V2_PARSED.txt")
full.seq.info.equimolar <- fread("~/Documents/FINAL_Equimolar_pool_SynthRNA_1152-1.csv")
full.seq.info.ratiometric <- fread("~/Documents/FINAL_Ratiometric_SynthA_and_SynthB-1.csv")
setnames(full.seq.info.equimolar, "Sequence", "sequence")
setnames(full.seq.info.ratiometric, c("Sequence ID", "Sequence"), c("ratio.seqID","sequence"))
setnames(full.seq.info.ratiometric, colnames(full.seq.info.ratiometric), make.names(colnames(full.seq.info.ratiometric)))

CROSS.U01.SYNTH.METADATA <- fread("20160922_METADATA_ALL_CROSS-U01_MANUSCRIPT_FILES_SYNTHETIC_POOLS.csv")
CROSS.U01.SYNTH.METADATA[, file.exists:=ifelse(CalibratorCountsPath=="NOT_FOUND", FALSE, file.exists(CalibratorCountsPath))]
setkey(CROSS.U01.SYNTH.METADATA, CalibratorCountsPath)

counts.dt <- do.call(rbind, sapply(CROSS.U01.SYNTH.METADATA[file.exists==TRUE]$CalibratorCountsPath, USE.NAMES=FALSE, simplify=FALSE, FUN=function(x){
  dt <- fread(x); # the FileName.Loc variable is passed to sapply as the variable "x". Read this in as a data.table
  setnames(dt, 1:2, c("count", "orig.seqID")) # Column names will be V1 and V2, so set them appropriately
  dt$CalibratorCountsPath <- x # add the FileName.Loc value as a new variable in the count data table. Will be used as a key
  setkey(dt, "CalibratorCountsPath") # Set the key
  dt <- merge(CROSS.U01.SYNTH.METADATA, dt) # Use the key to merge with the CROSS_U01 table, which contains the relevant sample info
  setkey(dt, orig.seqID)
  # Gets the parsed sequence information from a new table. If all names match the corrected names (new.seqID) use that, otherwise use the old ones
  if(all(dt$orig.seqID%in%sequence.annot$orig.seqID)){
    setkey(sequence.annot, orig.seqID)      
  } else if(all(dt$orig.seqID%in%sequence.annot$new.seqID)){
    setkey(sequence.annot, new.seqID)      
  }
  return(dt[sequence.annot][!is.na(CalibratorCountsPath)]) # Return it to sapply, which will add it to the next element in the output list, which will be concatenated using rbind
}))

save.image(file="20160922_SYNTHETIC_POOL_MANUSCRIPT_DATA.RData")

# not get plasma pool
setwd("/data/ANALYSIS/CROSS-U01/MANUSCRIPT_DATA_CROSS-U01_SEPT2016/MANUSCRIPT_PLASMA_POOL_FILES/")

CROSS.U01.PLASM.METADATA <- fread("/data/ANALYSIS/CROSS-U01/MANUSCRIPT_DATA_CROSS-U01_SEPT2016/20160922_METADATA_ALL_CROSS-U01_MANUSCRIPT_FILES_PLASMA_POOLS.csv")
CROSS.U01.PLASM.METADATA[, file.exists:=ifelse(miRNASensePath=="NOT_FOUND", FALSE, file.exists(miRNASensePath))]
setkey(CROSS.U01.PLASM.METADATA, miRNASensePath)
CROSS.U01.PLASM.METADATA <- unique(CROSS.U01.PLASM.METADATA)
mir.counts.dt <- do.call(rbind, 
                         sapply( CROSS.U01.PLASM.METADATA[file.exists==TRUE]$miRNASensePath, 
                                 USE.NAMES=FALSE,
                                 simplify=FALSE,
                                 FUN=function(x){
                                   dt <- fread(x)
                                   dt[, miRNASensePath:=x]
                                   dt[, c("mir.ID", "miR.Accession"):=tstrsplit(ReferenceID, split=":", fixed=TRUE)[1:2]]
                                   return(dt)
                                 }
                         )
)

save.image(file="20160922_PLASMA_POOL_MANUSCRIPT_DATA.RData")

setkey(mir.counts.dt, miRNASensePath)

mir.counts.dt <- merge(CROSS.U01.PLASM.METADATA, mir.counts.dt)
mature.counts.dt <- subset(mir.counts.dt, substr(miR.Accession, 1, 5)=="MIMAT")
mature.counts.dt[, mature.fam.id:=tstrsplit(ReferenceID, split=":", fixed=TRUE)[5]]

# cast for zero.counts
mature.counts.dt.cast <- dcast.data.table(mature.counts.dt, mir.ID + miR.Accession + mature.fam.id ~ Lab + lib.method + replicate, fill=0, value.var="multimapAdjustedReadCount", fun.aggregate=max)
mature.counts.dt.zerocounts <- melt(mature.counts.dt.cast, value.name="multimapAdjustedReadCount", variable.name="SAMPLE.ID", id.vars=c("mir.ID", "miR.Accession", "mature.fam.id"))
mature.counts.dt.zerocounts[, total.mir.adjusted.counts:=sum(multimapAdjustedReadCount), by=.(SAMPLE.ID)]
mature.counts.dt.zerocounts[, cpm.multimap.adjust:=(multimapAdjustedReadCount*10^6)/total.mir.adjusted.counts]
mature.counts.dt.zerocounts[, c("Lab", "lib.method", "replicate"):=tstrsplit(SAMPLE.ID, split="_")]
label.vars <- c("mir.ID", "miR.Accession", "mature.fam.id")

mature.counts.dt.cast.cpm <- dcast.data.table(mature.counts.dt.zerocounts, mir.ID + miR.Accession + mature.fam.id ~ Lab + lib.method + replicate, fill=0, value.var="cpm.multimap.adjust", fun.aggregate=max)
mature.counts.dt.cast.cpm.df <- data.frame(subset(mature.counts.dt.cast.cpm, select= !colnames(mature.counts.dt.cast)%in%label.vars), row.names=mature.counts.dt.cast$mature.fam.id)
library(pheatmap)
pheatmap(log10(mature.counts.dt.cast.cpm.df+1))


mature.counts.dt.zerocounts.by.sample <- mature.counts.dt.zerocounts[, .(mir.tot.adj.read.count=sum(multimapAdjustedReadCount), mean.cpm.multimap.adj=mean(cpm.multimap.adjust)), by=.(mir.ID, miR.Accession, mature.fam.id, Lab, lib.method)]
mature.counts.dt.zerocounts.by.sample[, SAMPLE:=paste(Lab, lib.method, sep=":")]
unique.sample.ids <- unique(mature.counts.dt.zerocounts.by.sample$SAMPLE)
levels <- c(grep("truseq", unique.sample.ids, ignore.case=T, value=TRUE), grep("NebNext", unique.sample.ids, ignore.case=T, value=TRUE), grep("CleanTag", unique.sample.ids, ignore.case=T, value=TRUE),  grep("4N", unique.sample.ids, ignore.case=T, value=TRUE))
mature.counts.dt.zerocounts.by.sample[, SAMPLE:=factor(SAMPLE, levels=levels)]


ggplot(mature.counts.dt.zerocounts.by.sample, aes(x=SAMPLE, fill=lib.method, y=mean.cpm.multimap.adj)) + geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) + scale_y_log10()







# get some info on detected mirs
min.adjusted.count=1
unique.mirs.detected <- mature.counts.dt.zerocounts[multimapAdjustedReadCount>=min.adjusted.count, .(n.mirs.detected=.N, total.mir.adj.counts=sum(multimapAdjustedReadCount)), by=.(SAMPLE.ID)]

mirs.n.samples.detected <- mature.counts.dt.zerocounts[multimapAdjustedReadCount>=min.adjusted.count, .(n.samples.detected=.N, median.adjusted.cpm=median(cpm.multimap.adjust)), by=.(mir.ID, miR.Accession, mature.fam.id)]
mirs.n.samples.detected.by.lib <- mature.counts.dt.zerocounts[multimapAdjustedReadCount>=min.adjusted.count, .(n.samples.detected=.N, median.adjusted.cpm=median(cpm.multimap.adjust)), by=.(mir.ID, miR.Accession, mature.fam.id, Lab , lib.method)]
mirs.n.samples.detected.by.lib[, SAMPLE:=paste(Lab, lib.method, sep=":")]
mirs.n.samples.detected.by.lib[, SAMPLE:=factor(SAMPLE, levels=levels)]

mirs.countsummary <- mature.counts.dt.zerocounts[multimapAdjustedReadCount>=min.adjusted.count, .(total.mir.adjusted.counts=max(total.mir.adjusted.counts)), by=.(SAMPLE.ID, Lab, lib.method, replicate)]

# plot # detected
ggplot(mirs.n.samples.detected.by.lib, aes(x=SAMPLE, fill=factor(n.samples.detected, levels=c(4,3,2,1)))) + 
  geom_bar(color="black") +
  theme_bw() +
  theme(
    axis.text=element_text(size=14, face="bold"),
    axis.text.x=element_text(angle=50, hjust=1),
    axis.title=element_text(face="bold", size=16),
    plot.title=element_text(face="bold", size=18),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=14, face="bold"),
    legend.title=element_text(size=15, face="bold"),
    legend.key.size=unit(1.2, "lines")
  ) 
ggsave("MANUSCRIPT_unique_miRs_detected_by_sample.jpeg")

mirs.countsummary[, SAMPLE:=paste(Lab, lib.method, sep=":")]
mirs.countsummary[, SAMPLE:=factor(SAMPLE, levels=levels)]
ggplot(mirs.countsummary, aes(x=SAMPLE, y=total.mir.adjusted.counts, fill=replicate)) + 
  geom_bar(stat="identity", color="black", pos="dodge") +
  theme_bw() +
  theme(
    axis.text=element_text(size=14, face="bold"),
    axis.text.x=element_text(angle=50, hjust=1),
    axis.title=element_text(face="bold", size=16),
    plot.title=element_text(face="bold", size=18),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=14, face="bold"),
    legend.title=element_text(size=15, face="bold"),
    legend.key.size=unit(1.2, "lines")
  ) 
ggsave("MANUSCRIPT_Total_miR_counts_by_sample_linearY.jpeg", width = 12, height=10)

ggplot(mirs.countsummary, aes(x=SAMPLE, y=total.mir.adjusted.counts, fill=replicate)) + 
  geom_bar(stat="identity", color="black", pos="dodge") +
  theme_bw() + scale_y_log10() +
  theme(
    axis.text=element_text(size=14, face="bold"),
    axis.text.x=element_text(angle=50, hjust=1),
    axis.title=element_text(face="bold", size=16),
    plot.title=element_text(face="bold", size=18),
    panel.grid.major = element_line(color=NA),
    panel.border = element_rect(size=1.5, color="black"),
    axis.ticks = element_line(color="black", size=1.5),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text=element_text(size=14, face="bold"),
    legend.title=element_text(size=15, face="bold"),
    legend.key.size=unit(1.2, "lines")
  ) 
ggsave("MANUSCRIPT_Total_miR_counts_by_sample_logY.jpeg", width = 12, height=10)

ggplot(mirs.n.samples.detected.by.lib, aes(x=SAMPLE, fill=factor(n.samples.detected, levels=c(4,3,2,1)))) + geom_bar() + theme_bw()


mir.max.nreplicates <- mirs.n.samples.detected.by.lib[, .(n.samples.detected.bin=.N), by=.(mir.ID, miR.Accession, mature.fam.id, n.samples.detected)]
mir.id.filter <- unique(subset(mir.max.nreplicates, n.samples.detected==4 & n.samples.detected.bin>1)$mature.fam.id)
pheatmap(log10(mature.counts.dt.cast.cpm.df[mir.id.filter,]+1))
pheatmap(cor(log10(mature.counts.dt.cast.cpm.df[mir.id.filter,]+1), method="spearman", use="pair")^2)
library(edgeR)
plotMDS(log10(mature.counts.dt.cast.cpm.df[mir.id.filter,]+1), col=)
save.image(file="20160922_PLASMA_POOL_MANUSCRIPT_DATA.RData")
