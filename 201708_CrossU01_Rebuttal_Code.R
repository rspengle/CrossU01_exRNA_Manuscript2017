# Load packages -------------------------------------------------------------------
library(limma)
library(data.table)
library(ggplot2)
library(Biostrings)
library(RColorBrewer)
library(scales)
library(tools)
library(xlsx)
library(stringr)
library(edgeR)
library(pheatmap)
library(ggbeeswarm)
library(Hmisc)
library(cowplot)
library(MASS)
library(DESeq2)
library(vegan)
library(dendsort)

# Set global variables -----
# Set min and max lenghts to use (based on expected size of sequences in pool)
min.seqLen=16
max.seqLen=25
min.adjusted.count=1 # minimum adjusted count for plasma pool filtering
# Annotation files -----
full.seq.info.ratiometric.file <- "data/ratiometric.sequence.annotations.txt"
full.seq.info.equimolar.file <- "data/equimolar.sequence.annotations.txt"
synth.fasta.exceRpt.calibrator.seq.file <- "data/20160915_SPIKE-IN_SYNTHETIC_POOL_FULL_CORRECTED.fa"
CROSS.U01.METADATA.File <- "data/20170830_CrossU01_Metadata_Cleaned.txt" 
# Subset of NEB labs that were done using identical protocols
NEB.subset <- c("NEBNext.Lab1", "NEBNext.Lab2", "NEBNext.Lab4", "NEBNext.Lab5")


# FUNCTIONS --------
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

FUNCTION.collapse.IDs <- function(id){
  multiID = unlist(strsplit(id,"\\|"));
  multiIDs = sapply(
    multiID, 
    function(idPart){
      unlist(strsplit(idPart,":"))[1]
    });   
  if(length(multiIDs) == 1){
    multiIDs 
  } else{
    paste(sort(multiIDs),collapse="|") 
  }}

# Collapses filenames, delimited by a semicolon. Used for GEO submisison tables
agg.fun <- function(x){
  paste(unique(basename(x)), collapse=";") 
}

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

# Import and wrangle study data -----
full.seq.info.ratiometric <- fread(full.seq.info.ratiometric.file)
full.seq.info.equimolar <- fread(full.seq.info.equimolar.file)
synth.fasta.exceRpt.calibrator.seqs <- readDNAStringSet(synth.fasta.exceRpt.calibrator.seq.file, use.names=TRUE)
sequence.annot <- data.table(data.frame(new.seqID=names(synth.fasta.exceRpt.calibrator.seqs), sequence=as.character(synth.fasta.exceRpt.calibrator.seqs), stringsAsFactors = FALSE))
sequence.annot[, c("equimolar.seqID", "ratio.seqID", "ratio.A", "ratio.B"):=tstrsplit(new.seqID, split=";|\\|")]
sequence.annot[, colnames(sequence.annot):=lapply(.SD, function(x) ifelse(x=="NA", NA, x))]
CROSS.U01.METADATA.ALL <- fread(CROSS.U01.METADATA.File)
setnames(CROSS.U01.METADATA.ALL, "pool.ID", "pool")
CROSS.U01.METADATA.ALL[, replicate:=as.character(replicate)]
# Add a couple more grouping variables based on lab, method, pool, replicate, etc
# Make sure they're valid names as well
CROSS.U01.METADATA.ALL[, `:=`(lab.libMethod.pool.replicate=paste(lib.method.detail, Lab, pool, replicate, sep="."),
                              lab.libMethod.pool=paste(lib.method.detail, Lab, pool, sep="."),
                              lab.libMethod.replicate=paste(lib.method.detail, Lab , replicate, sep="."),
                              lab.libMethod=paste(lib.method.detail, Lab, sep="."))]

synth.files <- c("SynthA", "SynthB", "SynthEQ")
CROSS.U01.SYNTH.METADATA <- dcast.data.table(subset(CROSS.U01.METADATA.ALL, pool%in%synth.files), ...~file.type, value.var="file.path", fun.aggregate=paste, sep=";", fill=NA)
setkey(CROSS.U01.SYNTH.METADATA, "CalibratorCounts")

counts.dt <- rbindlist(sapply(CROSS.U01.SYNTH.METADATA[file.exists(CalibratorCounts)]$CalibratorCounts, USE.NAMES=FALSE, simplify=FALSE, FUN=function(x){
  dt <- fread(x); # the FileName.Loc variable is passed to sapply as the variable "x". Read this in as a data.table
  setnames(dt, 1:2, c("count", "new.seqID")) # Column names will be V1 and V2, so set them appropriately
  dt$CalibratorCounts <- x # add the FileName.Loc value as a new variable in the count data table. Will be used as a key
  setkey(dt, "CalibratorCounts") # Set the key
  dt <- merge(CROSS.U01.SYNTH.METADATA, dt) # Use the key to merge with the CROSS_U01 table, which contains the relevant sample info
  setkey(dt, new.seqID)
  # Gets the parsed sequence information from a new table. If all names match the corrected names (new.seqID) use that, otherwise use the old ones
  setkey(sequence.annot, new.seqID)   
  dt.f <- dt[sequence.annot][!is.na(CalibratorCounts)]        
  return(dt.f) # Return it to sapply, which will add it to the next element in the output list, which will be concatenated using rbind
}
), use.names=TRUE)


# Re-ran the genboree pipeline, mapping the synthetic pool sequences all the way through the genboree pipeline without using the calibrator reads. This is for direct comparison of the equimolar pools with the plasma pools.
CROSS.U01.METADATA.ALL.CAST <- dcast.data.table(CROSS.U01.METADATA.ALL, ...~file.type, value.var="file.path", fun.aggregate=paste, sep=";", fill=NA)
setkey(CROSS.U01.METADATA.ALL.CAST, miRNASense)
counts.dt.endog.miRs.all <- rbindlist(sapply(CROSS.U01.METADATA.ALL.CAST[file.exists(miRNASense)]$miRNASense, USE.NAMES=FALSE, simplify=FALSE, FUN=function(x){
  dt <- fread(x)
  dt[, miRNASense:=x]
  setkey(dt, miRNASense)
  dt <- merge(CROSS.U01.METADATA.ALL.CAST, dt)
  dt[, miR.ID:=FUNCTION.collapse.IDs(ReferenceID), by=1:nrow(dt)]
  return(dt)
}), use.names=TRUE)
counts.dt.synth.endog.miRs <- subset(counts.dt.endog.miRs.all, pool%in%synth.files)

# get stats and qc files too
excerpt.stats.genome.dt <- rbindlist(sapply(CROSS.U01.METADATA.ALL.CAST[file.exists(StatsPath.genome)]$StatsPath.genome, USE.NAMES=FALSE, simplify=FALSE, FUN=function(x){
  dt <- fread(x)
  dt[, StatsPath.genome:=x]
  setkey(dt, StatsPath.genome)
  setkey(CROSS.U01.METADATA.ALL.CAST, StatsPath.genome)
  dt <- merge(CROSS.U01.METADATA.ALL.CAST, dt)  
  dt[, runType:="NoCalibrator"]
  return(dt)
}), use.names=TRUE)
excerpt.stats.calib.dt <- rbindlist(sapply(CROSS.U01.METADATA.ALL.CAST[file.exists( StatsPath.calibrator)]$ StatsPath.calibrator, USE.NAMES=FALSE, simplify=FALSE, FUN=function(x){
  dt <- fread(x)
  dt[,  StatsPath.calibrator:=x]
  setkey(dt,  StatsPath.calibrator)
  setkey(CROSS.U01.METADATA.ALL.CAST,  StatsPath.calibrator)
  dt <- merge(CROSS.U01.METADATA.ALL.CAST, dt)
  dt[, runType:="Calibrator"]
  return(dt)
}), use.names=TRUE)

EXCERPT.STATS.PLASMA.AND.SYNTHETIC <- dcast.data.table(rbind(excerpt.stats.calib.dt, excerpt.stats.genome.dt), ...~runType, fun.aggregate=sum, fill=NA, value.var="ReadCount")


# annotate number of mods for filtering later
full.seq.info.equimolar[, n.mods.total:=.N, by=sequence]
full.seq.info.equimolar[, seq.len:=nchar(sequence)]
seq.ids.equimolar.5p.only <- subset(full.seq.info.equimolar, n.mods.total==1 & modification == "5'-phosphorylation")$equimolar.seqID


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

# Set custom plot theme ------------
my_theme_bw <- theme_bw(base_size=8)
theme_set(my_theme_bw)
my_theme_bw <- theme_update(text=element_text(color="black"), axis.text=element_text(color="black", size=rel(1.0)))

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


# MAKE OUTPUT DIRECTORIES -----------------------------------------------------------
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



# GEO Submission tables -----
CROSS.U01.METADATA.ALL.GEO <- copy(CROSS.U01.METADATA.ALL)
CROSS.U01.METADATA.ALL.GEO.synth <- subset(CROSS.U01.METADATA.ALL.GEO, pool!="PlasmaPool")
CROSS.U01.METADATA.ALL.GEO.plasma <- subset(CROSS.U01.METADATA.ALL.GEO, pool=="PlasmaPool")
CROSS.U01.METADATA.ALL.GEO.synth[, `:=`(source.name=ifelse(pool=="SynthEQ", 
                                                                "Synthetic smallRNA; Equimolar Pool",
                                                                ifelse(pool=="SynthA",
                                                                       "Synthetic smallRNA; Ratiometric Pool A",
                                                                       ifelse(pool=="SynthB",
                                                                              "Synthetic smallRNA; Ratiometric Pool B",
                                                                              "Unknown"))),
                                             organism="Synthetic",
                                             molecule="total RNA",
                                             description="",
                                             file.checksum=md5sum(file.path))]
CROSS.U01.METADATA.ALL.GEO.plasma[, `:=`(source.name="Human Plasma Pool; 11 healthy individuals",
                                             organism="Homo sapiens",
                                             sex="Male",
                                             Age="21-45 years",
                                             molecule="total RNA",
                                             description="",
                                             file.checksum=md5sum(file.path))]


CROSS.U01.METADATA.ALL.GEO.EQ.metadata <- dcast.data.table(CROSS.U01.METADATA.ALL.GEO.synth[pool=="SynthEQ" & (file.type == "miRNASense" | file.type == "CalibratorCounts")], 
                                                                lab.libMethod.pool + lab.libMethod.pool + source.name + organism + Lab + lib.method.detail + lib.method.simple + pool + molecule + description ~  file.type, 
                                                           fun.aggregate = agg.fun, fill="", value.var = "file.path")
CROSS.U01.METADATA.ALL.GEO.RATIO.metadata <- dcast.data.table(CROSS.U01.METADATA.ALL.GEO.synth[pool!="SynthEQ" & file.type == "CalibratorCounts"], 
                                                           lab.libMethod.pool + lab.libMethod.pool + source.name + organism + Lab + lib.method.detail + lib.method.simple + pool + molecule + description ~  file.type, 
                                                           fun.aggregate = agg.fun, fill="", value.var = "file.path")
CROSS.U01.METADATA.ALL.GEO.PLASMA.metadata <- dcast.data.table(CROSS.U01.METADATA.ALL.GEO.plasma[file.type == "miRNASense"], 
                                                              lab.libMethod.pool + lab.libMethod.pool + source.name + organism + sex + Age + Lab + lib.method.detail + lib.method.simple + pool + molecule + description ~  file.type, 
                                                              fun.aggregate = agg.fun, fill="", value.var = "file.path")

CROSS.U01.METADATA.ALL.GEO.synth.procfiles <- subset(CROSS.U01.METADATA.ALL.GEO.synth, file.type%in%c("miRNASense", "CalibratorCounts"))
CROSS.U01.METADATA.ALL.GEO.plasma.procfiles <- subset(CROSS.U01.METADATA.ALL.GEO.plasma, file.type%in%c("miRNASense"))
CROSS.U01.METADATA.ALL.GEO.synth.procfiles[, file.type.geo:=ifelse(file.type=="miRNASense", 
                                                                   "Mature miRNA Read Counts",
                                                                   "Mapped Calibrator Counts")]
CROSS.U01.METADATA.ALL.GEO.plasma.procfiles[, file.type.geo:=ifelse(file.type=="miRNASense", 
                                                                   "Mature miRNA Read Counts",
                                                                   "Mapped Calibrator Counts")]

CROSS.U01.METADATA.ALL.GEO.EQ.procfiles <- subset(CROSS.U01.METADATA.ALL.GEO.synth.procfiles, pool=="SynthEQ", select=c("file.path", "file.type.geo", "file.checksum"))
CROSS.U01.METADATA.ALL.GEO.RATIO.procfiles <- subset(CROSS.U01.METADATA.ALL.GEO.synth.procfiles, pool!="SynthEQ" & file.type=="CalibratorCounts", select=c("file.path", "file.type.geo", "file.checksum"))
CROSS.U01.METADATA.ALL.GEO.PLASMA.procfiles <- subset(CROSS.U01.METADATA.ALL.GEO.plasma.procfiles, select=c("file.path", "file.type.geo", "file.checksum"))

# copy the files to our geo subdirectory
old.file.paths <- c(CROSS.U01.METADATA.ALL.GEO.EQ.procfiles$file.path, CROSS.U01.METADATA.ALL.GEO.RATIO.procfiles$file.path, CROSS.U01.METADATA.ALL.GEO.PLASMA.procfiles$file.path)
new.file.paths <- paste0(outdirs["geo_processed_data_files"], "/", basename(old.file.paths))
file.copy(old.file.paths, new.file.paths, overwrite=TRUE)
CROSS.U01.METADATA.ALL.GEO.EQ.procfiles[, file.path:=basename(file.path)]
CROSS.U01.METADATA.ALL.GEO.RATIO.procfiles[, file.path:=basename(file.path)]
CROSS.U01.METADATA.ALL.GEO.PLASMA.procfiles[, file.path:=basename(file.path)]

# Make the FASTQ file tabe that will be filled in later by the DMRR, since FASTQ files were located on their remote server
# ExceRpt generates sample names from the original file names. Use those to get the original fastq file names.
# Missing info will be filled in manually. Leave as X for now
CROSS.U01.METADATA.ALL.GEO.fastqs <- CROSS.U01.METADATA.ALL.GEO[, .(file.name=unique(sub("_fastq", ".fastq.gz", sub("sample_", "", FileBaseName))), file.type="fastq", file.checksum="X", instrument.model="X", read.length="X", single.or.paired.end="single"), by=.(FileBaseName, pool)]
this.outfile <- paste0(outdirs["geo_metadata_files"], "/geo_metadata_information.xlsx")
write.xlsx(CROSS.U01.METADATA.ALL.GEO.fastqs, this.outfile, sheetName = "ALL_raw_files", row.names = FALSE)
write.xlsx(CROSS.U01.METADATA.ALL.GEO.EQ.procfiles, this.outfile, sheetName = "EQ_processed_files", row.names = FALSE, append = TRUE)
write.xlsx(CROSS.U01.METADATA.ALL.GEO.RATIO.procfiles, this.outfile, sheetName = "RATIO_processed_files", row.names = FALSE, append = TRUE)
write.xlsx(CROSS.U01.METADATA.ALL.GEO.PLASMA.procfiles, this.outfile, sheetName = "PLASMA_processed_files", row.names = FALSE, append = TRUE)
write.xlsx(CROSS.U01.METADATA.ALL.GEO.EQ.metadata, this.outfile, sheetName = "EQ_SAMPLES", row.names = FALSE, append = TRUE)
write.xlsx(CROSS.U01.METADATA.ALL.GEO.RATIO.metadata, this.outfile, sheetName = "RATIO_SAMPLES", row.names = FALSE, append = TRUE)
write.xlsx(CROSS.U01.METADATA.ALL.GEO.PLASMA.metadata, this.outfile, sheetName = "PLASMA_SAMPLES", row.names = FALSE, append = TRUE)


# EQUIMOLAR WRANGLING -------
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

# ExceRpt QC Metrics ----
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

# MORE EQUIMOLAR WRANGLING -------
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

# KEY VARIABLES FROM EQUIMOLAR POOL WRANGLING ----
# equimolar.counts --> full results from equimolar pool. Unfiltered. Includes zerocounts. 
# equimolar.counts.5p.only --> Full equimolar pool. Filtered for 5'p only.
# FINAL.EQUIMOLAR.LONG --> Further filtered for size contraints used in paper. All replicates included. 
# FINAL.EQUIMOLAR.LONG.AVERAGE.REPS.FILT --> Averaged multiple things over replicates
# equimolar.counts.5p.only.SizeFilt.cast.df.counts --> data frame of count matrix, including all replicates
# equimolar.counts.5p.only.SizeFilt.cast.df.pseudo.cpm --> data frame of cpm, including all replicates
# sample.info.eq --> sample info for equimolar pool in data.frame format
# gene.info.eq --> gene info for equimolar pool in data.frame format

# FIG2 -----------------
# HEATMAP -- Unnormalized cpm +1 pseudo count
equimolar.counts.5p.only.SizeFilt.cast.df.pseudo.cpm.log2 <- log2(equimolar.counts.5p.only.SizeFilt.cast.df.pseudo.cpm)
colnames(equimolar.counts.5p.only.SizeFilt.cast.df.pseudo.cpm.log2) <- sub("^X", "", colnames(equimolar.counts.5p.only.SizeFilt.cast.df.pseudo.cpm.log2))
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
  ) ; print(f2b)
ggsave(this.outfile, width=4.6, height=2.75, units = "in")


## Now use a cutoff to get miRs > x fold from expected and see what % miRs in the different methods
n.seqs <- length(unique(FINAL.EQUIMOLAR.LONG$equimolar.seqID))
cutoff.from.exp = log2(10)
FINAL.EQUIMOLAR.LONG.count.simple.expected <- FINAL.EQUIMOLAR.LONG[, .(n.outside.cutoff=sum(ifelse(abs(logratio.countplus1.filt.vs.expected)>=cutoff.from.exp, 1, 0)), n.total=.N), by=.(lab.libMethod, lab.libMethod.replicate, lib.method.detail, lib.method.simple, Lab)]
FINAL.EQUIMOLAR.LONG.count.simple.expected[, percent.outside.cutoff:=n.outside.cutoff/n.total]
summary.from.expected.eq <- FINAL.EQUIMOLAR.LONG.count.simple.expected[, .(min=min(percent.outside.cutoff), max=max(percent.outside.cutoff), median=median(percent.outside.cutoff), mean=mean(percent.outside.cutoff)), .(by=lib.method.detail)]

# Table with median % sequences more than 10x from expectation for in-text values ----
median.percent.10x.expected.eq.by.method <- FINAL.EQUIMOLAR.LONG.count.simple.expected[, .(median.percent.10x.exp=median(percent.outside.cutoff)), by=lib.method.detail]
this.outfile <- paste0(outdirs["tables"], "/TextSummary_Fig2C_median_percent_10x_expected_eqPool.txt")
write.table(median.percent.10x.expected.eq.by.method, file = this.outfile, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

# FIG2C ---------
# NOTE: Had to do minor changes in Illustrator to horizontally shift point/range plots where the ranges overlapped eachother
this.outfile <- paste0(outdirs["FIG2"], "/FIG2C_SIMPLE_10x_expected_labSeparate.pdf")
g <- ggplot(FINAL.EQUIMOLAR.LONG.count.simple.expected, 
       aes(
         x=lib.method.detail,
         y=percent.outside.cutoff, 
         color=lib.method.simple,
         pos=lab.libMethod)) + 
  stat_summary(
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

# FIG S3 Spearman Correlation Heatmaps for Eq Pool -----------
# For visualization on screen
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
# Output without legend for correct sizing
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
# Output with legend so legend can be added to Illustrator file
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

# Equimolar Correlation calculations ----
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
n.miRs <- dim(cpm.rle)[1]
cor.matrix.all.eq.sample.info.cor.summary.details <- cor.matrix.all.eq.sample.info.unique[, { nr=length(unique(c(lab.libMethod.replicate.A, lab.libMethod.replicate.B))); FUNCTION.meanRhoFromCoefs(rho, ns=n.miRs, nr=nr) }, by=comparison.detail]
cor.matrix.all.eq.sample.info.cor.summary.simple <- cor.matrix.all.eq.sample.info.unique[, { nr=length(unique(c(lab.libMethod.replicate.A, lab.libMethod.replicate.B))); FUNCTION.meanRhoFromCoefs(rho, ns=n.miRs, nr=nr) }, by=comparison.simple]
setnames(cor.matrix.all.eq.sample.info.cor.summary.details, "comparison.detail", "comparison")
setnames(cor.matrix.all.eq.sample.info.cor.summary.simple, "comparison.simple", "comparison")
cor.matrix.all.eq.sample.info.cor.summary.all <- rbind(cor.matrix.all.eq.sample.info.cor.summary.simple, cor.matrix.all.eq.sample.info.cor.summary.details)

# meanRho Summary Equimolar: All replicates, but no comparison across technical replicates----
cor.matrix.all.eq.sample.info.unique[, `:=`(lab.libMethod.A=sub(".[0-9]$", "", lab.libMethod.replicate.A), lab.libMethod.B=sub(".[0-9]$", "", lab.libMethod.replicate.B)) ]
cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison <- cor.matrix.all.eq.sample.info.unique[lab.libMethod.A!=lab.libMethod.B]

cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.details <- cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison[, { nr=length(unique(c(lab.libMethod.replicate.A, lab.libMethod.replicate.B))); FUNCTION.meanRhoFromCoefs(rho, ns=n.miRs, nr=nr) }, by=comparison.detail]
cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.simple <- cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison[, { nr=length(unique(c(lab.libMethod.replicate.A, lab.libMethod.replicate.B))); FUNCTION.meanRhoFromCoefs(rho, ns=n.miRs, nr=nr) }, by=comparison.simple]
setnames(cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.details, "comparison.detail", "comparison")
setnames(cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.simple, "comparison.simple", "comparison")
cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.all <- unique(rbind(cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.details, cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.simple), by="comparison")

# Summary Table: Equimolar Pool Correlation Coefs ----
this.outfile <- paste0(outdirs["tables"], "/Spearman_correlation_equimolar_Pool_summaries_noTechRep_comparison.txt")
write.table(cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.all, this.outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

comparison.levels.order <- unique(as.character(melt(outer(lib.detail.levels, lib.detail.levels, FUN = paste, sep=".VS."))$value))
comparison.levels.order <- comparison.levels.order[comparison.levels.order%in%cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.simple$comparison]
cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.all[, comparison:=factor(comparison, levels=comparison.levels.order)]
setorder(cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.all, comparison)
# IN-TEXT STATS: INTRA-Protocol Equimolar Correlation----
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

# IN-TEXT STATS: INTER-Protocol Equimolar Correlation----
setkey(cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.all, comparison)
gsub(".VS.", " vs ", paste(cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.all[Cs(TruSeq.VS.NEBNext, TruSeq.VS.CleanTag, TruSeq.VS.4N_B, CleanTag.VS.NEBNext, NEBNext.VS.4N_B, CleanTag.VS.4N_B), paste0(comparison, ": ", round(simple.mean.r, digits = 2), " (", round(r02, digits = 2), ", ", round(r98, digits = 2), "; n=", nlibs,  ")")], collapse="; "))

# Table: Equimolar Correlation Coefficents ----
this.outfile <- paste0(outdirs["tables"], "/Spearman_correlation_equimolar_Pool_summaries.txt")
write.table(cor.matrix.all.eq.sample.info.unique.no.tech.rep.comparison.cor.summary.all, this.outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)


# FIG4A: Equimolar Intralab %CV ----------------
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


# FIG4A-- Equimolar Intralab QCD --------------------------------------------------------------
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


# Intralab CV Equimolar Summary Table -----------
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

# IN-TEXT STATS: Equimolar INTRAlab CV/QCD summary numbers----
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
# FIG4b--Equimolar INTERlab %CV --------------------------------------------------------------
MEAN.CV.WI.LABS.drop.NEB.1to2.dilutions <- subset(MEAN.CV.WI.LABS, lib.method.detail!="NEBNext" | lab.libMethod%in%NEB.subset)
MEAN.CV.ACROSS.LABS <- MEAN.CV.WI.LABS.drop.NEB.1to2.dilutions[(lib.method.detail=="4N_B" | lib.method.simple=="TruSeq" | lib.method.simple == "NEBNext"), .(cv.inter=sd(mean.cpm)/(sum(mean.cpm*n)/sum(n))), by=.(lib.method.simple, lib.method.detail, equimolar.seqID)]
this.outfile <- paste0(outdirs["FIG4"], "/FIG4B_EQUIMOLAR_PERCENT_CV_INTERLAB_TS_NEB_4N_B_NEB1to2DilutionsOnly.pdf")
ggplot(MEAN.CV.ACROSS.LABS,
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

# FIG4b--Equimolar INTERlab QCD --------------------------------------------------------------
QCD.ACROSS.LABS <- MEAN.CV.WI.LABS.drop.NEB.1to2.dilutions[(lib.method.detail=="4N_B" | lib.method.simple=="TruSeq" | lib.method.simple == "NEBNext"), .(q1=quantile(mean.cpm, 0.25), q3=quantile(mean.cpm, 0.75)), by=.(lib.method.simple, lib.method.detail, equimolar.seqID)]
QCD.ACROSS.LABS[, `:=`(iqr=(q3-q1)/2, midhinge=(q1+q3)/2)]
QCD.ACROSS.LABS[, qcd:=iqr/midhinge]
this.outfile <- paste0(outdirs["FIG4"], "/FIG4B_EQUIMOLAR_PERCENT_QCD_INTERLAB_TS_NEB_4N_B_NEB1to2DilutionsOnly.pdf")
ggplot(QCD.ACROSS.LABS,
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

# IN-TEXT STATS: INTER-lab equimolar CV/QCD summary numbers----
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
# "4N_B: 0.13 (0.03, 0.42; n=4); NEBNext_subset: 0.18 (0.04, 0.46; n=4); TruSeq: 0.18 (0.06, 0.45; n=8); NEBNext: 0.25 (0.06, 0.56; n=6); 4N: 0.3 (0.07, 0.82; n=7)"

# Interlab-CV Equimolar All 3V3 Subsets ----
# Reviewer requested we perform all combinations of 3 labs for inter-lab measures
comparisons.dt <- data.frame(MEAN.CV.WI.LABS.subgroups[n>2, .(lab.libMethod=unique(lab.libMethod)), by=c("lib.method.group", "lib.method.detail")])

samples.by.method <- split(as.character(comparisons.dt$lab.libMethod), f = comparisons.dt$lib.method.group, drop=TRUE)
all.combinations <- lapply(lapply(samples.by.method, combn, m=3), t)
selected.groups <- c("TruSeq", "NEBNext_subset", "4N_B")
selected.combinations <- all.combinations[selected.groups]
setkey(MEAN.CV.WI.LABS.subgroups, lab.libMethod)
selected.combin.full.inter.cv <- rbindlist(
                      sapply(1:length(selected.combinations),
                                  simplify=FALSE, 
                                  FUN=function(x){
                                    s.group <- names(selected.combinations)[x]
                                    all.comb.matr <- selected.combinations[[x]]
                                    dt.wi.labs <- subset(MEAN.CV.WI.LABS.subgroups, lib.method.group==s.group)
                                    setkey(dt.wi.labs, lab.libMethod)
                                    dt.betw <- rbindlist(
                                      sapply(1:dim(all.comb.matr)[1], simplify=FALSE, FUN=function(y){
                                        this.set <- all.comb.matr[y,]
                                        set.cat <- paste(this.set, collapse="|")
                                        dt.s <- dt.wi.labs[this.set, .(mean.cpm=mean(mean.cpm), sd.cpm=sd(mean.cpm), q1.cpm=quantile(mean.cpm, 0.25), q3.cpm=quantile(mean.cpm, 0.75), n=.N), by=.(equimolar.seqID, lib.method.group, lib.method.detail)]
                                      dt.s[, `:=`(interlab.cv=100*sd.cpm/mean.cpm, interlab.qcd=((q3.cpm-q1.cpm)/2)/((q3.cpm+q1.cpm)/2))]
                                      dt.s[, set.combin:=set.cat]
                                      if(y==1){
                                        set.cat <- dt.wi.labs[, paste(unique(lab.libMethod), collapse="|")]
                                        dt.full <- dt.wi.labs[, .(mean.cpm=mean(mean.cpm), sd.cpm=sd(mean.cpm), q1.cpm=quantile(mean.cpm, 0.25), q3.cpm=quantile(mean.cpm, 0.75), n=.N), by=.(equimolar.seqID, lib.method.group, lib.method.detail)]
                                        dt.full[, `:=`(interlab.cv=100*sd.cpm/mean.cpm, interlab.qcd=((q3.cpm-q1.cpm)/2)/((q3.cpm+q1.cpm)/2))]
                                        dt.full[, set.combin:=set.cat]
                                        dt.s <- rbind(dt.full, dt.s)
                                      }
                                      return(dt.s)
                                    }))
                                    
                                    
                                    
                                    return(dt.betw)
}))

med.cv.dt <- selected.combin.full.inter.cv[n==3, .(med.cv=median(interlab.cv)), by=.(set.combin)]
setorder(med.cv.dt, -med.cv)
full.set.ids <- selected.combin.full.inter.cv[n>3, unique(as.character(set.combin))]
#full.set.labs <- paste0("ALL:n=", nchar(full.set.ids)-nchar(gsub("Lab", "La", full.set.ids)))
full.set.labs <- paste0("Labs:", gsub("[A-Za-z]|\\.|_", "", gsub("4N", "", full.set.ids)))
cv.plot.order <- c(med.cv.dt$set.combin, full.set.ids)
cv.plot.labels <- c(paste0("Labs:", gsub("[A-Za-z]|\\.|_", "", gsub("4N", "", med.cv.dt$set.combin))), full.set.labs)
names(cv.plot.labels) <- cv.plot.order
selected.combin.full.inter.cv[, set.combin:=factor(set.combin, levels=cv.plot.order)]
# SUPP FIG S6a: Equimolar CV combinations -----
g <- ggplot(selected.combin.full.inter.cv, aes(x=set.combin, y=interlab.cv, fill=lib.method.detail)) + 
    geom_boxplot(outlier.alpha = 0.5, size=0.5) + 
  coord_flip() +
  scale_y_continuous("Interlab %CV", limits=c(0,280), breaks=seq(0,300,25)) +
  scale_x_discrete(labels=cv.plot.labels) + labs(x=NULL) +
  scale_fill_manual(values=ann_colors$lib.method.detail) +
  theme(
    panel.grid.major = element_line(color = NA),
    panel.grid.minor = element_line(color = NA),
    panel.border = element_rect(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  )
g1 <- g %+% selected.combin.full.inter.cv[lib.method.group=="TruSeq"]; g1
g2 <- g %+% selected.combin.full.inter.cv[lib.method.group=="NEBNext_subset"]; g2
g3 <- g %+% selected.combin.full.inter.cv[lib.method.detail=="4N_B"]; g3
g23 <- plot_grid(plotlist = list(g2, g3), ncol = 1, align = "v", rel_heights = c(1, 1))
g123 <- plot_grid(plotlist = list(g1, g23), ncol=2, align="h", rel_widths=c(1.25,1)); g123
this.outdir <- paste0(outdirs["supplemental_figures"], "/FIGS6A_Equimolar_CV_All3wayCombinations_Boxplot.pdf")
save_plot(this.outdir, g123, base_width = 7.5, base_height=8.5)

# SUPP FIG S6b: Equimolar CV combinations -----
g <- ggplot(selected.combin.full.inter.cv, aes(x=set.combin, y=interlab.qcd, fill=lib.method.detail)) + 
  geom_boxplot(outlier.alpha = 0.5, size=0.5) + 
  coord_flip() +
  scale_y_continuous("Interlab QCD", limits=c(0,1), breaks=seq(0, 1, 0.1)) +
  scale_x_discrete(labels=cv.plot.labels) + labs(x=NULL) +
  scale_fill_manual(values=ann_colors$lib.method.detail) +
  theme(
    panel.grid.major = element_line(color = NA),
    panel.grid.minor = element_line(color = NA),
    panel.border = element_rect(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.position = "none"
  )
g1 <- g %+% selected.combin.full.inter.cv[lib.method.group=="TruSeq"]; g1
g2 <- g %+% selected.combin.full.inter.cv[lib.method.group=="NEBNext_subset"]; g2
g3 <- g %+% selected.combin.full.inter.cv[lib.method.detail=="4N_B"]; g3
g23 <- plot_grid(plotlist = list(g2, g3), ncol = 1, align = "v", rel_heights = c(1, 1))
g123 <- plot_grid(plotlist = list(g1, g23), ncol=2, align="h", rel_widths=c(1.25,1)); g123
this.outdir <- paste0(outdirs["supplemental_figures"], "/FIGS6B_Equimolar_QCD_All3wayCombinations_Boxplot.pdf")
save_plot(this.outdir, g123, base_width = 7.5, base_height=8.5)

# SUP Table 4: Number of undetected sequences in the equimolar pool -----
total.eq.seqs <- length(unique(FINAL.EQUIMOLAR.LONG$equimolar.seqID))
FINAL.EQUIMOLAR.LONG[, `:=`(n.replicates.detected=sum(ifelse(count>0, 1, 0)), total.replicates=.N, mean.lib.size=exp(mean(log( count.total.by.sample.filt.lengths)))), by=.(lab.libMethod, Lab, lib.method.detail, lib.method.simple, equimolar.seqID, sequence, seq.len, gc.perc)]
FINAL.EQUIMOLAR.LONG[, mean.lib.size:=ifelse(n.replicates.detected==total.replicates, mean.lib.size, 0)]
FINAL.EQUIMOLAR.LONG[, mean.lib.size:=max(mean.lib.size), by=.(lab.libMethod)]
EQUIMOLAR.MISSING.MIR.SUMMARY.BY.MIR <- FINAL.EQUIMOLAR.LONG[, .(mean.lib.size=max(mean.lib.size)), by=.(lab.libMethod, Lab, lib.method.detail, lib.method.simple, equimolar.seqID, sequence, seq.len, gc.perc, n.replicates.detected, total.replicates)]
EQUIMOLAR.MISSING.MIR.SUMMARY <- EQUIMOLAR.MISSING.MIR.SUMMARY.BY.MIR[, .(n.miRs=.N), by=.(lab.libMethod, Lab, lib.method.detail, lib.method.simple, mean.lib.size, n.replicates.detected, total.replicates)]
EQUIMOLAR.MISSING.MIR.SUMMARY.drop.nonMissing <- EQUIMOLAR.MISSING.MIR.SUMMARY[n.miRs==total.eq.seqs | n.replicates.detected<total.replicates]
EQUIMOLAR.MISSING.MIR.SUMMARY.total.missing.in.any <- EQUIMOLAR.MISSING.MIR.SUMMARY[, .(n.miRs.missing=sum(ifelse(n.replicates.detected<total.replicates, n.miRs, 0))), by=.(lab.libMethod, Lab, lib.method.detail, lib.method.simple, mean.lib.size, total.replicates)]
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

FINAL.EQUIMOLAR.LONG.25nt.disp.replicates.top.10 <- subset(FINAL.EQUIMOLAR.LONG.AVERAGE.REPS.FILT,  equimolar.seqID%in%mirs.with.75.perc.agreement.in.top.10 & !lab.libMethod%in%unique(gsub(".[1-4]$|^X", "", drop.eq.samples)))
FINAL.EQUIMOLAR.LONG.25nt.disp.replicates.top.10.mean <- FINAL.EQUIMOLAR.LONG.25nt.disp.replicates.top.10[, .(lib.meanCPM=mean(mean.pseudo.cpm)), by=.(lib.method.simple, equimolar.seqID)]

# now, instead, divide by the expected value to get a difference from expected
FINAL.EQUIMOLAR.LONG.25nt.disp.replicates.top.10.mean[, `:=`(meanCPM.from.expected=(lib.meanCPM-expected.cpm)/expected.cpm, log.CPM.from.expected=log2(lib.meanCPM/expected.cpm))]

# SUPP FIG S1a: Over-represented miRs---------------------
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
    axis.ticks.y=element_line(color="black"),
    axis.line.x=element_line(color="black"),
    axis.line.y=element_line(color="black"),
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
  ) + geom_hline(yintercept = 0, size=1.5); g
ggsave(plot = g, filename = this.outfile, width=5, height=4, units="in")

# Bottom miRs by protocol
rank.matrix.all.down <- apply(cpm.lrt.filt, 2, rank, ties.method = "min")
in.bottom.ten <- ifelse(rank.matrix.all.down<=10, 1, 0)
n.samples <- ifelse(rank.matrix.all.down<=10, 1, 1)
in.bottom.ten.agreement <- sumTechReps(in.bottom.ten, ID = sample.info.eq[colnames(in.top.ten),]$lib.method.simple)
n.samples.summary <- sumTechReps(n.samples, ID = sample.info.eq[colnames(in.top.ten),]$lib.method.simple)
in.bottom.ten.agreement.perc <- in.bottom.ten.agreement/n.samples.summary
mirs.with.75.perc.agreement.in.bottom.10 <- names(which(apply(in.bottom.ten.agreement.perc, 1, max)>=0.75))

FINAL.EQUIMOLAR.LONG.25nt.disp.replicates.bottom.10 <- subset(FINAL.EQUIMOLAR.LONG.AVERAGE.REPS.FILT,  equimolar.seqID%in%mirs.with.75.perc.agreement.in.bottom.10 & !lab.libMethod%in%unique(gsub(".[1-4]$|^X", "", drop.eq.samples)))
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
    ) + geom_hline(yintercept = 0, size=1.5); g
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
FINAL.RATIOMETRIC.LONG.AVERAGE.REPS.CAST.AvsB[, ratio.AvsB.Measured:=ifelse(mean.cpm_SynthA==0 & mean.cpm_SynthB==0, 0, (mean.cpm_SynthA)/(mean.cpm_SynthB))]
FINAL.RATIOMETRIC.LONG.AVERAGE.REPS.CAST.AvsB[, RATIO.AvsB.Expected:=factor(paste(ratio.A, ratio.B, sep=":"), levels=ratio.order)]

# FIG3A -------------------------------------------------------------------
# Specify new order to fit into grid in a semi-resonable way
group.split <- split(new.order$ratiometric, strsplit2(new.order$ratiometric, split=".", fixed=TRUE)[,1])
new.order.ratio <- c(group.split[["TruSeq"]][1:2],
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


FINAL.RATIOMETRIC.LONG.AVERAGE.REPS.CAST.AvsB[, lab.libMethod2:=factor(lab.libMethod, levels=new.order.ratio)]
ratio.labels <- ratio.order.fract
names(ratio.labels) <- ratio.order
undetected.mirs <- FINAL.RATIOMETRIC.LONG.AVERAGE.REPS.CAST.AvsB[, .(n.missing=sum(ifelse(mean.cpm_SynthA==0 | mean.cpm_SynthB==0, 1, 0))), by=.(lab.libMethod2)]
g <- ggplot(FINAL.RATIOMETRIC.LONG.AVERAGE.REPS.CAST.AvsB[mean.cpm_SynthA>0 & mean.cpm_SynthB>0], 
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


# Differential Expression Ratio Pool --------------------------
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
all.ratio.topTables.both[, lab.libMethod2:=factor(lab.libMethod, levels=new.order.ratio)]
# SUPP FIG S4: Ratiometric Differential Expression Barplot ----------------------------------
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


# Table with ratiometric pool DE Calls
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
# TABLE S5: Ratiometric Limma/Voom Calls
this.outfile <- paste0(outdirs["tables"], "/TABLE_S5_RatioPool_limma_voom_calls_obs_vs_exp.xlsx")
write.xlsx(all.ratio.topTables.limma.summary, this.outfile, row.names=FALSE, col.names=TRUE)

# COR HEATMAPS FOR MAN

# FIG3B: Ratiometric Correlation HEatmap -------------------------------------------------------------------
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
                   breaks=seq(0,1, length.out = 101),
                   treeheight_row=10,
                   treeheight_col=10,
                   width=3,
                   height=3,
                   clustering_distance_rows = d.mean,
                   clustering_distance_cols = d.mean)$gtable
ratio.cor.heatmap <- plot_grid(plotlist = list(pA, pB, pRatio), nrow=1)
save_plot(this.outfile, ratio.cor.heatmap, ncol = 3, nrow = 1, base_height = 3, base_width = 2.5)
# Contains all the annotation/label info for adding later
this.outfile <- paste0(outdirs["FIG3"], "/FULL_ANNOT_FIG3B_SPEARMAN_COR_HEATMAP_CPM_unnormalized_with_labels.pdf")
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

# Ratiometric Correlation Summary Table ----
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

# INTERLAB CV/QCD Ratiometric Summary Table -------------------------------
INTRALAB.CV.RATIOMETRIC.subgroups <- subset(INTRALAB.CV.RATIOMETRIC, lib.method.detail!="NEBNext" | lab.libMethod%in%NEB.subset)
INTRALAB.CV.RATIOMETRIC.subgroups[, n:=.N, by=.(pool, lib.method.detail, ratio.seqID)]

MEAN.CV.ACROSS.LABS.RATIO.subgroups <- INTRALAB.CV.RATIOMETRIC.subgroups[n>2, .(mean.cpm=mean(mean.cpm), sd.cpm=sd(mean.cpm), q1.cpm=quantile(mean.cpm, 0.25), q3.cpm=quantile(mean.cpm, 0.75), n=.N), by=.(pool, lib.method.detail, ratio.seqID)]
MEAN.CV.ACROSS.LABS.RATIO.subgroups[, `:=`(interlab.cv=100*sd.cpm/mean.cpm, interlab.qcd=((q3.cpm-q1.cpm)/2)/((q3.cpm+q1.cpm)/2))]
MEAN.CV.ACROSS.LABS.RATIO.subgroups[, (sub("^0.", "interlab.cv", summary.quantiles)):=lapply(summary.quantiles, USE.NAMES=FALSE, SIMPLIFY=FALSE, quantile, x=interlab.cv), by=.(pool, lib.method.detail)]
MEAN.CV.ACROSS.LABS.RATIO.subgroups[, (sub("^0.", "interlab.qcd", summary.quantiles)):=lapply(summary.quantiles, USE.NAMES=FALSE, SIMPLIFY=FALSE, quantile, x=interlab.qcd), by=.(pool, lib.method.detail)]

INTERLAB.CV.RATIOMETRIC.SUMMARY <- subset(unique(MEAN.CV.ACROSS.LABS.RATIO.subgroups, by=c("lib.method.detail", "pool", "n")), select=c("lib.method.detail", "pool", "n", sub("^0.", "interlab.cv", summary.quantiles), sub("^0.", "interlab.qcd", summary.quantiles)))
this.outfile <- paste0(outdirs["tables"], "/Ratiometric_interlab_CV_QCD_summary.txt")
write.table(INTERLAB.CV.RATIOMETRIC.SUMMARY, this.outfile, row.names = FALSE, col.names=TRUE, sep="\t", quote=FALSE)

# IN-TEXT STATS: INTER-lab RATIOMETRIC CV/QCD summary numbers----

paste(
  INTERLAB.CV.RATIOMETRIC.SUMMARY[, 
                                paste0(lib.method.detail, " ", pool,
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
# "4N_B SynthA: 32.2 (7.89, 87.05; n=4); 4N_B SynthB: 33.67 (8.06, 85.24; n=4); NEBNext SynthA: 36.28 (11.69, 66.68; n=4); NEBNext SynthB: 44.33 (15.24, 73.92; n=4); TruSeq SynthA: 30.81 (14.48, 66.87; n=8); TruSeq SynthB: 32.09 (14.78, 84.25; n=8)"

paste(
  INTERLAB.CV.RATIOMETRIC.SUMMARY[, 
                                  paste0(lib.method.detail, " ", pool,
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
# 4N_B SynthA: 0.15 (0.04, 0.41; n=4); 4N_B SynthB: 0.16 (0.04, 0.45; n=4); NEBNext SynthA: 0.19 (0.06, 0.43; n=4); NEBNext SynthB: 0.22 (0.07, 0.51; n=4); TruSeq SynthA: 0.17 (0.05, 0.43; n=8); TruSeq SynthB: 0.18 (0.05, 0.43; n=8)


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

# IN-TEXT STATS: ratiometric correlation summary numbers----
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

# IN-TEXT STATS: inter-protocol comparisons
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
#g <- ggplot(cor.matrix.all.ratio.sample.info.cor.summary.all.dupl.filt, aes(x=Type)) +
#geom_pointrange(aes(y=simple.mean.r, ymin=r02, ymax=r98)) + facet_grid(MethodA~MethodB); g

# DOWNSAMPLING: Equimolar ------------------------------
dge.eq.counts <- dge.eq$counts
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

# FIG5A: Plasma Pool CPM HEatmap -------------------------------------------------------------------
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


# Plasma Intra- and Inter-lab variation ----
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

# IN-TEXT STATS: PLASMA Intralab %CV summary numbers----
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
# IN-TEXT STATS: PLASMA Intralab QCD summary numbers----
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

# FIG5B--%CV -------------------------------------------------------------------
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


# FIG5B--QCD --------------------------------------------------------------
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

# Plasma Cross-Lab %CV ----
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
MEAN.CV.ACROSS.LABS.plasma.subgroups <- MEAN.CV.WI.LABS.plasma.subgroups[, .(mean.cpm=mean(mean.cpm), sd.cpm=sd(mean.cpm), q1.cpm=quantile(mean.cpm, 0.25), q3.cpm=quantile(mean.cpm, 0.75), n=.N), by=.(lib.method.group, miR.ID)][n>2,]
MEAN.CV.ACROSS.LABS.plasma.subgroups[, `:=`(interlab.cv=100*sd.cpm/mean.cpm, interlab.qcd=((q3.cpm-q1.cpm)/2)/((q3.cpm+q1.cpm)/2))]
MEAN.CV.ACROSS.LABS.plasma.subgroups[, (sub("^0.", "interlab.cv", summary.quantiles)):=lapply(summary.quantiles, USE.NAMES=FALSE, SIMPLIFY=FALSE, quantile, x=interlab.cv), by=.(lib.method.group)]
MEAN.CV.ACROSS.LABS.plasma.subgroups[, (sub("^0.", "interlab.qcd", summary.quantiles)):=lapply(summary.quantiles, USE.NAMES=FALSE, SIMPLIFY=FALSE, quantile, x=interlab.qcd), by=.(lib.method.group)]

INTERLAB.CV.PLASMA.SUMMARY <- subset(unique(MEAN.CV.ACROSS.LABS.plasma.subgroups, by=c("lib.method.group", "n")), select=c("lib.method.group", "n", sub("^0.", "interlab.cv", summary.quantiles), sub("^0.", "interlab.qcd", summary.quantiles)))
this.outfile <- paste0(outdirs["tables"], "/Plasma_interlab_CV_QCD_summary.txt")
write.table(INTERLAB.CV.PLASMA.SUMMARY, this.outfile, row.names = FALSE, col.names=TRUE, sep="\t", quote=FALSE)

# IN-TEXT STATS: Plasma INTER-lab CV/QCD summary numbers----
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
est.prob.detected.n.miRs <- drr.all.plasma.filt[, .(n.detected=sum(ifelse(est.prob.detected>=min.prob.detected, 1, 0)), n.not.detected=sum(ifelse(est.prob.detected>=min.prob.detected, 0, 1)), n.total=.N), by=.(lab.libMethod.replicate, lab.libMethod, downsample.to, Lab, lib.method.detail, lib.method.simple)]

# Figure 5D: Plasma Pool miRs detected ------
this.outfile <- paste0(outdirs["FIG5"], "/FIG5D_miRs_detected_plasmaPool.pdf")
g <- ggplot(est.prob.detected.n.miRs,
       aes(
         x=lib.method.simple,
         y=n.detected, fill=lib.method.simple)) + 
  geom_boxplot() + 
  facet_wrap(~downsample.to, nrow = 1) +
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


# Figure 5D: Plasma Pool miRs detected points overlaid------
this.outfile <- paste0(outdirs["FIG5"], "/FIG5D_miRs_detected_plasmaPool_pointsOverlaid.pdf")
g <- ggplot(est.prob.detected.n.miRs,
            aes(
              x=lib.method.simple,
              y=n.detected, fill=lib.method.simple)) + 
  geom_boxplot(alpha=0.5, outlier.color = NA) + geom_jitter(position = position_quasirandom(), alpha=0.7) + 
  facet_wrap(~downsample.to, nrow = 1) +
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


# Figure 5D ALT (Not Used): Plasma Pool miRs detected Separate 4N ------
this.outfile <- paste0(outdirs["FIG5"], "/miRs_detected_plasmaPool_detail_NOT_USED.pdf")
g <- ggplot(est.prob.detected.n.miRs,
            aes(
              x=lib.method.detail,
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
ggsave(this.outfile, g, width = 7.5, height = 2); g

# Scaling factors ----
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

ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG[, c("from.method", "to.method"):=tstrsplit(comparison, split=".V.")]
ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset <- subset(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG, select=c("miR.id.full", "comparison", "from.method", "to.method", "logFC", "CI.L", "CI.R"))
setnames(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, c("CI.L", "CI.R", "comparison"), c("CI.95L", "CI.95R", "correction"))
ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset[, from.method:=sub("^X4N", "4N", from.method)]
ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset[, to.method:=sub("^X4N", "4N", to.method)]
ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset[, correction:=sub("^X4N", "4N", correction)]
this.outfile <- paste0(outdirs["tables"],"/TABLES6_INTER-METHOD_Correction_Factors_equimolar_pool_voom_limma.xlsx")
write.xlsx(ALL.COMPARISON.SCALING.FACTOR.DT.4N.SUBS.AND.ORIG.subset, this.outfile, row.names=FALSE, col.names=TRUE)

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
ann.colors.plasm = list(
  from.method = ann_colors$lib.method.detail,
  to.method = ann_colors$lib.method.detail,
  is.corrected = grey_pal()(2)
)

# SUPP FIG S9: Plot heatmaps Plasma Pool Corrected VS Uncorrected ----
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
  
  this.outfile <- paste0(outdirs["supplemental_figures"], "/FIG_S9_", this.comparison, "_Correction_Heatmap_PLASMAPOOL.pdf")
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
                                                                                
# Apply correction factors to EQUIMOLAR Pool ----
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
# SUPP FIG S8c-e: Plot heatmaps Equimolar Pool Corrected VS Uncorrected ----
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
  
  this.outfile <- paste0(outdirs["supplemental_figures"], "/FIG_S8_", this.comparison, "_Correction_Heatmap_EQUIMOLAR_POOL_WITH_LEGEND.pdf")
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
  this.outfile <- paste0(outdirs["supplemental_figures"], "/FIG_S8_", this.comparison, "_Correction_Heatmap_EQUIMOLAR_POOL_NOLEGEND.pdf")
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

# SUPP FIG S8C Bias_reduction violin plot equimolar----
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
this.outfile <- paste0(outdirs[["supplemental_figures"]], "/FIGS8C_Equimolar_Correction_Truseq_and_NEB_to_4N_biasReduce_violin.pdf")
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


# Supp Fig 9 Plasma Denstity Plot Correction Factors  -------------------------------------------------------------------
g.s9.all <- ggplot(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF, aes(x=abs.log.dif, color=comparison)) + geom_density(size=1, adjust=1, trim=TRUE) + facet_wrap(~group.id) + theme_bw() + 
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
g1 <- g.s9.all %+% MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF[group.id =="TruSeq.V.4N_B"]; g1
g2 <- g.s9.all %+% MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF[group.id =="TruSeq.V.NEBNext"]; g2
g3 <- g.s9.all %+% MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF[group.id =="NEBNext.V.4N_B"]; g3
this.outfile <- paste0(outdirs["supplemental_figures"], "/FigS9_DensityPlot_PlasmaPool_TS.V.4NB.pdf"); this.outfile
ggsave(filename = this.outfile, plot = g1, units = "in", width=4, height=3)
this.outfile <- paste0(outdirs["supplemental_figures"], "/FigS9_DensityPlot_PlasmaPool_TS.V.NEB.pdf"); this.outfile
ggsave(filename = this.outfile, plot = g2, units = "in", width=4, height=3)
this.outfile <- paste0(outdirs["supplemental_figures"], "/FigS9_DensityPlot_PlasmaPool_NEB.V.4NB.pdf"); this.outfile
ggsave(filename = this.outfile, plot = g3, units = "in", width=4, height=3)


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

g.s8.all.eq <- ggplot(MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF, aes(x=abs.log.dif, color=comparison)) + geom_density(adjust=1, trim=TRUE) + 
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
g1 <- g.s8.all.eq %+% MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF[group.id =="TruSeq.V.4N_B"]; g1
g2 <- g.s8.all.eq %+% MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF[group.id =="TruSeq.V.NEBNext"]; g2
g3 <- g.s8.all.eq %+% MANN.WHITNEY.SHIFTS.EQ.ABS.LOG.DIF[group.id =="NEBNext.V.4N_B"]; g3
g.all <- plot_grid(plotlist = list(g1, g2, g3), nrow = 1, labels = paste0("Correction: ", sub(".V.", " to ", show.groups)), label_size = 8); g.all
this.outfile <- paste0(outdirs["supplemental_figures"], "/FigS8_DensityPlot_EquimolarPool_ALL.pdf"); this.outfile
save_plot(filename = this.outfile, ncol = 3, nrow = 1, base_width = 2.5, base_height = 2, plot = g.all)


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

MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.plot <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB[group.id %in% show.groups]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.plot <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF[group.id %in% show.groups]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.plot[, `:=`(query.Lab="Aggregate", fisher.p.val=NA, fisher.df=NA)]

MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.Add.Agg <- rbind(MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.plot, MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.plot, fill=TRUE)
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.Add.Agg.STATS <- MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.Add.Agg[, .(abs.log.dif=max(abs.log.dif)), by=c("group.id", "from.method", "to.method", "query.Lab", "W", "p.value", "alternative", "estimate", "fisher.p.val", "fisher.df", "comparison")]
MANN.WHITNEY.SHIFTS.PLASMA.ABS.LOG.DIF.PLASMA.FISHER.BY.LAB.Add.Agg.STATS[, max.abs.log.dif.group:=max(abs.log.dif), by=group.id]

# SUPP FIG9: FigS9_BOXPLOT
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
  
  this.outfile <- paste0(outdirs["supplemental_figures"], "/FigS9_BOXPLOT_PLASMAPool_withAgg", gid, ".pdf")
  ggsave(filename = this.outfile, plot = g, units = "in", width=4, height=3)
})

