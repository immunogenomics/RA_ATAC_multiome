print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(ArchR)
    library(stringr)
    library(argparse)
})


parser <- ArgumentParser(description='extract TF motif info from ArchR project')
parser$add_argument('project_dir', metavar='project_dir', help='project directory')
parser$add_argument('--wk_dir', metavar='wk_dir', help='working directory. this directory must already exist or will default to one directory up from project_dir')
parser$add_argument('extract_dir', metavar='extract_dir', help='extract (aka output) directory')
parser$add_argument('motif_prefix', metavar='motif_prefix', help='motif set prefix; usually "" for cisBP and "J20" for JASPAR2020')
parser$add_argument('--topn', metavar='topn', default=20, type='integer', help='mlog10Padj of top enriched motifs per cluster; default=20')
parser$add_argument('--motif_fam', metavar='motif_fam', help='comma delimited list of motif families to grep for')
parser$add_argument('--extract_bed', action='store_true', help='if given, output BED file of motif positions.')
args <- parser$parse_args()

project_dir <- args$project_dir
wk_dir <- args$wk_dir
extract_dir <- args$extract_dir
motif_prefix <- args$motif_prefix
topn <- args$topn
motif_fam <- args$motif_fam
extract_bed <- args$extract_bed

print_time("Argument Checking")
if(!all(file.exists(project_dir))) stop("Input file(s) don't exist.")

if(substrRight(project_dir,1)!='/') project_dir<-paste(sep="",project_dir,'/')

if(substrRight(extract_dir,1)!='/') extract_dir<-paste(sep="",extract_dir,'/')
dir.create(extract_dir,showWarnings=FALSE)

cat("Arguments\n")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
}


print_time("Set wd")
if(length(wk_dir)==0){
    wk_dir <- paste(sep='',dirname(project_dir),'/')
} else {
    if(substrRight(wk_dir,1)!='/') wk_dir<-paste(sep="",wk_dir,'/')
    if(!file.exists(wk_dir)) wk_dir <- paste(sep='',dirname(project_dir),'/')
}
cat(paste('New working dir:',wk_dir,'\n'))
setwd(wk_dir)


print_time("motif matches bed")
if(extract_bed){
    motPos_file <- paste(sep='',project_dir,'Annotations/',motif_prefix,'Motif-Positions-In-Peaks.rds')
    if(file.exists(motPos_file)){
        motPos <- readRDS(motPos_file)

        motif_bed_df <- data.frame('chrom'=character(),'start'=integer(),'end'=integer(),
                                   'name'=character(),'score'=numeric(),stringsAsFactors=FALSE)
        for(idx in 1:length(motPos)){
            motif_name = names(motPos)[idx]
            motif_pos_GR = motPos[[idx]]

            motif_pos_df <- data.frame('chrom'=seqnames(motif_pos_GR),
                                       'start'=start(motif_pos_GR)-1,'end'=end(motif_pos_GR),
                                       'name'=c(rep(motif_name, length(motif_pos_GR))),
                                       'score'=motif_pos_GR$score,'strand'=strand(motif_pos_GR),
                                       stringsAsFactors=FALSE)

            motif_bed_df <- rbind(motif_bed_df,motif_pos_df)
        }
        print(dim(motif_bed_df))
        print(head(motif_bed_df))

        withr::with_options(c(scipen=99),
                            write.table(motif_bed_df,
                                        file=paste(sep='',extract_dir,basename(project_dir),'_',motif_prefix,'motifPos.bed'),
                                        sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE))
    } else {
        cat(paste('SKIPPING: motif positions file does not exist:',motPos_file,'\n.'))
    }
} else {
    cat('NOT REQUESTED: add flag extract_bed to script if BED file of motif matches is desired\n')
}


print_time("pxm")
summ_file <- paste(sep='',project_dir,'Annotations/',motif_prefix,'Motif-In-Peaks-Summary.rds')
if(file.exists(summ_file)){
    summ <- readRDS(summ_file)

    if('motifMatches' %in% names(summ)){
        if('matches' %in% names(assays(summ$motifMatches))){
            pxm <- assays(summ$motifMatches)$matches
            
            saveRDS(pxm,paste(sep='',extract_dir,basename(project_dir),'_',motif_prefix,
                              if(motif_prefix=='') {''} else {'_'},'pxm.rds'))
        } else {
            cat('SKIPPING: summ$motifMatches object does not have matches\n.')
        }
    } else {
        cat('SKIPPING: summ object does not have motifMatches\n.')
    }
    
} else {
    cat(paste('SKIPPING: motif summary file does not exist:',summ_file,'\n.'))
}


print_time('Enriched motifs')
eMot_file <- paste(sep='',extract_dir,basename(project_dir),'_enrich',motif_prefix,'Motifs.rds')
if(file.exists(eMot_file)){
    eMot <- readRDS(eMot_file)
    if('mlog10Padj' %in% names(assays(eMot))){
        eMot_df <- assays(eMot)$mlog10Padj
        saveRDS(eMot_df,paste(sep='',extract_dir,basename(project_dir),'_enrich',motif_prefix,'Motifs_df.rds'))
        
        for(idx in 1:ncol(eMot_df)){
            cat(paste('top motifs for',colnames(eMot_df)[idx],'\n'))
            print(head(eMot_df[order(eMot_df[,idx],decreasing=TRUE),],n=topn))
            cat('\n')
        }
        
        if(length(motif_fam)!=0){
            motif_fam_list <- unlist(str_split(motif_fam,','))
            cat(paste('motif enrichment log10Padj for motif families:',paste(motif_fam_list,collapse=', '),'\n'))
            print(eMot_df[grep(paste(motif_fam_list,collapse='|'),rownames(eMot_df),value=TRUE),])
            cat('\n')
        }
        
        heatmapEM_mat <- plotEnrichHeatmap(eMot, n = topn, transpose = TRUE, cutOff = 5, returnMatrix=TRUE)
        saveRDS(heatmapEM_mat,paste(sep='',extract_dir,basename(project_dir),'_enrich',motif_prefix,'Motifs_heatMat.rds'))
        
    } else {
        cat('SKIPPING: motif enrichment object does not have mlog10Padj\n.')
    }
} else {
    cat(paste('SKIPPING: motif enrichment file does not exist:',eMot_file,'\n.'))
}



print_time("chromVAR z")

proj <- loadArchRProject(path=project_dir,showLogo=FALSE)

colorBy <- paste(sep='',motif_prefix,"MotifMatrix")
allColorBy <-  c("colData", "cellColData", ArchR:::.availableArrays(head(getArrowFiles(proj), 2)))
colorBy <- allColorBy[match(tolower(colorBy), tolower(allColorBy))]
colorBy
if(is.na(colorBy)) stop("do not have colorBy matrix!")

units <- tryCatch({ArchR:::.h5read(getArrowFiles(proj)[1], paste0(colorBy, "/Info/Units"))[1]},error=function(e){"values"})
units

markerMotifs = getFeatures(proj, useMatrix = colorBy)
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
length(markerMotifs)
head(markerMotifs)

if(length(markerMotifs)>0){
    ll <- unlist(str_split(Sys.time(),' '))
    colorMat <- ArchR:::.getMatrixValues(ArchRProj = proj, name = markerMotifs, matrixName = colorBy, 
                                         log2Norm = FALSE, threads = 8,
                                         logFile = paste(sep='',dirname(project_dir),
                                                         '/ArchRLogs/ArchR-getMatrixValues-RtermKW-Date-',ll[1],
                                                         '_Time-',gsub(':','-',ll[2]),'.log'))

    saveRDS(colorMat,paste(sep='',extract_dir,basename(project_dir),'_',motif_prefix,'chromVARz.rds'))
} else {
    cat(paste('SKIPPING: no features selected in',colorBy,'\n'))
}


print_time('Done.')

