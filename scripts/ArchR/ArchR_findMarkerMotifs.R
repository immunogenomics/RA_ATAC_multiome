print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(ArchR)
    library(stringr)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(JASPAR2016)
    library(JASPAR2018)
    library(JASPAR2020)
    library(chromVARmotifs)
    library(argparse)
})


parser <- ArgumentParser(description='Add motifs to ArchR project')
parser$add_argument('project_dir', metavar='project_dir', help='project (aka output) directory - full path')
parser$add_argument('groupBy_col', metavar='groupBy_col', help='groupBy (aka cluster) column')
parser$add_argument('motif_set', metavar='motif_set', choices=c("JASPAR2016", "JASPAR2018", "JASPAR2020", "cisbp", "encode", "homer"), help='motif set')
parser$add_argument('motif_prefix', metavar='motif_prefix', help='motif set prefix')
parser$add_argument('--motif_width', metavar='motif_width', default=7, type='integer', help='motif search window width; default 7')
parser$add_argument('--motif_pval', metavar='motif_pval', default=5e-05, type='double', help='motif p-value; default 5e-05')
parser$add_argument('--motif_ver', metavar='motif_ver', default=2, type='integer', choices=c(1,2), help='motif version; default 2')
parser$add_argument('--motif_force', action='store_true', help='if given, force overwritting of all motif stuff in ArchR project')
parser$add_argument('--peak_force', action='store_true', help='if given, force overwritting of all peak stuff in ArchR project')
parser$add_argument('extract_dir', metavar='extract_dir', help='extract (aka output) directory')
parser$add_argument('--thread_num', metavar='thread_num', default=8, type='integer', help='number of threads to use; default 8')
parser$add_argument('--wk_dir', metavar='wk_dir', help='working directory. this directory must already exist or will default to one directory up from project_dir')
parser$add_argument('--seed', metavar='seed', default=1, type='integer', help='randomization seed; default 1')
args <- parser$parse_args()

project_dir <- args$project_dir
groupBy_col <- args$groupBy_col
motif_set <- args$motif_set
motif_prefix <- args$motif_prefix
motif_width <- args$motif_width
motif_pval <- args$motif_pval
motif_ver <- args$motif_ver
motif_force <- args$motif_force
peak_force <- args$peak_force
extract_dir <- args$extract_dir
thread_num <- args$thread_num
wk_dir <- args$wk_dir
seed <- args$seed

print_time("Argument Checking")
if(!all(file.exists(project_dir))) stop("Input file doesn't exist.")

if(substrRight(project_dir,1)!='/') project_dir<-paste(sep="",project_dir,'/')
if(substrRight(extract_dir,1)!='/') project_dir<-paste(sep="",extract_dir,'/')
if(!file.exists(extract_dir)) dir.create(extract_dir)

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


print_time('Set seed and threads')
set.seed(seed)
addArchRThreads(threads = thread_num)


print_time("Load Project")
proj <- loadArchRProject(path=project_dir, showLogo=FALSE)


print_time("Verify groupBy_col")
if(!(groupBy_col %in% colnames(getCellColData(ArchRProj = proj)))) stop('groupBy column does not exist in cellColData')


print_time('Add Motifs')
proj <- addMotifAnnotations(ArchRProj = proj,
                            motifSet = motif_set,
                            annoName=paste(sep='',motif_prefix,"Motif"),
                            width=motif_width,cutOff=motif_pval,version=motif_ver,force=motif_force) 

if(file.exists(paste(sep='',project_dir,'Background-Peaks.rds'))){
    if(peak_force){
        proj <- addBgdPeaks(ArchRProj = proj,force=TRUE)
    } else {
        cat('NOT redoing background peaks!\n')
    }
} else {
    proj <- addBgdPeaks(ArchRProj = proj)
}

proj <- addDeviationsMatrix(ArchRProj = proj, peakAnnotation=paste(sep='',motif_prefix,"Motif"), force=motif_force)

getAvailableMatrices(proj)


print_time('Extract Motif Matrix')
saveRDS(getVarDeviations(ArchRProj = proj, name = paste(sep='',motif_prefix,"MotifMatrix"), plot = FALSE),
        paste(sep='',extract_dir,basename(project_dir),'_',motif_prefix,'MotifMatrix.rds'))


print_time('Save')
proj <- saveArchRProject(ArchRProj = proj)


print_time('test matches')
mtchs <- getMatches(proj,paste(sep='',motif_prefix,"Motif"))


print_time('Find Marker Peaks')
markerPeaks_file <- paste(sep='',extract_dir,basename(project_dir),'_markerPeaks.rds')
if(!file.exists(markerPeaks_file) | peak_force){
    markerPeaks <- getMarkerFeatures(ArchRProj = proj,useMatrix = "PeakMatrix",groupBy = groupBy_col,
                                      bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
    saveRDS(markerPeaks,markerPeaks_file)
} else if(file.exists(markerPeaks_file)){
    markerPeaks <- readRDS(markerPeaks_file)
} else {
    cat('should not have gotten here...\n')
}


print_time('test matches')
mtchs <- getMatches(proj,paste(sep='',motif_prefix,"Motif"))


print_time('Find Enriched Motifs')
enrichMotifs <- peakAnnoEnrichment(seMarker = markerPeaks,ArchRProj = proj,peakAnnotation = paste(sep='',motif_prefix,"Motif"),
                                   cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
saveRDS(enrichMotifs,paste(sep='',extract_dir,basename(project_dir),'_enrich',motif_prefix,'Motifs.rds'))


print_time('Plot Enriched Motifs') 
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 5, transpose = TRUE, cutOff = 5)
plotPDF(heatmapEM, name = paste(sep='',basename(project_dir),'_enrich',motif_prefix,'MotifsHeatmap'),
        width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

print_time("Saving ArchR Project")
proj <- saveArchRProject(ArchRProj = proj)


print_time("Done.")

