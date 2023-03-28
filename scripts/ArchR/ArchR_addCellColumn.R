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


parser <- ArgumentParser(description='Add cell cluster metadata column to ArchR project')
parser$add_argument('project_dir', metavar='project_dir', help='project (aka output) directory - full path')
parser$add_argument('meta_file', metavar='meta_file', help='meta file')
parser$add_argument('meta_column', metavar='meta_column', help='cluster column in meta file')
parser$add_argument('extract_dir', metavar='extract_dir', help='extract (aka output) directory')
parser$add_argument('--thread_num', metavar='thread_num', default=4, type='integer', help='number of threads to use; default 4')
parser$add_argument('--wk_dir', metavar='wk_dir', help='working directory. this directory must already exist or will default to one directory up from project_dir')
parser$add_argument('--seed', metavar='seed', default=1, type='integer', help='randomization seed; default 1')
args <- parser$parse_args()

project_dir <- args$project_dir
meta_file <- args$meta_file
meta_column <- args$meta_column
extract_dir <- args$extract_dir
thread_num <- args$thread_num
wk_dir <- args$wk_dir
seed <- args$seed

print_time("Argument Checking")
if(!all(file.exists(project_dir,meta_file))) stop("Input file(s) don't exist.")

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

print_time("Load Meta")
meta <- readRDS(meta_file)
if(!(meta_column %in% colnames(meta))) stop('meta column not in meta file')
cellCol <- 'ArchR_cellNames'
while(cellCol %in% colnames(meta)) cellCol <- paste(sep='',cellCol,'1')
if(!all(str_detect(rownames(meta),'^BRI-[0-9]+#[ACGT]+-1$'))){
    if(all(str_detect(rownames(meta),'^BRI-[0-9]+_[ACGT]+$'))){
        splits <- str_split_fixed(rownames(meta),'_',2)
        meta[,cellCol] <- paste(sep='',splits[,1],'#',splits[,2],'-1')
    } else {
        stop('Cell name pattern unknown.')
    }
} else {
    meta[,cellCol] <- rownames(meta)
}
dim(meta)
head(meta,n=3)
if(!all(str_detect(meta[,cellCol],'^BRI-[0-9]+#[ACGT]+-1$'))) stop('ArchR cell format not present.')

print_time("Load Project")
proj <- loadArchRProject(path=project_dir, showLogo=FALSE)

print_time("Verify cells")
archr_cells <- getCellNames(ArchRProj = proj)
table(archr_cells %in% meta[,cellCol])
if(!identical(sort(meta[,cellCol]), sort(archr_cells))){
    cat('WARNING: Cells not the same between meta and ArchR project.\n')
}
if(all(archr_cells %in% meta[,cellCol])){
    cat(paste(sep='','WARNING: Subsetting meta (n=',nrow(meta),') by ArchR cells (n=',length(archr_cells),')\n'))
    meta <- meta[which(meta[,cellCol] %in% archr_cells),]
} else {
    stop('ERROR: Not all ArchR cells in meta file\n')
}

print_time("Add meta col")
ArchR_column <- paste(sep='','my_',meta_column)
proj <- addCellColData(ArchRProj = proj, data = meta[,meta_column], name=ArchR_column, cells=meta[,cellCol])

print_time('Extract cellData')
cellColData <- getCellColData(ArchRProj = proj)
saveRDS(cellColData,paste(sep='',extract_dir,basename(project_dir),'_cellColData_wClusterCol.rds'))

print_time('Check extraction')
them <- as.data.frame(cellColData)
me <- meta[,c(meta_column,cellCol)]
me$mrnaCN <- rownames(me)
rownames(me) <- me[,cellCol]
table(me[,meta_column])
table(them[,ArchR_column])
if(!identical(sort(rownames(me)),sort(rownames(them)))) stop('cellNames are not the same')
if(!identical(me[rownames(me),meta_column],them[rownames(me),ArchR_column])) stop('cluster column not the same')

print_time("Saving ArchR Project")
proj <- saveArchRProject(ArchRProj = proj)

print_time("Done.")

