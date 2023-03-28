print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(stringr)
    library(argparse)
})


listify_string <- function(x) {substring(x, 1:nchar(x), 1:nchar(x))}

get_grep_region <- function(region){
    library(stringr)
    
    if(!str_detect(region,'^chr[0-9XY]+:[0-9]+-[0-9]+$')) stop("region must be in ^chr[0-9XY]+:[0-9]+-[0-9]+$ format")
    
    split1 <- str_split_fixed(region,':',2)
    split2 <- str_split_fixed(split1[,2],'-',2)
    thisChr <- split1[,1]
    thisStart <- split2[,1]
    thisEnd <- split2[,2]
    
    start_list <- listify_string(thisStart)
    end_list <- listify_string(thisEnd)
    
    start_grep <- ''
    end_grep <- ''
    for(idx in 1:length(start_list)){
        if(start_list[idx]==end_list[idx]){
            start_grep <- paste(sep='',start_grep,start_list[idx])
            end_grep <- paste(sep='',end_grep,end_list[idx])
        } else {
            start_grep <- paste(sep='',start_grep,'[',start_list[idx],'-',end_list[idx],
                                '][0-9]{',length(start_list)-idx,'}')
            end_grep <- paste(sep='',end_grep,'[',start_list[idx],'-',end_list[idx],
                              '][0-9]{',length(start_list)-idx,'}')
            break
        }
    }
    
    return(paste(sep='\\t',thisChr,start_grep,end_grep))

    
}

parser <- ArgumentParser(description='subset enriched TF motifs by region')
parser$add_argument('locus', metavar='locus', help='locus to get motifs')
parser$add_argument('motif_pos_file', metavar='motif_pos_file', help='motif position file')
parser$add_argument('motif_enrich_file', metavar='motif_enrich_file', help='motif enrichment file')
parser$add_argument('--minE_cutoff', metavar='minE_cutoff', default=1.3, type='double', help='min enrichement cutoff; default=1.3')
parser$add_argument('--cellStates', metavar='cellStates', help='cell states to apply cutoffs; if none given, do across cell states')
parser$add_argument('outDir', metavar='outDir', help='output directory')
parser$add_argument('prefix', metavar='prefix', help='output prefix')
args <- parser$parse_args()

locus <- args$locus
motif_pos_file <- args$motif_pos_file
motif_enrich_file <- args$motif_enrich_file
minE_cutoff <- args$minE_cutoff
cellStates <- args$cellStates
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(motif_pos_file,motif_enrich_file))) stop("Input file(s) don't exist.")
if(!str_detect(locus,'^chr[0-9XY]+:[0-9]+-[0-9]+$')) stop("locus must be in ^chr[0-9XY]+:[0-9]+-[0-9]+$ format")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
}


print_time("Get enriched motifs")
motifs_df <- readRDS(motif_enrich_file)
nrow(motifs_df)
head(motifs_df,n=2)

motifs_df$maxE <- apply(motifs_df,1,max)

if(length(cellStates)!=0){
    cellStates <- unlist(strsplit(cellStates,","))
    if(!all(cellStates %in% colnames(motifs_df))) {cat('Not all cell states are in enriched motif file\n')}
    if(!any(cellStates %in% colnames(motifs_df))) {stop('cell states not valid')}
    
    motifs_minE <- c()
    for(cs in cellStates){
        motifs_minE <- c(motifs_minE,rownames(motifs_df[which(motifs_df[,cs]>minE_cutoff),]))
    }
    motifs_minE <- unique(motifs_minE)
} else {
    motifs_minE <- rownames(motifs_df[which(motifs_df$maxE>minE_cutoff),])
}
length(motifs_minE)
head(motifs_minE,n=2)


print_time('Get locus matches')
cmd <- paste(sep='',"grep -P '",get_grep_region(locus),"' ",motif_pos_file)
cat(paste(cmd,'\n'))
motif_pos_subset <- system(cmd,intern=TRUE)
motif_pos_subset <- as.data.frame(str_split_fixed(motif_pos_subset,'\t',6),stringsAsFactors=FALSE)
colnames(motif_pos_subset) <- c('chr','start','stop','motif','score','strand')
nrow(motif_pos_subset)
head(motif_pos_subset,n=2)
write.table(motif_pos_subset,file=paste(sep='',outDir,prefix,'_motifs.bed'),
            sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)


print_time('Subset locus matches')
motif_pos_enrich <- motif_pos_subset[which(motif_pos_subset$motif %in% motifs_minE),]
nrow(motif_pos_enrich)
write.table(motif_pos_enrich,file=paste(sep='',outDir,prefix,'_motifsEnr.bed'),
            sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)


print_time('\nEnrichments for all motifs in file\n')
motifs_df[unique(motif_pos_subset$motif),]


print_time('Done.')

