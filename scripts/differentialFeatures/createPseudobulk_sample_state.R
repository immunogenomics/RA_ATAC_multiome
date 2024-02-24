print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(stringr)
    library(argparse)
})

parse_columns <- function(list_str,delimiter=','){
    if(!is.null(list_str)){
        list_toUse <- unlist(strsplit(list_str,delimiter))
    } else{
        list_toUse <- c()
    }
    
    return(list_toUse)
}

aggMat_by2col <- function(thisPxC,thisMeta,col1,col2,fun){

    uniq1 <- sort(unique(thisMeta[,col1]))
    uniq2 <- sort(unique(thisMeta[,col2]))
    
    firstIter <- TRUE
    for(idx1 in 1:length(uniq1)){
        for(idx2 in 1:length(uniq2)){
            
            val1 <- uniq1[idx1]
            val2 <- uniq2[idx2]
            
            these_cells <- rownames(thisMeta[which(thisMeta[,col1]==val1 & thisMeta[,col2]==val2),])
            this_subset <- thisPxC[,these_cells]

            if(length(these_cells)==0){
                cat(paste('NO CELLS:',col1,'==',val1,'and',col2,'==',val2,'\n'))
                thisAgg <- NULL
                next
            }
            if(length(these_cells)==1){
                thisAgg <- this_subset
            } else {
                thisAgg <- as.matrix(apply(this_subset,1,fun))
            }

            if(firstIter){
                catAgg <- thisAgg
                firstIter <- FALSE
            } else {
                catAgg <- cbind(catAgg,thisAgg)
            }
            colnames(catAgg)[ncol(catAgg)] <- paste(sep='_',val1,val2)
        }
    }

    return(as(catAgg,'dgCMatrix'))

}

checkMat <- function(thisMat){
    rs <- rowSums(thisMat)
    cat(paste(sep='','Number of all zero peaks: ',length(rs[rs==0]),'\n'))

    cs <- colSums(thisMat)
    cat(paste(sep='','Number of all zero cell (types): ',length(cs[cs==0]),'\n'))
}


parser <- ArgumentParser(description='create sample/state pseudobulks')
parser$add_argument('data_file', metavar='data_file', help='Matrix file: non-binary peaks x cells')
parser$add_argument('meta_file', metavar='meta_file', help='Cell metadata file')
parser$add_argument('state_col', metavar = 'state_col', help='cell state column in meta file')
parser$add_argument('--exclude_state_val', metavar = 'exclude_state_val', help='if given, exclude these values from the state column; use a comma-delimited list')
parser$add_argument('sample_col', metavar='sample_col', help='sample/donor column in meta file')
parser$add_argument('--state_cutoff', metavar='state_cutoff', type='integer', default=130, help='state cell count cutoff; default: 130')
parser$add_argument('--sample_cutoff', metavar='sample_cutoff', type='integer', default=150, help='sample/donor cell count cutoff; default: 150')
parser$add_argument('--combo_cutoff', metavar='combo_cutoff', type='integer', default=10, help='sample/state combination cell count cutoff; default: 10')
parser$add_argument('--peak_cutoff', metavar='peak_cutoff', type='integer', default=4, help='per peak pseudobulk count cutoff; default: 4')
parser$add_argument('pma_file', metavar='pma_file', help='Peaks with Minimal Accessibility file')
parser$add_argument('class_state_file', metavar='class_state_file', help='File with at least class and state columns')
parser$add_argument('--rmHyp', action='store_true', help='if given, remove hyphens from sample and state values')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

data_file <- args$data_file
meta_file <- args$meta_file
state_col <- args$state_col
exclude_state_val <- args$exclude_state_val
sample_col <- args$sample_col
state_cutoff <- args$state_cutoff
sample_cutoff <- args$sample_cutoff
combo_cutoff <- args$combo_cutoff
peak_cutoff <- args$peak_cutoff
pma_file <- args$pma_file
class_state_file <- args$class_state_file
rmHyp <- args$rmHyp
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(data_file,meta_file,pma_file))) stop("Input files don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_args.txt")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
}


print_time("Load files")
data <- readRDS(data_file)
meta <- readRDS(meta_file)
pma <- readRDS(pma_file)
class_state_df <- readRDS(class_state_file)
if(!all(c(state_col,sample_col) %in% colnames(meta))) stop("Required column(s) not in meta file")
if(!identical(sort(rownames(meta)),sort(colnames(data)))) stop("Cell mismatch")
if(!identical(colnames(data),rownames(meta))) meta <- meta[colnames(data),]
if(!all(pma %in% rownames(data))) stop("PMA peak mismatch")
if(!all(c('state','class') %in% colnames(class_state_df))) stop("Class/State file does not have columns named class and state")
if(!any(class_state_df$state %in% unique(meta[,state_col]))) stop("There is no overlap between Class/State file states and meta file states")


if(!is.null(exclude_state_val)){
    print_time("Exclude cells with given state values")
    exclude_state_val <- parse_columns(exclude_state_val)
    if(any(exclude_state_val %in% unique(meta[,state_col]))){
        meta <- meta[which(!(meta[,state_col] %in% exclude_state_val)),]
        data <- data[,rownames(meta)]
        cat(paste("Subset to",nrow(meta),"cells.\n"))
    } else {
        cat("None of the state values given exist in the given state column to exclude them. Not subsetting.\n")
    }
}


print_time("Pseudobulk matrix")
ps_data <- aggMat_by2col(data,meta,state_col,sample_col,sum)


print_time("QC and PMA subsetting")
state_cell_count <- as.data.frame(table(meta[,state_col]),stringsAsFactors=FALSE)
colnames(state_cell_count) <- c(state_col,'cellCt')
exclude_states <- state_cell_count[which(state_cell_count[,'cellCt']<state_cutoff),state_col]
cat(paste('Excluding states:',paste(collapse=', ',exclude_states),'\n'))

sample_cell_count <- as.data.frame(table(meta[,sample_col]),stringsAsFactors=FALSE)
colnames(sample_cell_count) <- c(sample_col,'cellCt')
exclude_samples <- sample_cell_count[which(sample_cell_count[,'cellCt']<sample_cutoff),sample_col]
cat(paste('Excluding samples:',paste(collapse=', ',exclude_samples),'\n'))

sample_state_cell_count <- as.data.frame(table(meta[,c(state_col,sample_col)]),stringsAsFactors=FALSE)
sample_state_cell_count$clnm <- paste(sep='_',sample_state_cell_count[,state_col],sample_state_cell_count[,sample_col])

exclude_combos <- sample_state_cell_count[which(sample_state_cell_count$Freq<combo_cutoff | 
                                                   sample_state_cell_count[,state_col] %in% exclude_states |
                                                   sample_state_cell_count[,sample_col] %in% exclude_samples),'clnm']

keep_combos <- sample_state_cell_count[which((!(sample_state_cell_count$clnm %in% exclude_combos)) & 
                                             (sample_state_cell_count$clnm %in% colnames(ps_data))),'clnm']

if(!all(keep_combos %in% colnames(ps_data))) stop('pseudobulks mismatch')
if(!all(pma %in% rownames(ps_data))) stop('PMA peaks mismatch')

ps_data <- ps_data[pma,keep_combos]


print_time("Collapse metadata")
splits <- str_split_fixed(colnames(ps_data),'_',2)

cc <- unlist(lapply(sort(unique(meta[,state_col])),
             function(x){sum(colSums(ps_data[,grep(x,colnames(ps_data))]))}))
names(cc) <- sort(unique(meta[,state_col]))

ss <- unlist(lapply(sort(unique(meta[,sample_col])),
             function(x){sum(colSums(ps_data[,grep(x,colnames(ps_data))]))}))
names(ss) <- sort(unique(meta[,sample_col]))

ps_meta <- data.frame(row.names=colnames(ps_data),
                      'sample'=splits[,2],'CITE_abbr'=splits[,1],'sample_nFrag'=ss[splits[,2]],'CITE_nFrag'=cc[splits[,1]],
                      'sample_CITE_nFrag'=colSums(ps_data),stringsAsFactors=FALSE)
colnames(ps_meta) <- c(sample_col,state_col,paste(sep='',sample_col,'_nFrag'),paste(sep='',state_col,'_nFrag'),
                       paste(sep='_',sample_col,state_col,'nFrag'))


if(rmHyp){
    print_time("Remove hyphens")
    rownames(ps_meta) <- gsub('-','',rownames(ps_meta))
    colnames(ps_data) <- gsub('-','',colnames(ps_data))

    ps_meta[,sample_col] <- gsub('-','',ps_meta[,sample_col])
    ps_meta[,state_col] <- gsub('-','',ps_meta[,state_col])
    
    class_state_df$class <- gsub('-','',class_state_df$class)
    class_state_df$state <- gsub('-','',class_state_df$state)
}


print_time("Verify order")
if(!identical(rownames(ps_meta),colnames(ps_data))) stop("Pseudobulks mismatch")


print_time("Save full files")
saveRDS(ps_data,paste(sep='',outPrefix,'_allClass_matrix.rds'))
saveRDS(ps_meta,paste(sep='',outPrefix,'_allClass_meta.rds'))


print_time("Split by class")
for(class in unique(class_state_df$class)){
    states_inClass <- unique(class_state_df[which(class_state_df$class==class),'state'])
    if(length(unique(ps_meta[which(ps_meta[,state_col] %in% states_inClass),state_col]))<2){
        cat(paste('SKIPPING: class',class,'as there are less than 2 post-QC states.\n'))
        next
    }
    
    class_meta <- ps_meta[which(ps_meta[,state_col] %in% class_state_df[which(class_state_df$class==class),'state']),]
    
    class_data <- ps_data[,rownames(class_meta)]
    
    if(!identical(rownames(class_meta),colnames(class_data))) stop('Pseudobulks mismatch')
        
    rs <- rowSums(class_data)
    class_data <- class_data[names(rs[rs>peak_cutoff]),]
    
    cs <- colSums(class_data)
    if(length(cs[cs==0])>0) stop('some pseudobulks not present!')
    
    if(!identical(rownames(class_meta),colnames(class_data))) stop("Pseudobulks mismatch")
    
    saveRDS(class_data,paste(sep='',outPrefix,'_',class,'_matrix.rds'))
    saveRDS(class_meta,paste(sep='',outPrefix,'_',class,'_meta.rds'))
    
}


print_time('Done.')

