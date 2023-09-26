print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}


print_time("Library Loading and Defining Functions")

suppressMessages({
    library(GenomicRanges)
    library(Matrix)
    library(Matrix.utils)
    library(gtools)
    library(argparse)
})


iterFiles <- function(toDo_files, toDo_outFile){
    file_list <- read.table(toDo_files,stringsAsFactors=F)
    file_list <- file_list$V1
    
    print_time("Starting concatenation")
    for(i in 1:length(file_list)){
        if(!file.exists(file_list[i])) stop(paste("ERROR: File doesn't exist:",file_list[i]))
        
        if(i==1){
            allCat_df <- readRDS(file_list[i])
        } else {
            placeholder <- readRDS(file_list[i])
            if(identical(rownames(allCat_df),rownames(placeholder))){
                allCat_df <- cbind(allCat_df,placeholder)
            } else {
                cat(paste("File",i,"rownames are different\n"))
            }
            rm(placeholder)
        }
        cat(paste(sep="","After adding matrix ",i,", new # col= ",ncol(allCat_df),"\n"))
    }
    print_time("Stopping concatenation")
    
    cat(paste(sep="",dim(allCat_df),"\n"))
    
    saveRDS(allCat_df,toDo_outFile)
    
    return(allCat_df)
}


parser <- ArgumentParser(description='Concatenate matrices')
parser$add_argument('byCells_file', metavar='byCells_file', help='file with a list of features by cells files to concatenate')
parser$add_argument('byDonors_file', metavar='byDonors_file', help='file with a list of features by donors files to concatenate')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

byCells_file <- args$byCells_file
byDonors_file <- args$byDonors_file
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(c(byCells_file,byDonors_file)))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

cat("Arguments\n")
arg_list <- list("list of features by cells matrices"=byCells_file,
                 "list of features by donors matrices"=byDonors_file,
                 "output directory"=outDir,
                 "output prefix"=prefix)

outFile <- paste(sep="",outPrefix,"_args.txt")
for(i in 1:length(arg_list)){
    line <- paste(sep="",names(arg_list)[i]," = ",arg_list[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
}


byCells_df <- iterFiles(byCells_file, paste(sep="",outPrefix,"Xcells.rds"))

byDonors_df <- iterFiles(byDonors_file, paste(sep="",outPrefix,"Xdonors.rds"))
rs <- rowSums(byDonors_df)
table(rs)

print_time('Done.')

