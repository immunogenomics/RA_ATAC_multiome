print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

roundUpTo <- function(num,to){ceiling(num / to) * to}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(RcppCNPy)
    library(argparse)
})


parser <- ArgumentParser(description='add scrublet results to mRNA meta')
parser$add_argument('meta_file', metavar='meta_file', help='meta mRNA file')
parser$add_argument('scrublet_file', metavar='scrublet_file', help='scrublet npy file')
args <- parser$parse_args()

meta_file <- args$meta_file
scrublet_file <- args$scrublet_file

print_time("Argument Checking")
if(!all(file.exists(meta_file,scrublet_file))) stop("Input files don't exist.")

print_time('load files')
meta <- readRDS(meta_file)
doublet_vec <- npyLoad(scrublet_file)

print_time('check lengths')
if(nrow(meta)!=length(doublet_vec)) stop('lengths not the same between meta rows and doublet vector')

print_time("add doublet column")
db_col <- 'doublet'
while(db_col %in% colnames(meta)){db_col <- paste(sep='',db_col,'1')}
meta[,db_col] <- doublet_vec

print_time('resave')
saveRDS(meta,meta_file)

print_time('Done.')
