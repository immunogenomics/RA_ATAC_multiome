print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(Matrix)
    library(lme4)
    library(stats)
    library(MASS)
    library(ROCR)
    library(rcompanion)
    library(ggplot2)
    library(gtools)
    library(stringr)
    library(argparse)
})

#modified from gtools - to allow for the hyphen to be used as a delimiter or a negative number.
mixedorder_new <- function (x, decreasing = FALSE, na.last = TRUE, blank.last = FALSE,
    numeric.type = c("decimal", "roman"), roman.case = c("upper",
        "lower", "both"),keepNegative=FALSE)
{
    numeric.type <- match.arg(numeric.type)
    roman.case <- match.arg(roman.case)
    if (length(x) < 1)
        return(NULL)
    else if (length(x) == 1)
        return(1)
    if (!is.character(x))
        return(order(x, decreasing = decreasing, na.last = na.last))
    delim = "\\$\\@\\$"
    if (numeric.type == "decimal") {
        if(keepNegative)
            regex <- "((?:(?i)(?:[-+]?)(?:(?=[.]?[0123456789])(?:[0123456789]*)(?:(?:[.])(?:[0123456789]{0,}))?)(?:(?:[eE])(?:(?:[-+]?)(?:[0123456789]+))|)))"
        else
            regex <- "((?:(?i)(?:(?=[.]?[0123456789])(?:[0123456789]*)(?:(?:[.])(?:[0123456789]{0,}))?)))"
        numeric <- function(x) as.numeric(x)
    }
    else if (numeric.type == "roman") {
        regex <- switch(roman.case, both = "([IVXCLDMivxcldm]+)",
            upper = "([IVXCLDM]+)", lower = "([ivxcldm]+)")
        numeric <- function(x) roman2int(x)
    }
    else stop("Unknown value for numeric.type: ", numeric.type)
    nonnumeric <- function(x) {
        ifelse(is.na(numeric(x)), toupper(x), NA)
    }
    x <- as.character(x)
    which.nas <- which(is.na(x))
    which.blanks <- which(x == "")
    delimited <- gsub(regex, paste(delim, "\\1", delim, sep = ""),
        x, perl = TRUE)
    step1 <- strsplit(delimited, delim)
    step1 <- lapply(step1, function(x) {x[x > ""]})
    suppressWarnings(step1.numeric <- lapply(step1, numeric))
    suppressWarnings(step1.character <- lapply(step1, nonnumeric))
    maxelem <- max(sapply(step1, length))
    step1.numeric.t <- lapply(1:maxelem, function(i) {sapply(step1.numeric,
        function(x) {x[i]})})
    step1.character.t <- lapply(1:maxelem, function(i) {sapply(step1.character,
        function(x) {x[i]})})
    rank.numeric <- sapply(step1.numeric.t, rank)
    rank.character <- sapply(step1.character.t, function(x) {as.numeric(factor(x))})
    rank.numeric[!is.na(rank.character)] <- 0
    rank.character <- t(t(rank.character) + apply(matrix(rank.numeric),
        2, max, na.rm = TRUE))
    rank.overall <- ifelse(is.na(rank.character), rank.numeric,
        rank.character)
    order.frame <- as.data.frame(rank.overall)
    if (length(which.nas) > 0)
        if (is.na(na.last))
            order.frame[which.nas, ] <- NA
        else if (na.last)
            order.frame[which.nas, ] <- Inf
        else order.frame[which.nas, ] <- -Inf
    if (length(which.blanks) > 0)
        if (is.na(blank.last))
            order.frame[which.blanks, ] <- NA
        else if (blank.last)
            order.frame[which.blanks, ] <- 1e+99
        else order.frame[which.blanks, ] <- -1e+99
    order.frame <- as.list(order.frame)
    order.frame$decreasing <- decreasing
    order.frame$na.last <- NA
    retval <- do.call("order", order.frame)
    return(retval)
}

mixedsort_new <- function (x, decreasing = FALSE, na.last = TRUE, blank.last = FALSE,
    numeric.type = c("decimal", "roman"), roman.case = c("upper",
        "lower", "both"),keepNegative=FALSE)
{
    ord <- mixedorder_new(x, decreasing = decreasing, na.last = na.last,
        blank.last = blank.last, numeric.type = numeric.type,
        roman.case = roman.case, keepNegative = keepNegative)
    x[ord]
}

createDonorMat <- function(this_meta,this_sample_col){
    donors <- unique(this_meta[,this_sample_col])
    for(idx in 1:length(donors)){
        thisDonor <- donors[idx]
        thisDonorCol <- ifelse(this_meta[,this_sample_col]==thisDonor,1,0)

        if(idx==1){
            this_donor_df <- data.frame('donor'=thisDonorCol)
            rownames(this_donor_df) <- rownames(this_meta)
            colnames(this_donor_df) <- thisDonor
        } else {
            this_donor_df[,thisDonor] <- thisDonorCol
        }   
    }

    ll <- colSums(this_donor_df)
    rr <- table(this_meta[,this_sample_col])
    ss <- as.numeric(rr)
    names(ss) <- names(rr)
    theseNames <- sort(names(rr))
    if(!identical(ll[theseNames],ss[theseNames])){
        print(ll[theseNames])
        print(ss[theseNames])
        stop("donor issues")
    }
    
    return(this_donor_df)
}

parse_columns <- function(list_str){
    if(!is.null(list_str)){
        list_toUse <- unlist(strsplit(list_str,","))
    } else{
        list_toUse <- c()
    }
    
    return(list_toUse)
}



parser <- ArgumentParser(description='Pairwise prediction of transcriptional cell states via ATAC hPCs')
parser$add_argument('meta_file', metavar='meta_file', help='metadata file with both mrna and atac info')
parser$add_argument('clust_col', metavar='clust_col', help='cluster/cell type column in meta file that will be predicted - likely mrna fine grain')
parser$add_argument('--clust_val_toRm', metavar='clust_val_toRm', help='cluster/cell type values in meta file to be excluded from further analysis - like scATAC; comma-separated list')
parser$add_argument('--CT', metavar='CT', default='cellType', help='cell type name; only used in the output; default: cellType')
parser$add_argument('sample_col', metavar='sample_col', help='sample/donor column in meta file')
parser$add_argument('--exclude_sample_flag', action='store_true', help='exclude sample/donor column in model')
parser$add_argument('frag_col', metavar='frag_col', help='number of fragments column in meta file')
parser$add_argument('--exclude_frag_flag', action='store_true', help='exclude number of fragments column in model')
parser$add_argument('--scale_within_pair', action='store_true', help='scale fragments within pair')
parser$add_argument('hPC_file', metavar='hPC_file', help='harmonized PC matrix files')
parser$add_argument('--nPCs', metavar='nPCs', default=10, type='integer', help='nPCs to use to look for harmony PC file; default 10')
parser$add_argument('--seed', metavar='seed', default=1234567890, type='double', help='randomization seed; default 1234567890') 
parser$add_argument('--train_perc', metavar='train_perc', default=0.75, type='double', help='percentage of cells to train LDA model; default 0.75') 
parser$add_argument('--cell_count_cutoff', metavar='cell_count_cutoff', default=50, type='integer', help='predicted cluster must have this amount of cells to be used in pairwise comparison; default 50')
parser$add_argument('--plot_AUC', metavar='plot_AUC', choices=c('all','extremes','none'), default='none', help='plot AUC diagrams? options: all, extremes (AUC>0.99 or AUC>0.3 or AUC in [0.52,0.48]), none; default none')
parser$add_argument('--saveAUC_file', metavar='saveAUC_file', help='dataframe file listing c0/c1 pairs to save AUC plot dataframe RDS files')
parser$add_argument('outDir', metavar='outDir', help='Output directory')
parser$add_argument('prefix', metavar='prefix', help='Prefix for output files')
args <- parser$parse_args()

meta_file <- args$meta_file
clust_col <- args$clust_col
clust_val_toRm <- args$clust_val_toRm #optional
CT <- args$CT #optional
sample_col <- args$sample_col
exclude_sample_flag <- args$exclude_sample_flag #optional
frag_col <- args$frag_col
exclude_frag_flag <- args$exclude_frag_flag #optional
scale_within_pair <- args$scale_within_pair
hPC_file <- args$hPC_file
nPCs <- args$nPCs
seed <- args$seed
rm_other <- args$rm_other
train_perc <- args$train_perc
cell_count_cutoff <- args$cell_count_cutoff
plot_AUC <- args$plot_AUC
saveAUC_file <- args$saveAUC_file #optional
outDir <- args$outDir
prefix <- args$prefix

print_time("Argument Checking")
if(!all(file.exists(meta_file,hPC_file))) stop("Input file(s) don't exist.")

if(substrRight(outDir,1)!='/') outDir<-paste(sep="",outDir,'/')
dir.create(outDir,showWarnings=FALSE)
outPrefix <- paste(sep="",outDir,prefix)

if(plot_AUC %in% c('all','extremes') | length(saveAUC_file)!=0){
    plotDir <- paste(sep="",outDir,'AUC_plots/')
    dir.create(plotDir,showWarnings=FALSE)
    plotPrefix <- paste(sep="",plotDir,prefix)
}

saveAUC_flag <- FALSE
if(length(saveAUC_file)!=0){
    if(file.exists(saveAUC_file)){
        saveAUC_df <- readRDS(saveAUC_file)
        saveAUC_flag <- TRUE
    }
}

cat("Arguments\n")
outFile <- paste(sep="",outPrefix,"_args.txt")
for(i in 1:length(args)){
    line <- paste(sep="",names(args)[i]," = ",args[[i]],"\n")
    cat(line)
    if(i==1) cat(line,file=outFile) else cat(line,file=outFile,append=TRUE)
}


print_time("Load meta file")
meta <- readRDS(meta_file)
if(!all(c(clust_col,sample_col,frag_col) %in% colnames(meta))) stop('meta does not have all required columns.')


print_time("Load hPC file")
harmony_mat <- readRDS(hPC_file)
if(!all(rownames(harmony_mat) %in% rownames(meta))) stop('not all hPC mat cells in meta')
if(ncol(harmony_mat)<nPCs) stop('asked for more PCs than given')
if(ncol(harmony_mat)>nPCs) harmony_mat[,1:nPCs]
meta_subset <- meta[rownames(harmony_mat),]
if(!identical(rownames(harmony_mat),rownames(meta_subset))) stop('cell name mismatch')


if(length(clust_val_toRm)!=0){
    print_time(paste('Removing',clust_val_toRm,'values from',clust_col,'column'))
    clust_val_toRm_list <- parse_columns(clust_val_toRm)
    if(any(clust_val_toRm_list %in% unique(meta_subset[,clust_col]))){
        kept_cells <- rownames(meta_subset[which(!(meta_subset[,clust_col] %in% clust_val_toRm_list)),])
        cat(paste('Keeping',length(kept_cells),'cells\n'))
        meta_subset <- meta_subset[kept_cells,]
        harmony_mat <- harmony_mat[kept_cells,]
        if(!identical(rownames(harmony_mat),rownames(meta_subset))) stop('cell name mismatch')
    } else {
        cat(paste('Skipping:',clust_val_toRm,'values NOT FOUND in',clust_col,'column\n'))
    }
}


print_time("Donor matrix")
donor_mat <- as.matrix(createDonorMat(meta_subset,sample_col)) #only using if rownames match
colnames(donor_mat) <- str_replace(colnames(donor_mat),'-','')
colSums(donor_mat)
table(meta_subset[,sample_col])
#need to exclude one anyway, might as well be the smallest donor
fewest_donor <- names(sort(colSums(donor_mat))[1])
donor_mat <- donor_mat[,which(colnames(donor_mat)!=fewest_donor)] 
head(donor_mat)


if(!scale_within_pair){
    print_time("Scaling fragments across ALL cells")
    scale_frag_col <- 'scale_log10_frag'
    while(scale_frag_col %in% colnames(meta_subset)) scale_frag_col <- paste(sep='',scale_frag_col,1)
    meta_subset[,scale_frag_col] <- scale(log10(meta_subset[,frag_col]))
    head(meta_subset[,c(frag_col,scale_frag_col)])
}


print_time("big for loops")
set.seed(seed)

stats_df <- data.frame('cellType'=character(),'c0'=character(),'c1'=character(),'feature'=character(),
                       'value'=numeric(),stringsAsFactors=FALSE)


meta_subset$newClustNames <- meta_subset[,clust_col]

clusters_toPredict_toIter <- mixedsort_new(unique(meta_subset$newClustNames)) #sorting for upper triangle

first_flag <- TRUE
for(c0 in clusters_toPredict_toIter){
    for(c1 in clusters_toPredict_toIter){        
        
        #only need one triangle since pairwise!
        if(which(clusters_toPredict_toIter==c0)>=which(clusters_toPredict_toIter==c1)){ 
            next
        }
        
        
        #cells in this pair of clusters
        c0_cells <- rownames(meta_subset[which(meta_subset$newClustNames==c0),])
        c1_cells <- rownames(meta_subset[which(meta_subset$newClustNames==c1),])
        pair_cells <- c(c0_cells,c1_cells)
        if((length(c0_cells)<cell_count_cutoff) | (length(c1_cells)<cell_count_cutoff)){
            cat(paste('Skipping due to lack of cell counts:',c0,length(c0_cells),c1,length(c1_cells),'\n'))
            next
        }
        
        
        #subset data structures
        meta_pair <- meta_subset[pair_cells,]
        
        pred_clust <- ifelse(meta_pair$newClustNames==c0,0,1)
        
        harmony_mat_pair <- harmony_mat[pair_cells,]
        
        donor_col_pair <- meta_pair[,sample_col]
        
        donor_mat_pair <- donor_mat[pair_cells,]
        
        if(scale_within_pair){
            fragments_scale_log10_col_pair <- scale(log10(meta_pair[,frag_col]))
        } else {
            fragments_scale_log10_col_pair <- meta_pair[,scale_frag_col]
        }
        
        
        #combine data structures; explicitly naming single column covariates!
        data_wDonorMat_pair <- as.data.frame(harmony_mat_pair)
        if(!exclude_sample_flag) data_wDonorMat_pair <- cbind(data_wDonorMat_pair,donor_mat_pair) #excluded 1 column above!
        if(!exclude_frag_flag) data_wDonorMat_pair <- cbind(data_wDonorMat_pair,'frags'=fragments_scale_log10_col_pair)
        data_wDonorMat_pair <- cbind(data_wDonorMat_pair,pred_clust)
                
        
        #formula
        formula_full_wDonorMat <- as.formula(paste(colnames(data_wDonorMat_pair)[ncol(data_wDonorMat_pair)], "~", 
                                                   paste(colnames(data_wDonorMat_pair)[-ncol(data_wDonorMat_pair)],
                                                         collapse=" + ")))
                
        
        if(first_flag){
            cat('Formula\n')
            
            print(formula_full_wDonorMat)
            
            first_flag <- FALSE
        }
        
        
        #split train/test sets
        smp_size <- floor(train_perc * length(pair_cells))
        train_ind <- sample(length(pair_cells), size = smp_size)
        train_df <- as.data.frame(data_wDonorMat_pair[train_ind, ])
        test_df <- as.data.frame(data_wDonorMat_pair[-train_ind, ])
        
        
        #validate that both train and test have both values
        while(length(unique(train_df[,ncol(train_df)]))!=2 | length(unique(test_df[,ncol(test_df)]))!=2){
            cat(paste('re-picking train/test sets since both or either set missing both or either cluster', 
                      c0,c1,length(unique(train_df[,ncol(train_df)])),length(unique(test_df[,ncol(test_df)])),'\n'))

            train_ind <- sample(length(pair_cells), size = smp_size)
            train_df <- as.data.frame(data_wDonorMat_pair[train_ind, ])
            test_df <- as.data.frame(data_wDonorMat_pair[-train_ind, ])
        }
        
        
        #LDA
        resLDA <- tryCatch(withCallingHandlers({
                this.lda <- lda(formula_full_wDonorMat, data = train_df)
            }, error = function(err) {
                cat(paste(sep="\n",paste("ERROR: LDA",c0,c1),err))
            }, warning = function(warn) {
                cat(paste(trimws(warn),"LDA",c0,c1,'\n'))
            }),
        error = function(err){
            doNothing = 1
        })

        if(exists("this.lda")){
            this.lda.predict <- predict(this.lda, newdata = test_df)

            ct <- table(test_df$pred_clust, this.lda.predict$class)
            percCorr_byClust <- diag(prop.table(ct, 1)) # percent correct for each category 
            percCorr_tot <- sum(diag(prop.table(ct))) # total percent correct
            percCorr_c0 <- as.numeric(percCorr_byClust[which(names(percCorr_byClust)==0)])
            percCorr_c1 <- as.numeric(percCorr_byClust[which(names(percCorr_byClust)==1)])


            #AUC plot
            #Get the posteriors as a dataframe.
            this.lda.predict.posteriors <- as.data.frame(this.lda.predict$posterior)
            #Evaluate the model
            pred <- prediction(this.lda.predict.posteriors[,2], test_df$pred_clust)
            roc.perf = performance(pred, measure = "tpr", x.measure = "fpr")
            auc.train <- performance(pred, measure = "auc")
            auc.train <- auc.train@y.values
            auc_val <- auc.train[[1]]
            
            if(saveAUC_flag){
                if(nrow(saveAUC_df[which(saveAUC_df$c0==c0 & saveAUC_df$c1==c1),])>0){
                    c0_subname <- str_split_fixed(c0,':',2)[,1]
                    c1_subname <- str_split_fixed(c1,':',2)[,1]
                    saveRDS(pred,paste(sep='',plotPrefix,'_pred_',c0_subname,'_',c1_subname,'.rds'))
                }
            }

            if((plot_AUC=='all') | (plot_AUC=='extremes' & (auc_val>0.99 | auc_val<0.3 | (auc_val>0.47 & auc_val<0.53)))){
                c0_subname <- str_split_fixed(c0,':',2)[,1]
                c1_subname <- str_split_fixed(c1,':',2)[,1]
                png(filename=paste(sep='',plotPrefix,'_AUC_',c0_subname,'_',c1_subname,'.png'),units='in',height=5,width=5,
                    res=300)
                plot(roc.perf)
                abline(a=0, b= 1)
                title(paste(c0,c1))
                text(x = .85, y = .15 ,paste(sep='','c0: ',nrow(train_df[which(train_df[,ncol(train_df)]==0),]),'/',
                                             nrow(test_df[which(test_df[,ncol(test_df)]==0),]),
                                             '\nc1: ',nrow(train_df[which(train_df[,ncol(train_df)]==1),]),'/',
                                             nrow(test_df[which(test_df[,ncol(test_df)]==1),]),
                                             '\nAUC = ', round(auc_val,2)))
                dev.off()
            }
            
            rm(this.lda) #to reset block
        } else {
            percCorr_c0 <- NA
            percCorr_c1 <- NA
            percCorr_tot <- NA
            auc_val <- NA
        }
        
        
        numFeatures <- 10
        this_df <- data.frame('cellType'=rep(CT,numFeatures),
                              'c0'=rep(c0,numFeatures),
                              'c1'=rep(c1,numFeatures),
                              'feature'=c('nC0','nC1','nTrain','nTest','nTrain_c0','nTest_c0',
                                          'percCorr_c0','percCorr_c1','percCorr_tot','AUC'),
                              'value'=c(length(c0_cells),length(c1_cells),nrow(train_df),nrow(test_df),
                                        nrow(train_df[which(train_df[,ncol(train_df)]==0),]),
                                        nrow(test_df[which(test_df[,ncol(test_df)]==0),]),
                                        percCorr_c0,percCorr_c1,percCorr_tot,auc_val),
                              stringsAsFactors=FALSE)
        
        stats_df <- rbind(stats_df,this_df)

    }
}


print_time("Saving")
saveRDS(stats_df,paste(sep='',outPrefix,'_statsDF.rds')) 


print_time("Done.")

