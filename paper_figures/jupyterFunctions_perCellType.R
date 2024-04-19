suppressMessages({
    library(Matrix)
    library(gtools)
    library(Rmisc)
    library(ROCR)
    library(presto)
    library(ggpubr)
    
    library(plyr)
    library(stringr)
    library(tidyr)
    
    library(ggplot2)
    library(ggrastr)
    library(ggrepel)
    library(viridis)
    library(scales)
    library(RColorBrewer)
    library(grid)
    library(gridExtra)
    library(repr)
})


##UTILITY FUNCTIONS

#SORTING - to allow for the hyphen to be used as a delimiter or a negative number.
#modified from gtools
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

#SORTING
#modified from gtools
mixedsort_new <- function (x, decreasing = FALSE, na.last = TRUE, blank.last = FALSE,
    numeric.type = c("decimal", "roman"), roman.case = c("upper",
        "lower", "both"),keepNegative=FALSE)
{
    ord <- mixedorder_new(x, decreasing = decreasing, na.last = na.last,
        blank.last = blank.last, numeric.type = numeric.type,
        roman.case = roman.case, keepNegative = keepNegative)
    x[ord]
}

gene_order <- function(mat){
    ll <- as.data.frame(apply(mat,1,which.max),stringsAsFactors=F)
    colnames(ll) <- c('whichMax')
    ll$gene <- rownames(ll)
    ll$scaleExp <- NA
    for(idx in 1:nrow(ll)){
        gg <- ll[idx,'gene']
        cc <- as.numeric(ll[idx,'whichMax'])
        ll[idx,'scaleExp'] <- mat[gg,cc]
    }
    ll <- ll[order(ll$whichMax,ll$scaleExp,decreasing=T),]
    return(rownames(ll))
}

full_order_from_abbr <- function(df,abbr_col,full_col,abbr_order){
    ll <- unique(df[,c(abbr_col,full_col)])
    ll[,abbr_col] <- factor(ll[,abbr_col],levels=abbr_order)
    ll <- ll[order(ll[,abbr_col]),]
    full_order <- ll[,full_col]
    return(full_order)
}

reorder_peaks <- function(df,cutHeight,cutOrder,xCol='peak',yCol='CITE_abbr',valCol='value',
                          distMethod='euclidean',clustMethod='ward.D2'){
    
    if(!all(c(xCol,yCol,valCol) %in% colnames(df))) stop('not all columns in df')
    
    toMat <- spread(df[,c(xCol,yCol,valCol)],key=yCol,value=valCol)
    thisMat <- as.matrix(toMat[,2:ncol(toMat)])
    rownames(thisMat) <- toMat[,1]
    
    row_clust <- hclust(dist(thisMat, method = distMethod), method = clustMethod)
    cuttree <- cutree(row_clust,h=cutHeight)

    heatmap_peak_order <- rev(row_clust$labels[row_clust$order])
    heatmap_peak_reorder <- c()
    for(idx in rev(cutOrder)){ 
        thesePeaks <- names(cuttree[which(cuttree==idx)])
        thesePeaks <- factor(thesePeaks,levels=heatmap_peak_order)
        thesePeaks <- thesePeaks[order(thesePeaks)]
        heatmap_peak_reorder <- c(heatmap_peak_reorder,as.character(thesePeaks))
    }
    if(!identical(sort(heatmap_peak_order),sort(heatmap_peak_reorder))) stop('peak order issue')
    
    return(heatmap_peak_reorder)
}

fix_hyphen_names <- function(s){
    if(str_detect(s,'[A-Za-z]{2}[0-9]+')){paste(sep='',substr(s,1,2),'-',substr(s,3,nchar(s)))}
    else if(str_detect(s,'[A-Za-z][0-9]+')){paste(sep='',substr(s,1,1),'-',substr(s,2,nchar(s)))}
    else{s}
}


##CELL COUNTS

#CELL COUNTS BY SAMPLE BAR PLOT
cellCount_bySample_barPlot <- function(thisMeta,thisCol,res_toPlot,fLab,cluster_colors,assayCol='assay',colorCol='disease',
                                       yLab_marker='>',res_order=NA,
                                       colOrder=c('BRI-449','BRI-450',
                                                  'BRI-542','BRI-544','BRI-545','BRI-546','BRI-547',
                                                  'BRI-639','BRI-640','BRI-641','BRI-643','BRI-644','BRI-646','BRI-647',
                                                  'BRI-448','BRI-543','BRI-642','BRI-645',
                                                  'BRI-1281','BRI-1385','BRI-1387','BRI-1389',
                                                  'BRI-1586','BRI-1588','BRI-1590','BRI-1592','BRI-1594','BRI-1607',
                                                  'BRI-1609','BRI-1611'),
                                       OA_samples=c('BRI-448', 'BRI-543', 'BRI-642', 'BRI-645', 'BRI-1281')){

    if(!any(colOrder %in% thisMeta[,thisCol])){
        colOrder <- mixedsort_new(unique(thisMeta[,thisCol]))
        modalityLine <- length(unique(thisMeta[which(thisMeta[,assayCol]=='snATAC'),thisCol]))+0.5
    } else {
        modalityLine <- 12.5
    }
    thisMeta[,thisCol] <- factor(thisMeta[,thisCol],levels=colOrder)
    
    toPlot <- as.data.frame(table(thisMeta[,c(thisCol,res_toPlot)]),stringsAsFactors=FALSE)
    if(!all(unique(toPlot[,res_toPlot]) %in% res_order)) res_order <- sort(unique(toPlot[,res_toPlot]))
    toPlot[,res_toPlot] <- factor(toPlot[,res_toPlot],levels=res_order)
    toPlot[,thisCol] <- factor(toPlot[,thisCol],levels=rev(colOrder))
    if(colorCol=='disease' & thisCol=='sample'){
        sample_yLab <- rep('',length(colOrder))
        names(sample_yLab) <- colOrder
        sample_yLab[which(names(sample_yLab) %in% OA_samples)] <- yLab_marker
        
        conv_df <- data.frame('sample'=colOrder,'disease'='RA',stringsAsFactors=FALSE)
        conv_df[which(conv_df$sample %in% OA_samples),'disease'] <- 'OA'
        toPlot[,colorCol] <- mapvalues(toPlot[,thisCol],from=conv_df$sample,to=conv_df$disease)
    }

    g <- ggplot(toPlot,aes_string(fill=res_toPlot,x='Freq',y=thisCol,color=colorCol)) + 
            geom_bar(position="stack", stat="identity") + theme_classic(base_size=25) + 
            labs(x='Cell Counts',y='Sample\nmultimodal                   unimodal',fill=fLab) + 
            guides(fill=guide_legend(order=1),color=guide_legend(order=2)) + 
            geom_hline(yintercept=modalityLine,linetype='dashed') + guides(color=guide_legend(override.aes=list(fill='white')))
    
    if(colorCol=='disease' & thisCol=='sample'){
        g <- g + scale_y_discrete(labels=sample_yLab) + labs(color='Disease') + 
                    scale_color_manual(values=c(NA,'grey15'),labels=c('RA',paste(sep='',yLab_marker,'OA'))) + 
                    theme(axis.ticks.y=element_blank())
        toSave <- toPlot[,c(thisCol,'disease',res_toPlot,'Freq')]
        colnames(toSave) <- c('Sample','Disease',gsub('\n',' ',fLab),'Cell Counts')
    } else {
        toSave <- toPlot[,c(thisCol,res_toPlot,'Freq')]
        colnames(toSave) <- c('Sample',gsub('\n',' ',fLab),'Cell Counts')
    }

    if(all(unique(toPlot[,res_toPlot]) %in% names(cluster_colors))){g <- g + scale_fill_manual(values=cluster_colors)}
    
    return(list('plot'=g,'data'=toSave))

}


##MARKER PEAK/GENE FUNCTIONS

#features by cell on UMAP
plot_markerPeaks_norm_v2 <- function (thisMeta, thisGxC, xCol, yCol, plotCol = 5, titleSize = 30, ptSize = 0.1, cutCap=0.01, 
                                   plot_genes = c("CD3D", "GNLY", "MS4A1", "XBP1", "C1QA",
                                                  "VWF", "PDPN", "PRG4", "THY1", "NOTCH3"), 
                                   orderValues = FALSE, colorOpt = "viridis",exclScale=TRUE,scaleLab=NA,titleFace='italic') 
{
    if (!identical(rownames(thisMeta), colnames(thisGxC))) 
        stop("cells have to be ordered the same.")
    toPlot <- thisMeta[1:nrow(thisMeta), 1:ncol(thisMeta)]
    myplots <- list()
    for (i in 1:length(plot_genes)) {
        gene <- plot_genes[i]
        
        if(cutCap!=0){
            max.cutoff = quantile(thisGxC[gene, ], 1-cutCap)
            min.cutoff = quantile(thisGxC[gene, ], cutCap)
            tmp <- sapply(X = thisGxC[gene, ], FUN = function(x) {
                return(ifelse(test = x > max.cutoff, yes = max.cutoff, 
                    no = x))
            })
            tmp <- sapply(X = tmp, FUN = function(x) {
                return(ifelse(test = x < min.cutoff, yes = min.cutoff, 
                    no = x))
            })
            toPlot$gene <- as.numeric(tmp)
        } else {
            toPlot$gene <- thisGxC[gene, ]
        }
        
        ind <- paste("p", i, sep = "")
        if (orderValues) {
            ind <- ggplot(data = toPlot[order(toPlot$gene), ], 
                aes_string(x = xCol, y = yCol))
        }
        else {
            ind <- ggplot(data = toPlot[sample(nrow(toPlot)), 
                ], aes_string(x = xCol, y = yCol))
        }
        ind <- ind + geom_point(mapping = aes(color = gene), 
            size = ptSize) + scale_color_viridis(option = colorOpt, 
            end = 0.9) + labs(x = "", y = "") + theme_bw(base_size = 15) + 
            theme(axis.text = element_blank(), axis.ticks = element_blank(), 
                panel.grid = element_blank(), plot.title = element_text(color = "black", 
                  size = titleSize, face = titleFace)) + labs(title = gene)
        if(exclScale){
            ind <- ind + theme(legend.position = "none")
        } else{
            if(!is.na(scaleLab)) ind <- ind + labs(color=scaleLab)
        }
            
        myplots[[i]] <- ind
    }
    
    p <- arrangeGrob(grobs = myplots, ncol=plotCol)
    
    return(p)
}

#features by hex on UMAP
plot_markerPeaks_norm_hex_v2 <- function(thisMeta, thisGxC, xCol, yCol, plotCol=5, titleSize=30, cutCap=0.01,
                                      colorOpt='viridis',hex_bins=30,exclScale=TRUE,scaleLab=NA,
                                      plot_genes=c('CD3D','GNLY','MS4A1','XBP1','C1QA',
                                                   'VWF','PDPN','PRG4','THY1','NOTCH3'), titleFace='italic'){

    if(!identical(rownames(thisMeta),colnames(thisGxC))) stop("cells have to be ordered the same.")

    toPlot <- thisMeta[1:nrow(thisMeta),1:ncol(thisMeta)]

    myplots <- list()
    for (i in 1:length(plot_genes)) {
        gene <- plot_genes[i]
        
        if(cutCap!=0){
            max.cutoff = quantile(thisGxC[gene, ], 1-cutCap)
            min.cutoff = quantile(thisGxC[gene, ], cutCap)
            tmp <- sapply(X = thisGxC[gene, ], FUN = function(x) {
                return(ifelse(test = x > max.cutoff, yes = max.cutoff, 
                    no = x))
            })
            tmp <- sapply(X = tmp, FUN = function(x) {
                return(ifelse(test = x < min.cutoff, yes = min.cutoff, 
                    no = x))
            })
            toPlot$gene <- as.numeric(tmp)
        } else {
            toPlot$gene <- thisGxC[gene, ]
        }

        ind <- paste("p", i, sep = "")
        ind <- ggplot(data = toPlot,aes_string(x = xCol, y = yCol, z = "gene")) +
          stat_summary_hex(bins=hex_bins) + scale_fill_viridis(option = colorOpt) +
          labs(x="", y="") + theme_bw(base_size = 15) +
          theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank(),
            plot.title = element_text(color="black", size=titleSize, face = titleFace)) + labs(title = gene)
        if(exclScale){
            ind <- ind + theme(legend.position = "none")
        } else{
            if(!is.na(scaleLab)) ind <- ind + labs(fill=scaleLab)
        }
        myplots[[i]] <- ind

    }

    p <- arrangeGrob(grobs = myplots, ncol=plotCol)
    
    return(p)

}

#Process normalized matrices for heatmaps
scaleFeat_forHeatmap <- function(gOrd,ctOrd,gp_map,gCTnorm,pCTnorm,geneCTcutoff=0.15,peakCTcutoff=0.05){
    
    if(!all(gOrd %in% rownames(gCTnorm))) stop('genes from order must be in normalized gene expression')
    if(!all(ctOrd %in% colnames(gCTnorm))) stop('cell types from order must be in normalized gene expression')
    if(!all(unname(gp_map) %in% rownames(pCTnorm))) stop('peaks from gene/peak relationship must be in normalized peak accessibility')
    if(!all(ctOrd %in% colnames(pCTnorm))) stop('cell types from order must be in normalized peak accessibility')
    
    ll <- apply(gCTnorm[gOrd,],1,max)
    gOrd <- names(ll[ll>=geneCTcutoff])
    rr <- apply(pCTnorm[unname(gp_map[which(names(gp_map) %in% gOrd)]),],1,max)
    names(rr) <- mapvalues(names(rr),from=unname(gp_map),to=names(gp_map),warn_missing=FALSE)

    peaks_toAdd <- gOrd[which(!(gOrd %in% names(gp_map)))]
    peaks_toWhite <- c(peaks_toAdd,names(rr[rr<peakCTcutoff]))
    
    gCTnorm_subset <- gCTnorm[rev(gOrd),rev(ctOrd)]

    ll <- gp_map[which(names(gp_map) %in% gOrd)]
    pCTnorm_subset <- pCTnorm[unname(ll),rev(ctOrd)]
    rownames(pCTnorm_subset) <- names(ll)
    
    if(length(peaks_toAdd)>0){
        toAdd <- matrix(0,nrow=length(peaks_toAdd),ncol=ncol(pCTnorm_subset))
        rownames(toAdd) <- peaks_toAdd
        colnames(toAdd) <- colnames(pCTnorm_subset)

        pCTnorm_subset <- rbind(pCTnorm_subset,toAdd)
    }
    
    if(!identical(sort(rownames(gCTnorm_subset)),sort(rownames(pCTnorm_subset)))) stop('rowname issue')
    pCTnorm_subset <- pCTnorm_subset[rownames(gCTnorm_subset),]
    if(length(peaks_toWhite)>0) pCTnorm_subset[peaks_toWhite,] <- 0
    
    gCTnorm_subset_scaled <- t(scale(t(gCTnorm_subset)))
    pCTnorm_subset_scaled <- t(scale(t(pCTnorm_subset)))
    
    gCTnorm_subset_scaled_toGather <- as.data.frame(gCTnorm_subset_scaled)
    gCTnorm_subset_scaled_toGather$gene <- rownames(gCTnorm_subset_scaled_toGather)
    gCTnorm_subset_scaled_toGather <- gather(gCTnorm_subset_scaled_toGather,'cluster_abbr','gene_norm_scale',
                                             all_of(colnames(gCTnorm_subset_scaled)))

    pCTnorm_subset_scaled_toGather <- as.data.frame(pCTnorm_subset_scaled)
    pCTnorm_subset_scaled_toGather$gene <- rownames(pCTnorm_subset_scaled_toGather)
    pCTnorm_subset_scaled_toGather <- gather(pCTnorm_subset_scaled_toGather,'cluster_abbr','peak_norm_scale',
                                             all_of(colnames(pCTnorm_subset_scaled)))

    if(!identical(gCTnorm_subset_scaled_toGather$gene,pCTnorm_subset_scaled_toGather$gene)) stop('genes do not match')
    if(!identical(gCTnorm_subset_scaled_toGather$cluster_abbr,pCTnorm_subset_scaled_toGather$cluster_abbr)) stop('clusters do not match')
    
    fxCT_norm_subset_scaled_gathered <- cbind(gCTnorm_subset_scaled_toGather,
                                              'peak_norm_scale'=pCTnorm_subset_scaled_toGather$peak_norm_scale)

    
    return(list('gxCT_norm_subset_scaled'=gCTnorm_subset_scaled,
                'pxCT_norm_subset_scaled'=pCTnorm_subset_scaled,
                'fxCT_norm_subset_scaled'=fxCT_norm_subset_scaled_gathered))
    
}

#feature heatmaps
#assumes that the given row/col is the order you want!
pseudobulk_scaled_heatmap <- function(toPlot,xlab,ylab,fillLab,plotTit=NULL,scale_lim=NA,clustColors=NA,includeZero=FALSE,
                                      colorLow='blue',colorMid='white',colorHigh='red'){
    
    toGather <- as.data.frame(toPlot)
    cols_toGather <- colnames(toGather)
    toGather$peak <- rownames(toGather)
    
    gathered <- gather(toGather,'CT','agg',all_of(cols_toGather))
    gathered$peak <- factor(gathered$peak,levels=rev(rownames(toGather)))
    gathered$CT <- factor(gathered$CT,levels=cols_toGather)
    
    if(!is.na(scale_lim)){
        scale_limits <- c(-scale_lim,scale_lim)
    } else {
        if(min(gathered$agg)>0 & includeZero){
            ll <- 0
        } else {
            ll <- min(gathered$agg)
        }
        scale_limits <- c(ll,max(gathered$agg))
    }
    
    g <- ggplot(gathered,aes_string(x='peak',y='CT',fill='agg')) + geom_tile() +
            theme_classic(base_size=20) + labs(x=xlab,y=ylab,fill=fillLab) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5)) + 
            scale_fill_gradient2(high=colorHigh,mid=colorMid,low=colorLow,midpoint=0,limits=scale_limits)
    if(all(levels(gathered$CT) %in% names(clustColors))){
        suppressWarnings(g <- g + theme(axis.text.y = element_text(color=clustColors[levels(gathered$CT)],face='bold',size=20)))
    }
    if(!is.null(plotTit)) {g <- g + ggtitle(plotTit)}
    
    return(g)
    
}

get_differential_features <- function(dP_df,dG_df,name_abbr_df,gxCT_nss,withinE=0.75,top_n=5,
                                      log2FC_cutoff=0.5,logFC_cutoff=0.5,log10padj_cutoff=5){
    
    #within max for class markers; some markers are useful for multiple classes
    this_gxCT_scaled <- as.data.frame(gxCT_nss)
    this_gxCT_scaled$maxE2 <- apply(this_gxCT_scaled,1,function(x){sort(x,decreasing=TRUE)[2]})
    this_gxCT_scaled$maxE <- apply(this_gxCT_scaled,1,max)
    this_gxCT_scaled$maxE2_cutoff <- this_gxCT_scaled$maxE*withinE
    this_gxCT_scaled$withinMaxE <- this_gxCT_scaled$maxE2>=this_gxCT_scaled$maxE2_cutoff
    
    dF_df <- data.frame('class'=character(),'peak'=character(),'gene'=character(),
                                 'peakFC'=numeric(),'peakPadj'=numeric(),'geneFC'=numeric(),'genePadj'=numeric(),
                                 'topPeak'=logical(),'topGene'=logical(),'marker'=logical(),stringsAsFactors=FALSE)

    for(cl in unique(dP_df$cellType)){
        #subset markers by class
        full_name <- name_abbr_df[which(name_abbr_df$cluster_abbr==cl),'cluster_name']
        markGenes <- unique(c(rownames(this_gxCT_scaled[which(this_gxCT_scaled[,cl]>=this_gxCT_scaled$maxE2_cutoff),]), 
                              sub('\\+','',grep('[A-Za-z0-9]+\\+',unlist(str_split(full_name,' ')),value=TRUE)),
                              sub('hi','A',grep('[A-Za-z0-9]+hi',unlist(str_split(full_name,' ')),value=TRUE))))
        markPeaks <- chosenPeaks[markGenes]

        #chosen markers
        mp <- dP_df[which(dP_df$cellType==cl & dP_df$peak %in% unname(markPeaks)),]
        #in case a peak overlaps multiple TSS
        mp$gene <- mapvalues(mp$peak,from=unname(markPeaks),names(markPeaks),warn_missing=FALSE)
        mg <- dG_df[which(dG_df$group==cl & dG_df$feature %in% markGenes),]
        mf <- merge(mg,mp,by.x='feature',by.y='gene',all=TRUE)
        #making a new dataframe for this class
        toAdd <- data.frame('class'=rep(cl,nrow(mf)),'peak'=mf$peak,'gene'=mf$feature,
                            'peakFC'=mf$log2FC,'peakPadj'=mf$log10padj.y,
                            'geneFC'=mf$logFC,'genePadj'=mf$log10padj.x,
                            'topPeak'=rep(FALSE,nrow(mf)),'topGene'=rep(FALSE,nrow(mf)),'marker'=rep(TRUE,nrow(mf)),
                            stringsAsFactors=FALSE)

        #top peaks
        tp <- dP_df[which(dP_df$cellType==cl & dP_df$log2FC>log2FC_cutoff & dP_df$log10padj>log10padj_cutoff),]
        tp <- head(tp[order(tp$log10padj,decreasing=TRUE),],n=top_n)
        #don't need to readd peaks accounted for in markers
        toAdd[which(toAdd$peak %in% tp$peak),'topPeak'] <- TRUE
        tp <- tp[which(!(tp$peak %in% toAdd$peak)),]
        if(nrow(tp)>0){
            #get gene info
            tp_genes <- dG_df[which(dG_df$group==cl & dG_df$feature %in% tp$gene),]
            for(gn in tp[which(!(tp$gene %in% tp_genes$feature)),'gene']) {tp_genes <- rbind(tp_genes,c(gn,cl,rep(NA,9)))}
            rownames(tp_genes) <- tp_genes$feature
            tp_genes <- tp_genes[tp$gene,]
            toAdd <- rbind(toAdd,data.frame('class'=rep(cl,nrow(tp)),'peak'=tp$peak,'gene'=tp$gene,
                                            'peakFC'=tp$log2FC,'peakPadj'=tp$log10padj,
                                            'geneFC'=tp_genes$logFC,'genePadj'=tp_genes$log10padj,
                                            'topPeak'=rep(TRUE,nrow(tp)),'topGene'=rep(FALSE,nrow(tp)),
                                            'marker'=rep(FALSE,nrow(tp)),stringsAsFactors=FALSE))
        }

        #top genes
        tg <- dG_df[which(dG_df$group==cl & dG_df$logFC>logFC_cutoff & dG_df$log10padj>log10padj_cutoff),]
        tg <- head(tg[order(tg$log10padj,tg$logFC,decreasing=TRUE),],n=top_n)
        #don't need to readd genes already accounted for
        toAdd[which(toAdd$gene %in% tg$feature),'topGene'] <- TRUE
        tg <- tg[which(!(tg$feature %in% toAdd$gene)),]
        if(nrow(tg)>0){
            #assign peaks to genes
            tg$peak <- NA; tg$peakFC <- NA; tg$peakPadj <- NA
            for(gn in tg$feature){
                tg_peaks <- dP_df[which(dP_df$cellType==cl & dP_df$gene==gn),]
                if(nrow(tg_peaks)==0){
                    #if no genes, leave NA
                    next
                } else {
                    #choose most significant - not including log2FC check here
                    tg_peaks <- tg_peaks[order(tg_peaks$log10padj,decreasing=TRUE),]
                    tg[which(tg$feature==gn),'peak'] <- tg_peaks[1,'peak']
                    tg[which(tg$feature==gn),'peakFC'] <- tg_peaks[1,'log2FC']
                    tg[which(tg$feature==gn),'peakPadj'] <- tg_peaks[1,'log10padj']
                }
            }
            toAdd <- rbind(toAdd,data.frame('class'=rep(cl,nrow(tg)),'peak'=tg$peak,'gene'=tg$feature,
                                            'peakFC'=tg$peakFC,'peakPadj'=tg$peakPadj,
                                            'geneFC'=tg$logFC,'genePadj'=tg$log10padj,
                                            'topPeak'=rep(FALSE,nrow(tg)),'topGene'=rep(TRUE,nrow(tg)),
                                            'marker'=rep(FALSE,nrow(tg)),stringsAsFactors=FALSE))
        }

        dF_df <- rbind(dF_df,toAdd)
    }
    
    return(dF_df)
}


##SYMPHONY RELATED FUNCTIONS

#proportions
symp_prop_df <- function(queryDF,refDF,xLab,yLab,cLab,tLab='',
                         qSamp_col='sample',qCT_col='CITE_abbr',rSamp_col='sample',rCT_col='cluster_abbr',
                         clustColors=NA,expand_flag=FALSE){
    
    query_CT_sample_freq_df <- as.data.frame(table(queryDF[,c(qCT_col,qSamp_col)]),stringsAsFactors=FALSE)
    query_sample_freq_df <- as.data.frame(table(queryDF[,qSamp_col]),stringsAsFactors=FALSE)
    query_CT_sample_freq_df$sample_tot <- as.numeric(mapvalues(query_CT_sample_freq_df[,qSamp_col],
                                                               from=query_sample_freq_df$Var1,to=query_sample_freq_df$Freq))
    query_CT_sample_freq_df$perc <- query_CT_sample_freq_df$Freq/query_CT_sample_freq_df$sample_tot
    if(!all(for(idx in unique(queryDF[,qSamp_col])){
        sum(query_CT_sample_freq_df[which(query_CT_sample_freq_df[,qSamp_col]==idx),'perc'])==1})) stop('not all query fractions sum to 1')
    query_CT_freq_df <- as.data.frame(table(queryDF[,qCT_col]),stringsAsFactors=FALSE)
    query_CT_freq_df$perc <- query_CT_freq_df$Freq/nrow(queryDF)
    if(sum(query_CT_freq_df$perc)!=1) stop('query percentages do not sum to 1')
    
    query_CT_bySample_stats <- summarySE(query_CT_sample_freq_df, 'perc', groupvars = qCT_col)
    if(sum(query_CT_bySample_stats$perc)!=1) stop('query percentages in stats do not sum to 1')
    
    
    ref_CT_sample_freq_df <- as.data.frame(table(refDF[,c(rCT_col,rSamp_col)]),stringsAsFactors=FALSE)
    ref_sample_freq_df <- as.data.frame(table(refDF[,rSamp_col]),stringsAsFactors=FALSE)
    ref_CT_sample_freq_df$sample_tot <- as.numeric(mapvalues(ref_CT_sample_freq_df[,rSamp_col],
                                                               from=ref_sample_freq_df$Var1,to=ref_sample_freq_df$Freq))
    ref_CT_sample_freq_df$perc <- ref_CT_sample_freq_df$Freq/ref_CT_sample_freq_df$sample_tot
    if(!all(for(idx in unique(refDF[,rSamp_col])){
        sum(ref_CT_sample_freq_df[which(ref_CT_sample_freq_df[,rSamp_col]==idx),'perc'])==1})) stop('not all ref fractions sum to 1')
    ref_CT_freq_df <- as.data.frame(table(refDF[,rCT_col]),stringsAsFactors=FALSE)
    ref_CT_freq_df$perc <- ref_CT_freq_df$Freq/nrow(refDF)
    if(sum(ref_CT_freq_df$perc)!=1) stop('ref percentages do not sum to 1')
    
    ref_CT_bySample_stats <- summarySE(ref_CT_sample_freq_df, 'perc', groupvars = rCT_col)
    if(sum(ref_CT_bySample_stats$perc)!=1) stop('ref percentages in stats do not sum to 1')
    
    if(!identical(ref_CT_bySample_stats[,rCT_col],query_CT_bySample_stats[,qCT_col])) stop('cell states do not match between stats')
    
    both_CT_bySample_stats <- ref_CT_bySample_stats[,c(rCT_col,'perc','sd','se','ci')]
    colnames(both_CT_bySample_stats) <- c('cellType','ref_mean','ref_sd','ref_se','ref_ci')
    both_CT_bySample_stats <- cbind(both_CT_bySample_stats,query_CT_bySample_stats[,c('perc','sd','se','ci')])
    colnames(both_CT_bySample_stats) <- c('cellType','ref_mean','ref_sd','ref_se','ref_ci',
                                          'qry_mean','qry_sd','qry_se','qry_ci')
    
    
    ll <- cor.test(both_CT_bySample_stats$ref_mean,both_CT_bySample_stats$qry_mean)

    g <- ggplot(both_CT_bySample_stats,aes_string(x='qry_mean',y='ref_mean',color='cellType')) + geom_point(size=3) + 
            theme_bw(base_size=22) + 
            labs(x=xLab,y=yLab,color=cLab) +
            ggtitle(paste(sep='',tLab,'R=',round(ll$estimate,2),' P=',signif(ll$p.value,3))) + 
            theme(plot.title=element_text(hjust = 0.5)) +
            geom_abline(slope=1,intercept=0, linetype='dashed') + 
            expand_limits(x=0) + expand_limits(y=0) + 
            geom_errorbarh(aes(xmin = qry_mean-qry_se,xmax = qry_mean+qry_se)) + 
            geom_errorbar(aes(ymin = ref_mean-ref_se,ymax = ref_mean+ref_se))
    if(all(unique(both_CT_bySample_stats$cellType) %in% names(clustColors))){g <- g + scale_color_manual(values=clustColors)}
    if(expand_flag){
        expand_lim <- max(both_CT_bySample_stats$qry_mean+both_CT_bySample_stats$qry_se,
                          both_CT_bySample_stats$ref_mean+both_CT_bySample_stats$ref_se)
        expand_lim <- round_any(expand_lim,0.05,f = ceiling)
        
        g <- g + expand_limits(x=expand_lim) + expand_limits(y=expand_lim)
    }
    
    toSave <- both_CT_bySample_stats[,c('cellType','qry_mean','ref_mean')]
    toSave$xMin_err <- both_CT_bySample_stats$qry_mean-both_CT_bySample_stats$qry_se
    toSave$xMax_err <- both_CT_bySample_stats$qry_mean+both_CT_bySample_stats$qry_se
    toSave$yMin_err <- both_CT_bySample_stats$ref_mean-both_CT_bySample_stats$ref_se
    toSave$yMax_err <- both_CT_bySample_stats$ref_mean+both_CT_bySample_stats$ref_se
    qry_prefix <- str_split_fixed(xLab,' ',2)[,1]
    ref_prefix <- str_split_fixed(yLab,' ',2)[,1]
    colnames(toSave) <- c(gsub('\n',' ',cLab),gsub('\n',' ',xLab),gsub('\n',' ',yLab),
                          paste(qry_prefix,'mean-se'),paste(qry_prefix,'mean+se'),
                          paste(ref_prefix,'mean-se'),paste(ref_prefix,'mean+se'))
    return(list('plot'=g,'data'=toSave))
        
}

#accuracy of mapping score
symphony_mapScore_byCT <- function(inDir,cs_suffix='_class_state_df.rds',meta_suffix='_ATAC_meta.rds',
                                   CTs=c('Bplasma','endothelial','stromal','myeloid','Tcell'),
                                   assay_col='assay',assayVal='snATAC',
                                   val_col='CITE_prob',outVal='out',inVal='in',
                                   classA_col='cluster_abbr',stateA_col='CITE_abbr',
                                   classB_col='class',stateB_col='state'){
    
    if(!file.exists(inDir)) stop('Input directory does not exist.')
    
    df <- data.frame('cellType'=character(),'group'=character(),'perc'=numeric(),stringsAsFactors=FALSE)
    for(CT in CTs){
        meta_file <- paste(sep='',inDir,CT,'/',CT,meta_suffix)
        class_state_file <- paste(sep='',inDir,CT,'/',CT,cs_suffix)
        if(!all(file.exists(c(meta_file,class_state_file)))) stop('Input file(s) do not exist.')
        
        meta <- readRDS(meta_file)
        class_state_df <- readRDS(class_state_file)
        
        if(!all(c(assay_col,val_col,stateA_col,classA_col) %in% colnames(meta))) stop('meta column(s) missing')
        if(!all(c(stateB_col,classB_col) %in% colnames(class_state_df))) stop('class/state column(s) missing')

        meta <- meta[which(meta[,assay_col]==assayVal),]

        meta$group <- outVal
        for(idx in 1:nrow(class_state_df)){
            meta[which(meta[,stateA_col]==class_state_df[idx,stateB_col] & 
                       meta[,classA_col]==class_state_df[idx,classB_col]),'group'] <- inVal
        }
        
        bin_col <- paste(sep='_',val_col,'binary')
        meta[,bin_col] <- ifelse(meta[,val_col]==1,1,0)

        in_group_perc <- nrow(meta[which(meta[,bin_col]==1 & meta$group==inVal),])/nrow(meta[which(meta$group==inVal),])
        out_group_perc <- nrow(meta[which(meta[,bin_col]==1 & meta$group==outVal),])/nrow(meta[which(meta$group==outVal),])

        df <- rbind(df,data.frame('cellType'=c(CT,CT),'group'=c(inVal,outVal),'perc'=c(in_group_perc,out_group_perc),
                                  stringsAsFactors=FALSE))

    }
    
    return(df)
}


##ODDS RATIO FUNCTIONS

#OR PLOT ORDER - plot classes/states in order
plot_order <- function(uniq_val){
    if(all(uniq_val %in% grep('[a-zA-Z]+-[0-9]+:',uniq_val,value=TRUE)))
        ret <- mixedsort_new(uniq_val)
    else if(!any(is.na(as.numeric(uniq_val))))
        ret <- as.character(sort(as.numeric(uniq_val)))
    else
        ret <- sort(uniq_val)
    
    return(ret)
}

#OR PLOT ORDER - plot classes/states such that the diagonal has the highest OR
reorder_col_diag_plotOR <- function(fish_df,xCol,yCol,yOrd=NA,mCol='lnOR',op='max'){
    if(!(op %in% c('min','max'))) op='max' #defaulting to max.
    
    if(!any(is.na(yOrd))){fish_df[,yCol] <- factor(fish_df[,yCol],levels=rev(yOrd))}
    
    row_col_maxOR_df <- data.frame('row'=character(),'col'=character(),'val'=numeric(),stringsAsFactors=FALSE)
    for(thisX in unique(fish_df[,xCol])){
        this_subset <- fish_df[which(fish_df[,xCol]==thisX),]
        if(op=='min'){
            this_idx <- which.min(this_subset[,mCol])
        } else {
            this_idx <- which.max(this_subset[,mCol])
        }
        row_col_maxOR_df <- rbind(row_col_maxOR_df,data.frame('row'=this_subset[this_idx,yCol],
                                                              'col'=thisX,'val'=this_subset[this_idx,mCol],
                                                              stringsAsFactors=FALSE))
    }
    #reverse row order since that's what the plot does!
    row_col_maxOR_df <- row_col_maxOR_df[order(row_col_maxOR_df$row,decreasing=TRUE),]
    new_xOrd <- row_col_maxOR_df$col
    return(new_xOrd)
}

#For plotting purposes
label_spacing <- function(s,wiggle=1){
    if(substr(s,nchar(s)-wiggle,nchar(s)-wiggle)=='\n'){substr(s,nchar(s)-wiggle,nchar(s)-wiggle) <- ' '}
    return(s)
}

#For plotting purposes
fix_infinite <- function(this_vec,addVal=1,multVal=1){
    nonInf <- this_vec[which(!is.infinite(this_vec))]
    limit <- ceiling(max(nonInf,abs(min(nonInf)))*multVal)+addVal
    this_vec[which(is.infinite(this_vec) & this_vec>0)] <- limit
    this_vec[which(is.infinite(this_vec) & this_vec<0)] <- -limit
    
    return(this_vec)
}

calc_OR <- function(toPlot, xCol, yCol, pLim = 0.05,includeNomP = FALSE){
    if(!all(c(xCol,yCol) %in% colnames(toPlot))) stop('xCol and yCol need to be in toPlot columns')
    
    fisher_res_df <- data.frame('xCol'=character(),'yCol'=character(),
                                'inX_inY'=integer(),'inX_noY'=integer(),'noX_inY'=integer(),'noX_noY'=integer(),
                                'pval'=numeric(),'OR'=numeric(),'nullOR'=numeric(),
                                'CI_low'=numeric(),'CI_high'=numeric(),stringsAsFactors=FALSE)
    colnames(fisher_res_df)[1] <- xCol
    colnames(fisher_res_df)[2] <- yCol
    for(xVal in sort(unique(toPlot[,xCol]))){
        for(yVal in sort(unique(toPlot[,yCol]))){
            inX_inY <- nrow(toPlot[which(toPlot[,xCol]==xVal & toPlot[,yCol]==yVal),])
            inX_noY <- nrow(toPlot[which(toPlot[,xCol]==xVal & toPlot[,yCol]!=yVal),])
            noX_inY <- nrow(toPlot[which(toPlot[,xCol]!=xVal & toPlot[,yCol]==yVal),])
            noX_noY <- nrow(toPlot[which(toPlot[,xCol]!=xVal & toPlot[,yCol]!=yVal),])

            fish_df <- data.frame('inX'=c(inX_inY,inX_noY),'noX'=c(noX_inY,noX_noY),
                                  row.names=c("inY", "noY"),
                                  stringsAsFactors = FALSE)

            fish_test <- fisher.test(fish_df)

            placeholder <- data.frame('xCol'=xVal,'yCol'=yVal,
                                      'inX_inY'=inX_inY,'inX_noY'=inX_noY,
                                      'noX_inY'=noX_inY,'noX_noY'=noX_noY,
                                      'pval'=fish_test$p.value,'OR'=fish_test$estimate,
                                      'nullOR'=fish_test$null.value,'CI_low'=fish_test$conf.int[1],
                                      'CI_high'=fish_test$conf.int[2],stringsAsFactors=FALSE)
            colnames(placeholder)[1] <- xCol
            colnames(placeholder)[2] <- yCol

            fisher_res_df <- rbind(fisher_res_df,placeholder)
        }
    }
    rownames(fisher_res_df) <- NULL
    head(fisher_res_df)

    fisher_res_df$sig <- ifelse(fisher_res_df$pval<=pLim,TRUE,FALSE)
    fisher_res_df$padj <- p.adjust(fisher_res_df$pval,method='BH')
    fisher_res_df$sigBH <- ifelse(fisher_res_df$padj<=pLim,TRUE,FALSE)
    fisher_res_df$signif <- 'not'
    if(includeNomP){
        fisher_res_df[which(fisher_res_df$sigBH==FALSE & fisher_res_df$sig==TRUE),'signif'] <- 'nominal'
    }
    fisher_res_df[which(fisher_res_df$sigBH==TRUE),'signif'] <- 'BH'

    fisher_res_df$lnOR <- fix_infinite(log(fisher_res_df$OR))
    
    return(fisher_res_df)

}

#OR - calculate and plot OR values between two columns of a dataframe
plot_OR <- function(toPlot, xCol, yCol, xLab, yLab, xOrder, yOrder, xAxisVert = FALSE, fCol='lnOR', flab='ln(OR)',clustColors=NA,
                    wiggle=1,yLab_charLim=50){
    
    if(!all(c(xCol,yCol,fCol) %in% colnames(toPlot))) stop('xCol and yCol need to be in toPlot columns')
    
    if(!all(unique(toPlot[,xCol]) %in% xOrder)){
        print('xOrder not sufficient, just sorting')
        xOrder <- sort(unique(toPlot[,xCol]))
    } else {
        xOrder <- xOrder[which(xOrder %in% unique(toPlot[,xCol]))]
    }
    
    if(!all(unique(toPlot[,yCol]) %in% yOrder)){
        print('yOrder not sufficient, just sorting')
        yOrder <- sort(unique(toPlot[,yCol]))
    } else {
        yOrder <- yOrder[which(yOrder %in% unique(toPlot[,yCol]))]
    }
    
    toPlot[,yCol] <- factor(toPlot[,yCol],levels=rev(yOrder))
    toPlot[,xCol] <- factor(toPlot[,xCol],levels=xOrder)
    
    xLabels_wSpaces <- str_replace_all(unique(toPlot[,xCol]),' ','\n')
    xLabels_wSpaces <- lapply(xLabels_wSpaces,FUN=label_spacing,wiggle=wiggle)
    names(xLabels_wSpaces) <- unique(toPlot[,xCol])
    
    yLabels_wSpaces <- c()
    for(yy in sort(unique(toPlot[,yCol]))){
        if(nchar(yy)>yLab_charLim){
            yLabels_wSpaces <- c(yLabels_wSpaces,replace_space_newline_afterHalf(yy,wiggle=wiggle))
        } else {
            yLabels_wSpaces <- c(yLabels_wSpaces,yy)
        }
    }
    names(yLabels_wSpaces) <- sort(unique(toPlot[,yCol]))
    
    
    g <- ggplot(toPlot,aes_string(x=xCol,y=yCol,fill=fCol)) + 
            geom_tile() + theme_classic(base_size=25) +
            scale_fill_gradient2(low = "grey85", mid = "white", high = "red", midpoint = 0, na.value='grey85') +
            labs(x=xLab,y=yLab,fill=flab) +
            geom_tile(data=toPlot[which(toPlot$signif=='not'),],fill='white') + 
            scale_x_discrete(labels=xLabels_wSpaces[levels(toPlot[,xCol])]) +
            scale_y_discrete(labels=yLabels_wSpaces[levels(toPlot[,yCol])])
    if(all(c(levels(toPlot[,xCol]),levels(toPlot[,yCol])) %in% names(clustColors))){
        suppressWarnings(g <- g + theme(axis.text.x = element_text(color=clustColors[levels(toPlot[,xCol])],
                                                                   face='bold',size=17)) + 
                                  theme(axis.text.y = element_text(color=clustColors[levels(toPlot[,yCol])],
                                                                   face='bold',size=17)))
    }
    
    return(g)
}


##LDA PLOTTING FUNCTIONS

#plot LDA results
LDA_plots <- function(predMod_df,thisCT,c01Lab,c0_col='c0',c1_col='c1',thisFeat='AUC',
                      class_state_df=NA,class_col='class',state_col='state',ctOrd_col=NA,
                      ctOrd=NA,toColor=NA,clustColors=NA,plot_square=TRUE){
    
    if(!all(c('cellType',c0_col,c1_col,'feature','value') %in% colnames(predMod_df))) stop('predMod_df columns: cellType, c0, c1, feature, value')
    
    all_LDA_states <- unique(c(predMod_df[,c0_col],predMod_df[,c1_col]))
    
    if(!any(is.na(class_state_df))){
        if(!all(c(class_col,state_col) %in% colnames(class_state_df))) stop('class_state_df must have class and state columns')
        if(length(class_state_df[,state_col])!=length(unique(class_state_df[,state_col]))) stop('state must be unique in class_state_df')
        if(!all(all_LDA_states %in% class_state_df[,state_col])) stop('all LDA states must between be in class_state_df')
        if(!all(ctOrd %in% unique(class_state_df[,class_col]))) stop('not all classes in ctOrd in class_state_df; if class_state_df is given, assumed that classes are given.')
        
        class_state_df <- class_state_df[which(class_state_df$state %in% all_LDA_states),]

        if(!(ctOrd_col %in% colnames(class_state_df))){
            print('no ctOrd_col given, so assuming the order given is desired')
            ctOrd_col <- 'intOrd'
            class_state_df[,ctOrd_col] <- 1:nrow(class_state_df)
        }
        
        class_state_df[,class_col] <- factor(class_state_df[,class_col],levels=ctOrd)
        ctOrd <- class_state_df[order(class_state_df[,class_col],class_state_df$intOrd),state_col]

    } else {
        if(!all(all_LDA_states %in% ctOrd)){
            ctOrd <- plot_order(all_LDA_states)
        }
    }
    
    if(any(is.na(toColor))){
        color_bins <- 6
        ramp <- colorRamp(c("darkred","red","grey90"))
        toColor <- rgb( ramp(seq(0, 1, length = color_bins)), max = 255)
    }
    
    toPlot <- predMod_df[which(predMod_df$cellType==thisCT & predMod_df$feature==thisFeat),
                         c('cellType',c0_col,c1_col,'feature','value')]

    if(plot_square){
        toAdd <- toPlot[,c('cellType',c1_col,c0_col,'feature','value')]
        colnames(toAdd) <- c('cellType',c0_col,c1_col,'feature','value')
        toPlot <- rbind(toPlot,toAdd)

        for(x in unique(toPlot[,c0_col])){
            ll <- data.frame('cellType'=CT,'c0_abbr'=x,'c1_abbr'=x,'feature'=thisFeat,'value'=NA)
            colnames(ll) <- c('cellType',c0_col,c1_col,'feature','value')
            toPlot <- rbind(toPlot,ll)
        }
    }

    toPlot[,c0_col] <- factor(toPlot[,c0_col], levels=ctOrd)
    toPlot[,c1_col] <- factor(toPlot[,c1_col], levels=rev(ctOrd))
    g <- ggplot(toPlot, aes_string(x=c0_col,y=c1_col,fill='value')) + geom_tile() + theme_classic(base_size=30) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
            labs(x=c01Lab,y=c01Lab,fill=thisFeat) +
            scale_fill_gradientn(colors=toColor,limit=c(0.5,1)) 
    if(all(unique(c(levels(toPlot[,c0_col]),levels(toPlot[,c1_col]))) %in% names(clustColors))){
        suppressWarnings(g <- g + theme(axis.text.x = element_text(color=clustColors[levels(toPlot[,c0_col])],face='bold')) +
                                  theme(axis.text.y = element_text(color=clustColors[levels(toPlot[,c1_col])],face='bold')))
    }
    
    if(plot_square & !any(is.na(class_state_df)) & all(unique(as.character(class_state_df$class)) %in% names(clustColors))){
        yMax <- nrow(class_state_df)+0.5
        xMin <- 0.5

        for(cc in sort(unique(class_state_df$class))){
            numStates <- nrow(class_state_df[which(class_state_df$class==cc),])
            yMin <- yMax-numStates
            xMax <- xMin+numStates
            g <- g + geom_rect(xmin=xMin, xmax=xMax, ymin=yMin, ymax=yMax,alpha=0,
                               color=unname(clustColors[as.character(cc)]),size=2.5) 
            yMax <- yMin
            xMin <- xMax
        }

    }
    
    return(g)
    
    
}

#state within/out class AUCs
LDA_class_state_AUC <- function(class_state_df,predMod_df,thisCT,c0_col='c0',c1_col='c1',simplifyCT=TRUE,thisFeat='AUC'){
    
    if(!all(c('cellType',c0_col,c1_col,'feature','value') %in% colnames(predMod_df))) stop('predMod_df columns: cellType, c0, c1, feature, value')
    if(!all(c('class','state') %in% colnames(class_state_df))) stop('class_state_df must have class and state columns')
    
    if(simplifyCT){
        predMod_df$c0_abbr <- str_split_fixed(predMod_df[,c0_col],":",2)[,1]
        predMod_df$c1_abbr <- str_split_fixed(predMod_df[,c1_col],":",2)[,1]
        
        c0_col <- 'c0_abbr'
        c1_col <- 'c1_abbr'
    }
    
    toPlot <- predMod_df[which(predMod_df$cellType==thisCT & predMod_df$feature==thisFeat),
                         c('cellType',c0_col,c1_col,'feature','value')]
    toPlot$withinCC <- FALSE

    for(cc in unique(class_state_df$class)){
        theseStates <- class_state_df[which(class_state_df$class==cc),'state']
        toPlot[which((toPlot[,c0_col] %in% theseStates) & (toPlot[,c1_col] %in% theseStates)),'withinCC'] <- TRUE
    }

    inCC <- toPlot[which(toPlot$withinCC==TRUE & !is.na(toPlot$value)),'value']
    outCC <- toPlot[which(toPlot$withinCC==FALSE & !is.na(toPlot$value)),'value']

    ll <- data.frame('relationship'=c('within','without'),
                     'numCombos'=c(length(inCC),length(outCC)),
                     'mean'=c(mean(inCC),mean(outCC)),
                     'median'=c(median(inCC),median(outCC)),
                     stringsAsFactors=FALSE)
    return(ll)
}

AUROC_plots <- function(AUC_pred,tLab,tSize=20,aXpos=0.75,aYpos=0.15,aSize=8){
    roc.perf = performance(AUC_pred, measure = "tpr", x.measure = "fpr")
    auc.train <- performance(AUC_pred, measure = "auc")
    auc.train <- auc.train@y.values
    auc_val <- auc.train[[1]]
    auc_val

    toPlot <- data.frame('FPR'=unlist(roc.perf@x.values),'TPR'=unlist(roc.perf@y.values),stringsAsFactors=FALSE)
    g <- ggplot(toPlot,aes_string(x='FPR',y='TPR')) + geom_line(size=0.8) + theme_bw(base_size=22) +
            geom_abline(slope=1,intercept=0,linetype='dashed',size=0.8) + 
            ggtitle(tLab) +
            theme(plot.title=element_text(size=tSize,hjust = 0.5),panel.grid = element_blank()) + 
            annotate("text",x=aXpos,y=aYpos,label=paste('AUC =', round(auc_val,2)),size=aSize)
    
    toSave <- toPlot[,c('FPR','TPR')]
    return(list('plot'=g,'data'=toSave))
}


##CITE-SEQ COMPARISONS

#STRING MANIPULATIONS
replace_space_newline_afterHalf <- function(s,wiggle=0){
    if(!grepl(' ', s, fixed = TRUE)) return(s)
    ll <- gregexpr(' ', s)[[1]]
    idx <- ll[ll > (nchar(s)/2-wiggle)][1]
    if(is.na(idx)) {idx <- ll[length(ll)]}
    substr(s, idx, idx) <- "\n"
    return(s)
}

#DONOR PROPORTIONS
donor_prop_comp_plot <- function(ID_conv_df,scATAC_meta,CITE_meta,
                                 cMrna_col='mrna_run',cAtac_col='atac_run',cID_col='subject_id',
                                 aSamp_col='sample',aCT_col='cluster_name',
                                 mSamp_col='sample',mCT_col='ATAC_cluster_name',
                                 cellCutoff=200,
                                 aLab='Class proportion\nin scATAC samples',
                                 mLab='Class proportion\nin CITE samples',tSize=8,tVec=NA,clustColors=NA,cs_order=NA){
    
    if(!all(c(cMrna_col,cAtac_col,cID_col) %in% colnames(ID_conv_df))) stop('ID_conv_df does not have all 3 required columns')
    if(!all(c(aSamp_col,aCT_col) %in% colnames(scATAC_meta))) stop('scATAC_meta does not have all required columns')
    if(!all(c(mSamp_col,mCT_col) %in% colnames(CITE_meta))) stop('CITE_meta does not have all required columns')
    if(!any(ID_conv_df[,cAtac_col] %in% unique(scATAC_meta[,aSamp_col]))) stop('atac mismatch')
    if(!any(ID_conv_df[,cMrna_col] %in% unique(CITE_meta[,mSamp_col]))) stop('mrna mismatch')
    
    if(!all(unique(scATAC_meta[,aCT_col]) %in% names(tVec))){
        tVec <- sort(unique(scATAC_meta[,aCT_col]))
        names(tVec) <- sort(unique(scATAC_meta[,aCT_col]))
    }
    
    ll <- table(scATAC_meta[,aSamp_col])
    atac_ID_tooSmall <- c(names(ll[ll<cellCutoff]),
                          ID_conv_df[which(!(ID_conv_df[,cAtac_col] %in% unique(scATAC_meta[,aSamp_col]))),cAtac_col])
    
    ll <- table(CITE_meta[,mSamp_col])
    mrna_ID_tooSmall <- c(names(ll[ll<cellCutoff]),
                          ID_conv_df[which(!(ID_conv_df[,cMrna_col] %in% unique(CITE_meta[,mSamp_col]))),cMrna_col])
    
    
    ID_conv_df <- ID_conv_df[which(!(ID_conv_df[,cAtac_col] %in% atac_ID_tooSmall) & 
                                   !(ID_conv_df[,cMrna_col] %in% mrna_ID_tooSmall)),]
    cat(paste('Using',nrow(ID_conv_df),'donors.\n'))
    
    CITE_CT_sample_freq_df <- as.data.frame(table(CITE_meta[which(CITE_meta[,mSamp_col] %in% ID_conv_df[,cMrna_col]),
                                                                  c(mCT_col,mSamp_col)]),stringsAsFactors=FALSE)
    colnames(CITE_CT_sample_freq_df) <- c('scsnATAC','sample','Freq')
    CITE_CT_sample_freq_df[,cID_col] <- mapvalues(CITE_CT_sample_freq_df[,mSamp_col],
                                                   from=ID_conv_df[,cMrna_col],to=ID_conv_df[,cID_col])
    CITE_CT_sample_freq_df <- CITE_CT_sample_freq_df[,c('scsnATAC','sample',cID_col,'Freq')]
    CITE_sample_freq_df <- as.data.frame(table(CITE_meta[which(CITE_meta[,mSamp_col] %in% ID_conv_df[,cMrna_col]),
                                                               mSamp_col]),stringsAsFactors=FALSE)
    CITE_CT_sample_freq_df$sample_tot <- as.numeric(mapvalues(CITE_CT_sample_freq_df[,mSamp_col],
                                                               from=CITE_sample_freq_df$Var1,to=CITE_sample_freq_df$Freq))
    CITE_CT_sample_freq_df$perc <- CITE_CT_sample_freq_df$Freq/CITE_CT_sample_freq_df$sample_tot
    rownames(CITE_CT_sample_freq_df) <- paste(sep='_',str_split_fixed(CITE_CT_sample_freq_df$scsnATAC,":",2)[,1],
                                              CITE_CT_sample_freq_df[,cID_col])
    if(!all(for(idx in unique(ID_conv_df[,cMrna_col])){
        sum(CITE_CT_sample_freq_df[which(CITE_CT_sample_freq_df[,mSamp_col]==idx),'perc'])!=1})) stop('not all CITE fractions sum to 1')
    

    scATAC_CT_sample_freq_df <- as.data.frame(table(scATAC_meta[which(scATAC_meta[,aSamp_col] %in% ID_conv_df[,cAtac_col]),
                                                            c(aCT_col,aSamp_col)]),stringsAsFactors=FALSE)
    colnames(scATAC_CT_sample_freq_df) <- c('scsnATAC','sample','Freq')
    scATAC_CT_sample_freq_df[,cID_col] <- mapvalues(scATAC_CT_sample_freq_df[,aSamp_col],
                                                     from=ID_conv_df[,cAtac_col],to=ID_conv_df[,cID_col])
    scATAC_CT_sample_freq_df <- scATAC_CT_sample_freq_df[,c('scsnATAC','sample',cID_col,'Freq')]
    scATAC_sample_freq_df <- as.data.frame(table(scATAC_meta[which(scATAC_meta[,aSamp_col] %in% ID_conv_df[,cAtac_col]),
                                                         aSamp_col]),stringsAsFactors=FALSE)
    scATAC_CT_sample_freq_df$sample_tot <- as.numeric(mapvalues(scATAC_CT_sample_freq_df[,aSamp_col],
                                                               from=scATAC_sample_freq_df$Var1,to=scATAC_sample_freq_df$Freq))
    scATAC_CT_sample_freq_df$perc <- scATAC_CT_sample_freq_df$Freq/scATAC_CT_sample_freq_df$sample_tot
    rownames(scATAC_CT_sample_freq_df) <- paste(sep='_',str_split_fixed(scATAC_CT_sample_freq_df$scsnATAC,":",2)[,1],
                                                scATAC_CT_sample_freq_df[,cID_col])
    if(!all(for(idx in unique(ID_conv_df[,cAtac_col])){
        sum(scATAC_CT_sample_freq_df[which(scATAC_CT_sample_freq_df[,aSamp_col]==idx),'perc'])==1})) stop('not all ATAC fractions sum to 1')
    
    if(!identical(sort(rownames(CITE_CT_sample_freq_df)),sort(rownames(scATAC_CT_sample_freq_df)))) stop('CITE/ATAC mismatched names')
    
    comp_scsnATAC_sID_df <- scATAC_CT_sample_freq_df[,c('scsnATAC',cID_col,'perc')]
    colnames(comp_scsnATAC_sID_df) <- c('scsnATAC','subject_id','scATAC_percent')
    comp_scsnATAC_sID_df$CITE_perc <- CITE_CT_sample_freq_df[rownames(comp_scsnATAC_sID_df),'perc']
    colnames(comp_scsnATAC_sID_df) <- c('scsnATAC','subject_id','scATAC_percent','CITE_percent')
    head(comp_scsnATAC_sID_df)
    
    if(!all(unique(comp_scsnATAC_sID_df$scsnATAC) %in% cs_order)) cs_order <- sort(unique(comp_scsnATAC_sID_df$scsnATAC))
    
    myplots <- list()
    i=1
    for(cs in cs_order){
        toPlot <- comp_scsnATAC_sID_df[which(comp_scsnATAC_sID_df$scsnATAC==cs),]

        ll <- cor.test(toPlot$scATAC_percent,toPlot$CITE_percent,method='pearson')

        expand_val <- round(max(toPlot$scATAC_percent,toPlot$CITE_percent),2)

        g <- ggplot(toPlot,aes_string(x='scATAC_percent',y='CITE_percent',color='scsnATAC')) + 
                geom_point(size=3) + theme_bw(base_size=25) + 
                geom_abline(slope=1,intercept=0,linetype='dashed') + 
                expand_limits(x=0) + expand_limits(y=0) + expand_limits(x=expand_val) + expand_limits(y=expand_val) +
                labs(x=aLab,y=mLab,title=paste(sep='',tVec[cs],'\nR=',round(ll$estimate,2),'\nP=',signif(ll$p.value,3))) +
                theme(plot.title = element_text(hjust = 0.5,size=tSize)) + theme(legend.position="none")
        if(all(unique(toPlot$scsnATAC) %in% names(clustColors))){
            g <- g + scale_color_manual(values=clustColors)
        }
        myplots[[i]] <- g
        i=i+1

    }
    
    p <- arrangeGrob(grobs = myplots, ncol=length(unique(comp_scsnATAC_sID_df$scsnATAC)))
    
    colnames(comp_scsnATAC_sID_df) <- c('Class','Subject ID','scATAC sample percent','CITE sample percent')
    return(list('plot'=p,'data'=comp_scsnATAC_sID_df))


}


##CNA FUNCTIONS

#CNA ADD COLUMN
CNA_add_col <- function(CITE_meta, cna_res, cna_col, ncorr_col='res_ncorrs'){

    if(!(ncorr_col %in% colnames(cna_res))) stop('cna_res missing column')
    if(!all(rownames(cna_res) %in% rownames(CITE_meta))) stop('not all cna_res cells in CITE_meta rownames')

    #merge & subset
    CITE_meta[,cna_col] <- NA
    CITE_meta[rownames(cna_res),cna_col] <- cna_res[rownames(cna_res),ncorr_col]
    CITE_meta <- CITE_meta[!is.na(CITE_meta[,cna_col]),] 

    #verify
    #won't be exactly the same because of the state=NA filter above.
    if(!all.equal(summary(CITE_meta[,cna_col]),summary(cna_res[,ncorr_col]),tolerance=1e-03)) stop('cna transfer issue')
    
    return(CITE_meta)
}

#CNA VIOLIN PLOTS - group CNA correlations by class/state
CNA_violin_plots <- function(toPlot,xCol,yCol,fdr_thresh,toColor=NA,reord=TRUE,
                             thisTitle='',xLab='',yLab='Neighborhood\ncorrelation',clustColors=NA){
    
    if(any(is.na(toColor))){
        myPalette <- colorRampPalette(c(brewer.pal(11, "RdBu")[11:7], "#DCDCDC", 
                                        brewer.pal(11, "RdBu")[5:1]))
        toColor <- myPalette(100)
    }
        
    #plot classes/states by decreasing median values
    if(reord){
        ll <- aggregate(toPlot[,yCol] ~ toPlot[,xCol], FUN=median)
        thisOrder <- ll[order(ll[,2],decreasing=TRUE),1]
        toPlot[,xCol] <- factor(toPlot[,xCol],levels=rev(thisOrder))
    }
    
    g <- ggplot(toPlot,aes_string(x = xCol, y = yCol)) + coord_flip() +
      geom_quasirandom_rast(aes_string(color=yCol), width = 0.3, size = 0.001) +
      stat_summary(aes_string(group = xCol), fun = median, fun.min = median, fun.max = median, 
                       geom = "crossbar", color = "grey10", width = 0.5, lwd = 0.2, linetype = 2) +
      scale_colour_gradientn(na.value="#DCDCDC",colours = toColor,
                             limits = c(-max(abs(toPlot[,yCol]),na.rm=TRUE), max(abs(toPlot[,yCol]),na.rm=TRUE))) +
      geom_hline(yintercept = fdr_thresh, linetype = "longdash", color = "black") + 
      geom_hline(yintercept = -fdr_thresh, linetype = "longdash", color = "black") +
      geom_hline(yintercept = 0, linetype = "longdash", color = "darkgrey") +
      labs(title = thisTitle, x= xLab, y = yLab) + 
      theme_classic(base_size = 25) +
      theme(plot.title = element_text(hjust = 0.5, color="black", size = 20),
            legend.position = "none",panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
            axis.text.y = element_text(color = "black"),axis.title.x=element_text(size=22))
    
    if(all(levels(toPlot[,xCol]) %in% names(clustColors))){
        suppressWarnings(g <- g + theme(axis.text.y = element_text(color=clustColors[levels(toPlot[,xCol])],
                                                                   face='bold',size=18)))
    }
    
    return(g)

}

#CNA UMAP PLOTS - show CNA correlations by cell on UMAP
CNA_umap_plots <- function(toPlot,xCol,yCol,corrCol,thisTitle='',xLab='UMAP1',yLab='UMAP2',colorLab='r',smallPt=FALSE,
                           fdr_thresh=NA,fdr_color='#DCDCDC',raster_dpi=300){

    myPalette <- colorRampPalette(c(brewer.pal(11, "RdBu")[11:7], "#DCDCDC", 
                                    brewer.pal(11, "RdBu")[5:1]))

    #Any cells that do not pass the global FDR p-value should not be shown
    if(!is.na(fdr_thresh)) {toPlot[which(abs(toPlot[,corrCol])<fdr_thresh),corrCol] <- NA}
    
    g <- ggplot(toPlot,aes_string(x=xCol,y=yCol,color=corrCol)) 
    if(smallPt) {g <- g + rasterise(geom_point(size=0.001,alpha=0.5),dpi=raster_dpi)} else {g <- g + rasterise(geom_point(size=1,alpha=0.5),dpi=raster_dpi)}
    g <- g + theme_bw(base_size=20) + labs(color=colorLab,title=thisTitle,x=xLab,y=yLab) + 
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid = element_blank()) +
            scale_colour_gradientn(na.value="#DCDCDC",colours = myPalette(100), 
                                   limits = c(-max(abs(toPlot[,corrCol]),na.rm=TRUE),
                                              max(abs(toPlot[,corrCol]),na.rm=TRUE)))
    
    
    toSave <- toPlot[,c(xCol,yCol,corrCol)]
    colnames(toSave) <- c(gsub('\n',' ',xLab),gsub('\n',' ',yLab),
                          str_split_fixed(thisTitle,'\n',2)[,1])
    return(list('plot'=g,'data'=toSave))
}


##ARCHR PLOTTING FUNCTIONS

getTFfromMotif <- function(x){
    split1 <- str_split_fixed(x,'_',2)[,1]
    if(grepl('\\.\\.',split1)){
        mot_list <- unlist(as.list(str_split(split1,'\\.\\.')))
    } else {
        mot_list <- c(split1)
    }
    
    if(grepl('NKX',paste(collapse=',',mot_list))){
        mot_list <- sub('\\.','-',mot_list)
    }
    if(grepl('\\.',paste(collapse=',',mot_list))){
        mot_list <- str_split_fixed(mot_list,'\\.',2)[,1]
    }
    
    return(mot_list)
}

#get gene expression from JASPAR2020 set
getGEfromJ20 <- function(x,gene_vec){
    
    mot_list <- getTFfromMotif(x)
    
    if(!all(mot_list %in% names(gene_vec))) stop(paste('missing gene:',paste(mot_list,collapse=', ')))
        
    return(min(gene_vec[mot_list]))
    
}

pickTopMotifs <- function(this_mxCT, this_gxCT, minE=5, num_mot=7, minGE=0.05, withinE=0.95){
    
    if(!identical(sort(colnames(this_mxCT)),sort(colnames(this_gxCT)))) stop('cell states (colnames) different!')
    
    #added to allow for similar clusters to have similar top motifs (within reason)
    original_columns <- colnames(this_mxCT)
    this_mxCT$maxE2 <- apply(this_mxCT,1,function(x){sort(x,decreasing=TRUE)[2]})
    this_mxCT$maxE <- apply(this_mxCT,1,max)
    this_mxCT$maxE2_cutoff <- this_mxCT$maxE*withinE
    this_mxCT$withinMaxE <- this_mxCT$maxE2>=this_mxCT$maxE2_cutoff
    
    these_all_mot <- c()
    this_chosen_mot_df <- data.frame('cluster_abbr'=character(),'motif'=character(),'int_order'=numeric(),
                                     'maxE'=numeric(),'minGE'=numeric(),stringsAsFactors=FALSE)
    for(cs in names(sort(apply(this_mxCT[,original_columns],2,max),decreasing=TRUE))){
        ll <- this_mxCT[which(!(rownames(this_mxCT) %in% these_all_mot)),]
        ll$minGE <- unname(sapply(rownames(ll),FUN=getGEfromJ20,gene_vec=this_gxCT[,cs]))
        ll <- ll[which(ll$maxE>=minE & ll$minGE>=minGE & (ll$maxE==ll[,cs] | (ll$maxE2==ll[,cs] & ll$withinMaxE==TRUE))),]
        ll <- head(ll[order(ll[,cs],decreasing=TRUE),],n=num_mot)
        if(nrow(ll)!=0){
            these_all_mot <- c(these_all_mot,rownames(ll))
            this_chosen_mot_df <- rbind(this_chosen_mot_df,data.frame('cluster_abbr'=rep(cs,n=nrow(ll)),'motif'=rownames(ll),
                                                                      'int_order'=1:nrow(ll),'maxE'=round(ll$maxE),
                                                                      'minGE'=ll$minGE,stringsAsFactors=FALSE))
        }
    }
    this_chosen_mot_df$motif_colnames <- paste(sep='',this_chosen_mot_df$motif,' (',this_chosen_mot_df$maxE,')')
    
    return(this_chosen_mot_df)
}

plotMatMotif <- function(toPlot,
                         xOrd=NULL,yOrd=NULL,xColors=NULL,yColors=NULL,
                         fillGrad=colorRampPalette(c("1"="#E6E7E8","2"="#3A97FF","3"="#8816A7","4"="black"))(100),
                         xLab='',yLab='',fLab='Norm. Enrich.\n-log10(Padj)\n[0-Max]'){
    
    heat_df <- as.data.frame(toPlot)
    heat_df$yy <- rownames(heat_df)
    heat_df <- gather(heat_df,'xx','ff',all_of(colnames(toPlot)))
    if(is.null(yOrd)) yOrd <- sort(unique(heat_df$yy))
    heat_df$yy <- factor(heat_df$yy,levels=rev(yOrd))
    if(is.null(xOrd)) xOrd <- sort(unique(heat_df$xx))
    heat_df$xx <- factor(heat_df$xx,levels=xOrd)
    
    g <- ggplot(heat_df,aes_string(x='xx',y='yy',fill='ff')) + geom_tile() + 
            theme_classic(base_size=18) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
            scale_fill_gradientn(colors = fillGrad) + 
            labs(x=xLab,y=yLab,fill=fLab)
        
    #fix xColors to include the (max val)
    if(!is.null(xColors) & !any(levels(heat_df$xx) %in% names(xColors))){
        names(xColors) <- lapply(names(xColors),FUN=function(x){levels(heat_df$xx)[charmatch(x,levels(heat_df$xx))]})
        xColors <- xColors[!is.na(names(xColors))]
    }
    
    if(all(levels(heat_df$xx) %in% names(xColors))){
        g <- g + theme(axis.text.x = element_text(color=xColors[levels(heat_df$xx)]))
    }
    if(all(levels(heat_df$yy) %in% names(yColors))){
        suppressWarnings(g <- g + theme(axis.text.y = element_text(color=yColors[levels(heat_df$yy)],face='bold',size=20)))
    }

    return(g)
    
}

ArchR_topMotifs_KWspin <- function(this_mxCT, this_gxCT, minE=5, num_mot=7, minGE=0.05,withinE=0.95,
                                   mOrd=NULL,cOrd=NULL,mColors=NULL,cColors=NULL,
                                   fillGrad=colorRampPalette(c("1"="#E6E7E8","2"="#3A97FF","3"="#8816A7","4"="black"))(100),
                                   mLab='',cLab='',eLab='Normalized\nEnrichment\n-log10(Padj)\n[0-Max]'){
    
    this_mot_df <- pickTopMotifs(this_mxCT, this_gxCT, minE=minE, num_mot=num_mot, minGE=minGE, withinE=withinE)
    
    this_mot_df$cluster_abbr <- factor(this_mot_df$cluster_abbr,levels=cOrd)
    this_mot_df <- this_mot_df[order(this_mot_df$cluster_abbr),]
    
    if(!is.null(mOrd) & all(mOrd %in% rownames(this_mxCT))){
        this_heatMat <- t(this_mxCT[mOrd,cOrd])
        
        this_xOrd <- mOrd
    } else {
        this_heatMat <- t(this_mxCT[this_mot_df$motif,cOrd])
        colnames(this_heatMat) <- paste(sep='',colnames(this_heatMat),' (',this_mot_df$maxE,')')
        
        this_xOrd <- this_mot_df$motif_colnames
    }
    this_heatMat <- apply(this_heatMat,2,function(x){x <- x/max(x)*100}) #normalize to max
    
    g <- plotMatMotif(this_heatMat,xOrd=this_xOrd,yOrd=cOrd,xColors=mColors,yColors=cColors,
                      fillGrad=fillGrad,xLab=mLab,yLab=cLab,fLab=eLab)
    
    #for top TFs df
    ll <- unlist(sapply(this_mot_df$motif,FUN=getTFfromMotif))
    names(ll) <- unname(sapply(names(ll),function(x){ifelse(x %in% this_mot_df$motif,x,substr(x,1,nchar(x)-1))}))

    this_TFs_df <- data.frame('cluster_abbr'=mapvalues(names(ll),
                                                       from=this_mot_df$motif,to=as.character(this_mot_df$cluster_abbr)),
                              'motif'=names(ll),'tf'=unname(ll),stringsAsFactors=FALSE)
    this_TFs_df$class_tf <- paste(sep='_',this_TFs_df$cluster_abbr,this_TFs_df$tf)
    this_TFs_df$motif <- factor(this_TFs_df$motif,levels=unique(this_mot_df$motif))
    this_TFs_df$cluster_abbr <- factor(this_TFs_df$cluster_abbr,levels=cOrd)
    
    return(list('motE'=g,'df'=this_mot_df,'TFdf'=this_TFs_df,'mat'=this_heatMat))
}

TF_exp_wilcox_cells_byClass <- function(aMeta,this_gxc,this_TFs_df,yCol,yLab,cOrd=NULL,xLab='Transcription Factor Gene',
                                        fLab='-log10\nWilcoxon\nFDR',padj_method='BH',pThresh=0.05,
                                        clustColors=NA,motColors=NULL){
    
    tf_exp_wilcox_df <- data.frame('tf'=character(),'class'=character(),'meanIn'=numeric(),'meanOut'=numeric(),
                                   'wilcoxW'=numeric(),'wilcoxP'=numeric(),stringsAsFactors=FALSE)
    for(cl in unique(aMeta[,yCol])){
        if(cl=='scATAC') next
        cl_cells <- rownames(aMeta[which(aMeta$assay=='snATAC' & aMeta[,yCol]==cl),])
        non_cl_cells <- rownames(aMeta[which(aMeta$assay=='snATAC' & aMeta[,yCol]!=cl),])
        for(tf in this_TFs_df[which(this_TFs_df[,yCol]==cl),'tf']){
            ll <- wilcox.test(this_gxc[tf,cl_cells],this_gxc[tf,non_cl_cells],alternative='greater')
            tf_exp_wilcox_df <- rbind(tf_exp_wilcox_df,data.frame('tf'=tf,'class'=cl,
                                                                  'meanIn'=mean(this_gxc[tf,cl_cells]),
                                                                  'meanOut'=mean(this_gxc[tf,non_cl_cells]),
                                                                  'wilcoxW'=unname(ll$statistic),'wilcoxP'=ll$p.value,
                                                                  stringsAsFactors=FALSE))
        }
    }
    tf_exp_wilcox_df$wilcoxPadj <- p.adjust(tf_exp_wilcox_df$wilcoxP,method = padj_method)
    tf_exp_wilcox_df$log10_wilcoxPadj <- -log10(tf_exp_wilcox_df$wilcoxPadj)
    tf_exp_wilcox_df$wilcoxPadj_binary <- tf_exp_wilcox_df$wilcoxPadj<pThresh
    tf_exp_wilcox_df$log10_wilcoxPadj_wNA <- tf_exp_wilcox_df$log10_wilcoxPadj
    tf_exp_wilcox_df[which(tf_exp_wilcox_df$wilcoxPadj>=pThresh),'log10_wilcoxPadj_wNA'] <- NA
    
    top_TFs <- unique(this_TFs_df[order(this_TFs_df$motif),'tf'])
    tf_exp_wilcox_df$tf_maxVal <- NA
    for(tf in top_TFs){ 
        tf_exp_wilcox_df[which(tf_exp_wilcox_df$tf==tf),'tf_maxVal'] <- max(tf_exp_wilcox_df[which(tf_exp_wilcox_df$tf==tf),
                                                                                             'log10_wilcoxPadj'])
    }
    tf_exp_wilcox_df$tf_scaleVal <- (tf_exp_wilcox_df$log10_wilcoxPadj/tf_exp_wilcox_df$tf_maxVal)*100
    tf_exp_wilcox_df$class <- factor(tf_exp_wilcox_df$class,levels=rev(cOrd))
    tf_exp_wilcox_df$tf <- factor(tf_exp_wilcox_df$tf,levels=top_TFs)
    
    yCol <- 'class'
    if(!is.null(cOrd)) tf_exp_wilcox_df[,yCol] <- factor(tf_exp_wilcox_df[,yCol],levels=cOrd)
    tf_exp_wilcox_df <- tf_exp_wilcox_df[order(tf_exp_wilcox_df[,yCol],tf_exp_wilcox_df$tf),
                           c(yCol,'tf','meanIn','meanOut','wilcoxW','wilcoxP','wilcoxPadj','log10_wilcoxPadj')]
    
    return(tf_exp_wilcox_df)

}


##PER CELL SCORES

perCell_score <- function(mat,ps1,ps2){
    if(!all(c(ps1,ps2) %in% rownames(mat))) stop('all peaks in sets 1 and 2 must be in matrix rownames')
    
    mat_subset <- rbind(mat[ps1,]/length(ps1),mat[ps2,]*-1/length(ps2))
    
    mat_score_byCell <- colSums(mat_subset)
    
    return(mat_score_byCell)
}


##SUPERSTATE EXPERIMENT FUNCTIONS

unique_p2g_function <- function(this_p2g,peakRank_vec){
    if(!all(c('peak','gene') %in% colnames(this_p2g))) stop('need peak and gene columns')
    
    ll <- as.data.frame(table(this_p2g$gene),stringsAsFactors=FALSE)
    ll <- ll[which(ll$Freq!=1),]
    this_p2g_all1 <- this_p2g[which(!(this_p2g$gene %in% ll$Var1)),]
    rr <- unlist(lapply(ll$Var1,function(x){xx <- peakRank_vec[this_p2g[which(this_p2g$gene==x),'peak']]; return(names(xx[which.max(xx)]))}))
    this_p2g_non1 <- data.frame('peak'=rr,'gene'=ll$Var1,stringsAsFactors=FALSE)
    toRet <- rbind(this_p2g_all1,this_p2g_non1)
    if(!identical(sort(unique(this_p2g$gene)),sort(toRet$gene))) stop('not the same genes')
    return(toRet)
}

mrna_atac_differential_function <- function(dP_df,dG_df,p2g,label_peakBound=1e6, label_geneBound=1e6,
                                            plot_full=TRUE,use_BH=TRUE,color_vec=NA,this_title=NA){
    if(!all(c('peak','gene','cellType','log2FC','pval') %in% colnames(dP_df))) stop('missing colnames in diffPeaks')
    if(!all(c('feature','group','logFC','pval') %in% colnames(dG_df))) stop('missing colnames in diffGenes')
    if(!all(c('peak','gene') %in% colnames(p2g))) stop('missing colnames in p2g')
    
    #directionality does not change significance, but it does change foldchanges.
    #Only want one statistic for each feature going forward
    if(any(dP_df$log2FC==0)) stop('diffPeaks has log2FC of 0')
    dP_df_pos <- dP_df[which(dP_df$log2FC>0),]
    rownames(dP_df_pos) <- dP_df_pos$peak
    dG_group <- unique(dG_df$group)[1]
    dG_df_pos <- dG_df[which((dG_df$logFC>0) | (dG_df$logFC==0 & dG_df$group==dG_group)),]
    rownames(dG_df_pos) <- dG_df_pos$feature
    
    toPlot <- data.frame('peak'=p2g$peak,'gene'=p2g$gene,
                         'peakCT'=dP_df_pos[p2g$peak,'cellType'],
                         'geneCT'=dG_df_pos[p2g$gene,'group'],
                         'peakLog2FC'=dP_df_pos[p2g$peak,'log2FC'],
                         'geneLogFC'=dG_df_pos[p2g$gene,'logFC'],
                         'peakPval'=dP_df_pos[p2g$peak,'pval'],
                         'genePval'=dG_df_pos[p2g$gene,'pval'],
                         stringsAsFactors=FALSE)
    #everything is the same if gene logFC==0, so help the cell types agree
    toPlot[which(toPlot$geneLogFC==0),'geneCT'] <- toPlot[which(toPlot$geneLogFC==0),'peakCT']
    toPlot$sameDir <- toPlot$peakCT==toPlot$geneCT
    
    if(use_BH){
        toPlot$genePadj <- p.adjust(toPlot$genePval,method='BH')
        toPlot$peakPadj <- p.adjust(toPlot$peakPval,method='BH')
        toPlot$log10genePadj <- -log10(toPlot$genePadj)
        toPlot$log10peakPadj <- -log10(toPlot$peakPadj)
        toPlot$log10genePadjPlot <- fix_infinite(toPlot$log10genePadj)
        toPlot$log10peakPadjPlot <- fix_infinite(toPlot$log10peakPadj)
        
        gene_plotCol <- 'log10genePadjPlot'
        peak_plotCol <- 'log10peakPadjPlot'
        line_intercept <- -log10(0.1)
        gene_lab <- 'Differential Gene -log10(FDR)'
        peak_lab <- 'Differential Promoter Peak\n-log10(FDR)'
    } else {
        toPlot$log10genePval <- -log10(toPlot$genePval)
        toPlot$log10peakPval <- -log10(toPlot$peakPval)
        toPlot$log10genePvalplot <- fix_infinite(toPlot$log10genePval)
        toPlot$log10peakPvalplot <- fix_infinite(toPlot$log10peakPval)
        
        gene_plotCol <- 'log10genePvalplot'
        peak_plotCol <- 'log10peakPvalplot'
        line_intercept <- -log10(0.05)
        gene_lab <- 'Differential Gene -log10(P-value)'
        peak_lab <- 'Differential Promoter Peak\n-log10(P-value)'
    }

    bound <- ceiling(max(toPlot[,gene_plotCol],toPlot[,peak_plotCol])*1.05)

    if(plot_full){
        g <- ggplot(toPlot,aes_string(x=gene_plotCol,y=peak_plotCol,color='geneCT',shape='sameDir',label='gene'))
        toSave <- toPlot[,c('gene','geneCT','sameDir',gene_plotCol,peak_plotCol)]
    } else {
        g <- ggplot(toPlot[which(toPlot[,gene_plotCol]>=line_intercept | toPlot[,peak_plotCol]>=line_intercept),],
                    aes_string(x=gene_plotCol,y=peak_plotCol,color='geneCT',shape='sameDir',label='gene'))
        toSave <- toPlot[which(toPlot[,gene_plotCol]>=line_intercept | toPlot[,peak_plotCol]>=line_intercept),c('gene','geneCT','sameDir',gene_plotCol,peak_plotCol)]
    }
    
    g <- g +
            geom_point(size=3,alpha=0.5) + theme_bw(base_size=25) +  
            labs(x=gene_lab,y=peak_lab,color='Greater Cell State\nIn Gene Exp',shape='Greater Cell State\nPeak Agreement') + 
            geom_vline(xintercept=line_intercept,linetype='dashed') + 
            geom_hline(yintercept=line_intercept,linetype='dashed') +
            expand_limits(x=0) + expand_limits(y=0) + expand_limits(x=bound) + expand_limits(y=bound) +
            guides(color=guide_legend(order=1),shape=guide_legend(order=2)) + 
            geom_text_repel(data=toPlot[which(toPlot[,gene_plotCol]>label_geneBound | 
                                              toPlot[,peak_plotCol]>label_peakBound),],
                            seed=0,size=6,box.padding=1,ylim=c(0.25,20))
    if(all(toPlot$label %in% names(color_vec))){ g <- g + scale_color_manual(values=color_vec)}
    if(!is.na(this_title)) {g <- g + ggtitle(this_title) + theme(plot.title = element_text(hjust = 0.5))}
    
    colnames(toSave) <- c('Gene','Greater Cell State In Gene Exp',
                          'Greater Cell State Peak Agreement',
                          gsub('\n',' ',gene_lab),gsub('\n',' ',peak_lab))
    return(list('plot'=g,'data'=toSave))
    
}

