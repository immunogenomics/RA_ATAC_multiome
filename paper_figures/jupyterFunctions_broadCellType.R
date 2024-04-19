suppressMessages({
    library(ggplot2)
    library(Matrix)
    library(tidyr)
})

QC_steps_barplot <- function(cellCount_df,this_meta,CT_col,CT_colors,CT_order,
                             xLab='Cell Count',yLab='QC step',fLab='Broad\nCell Type',tLab='',
                             allCell_name='allCells',allCell_color='grey'){
    if(!(CT_col %in% colnames(this_meta))) stop('CT_col must be in this_meta')
    if(!all(unique(this_meta[,CT_col]) %in% names(CT_colors))) stop('all CT_col values must be in CT_colors names')
    if(!all(unique(this_meta[,CT_col]) %in% CT_order)) stop('all CT_col values must be in CT_order')
    
    toPlot <- as.data.frame(colSums(cellCount_df))
    colnames(toPlot) <- 'cellCt'
    toPlot$QCstep <- factor(rownames(toPlot),levels=rev(rownames(toPlot)))
    
    if(toPlot[nrow(toPlot),'cellCt']!=nrow(this_meta)) stop('cell counts do not match between last step and this_meta')

    lastStep <- rownames(toPlot[nrow(toPlot),])
    rownames(toPlot) <- NULL
    toPlot$cellType <- allCell_name
    toPlot <- toPlot[which(toPlot$QCstep!=lastStep),]
    ll <- as.data.frame(table(this_meta[,CT_col]),stringsAsFactors=FALSE)
    colnames(ll) <- c('cellType','cellCt')
    ll$QCstep <- lastStep
    ll <- ll[,c('cellCt','QCstep','cellType')]
    toPlot <- rbind(toPlot,ll)
    toPlot$cellType <- factor(toPlot$cellType,levels=c(allCell_name,CT_order))
    
    
    toPlot_colors <- c(unname(CT_colors),allCell_color)
    names(toPlot_colors) <- c(names(CT_colors),allCell_name)
    
    g <- ggplot(toPlot,aes_string(x='cellCt',y='QCstep',fill='cellType')) + geom_bar(stat='identity',position='stack') +
            theme_bw(base_size=25) + labs(x=xLab,y=yLab,title=tLab,fill=fLab) + 
            theme(plot.title = element_text(hjust = 0.5)) + 
            scale_fill_manual(values=toPlot_colors)
    
    toSave <- toPlot[,c('cellType','QCstep','cellCt')]
    colnames(toSave) <- c(gsub('\n',' ',fLab),gsub('\n',' ',yLab),gsub('\n',' ',xLab))
    return(list('plot'=g,'data'=toSave))
}

readCount_byCT <- function(df,CT_col,fun='mean',
                           cols_vec=c('nReads_noDedup','nReads_dedup',
                                      'nReads_peak','nReads_promoter','nReads_chrM','nReads_blacklist')){
    if(!all(c(cols_vec,CT_col) %in% colnames(df))) stop('colname issue')
    
    for(cols in cols_vec){
        if(cols==cols_vec[1]){
            means_df <- aggregate(df[,cols] ~ df[,CT_col], FUN=fun)
            colnames(means_df) <- c('cellType',cols)
        } else {
            ll <- aggregate(df[,cols] ~ df[,CT_col], FUN=fun)
            if(!identical(ll[,1],means_df[,CT_col])) stop('cell type out of order')
            means_df[,cols] <- ll[,2]
        }
    }
    return(means_df)
}

scalePeak_forHeatmap <- function(gOrd,ctOrd,gp_map,pCTnorm,peakCTcutoff=0.05){
    if(!all(gOrd %in% names(gp_map))) stop('all genes from order need to be in gp_map')
    if(!all(unname(gp_map) %in% rownames(pCTnorm))) stop('peaks from gene/peak relationship must be in normalized peak accessibility')
    if(!all(ctOrd %in% colnames(pCTnorm))) stop('cell types from order must be in normalized peak accessibility')
    if(!all(apply(pCTnorm[unname(gp_map),],1,max)>peakCTcutoff)) stop('Minimal accessibility cutoff not met')
    
    gp_map <- gp_map[gOrd]
    
    pCTnorm_subset <- pCTnorm[rev(unname(gp_map)),rev(ctOrd)]
    rownames(pCTnorm_subset) <- rev(names(gp_map))

    pCTnorm_subset_scaled <- t(scale(t(pCTnorm_subset)))
    
    return(pCTnorm_subset_scaled)
}

gather_SNPs <- function(toPlot,SNP_peak_vec,class_broadCT_vec,
                        gather_key='class',gather_val='norm_pxCTmean_scale',gather_cols=NA,rowname_col='rsID'){
    if(!all(rownames(toPlot) %in% names(SNP_peak_vec))) stop('all peaks in toPlot must be in SNP_peak_vec')
    if(!all(unique(substr(colnames(toPlot),1,1)) %in% 
            names(class_broadCT_vec))) stop('all first letters of colnames in toPlot must be in class_broadCT_vec')
    
    if(any(is.na(gather_cols)) | !all(gather_cols %in% colnames(toPlot))) gather_cols <- colnames(toPlot)
    
    toGather_df <- as.data.frame(toPlot,stringsAsFactors=FALSE)
    if(!any(is.na(rowname_col))) toGather_df[,rowname_col] <- rownames(toGather_df)
    toGather_df$peak <- unname(SNP_peak_vec[rownames(toGather_df)])
    gathered <- gather(toGather_df,gather_key,gather_val,all_of(gather_cols))
    colnames(gathered)[ncol(gathered)-1] <- gather_key
    colnames(gathered)[ncol(gathered)] <- gather_val
    
    gathered$cellType <- unname(class_broadCT_vec[substr(gathered[,gather_key],1,1)])
    
    return(gathered)
}

