print_time <- function(msg) {
    cat(paste(sep='',"\n",Sys.time()," - ",msg,"\n"))
}
substrRight <- function(x, n){substr(x, nchar(x)-n+1, nchar(x))}

print_time("Library Loading and Defining Functions")

suppressMessages({
    library(ArchR)
    library(tidyr)
    library(stringr)
    library(ggplot2)
    library(argparse)
})

getGEfromJ20 <- function(x,gene_vec){ #gene vec would be gxCT_norm[cs]
    suppressMessages(library(stringr))
    
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
        
    if(!all(mot_list %in% names(gene_vec))) stop(paste('missing gene:',paste(mot_list,collapse=', ')))
        
    return(min(gene_vec[mot_list]))
    
}

pickTopMotifs <- function(this_mxCT, this_gxCT, minE=5, num_mot=7, minGE=0.05, withinE=0.95){
    
    suppressMessages(library(Matrix))
    
    if(!identical(sort(colnames(this_mxCT)),sort(colnames(this_gxCT)))) stop('cell states (colnames) different!')
    
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
        ll$minGE <- unname(sapply(rownames(ll),FUN=getGEfromJ20,gene_vec=gxCT_norm[,cs]))
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
    
    suppressMessages({
        library(tidyr)
        library(ggplot2)
    })
    
    heat_df <- as.data.frame(toPlot)
    heat_df$yy <- rownames(heat_df)
    heat_df <- gather(heat_df,'xx','ff',all_of(colnames(toPlot)))
    if(is.null(yOrd)) yOrd <- sort(unique(heat_df$yy))
    heat_df$yy <- factor(heat_df$yy,levels=rev(yOrd))
    if(is.null(xOrd)) xOrd <- sort(unique(heat_df$xx))
    heat_df$xx <- factor(heat_df$xx,levels=xOrd)
    
    g <- ggplot(heat_df,aes_string(x='xx',y='yy',fill='ff')) + geom_tile() + 
            theme_classic(base_size=15) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
            scale_fill_gradientn(colors = fillGrad) + 
            labs(x=xLab,y=yLab,fill=fLab)
    if(all(levels(heat_df$xx) %in% names(xColors))){
        g <- g + theme(axis.text.x = element_text(color=xColors[levels(heat_df$xx)]))
    }
    if(all(levels(heat_df$yy) %in% names(yColors))){
        g <- g + theme(axis.text.y = element_text(color=yColors[levels(heat_df$yy)],face='bold',size=15))
    }

    return(g)
    
}

ArchR_topMotifs_KWspin <- function(this_mxCT, this_gxCT, minE=5, num_mot=7, minGE=0.05,withinE=0.95,
                                   mOrd=NULL,cOrd=NULL,mColors=NULL,cColors=NULL,
                                   fillGrad=colorRampPalette(c("1"="#E6E7E8","2"="#3A97FF","3"="#8816A7","4"="black"))(100),
                                   mLab='',cLab='',eLab='Norm. Enrich.\n-log10(Padj)\n[0-Max]'){
    
    suppressMessages({
        library(Matrix)
        library(tidyr)
        library(stringr)
        library(ggplot2)
    })
    
    this_mot_df <- pickTopMotifs(this_mxCT, this_gxCT, minE=minE, num_mot=num_mot, minGE=minGE, withinE=withinE)
    
    this_mot_df$cluster_abbr <- factor(this_mot_df$cluster_abbr,levels=cOrd)
    this_mot_df <- this_mot_df[order(this_mot_df$cluster_abbr),]

    this_heatMat <- t(this_mxCT[this_mot_df$motif,cOrd])
    colnames(this_heatMat) <- paste(sep='',colnames(this_heatMat),' (',this_mot_df$maxE,')')
    this_heatMat <- apply(this_heatMat,2,function(x){x <- x/max(x)*100}) #normalize to max
    
    g <- plotMatMotif(this_heatMat,xOrd=this_mot_df$motif_colnames,yOrd=cOrd,xColors=mColors,yColors=cColors,
                      fillGrad=fillGrad,xLab=mLab,yLab=cLab,fLab=eLab)
    
    return(g)
}



parser <- ArgumentParser(description='remake TF motif heatmap from ArchR project')
parser$add_argument('project_dir', metavar='project_dir', help='project directory')
parser$add_argument('--wk_dir', metavar='wk_dir', help='working directory. this directory must already exist or will default to one directory up from project_dir')
parser$add_argument('extract_dir', metavar='extract_dir', help='extract (aka output) directory')
parser$add_argument('motif_prefix', metavar='motif_prefix', help='motif set prefix')
parser$add_argument('gxCT_norm_file', metavar='gxCT_norm_file', help='gene x cell type normalized file')
parser$add_argument('--topn', metavar='topn', default=7, type='integer', help='mlog10Padj of top enriched motifs per cluster; default=7')
parser$add_argument('--minE', metavar='minE', default=5, type='double', help='mlog10 minimum enrichment; default=5')
parser$add_argument('--minGE', metavar='minGE', default=0.05, type='double', help='minimum gxCT_norm value in cluster; default=0.05')
parser$add_argument('--withinE', metavar='withinE', default=0.95, type='double', help='motif in cluster must be within X percent of max across cluster to be in top motifs - to invalidate, use 1; default 0.95')
parser$add_argument('--color_file', metavar='color_file', help='chromatin class color file - must have at least cluster_abbr to be used!')
parser$add_argument('--CT_order', metavar='CT_order', help='chromatin class order in comma-delimited list string; e.g., BA-1,BA-3,BA-2')
parser$add_argument('--no_ArchR_plots', action='store_true', help='if given, do not remake plotEnrichHeatmap with these topn and minE requirements')
args <- parser$parse_args()

project_dir <- args$project_dir
wk_dir <- args$wk_dir
extract_dir <- args$extract_dir
motif_prefix <- args$motif_prefix
gxCT_norm_file <- args$gxCT_norm_file
topn <- args$topn
minE <- args$minE
minGE <- args$minGE
withinE <- args$withinE
color_file <- args$color_file
CT_order <- args$CT_order
no_ArchR_plots <- args$no_ArchR_plots

print_time("Argument Checking")
if(!all(file.exists(project_dir,gxCT_norm_file))) stop("Input file(s) don't exist.")

if(substrRight(project_dir,1)!='/') project_dir<-paste(sep="",project_dir,'/')

if(substrRight(extract_dir,1)!='/') extract_dir<-paste(sep="",extract_dir,'/')
dir.create(extract_dir,showWarnings=FALSE)

eMot_padj_file <- paste(sep='',extract_dir,basename(project_dir),'_enrich',motif_prefix,'Motifs_df.rds')
if(!file.exists(eMot_padj_file)) stop(paste('ERROR: motif enrichment df file does not exist:',eMot_padj_file,'\n.'))

if(!no_ArchR_plots){
    eMot_file <- paste(sep='',extract_dir,basename(project_dir),'_enrich',motif_prefix,'Motifs.rds')
    if(!file.exists(eMot_file)) stop(paste('ERROR: motif enrichment sumExp file does not exist:',eMot_file,'\n.'))
}

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


print_time('Enriched motifs - my ordering')
eMot_padj <- readRDS(eMot_padj_file)
gxCT_norm <- readRDS(gxCT_norm_file)
if(!identical(sort(colnames(eMot_padj)),sort(colnames(gxCT_norm))) & 
   all(str_detect(colnames(eMot_padj),'^[a-zA-Z]{2}[0-9]+$'))){
    colnames(eMot_padj) <- lapply(colnames(eMot_padj),FUN=function(s){paste(sep='',substr(s,1,2),'-',substr(s,3,nchar(s)))})
}
if(!identical(sort(colnames(eMot_padj)),sort(colnames(gxCT_norm)))) stop('mxCT and gxCT matrices do not have same CT.')

cluster_colors <- NULL
if(length(color_file)!=0){ if(file.exists(color_file)){ cluster_colors <- readRDS(color_file)}}

CTord <- NULL
if(length(CT_order)!=0){
    CT_order <- unlist(str_split(CT_order,','))
    if(all(CT_order %in% colnames(eMot_padj))) CTord <- CT_order
}

g <- ArchR_topMotifs_KWspin(eMot_padj,gxCT_norm,cOrd=CTord,cColors=cluster_colors,minE=minE,num_mot=topn,minGE=minGE,
                            withinE=withinE)
ggsave(file=paste(sep='',extract_dir,basename(project_dir),'_enrich',
                  motif_prefix,'Motifs_KWheatMat_top-',topn,'_withinE-',withinE*100,
                  '_minE-',sub('\\.','-',minE),'_minGE-',sub('\\.','-',minGE),'.png'),
       plot=g,units='in',height=5,width=10)

g <- ArchR_topMotifs_KWspin(eMot_padj,gxCT_norm,cOrd=CTord,cColors=cluster_colors,minE=minE,num_mot=topn,minGE=0,
                            withinE=withinE)
ggsave(file=paste(sep='',extract_dir,basename(project_dir),'_enrich',
                  motif_prefix,'Motifs_KWheatMat_top-',topn,'_withinE-',withinE*100,
                  '_minE-',sub('\\.','-',minE),'_minGE-0.png'),
       plot=g,units='in',height=5,width=10)


print_time('Enriched motifs - ArchR ordering')
if(!no_ArchR_plots){
    eMot <- readRDS(eMot_file)
    heatmapEM_mat <- plotEnrichHeatmap(eMot, n = topn, transpose = TRUE, cutOff = minE, returnMatrix=TRUE) #could lower cutoff
    saveRDS(heatmapEM_mat,paste(sep='',extract_dir,basename(project_dir),'_enrich',
                                motif_prefix,'Motifs_heatMat_top-',topn,'_minE-',sub('\\.','-',minE),'.rds'))

    heatmapEM <- plotEnrichHeatmap(eMot, n = topn, transpose = TRUE, cutOff = minE)
    plotPDF(heatmapEM, 
            name = paste(sep='',basename(project_dir),'_enrich',motif_prefix,'Motifs_heatMat_top-',topn,
                         '_minE-',sub('\\.','-',minE)),
            width = 8, height = 6, ArchRProj = NULL, addDOC = FALSE)
    oldFile <- paste(sep='',wk_dir,'Plots/',basename(project_dir),
                     '_enrich',motif_prefix,'Motifs_heatMat_top-',topn,'_minE-',sub('\\.','-',minE),'.pdf')
    newFile <- paste(sep='',extract_dir,basename(project_dir),
                     '_enrich',motif_prefix,'Motifs_heatMat_top-',topn,'_minE-',sub('\\.','-',minE),'.pdf')
    if(file.exists(oldFile)){
        cat('moving ArchR PDF file. might want to delete empty Plots/ in wk_dir\n')
        file.rename(from=oldFile,to=newFile)
    } else {
        cat(paste(sep='','ArchR filepath wrong: ',oldFile,'\n'))
    }

    g <- plotMatMotif(heatmapEM_mat,yColors=cluster_colors)
    ggsave(file=paste(sep='',extract_dir,basename(project_dir),
                      '_enrich',motif_prefix,'Motifs_heatMat_top-',topn,'_minE-',sub('\\.','-',minE),'.png'),
           plot=g,units='in',height=5,width=10)
} else {
    cat('SKIPPING: plotEnrichHeatmap with these topn and minE requirements.\n')
}


print_time('Done.')

