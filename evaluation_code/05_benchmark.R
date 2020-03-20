rm(list=ls())
source("/Users/jackfu/Documents/00_Postdoc/18-10-gCNV_deploy/00_functions.R")
source("/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/scripts/03_compile_results.R")
# source("/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/scripts/04_inheritance.R")
    
    # ### Select for comparison with XHMM, CoNIFER etc
    # # ssc <- m[!is.na(m$SFARI_ID),]
    # # mat <- match(ssc$SFARI_ID, gcnv_nodup$sample)
    # # ssc$batch <- gcnv_nodup$batch[mat]
    # # sort(table(ssc$batch[str_detect(ssc$batch, "COHORT")]))
    # mem <- read.table('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/data/manifest/mem_11_18_2019.txt', header=TRUE, stringsAsF=F)
    # cohorts <- c('COHORT_4_9', 'COHORT_4_7', 'COHORT_4_10',  'COHORT_4_5',  'COHORT_4_8')
    # mem <- mem[mem[,1] %in% cohorts,]
    #     colnames(mem)[1] <- "membership:sample_set_id"
    # write.table(mem[,1:2], file="/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/data/manifest/mem_comparison.txt", quote=F, row.names=F, sep="\t")

################################################################################################################################################
## WGS calls (overall)
################################################################################################################################################
meta_sscid <- m[!is.na(m$SSC_ID) & m$filter_raw_count==0 & m$filter_dup==0,]
gcnv_bm <- gcnv_nodup[gcnv_nodup$sample %in% meta_sscid$SFARI_ID & gcnv_nodup$filter_dup==FALSE & gcnv_nodup$filter_raw_count==FALSE,]
    mat <- match(gcnv_bm$sample, meta_sscid$SFARI_ID)
    gcnv_bm$sample <- meta_sscid$SSC_ID[mat]

################################################################################################################################################
################################################################################################################################################
samples <- as.character(unique(gcnv_bm$sample))
truth_dup <- processTruth("/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/data/truth/10027.pass.DUP.pass.exome.bed", 
    "/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/data/bins/calling/hg38_gencode_unsplit.bed", 
    samples)

truth_del <- processTruth("/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/data/truth/10027.pass.DEL.pass.exome.bed", 
    "/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/data/bins/calling/hg38_gencode_unsplit.bed", 
    samples)

truth_mcnv <- processTruth('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/data/truth/10027.pass.MCNV.pass.bed',
    "/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/data/bins/calling/hg38_gencode_unsplit.bed", 
    samples)

truth <- rbind(truth_del, truth_dup, truth_mcnv)
    truth[,7] <- unlist(truth[,7])

gr_truth <- GRanges(paste0(unlist(truth[,7]), "-", unlist(truth[,1])), IRanges(unlist(truth[,2]), unlist(truth[,3])), strand=unlist(truth[,5]), sample=unlist(truth[,7]), ac=unlist(truth[,9]), var=unlist(truth[,8]))
    gr_truth$af <- gr_truth$ac/length(unique(truth[,7]))
    
########################################################################
### Calculate performance
########################################################################
gcnv_clean <- gcnv_bm
    gcnv_clean$vaf <- gcnv_clean$sf
    gcnv_clean$vac <- gcnv_clean$sc
    gcnv_clean$CN <- as.numeric(as.character(gcnv_clean$CN))
  
########################################################################  
### First round sensitivity measurement
results_werling_sens <- compileResultsSens(samples, gcnv_clean, gr_truth, QS_dup=50, QS_del=100)
    sens <- do.call(c, results_werling_sens[[2]])
    sens_rare <- sens[sens$af<.01]
    
    ### Exclude outlier samples
    tbl <- sort(table(sens_rare$sample))
        samps_exclude_sens <- names(tbl[tbl>20])
        sens_rare <- sens_rare[sens_rare$sample %in% samps_exclude_sens==F]

    ### Sens Statistics
    info_sens_hit <- by(sens_rare$cov, paste0(as.character(strand(sens_rare)), "::", sens_rare$var), function(x) sum(x>=.3))
    info_sens_target <- by(sens_rare$cov, paste0(as.character(strand(sens_rare)), "::", sens_rare$var), function(x) length(x))
    info_sens_wid <- by(sens_rare$wid, paste0(as.character(strand(sens_rare)), "::", sens_rare$var), min)
    info_sens <- cbind(info_sens_hit, info_sens_target, info_sens_wid)
        info_sens <- data.frame(info_sens)
        info_sens$type <- "DUP"; info_sens$type[str_detect(rownames(info_sens), "^-")] <- "DEL"
        info_sens$success <- info_sens$info_sens_hit/info_sens$info_sens_target>=.75

    sen_mar_del <- by(info_sens$success[info_sens$type=="DEL"], info_sens$info_sens_wid[info_sens$type=="DEL"], mean)
    sen_mar_dup <- by(info_sens$success[info_sens$type=="DUP"], info_sens$info_sens_wid[info_sens$type=="DUP"], mean)
    sen_cum_del <- sapply(1:50, function(x) mean(info_sens$success[info_sens$info_sens_wid>=x & info_sens$type=="DEL"]))
    sen_cum_dup <- sapply(1:50, function(x) mean(info_sens$success[info_sens$info_sens_wid>=x & info_sens$type=="DUP"]))

    window <- 20
    par(mfrow=c(2,2))
    plot(sen_mar_del, xlim=c(1, window), ylim=c(0, 1))    
    plot(sen_mar_dup, xlim=c(1, window), ylim=c(0, 1))    
    plot(sen_cum_del, xlim=c(1, window), ylim=c(0, 1))    
    plot(sen_cum_dup, xlim=c(1, window), ylim=c(0, 1))    

########################################################################  
### First round ppv measurement
gr_truth_red <- reduce(gr_truth[strand(gr_truth)!="*"])
    gr_truth_red$sample <- str_replace(seqnames(gr_truth_red), "-.*", "")
    gcnv_ppv <- gcnv_clean[gcnv_clean$sample %in% gr_truth_red$sample,]

results_werling_ppv <- compileResultsPPV(unique(gcnv_ppv$sample), gcnv_ppv, gr_truth_red, QS_dup=50, QS_del=100)
    ppv <- do.call(c, results_werling_ppv[[2]])
    ppv_rare <- ppv[ppv$vaf<.01]

    info_ppv_hit <- by(ppv_rare$cov, paste0(as.character(strand(ppv_rare)), "::", ppv_rare$var_name), function(x) sum(x>=.3))
    info_ppv_target <- by(ppv_rare$cov, paste0(as.character(strand(ppv_rare)), "::", ppv_rare$var_name), function(x) length(x))
    info_ppv_wid <- by(ppv_rare$wid, paste0(as.character(strand(ppv_rare)), "::", ppv_rare$var_name), min)
    info_ppv <- cbind(info_ppv_hit, info_ppv_target, info_ppv_wid)
        info_ppv <- data.frame(info_ppv)
        info_ppv$type <- "DUP"; info_ppv$type[str_detect(rownames(info_ppv), "^-")] <- "DEL"
        info_ppv$success <- info_ppv$info_ppv_hit/info_ppv$info_ppv_target>=.75

    ppv_mar_del <- by(info_ppv$success[info_ppv$type=="DEL"], info_ppv$info_ppv_wid[info_ppv$type=="DEL"], mean)
    ppv_mar_dup <- by(info_ppv$success[info_ppv$type=="DUP"], info_ppv$info_ppv_wid[info_ppv$type=="DUP"], mean)
    ppv_cum_del <- sapply(1:50, function(x) mean(info_ppv$success[info_ppv$info_ppv_wid>=x & info_ppv$type=="DEL"]))
    ppv_cum_dup <- sapply(1:50, function(x) mean(info_ppv$success[info_ppv$info_ppv_wid>=x & info_ppv$type=="DUP"]))

    window <- 20
    par(mfrow=c(2,2))
    plot(ppv_mar_del, xlim=c(1, window), ylim=c(0, 1))    
    plot(ppv_mar_dup, xlim=c(1, window), ylim=c(0, 1))    
    plot(ppv_cum_del, xlim=c(1, window), ylim=c(0, 1))    
    plot(ppv_cum_dup, xlim=c(1, window), ylim=c(0, 1))    

########################################################################  
### Clean up genome calls
########################################################################  
### Sensitivity
inspect <- sens_rare[sens_rare$cov<.3 & sens_rare$wid>=3]

svtype <- rep("DEL", length(inspect)); svtype[which(strand(inspect)=="+")] <- "DUP"
    bed_geno <- cbind(str_replace(as.character(seqnames(inspect)), ".*-", ""), start(inspect), end(inspect), paste0("CNV_", svtype), as.character(inspect$sample), svtype)
    bed_geno <- bed_geno[order(bed_geno[,1], bed_geno[,2]),]

for(i in 1:length(inspect)){ 
    inspectPlot(i, bed_geno, '/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/plots/', sens=TRUE)
}

########################################################################  
### PPV
inspect <- ppv_rare[ppv_rare$wid>=2 & ppv_rare$cov<.3]
inspect <- ppv_rare[ppv_rare$wid==2 & ppv_rare$cov<.3]

svtype <- rep("DEL", length(inspect)); svtype[which(strand(inspect)=="+")] <- "DUP"
    bed_geno <- cbind(str_replace(as.character(seqnames(inspect)), ".*-", ""), start(inspect), end(inspect), paste0("CNV_", svtype), as.character(inspect$sample), svtype)
    bed_geno <- bed_geno[order(bed_geno[,1], bed_geno[,2]),]

for(i in 1:length(inspect)){ 
    inspectPlot(i, bed_geno, '/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/plots/', ppv=TRUE)
}

########################################################################  
### Performance, re-tuned
########################################################################  
### Second round sensitivity measurement, after re-sizing genome call boundaries
gr_sens <- gr_truth
torm <- str_replace(list.files('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/plots/sens'), ".pdf", "")
    torm <- str_split_fixed(torm, ":|-", n=4)
    mat <- match(paste0(torm[,4], "-", torm[,1], ":", torm[,2], "-", torm[,3]), str_replace(as.character(gr_sens), "((:\\+)$)|((:-)$)", ""))
    gr_sens <- gr_sens[-mat]

    ### Add in re-sized sens calls
    add <- str_replace(list.files('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/plots/sens_resized'), ".pdf", "")
        add <- str_split_fixed(add, ":|-", n=5)
        gr_add <- GRanges(paste0(add[,4], "-", add[,1]), IRanges(as.numeric(add[,2]), as.numeric(add[,3])), strand="-")
            gr_add$sample <- add[,4]
            gr_add$ac <- 0
                add_uniq <- unique(add[,1:3])
                mat_uiq <- match(paste0(start(gr_add), "-", end(gr_add)), paste0(add_uniq[,2], "-", add_uniq[,3]))
            gr_add$var = paste0("resized_", mat_uiq)
            gr_add$af <- 0
            strand(gr_add)[add[,5]=="DUP"] <- "+"
    gr_sens <- c(gr_sens, gr_add)

    ### Add in re-sized ppv calls
    add <- str_replace(list.files('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/plots/ppv_resized'), ".pdf", "")
        add <- str_split_fixed(add, ":|-", n=5)
        gr_add <- GRanges(paste0(add[,4], "-", add[,1]), IRanges(as.numeric(add[,2]), as.numeric(add[,3])), strand="-")
            gr_add$sample <- add[,4]
            gr_add$ac <- -1
                add_uniq <- unique(add[,1:3])
                mat_uiq <- match(paste0(start(gr_add), "-", end(gr_add)), paste0(add_uniq[,2], "-", add_uniq[,3]))
            gr_add$var = paste0("resized_", mat_uiq)
            gr_add$af <- -1
            strand(gr_add)[add[,5]=="DUP"] <- "+"
    gr_sens <- c(gr_sens, gr_add)

### Resized sensitivity
results_werling_sens <- compileResultsSens(samples, gcnv_clean, gr_sens, QS_dup=50, QS_del=100)
    sens <- do.call(c, results_werling_sens[[2]])
    sens_rare <- sens[sens$af<.01]
    ### Exclude outlier samples (previously computed list)
        sens_rare <- sens_rare[sens_rare$sample %in% samps_exclude_sens==F]

    ### Sens Statistics
    info_sens_hit <- by(sens_rare$cov, paste0(as.character(strand(sens_rare)), "::", sens_rare$var), function(x) sum(x>=.3))
    info_sens_target <- by(sens_rare$cov, paste0(as.character(strand(sens_rare)), "::", sens_rare$var), function(x) length(x))
    info_sens_wid <- by(sens_rare$wid, paste0(as.character(strand(sens_rare)), "::", sens_rare$var), min)
    info_sens <- cbind(info_sens_hit, info_sens_target, info_sens_wid)
        info_sens <- data.frame(info_sens)
        info_sens$type <- "DUP"; info_sens$type[str_detect(rownames(info_sens), "^-")] <- "DEL"
        info_sens$success <- info_sens$info_sens_hit/info_sens$info_sens_target>=.75

    sen_mar_del <- by(info_sens$success[info_sens$type=="DEL"], info_sens$info_sens_wid[info_sens$type=="DEL"], mean)
    sen_mar_dup <- by(info_sens$success[info_sens$type=="DUP"], info_sens$info_sens_wid[info_sens$type=="DUP"], mean)
    sen_cum_del <- sapply(1:1000, function(x) mean(info_sens$success[info_sens$info_sens_wid>=x & info_sens$type=="DEL"]))
    sen_cum_dup <- sapply(1:1000, function(x) mean(info_sens$success[info_sens$info_sens_wid>=x & info_sens$type=="DUP"]))

    window <- 20
    par(mfrow=c(2,2))
    plot(sen_mar_del, xlim=c(1, window), ylim=c(0, 1))    
    plot(sen_mar_dup, xlim=c(1, window), ylim=c(0, 1))    
    plot(sen_cum_del, xlim=c(1, window), ylim=c(0, 1))    
    plot(sen_cum_dup, xlim=c(1, window), ylim=c(0, 1))    

### Resized PPV
gr_truth_red <- reduce(gr_sens[strand(gr_sens)!="*"])
    gr_truth_red$sample <- str_replace(seqnames(gr_truth_red), "-.*", "")
    gcnv_ppv <- gcnv_clean[gcnv_clean$sample %in% gr_truth_red$sample,]

results_werling_ppv <- compileResultsPPV(unique(gcnv_ppv$sample), gcnv_ppv, gr_truth_red, QS_dup=50, QS_del=100)
    ppv <- do.call(c, results_werling_ppv[[2]])
    ppv_rare <- ppv[as.vector(ppv$sf)<.01]

    info_ppv_hit <- by(ppv_rare$cov, paste0(as.character(strand(ppv_rare)), "::", ppv_rare$var_name), function(x) sum(x>=.3))
    info_ppv_target <- by(ppv_rare$cov, paste0(as.character(strand(ppv_rare)), "::", ppv_rare$var_name), function(x) length(x))
    info_ppv_wid <- by(ppv_rare$wid, paste0(as.character(strand(ppv_rare)), "::", ppv_rare$var_name), min)
    info_ppv <- cbind(info_ppv_hit, info_ppv_target, info_ppv_wid)
        info_ppv <- data.frame(info_ppv)
        info_ppv$type <- "DUP"; info_ppv$type[str_detect(rownames(info_ppv), "^-")] <- "DEL"
        info_ppv$success <- info_ppv$info_ppv_hit/info_ppv$info_ppv_target>=.75

    ppv_mar_del <- by(info_ppv$success[info_ppv$type=="DEL"], info_ppv$info_ppv_wid[info_ppv$type=="DEL"], mean)
    ppv_mar_dup <- by(info_ppv$success[info_ppv$type=="DUP"], info_ppv$info_ppv_wid[info_ppv$type=="DUP"], mean)
    ppv_cum_del <- sapply(1:200, function(x) mean(info_ppv$success[info_ppv$info_ppv_wid>=x & info_ppv$type=="DEL"]))
    ppv_cum_dup <- sapply(1:200, function(x) mean(info_ppv$success[info_ppv$info_ppv_wid>=x & info_ppv$type=="DUP"]))

    window <- 20
    par(mfrow=c(2,2))
    plot(ppv_mar_del, xlim=c(1, window), ylim=c(0, 1))    
    plot(ppv_mar_dup, xlim=c(1, window), ylim=c(0, 1))    
    plot(ppv_cum_del, xlim=c(1, window), ylim=c(0, 1))    
    plot(ppv_cum_dup, xlim=c(1, window), ylim=c(0, 1))     

    # ### Add in missed 
    # add <- str_replace(list.files('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/plots/add'), ".pdf", "")
    #     add <- str_split_fixed(add, ":|-", n=4)
    #     gr_add <- GRanges(paste0(add[,4], "-", add[,1]), IRanges(as.numeric(add[,2]), as.numeric(add[,3])))
    #         gr_add$sample <- add[,4]
    #         gr_add$ac <- 1
    #         gr_add$var = paste0("add_", 1:length(gr_add))
    #         gr_add$af <- 0.005
    # gr_truth <- c(gr_truth, gr_add)


################################################
### Plot
################################################
    window <- 100
    par(mfrow=c(2,2))
    plot(as.numeric(sen_mar_del)~as.numeric(names(sen_mar_del)), ty="b", pch=19, col=2, xlim=c(0, window), ylim=c(0, 1), main="Marginal Sensitivity vs WGS", ylab="Sensitivity", xlab="# exons")
        lines(as.numeric(sen_mar_dup)~as.numeric(names(sen_mar_dup)), ty="b", pch=19, col=4)
    plot(as.numeric(ppv_mar_del)~as.numeric(names(ppv_mar_del)), ty="b", pch=19, col=2, xlim=c(0, window), ylim=c(0, 1), main="Marginal PPV vs WGS", ylab="PPV", xlab="# exons")
        lines(as.numeric(ppv_mar_dup)~as.numeric(names(ppv_mar_dup)), ty="b", pch=19, col=4)

    plot(y=sen_cum_del, x=1:length(sen_cum_del), ty="b", pch=19, col=2, xlim=c(0, window), ylim=c(0, 1), main="Cumulative Sensitivity vs WGS", ylab="Sensitivity", xlab=">= # exons")
        lines(y=sen_cum_dup, x=1:length(sen_cum_dup), ty="b", pch=19, col=4)
    plot(y=ppv_cum_del, x=1:length(ppv_cum_del), ty="b", pch=19, col=2, xlim=c(0, window), ylim=c(0, 1), main="Cumulative PPV vs WGS", ylab="PPV", xlab=">= # exons")
        lines(y=ppv_cum_dup, x=1:length(ppv_cum_dup), ty="b", pch=19, col=4)

    window <- 20
    par(mfrow=c(2,2))
    plot(as.numeric(sen_mar_del)~as.numeric(names(sen_mar_del)), ty="b", pch=19, col=2, xlim=c(0, window), ylim=c(0, 1), main="Marginal Sensitivity vs WGS", ylab="Sensitivity", xlab="# exons")
        lines(as.numeric(sen_mar_dup)~as.numeric(names(sen_mar_dup)), ty="b", pch=19, col=4)
    plot(as.numeric(ppv_mar_del)~as.numeric(names(ppv_mar_del)), ty="b", pch=19, col=2, xlim=c(0, window), ylim=c(0, 1), main="Marginal PPV vs WGS", ylab="PPV", xlab="# exons")
        lines(as.numeric(ppv_mar_dup)~as.numeric(names(ppv_mar_dup)), ty="b", pch=19, col=4)

    plot(y=sen_cum_del, x=1:length(sen_cum_del), ty="b", pch=19, col=2, xlim=c(0, window), ylim=c(0, 1), main="Cumulative Sensitivity vs WGS", ylab="Sensitivity", xlab=">= # exons")
        lines(y=sen_cum_dup, x=1:length(sen_cum_dup), ty="b", pch=19, col=4)
    plot(y=ppv_cum_del, x=1:length(ppv_cum_del), ty="b", pch=19, col=2, xlim=c(0, window), ylim=c(0, 1), main="Cumulative PPV vs WGS", ylab="PPV", xlab=">= # exons")
        lines(y=ppv_cum_dup, x=1:length(ppv_cum_dup), ty="b", pch=19, col=4)

tmp2 <- tmp2[order(tmp2$var_name)]
inspect <- tmp2
for(i in 19:33){
    inspectPlot2(tmp2[i], '/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/plots/SSC_vis/')
}

# ################################################
# ### Filtering refinement
# ################################################
# sc_clust <- by(gcnv_nodup$var_name, str_replace(gcnv_nodup$cluster, "_.*", ""), table)
#     sc_clust <- do.call(cbind, sc_clust)
#     denom_clust <- by(gcnv_nodup$sample, str_replace(gcnv_nodup$cluster, "_.*", ""), function(x) length(unique(x)))
#     sf_clust <- t(t(sc_clust)/as.numeric(denom_clust))

#     alarm <- apply(sf_clust, 1, function(x) sum(x==0)==4 & max(x)>=.005)
#     # alarm2 <- apply(sf_clust, 1, function(x) sum(x==0)==3 & max(x)>=.02)

# # sf <- table(gcnv_nodup$var_name)/sum(denom_clust)
# #     maxdif <- apply(sf_clust, 1, max)-apply(sf_clust, 1, min)
# #     plot(as.numeric(maxdif)~as.numeric(sf), xlim=c(0, .02))
# #     ppv$cov[ppv$var_name %in% rownames(sf_clust)[apply(sf_clust, 1, max)>.05 & sf<.01]]
#     ppv[ppv$var_name %in% rownames(sf_clust)[alarm]]
#     ppv[ppv$var_name %in% rownames(sf_clust)[alarm2]]

# ################################################################################################################################################
# ### Harold DN CNV vcf [Sensitivity]
# ################################################################################################################################################
# # samples <- unique(unlist(truth[,7]))
# samples <- unlist(read.table('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/data/truth/Talkowski.denovo.vcf.gz', comment="@", skip=60, nrows=1, stringsAsF=F)[1,])
# # vcf <- read.table('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/data/truth/Talkowski.denovo.vcf.gz', header=FALSE, stringsAsF=F)
# # svtk vcf2bed --split-cpx '/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/data/truth/Talkowski.denovo.vcf.gz' /Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/data/truth/Talkowski.denovo.bed

# dnwgs <- read.table('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/data/truth/Talkowski.denovo.bed', header=FALSE, stringsAsF=F) # 929
#     dnwgs[,1] <- paste0("chr", str_replace(dnwgs[,1], "chr", ""))
#     dnwgs <- dnwgs[dnwgs[,6] %in% trio$SSC_ID[!is.na(trio$SSC_ID)],] # 683
#     dnwgs <- dnwgs[dnwgs[,5] %in% c("DEL", "DUP"),] # 546
#     dnwgs <- dnwgs[!str_detect(dnwgs[,4], "mosaic"),] # 492

#     gr_dnwgs <- GRanges(paste0(dnwgs[,6], "::::", dnwgs[,1]), IRanges(dnwgs[,2], dnwgs[,3]))
#         gr_dnwgs$CNV_ID <- trio$CNV_ID[match(dnwgs[,6], trio$SSC_ID)]
#         gr_dnwgs$batch <- gcnv$batch[match(gr_dnwgs$CNV_ID, gcnv$sample)]
#         gr_dnwgs$wid_raw <- sapply(1:length(gr_dnwgs), function(x) length(toBinSpace(gr_dnwgs[x], list_bins[[gr_dnwgs$batch[x]]], pattern="::::")))
#         gr_dnwgs$wid <- 0
#         gr_dnwgs$wid[gr_dnwgs$wid_raw>0] <- sapply(which(gr_dnwgs$wid_raw>0), function(x) width(toBinSpace(gr_dnwgs[x], list_bins[[gr_dnwgs$batch[x]]], pattern="::::")))
#         gr_dnwgs$var_name <- dnwgs[,4]
#         table(str_replace(gr_dnwgs$var_name, "_.*", ""))

#     # gCNV calls to same subset of samples
#     calls_gcnv_ssc <- dn[dn$sample %in% trio$USE_ID[trio$SSC_ID %in% samples],]
#         calls_gcnv_ssc$sample <- trio$SSC_ID[match(calls_gcnv_ssc$sample, trio$USE_ID)]
#         gr_gcnv_ssc <- GRanges(paste0(calls_gcnv_ssc$sample, "::::", calls_gcnv_ssc$chr), IRanges(calls_gcnv_ssc$start, calls_gcnv_ssc$end))


#     ################################################################################################################################################
#     ### Compare exome to WGS
#     ################################################################################################################################################
#     ol <- findOverlaps(gr_dnwgs, gr_gcnv_ssc)
#     gr_dnwgs$ol <- 0
#         gr_dnwgs$ol[queryHits(ol)] <- subjectHits(ol)
#         gr_dnwgs$inheritance_exome <- "missing"
#         gr_dnwgs$inheritance_exome[queryHits(ol)] <- calls_gcnv_ssc$inheritance[subjectHits(ol)]

#     calls_gcnv_ssc[gr_dnwgs$ol,]

#     ################################################################################################################################################
#     ### Visualize
#     ind <- which(gr_dnwgs$wid>0 & gr_dnwgs$inheritance_exome != "denovo")
#     g_c3 <- dnwgs[ind,]
#         colnames(g_c3) <- c("chr", "start", "end", "name", "call", "sample")
#         g_c3 <- cbind(g_c3, values(gr_dnwgs[ind]))
#     files <- paste0('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/dn_mat/', g_c3$sample, "_dn_", g_c3$chr, ":", g_c3$start, "-", g_c3$end, ".txt")
#     for(i in which(file.exists(files)==F)){
#         message(i)
#         inspect <- g_c3[i,]
#         mat <- match(inspect$sample, trio$SSC_ID)
#         inspect$idp <- trio$FATHER_ID[mat]
#         inspect$idm <- trio$MOTHER_ID[mat]
#         inspect$idp <- m$SSC_ID[match(inspect$idp, m$SFARI_ID)]
#         inspect$idm <- m$SSC_ID[match(inspect$idm, m$SFARI_ID)]
#         inspectPlotDNWGS(inspect, '/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/plots/dn_ssc_missed/')
#     }
   
#     gr_dnwgs_view <- gr_dnwgs[gr_dnwgs$wid>0]
#         ### Variants to remove
#         ## SSC00331 called as de novo, but parental site frequency>1% [remove - denovo_435]
#         ## SSC10143 called, but genotyped as inherited because mother has smaller cover over [leave - denovo_639]
#         ## SSC10336, spans less exons than genotyped, exons involved are outside of the exons [remove - denovo_315]
#         ## SSC12701, overlap of complex event - no actual depth signal [remove - denovo_276]
#         ## SSC02077, overlap of complex event - no actual depth signal [remove - denovo_515]
#         gr_dnwgs_view <- gr_dnwgs[gr_dnwgs$wid>0]
#         gr_dnwgs_view <- gr_dnwgs_view[-which(gr_dnwgs_view$var_name %in% paste0("denovo_", c(435, 315, 276, 515)))]

#         par(mfrow=c(1,2))
#         plot(sapply(1:100, function(x) mean(gr_dnwgs_view$inheritance_exome[gr_dnwgs_view$wid>=x]=="denovo")), ty="b", pch=19, xlim=c(0, 25), ylim=c(0, 1))
#         plot(sapply(1:100, function(x) mean(gr_dnwgs_view$inheritance_exome[gr_dnwgs_view$wid==x]=="denovo")), ty="b", pch=19, xlim=c(0, 25), ylim=c(0, 1))

# ################################################################################################################################################
# ### Check all DN SSC calls are calls in raw callset [ppv]
# ################################################################################################################################################
# gr_truth_red <- reduce(gr_truth[strand(gr_truth)!="*"])
#     gr_truth_red$sample <- str_replace(seqnames(gr_truth_red), "-.*", "")

# results_werling_ppv <- compileResultsPPV(unique(calls_gcnv_ssc$sample), calls_gcnv_ssc, gr_truth_red, QS_dup=50, QS_del=100)
#     ppv <- do.call(c, results_werling_ppv[[2]])
#     view <- ppv[ppv$inheritance=="denovo" & ppv$cov<.3 & ppv$var_name %in% c('var_13597', 'var_323')==F]

#     view <- view[order(view$wid, decreasing=TRUE)]
#         view_mat <- data.frame(chr=str_extract(as.character(seqnames(view)), "chr.*"), start=start(view), end=end(view), values(view))
        
#     for(i in 1:dim(view_mat)[1]){
#         message(i)
#         inspect <- view_mat[i,]
#         mat <- match(inspect$sample, trio$SSC_ID)
#         inspect$idp <- trio$FATHER_ID[mat]
#         inspect$idm <- trio$MOTHER_ID[mat]
#         inspect$idp <- m$SSC_ID[match(inspect$idp, m$SFARI_ID)]
#         inspect$idm <- m$SSC_ID[match(inspect$idm, m$SFARI_ID)]
#         inspectPlotDNWGS(inspect, '/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/plots/dn_ssc_fp/')
#     }