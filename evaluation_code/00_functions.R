#################################################################################
### Function pool
#################################################################################
percBL <- function(gr, bl){
    gr <- GRanges(str_replace(gr, ".*?-", ""))
    ol_segdups <- findOverlaps(gr, bl)
    pint_segdups <- pintersect(gr[queryHits(ol_segdups)], bl[subjectHits(ol_segdups)])
        cov_segdups <- by(width(pint_segdups), queryHits(ol_segdups), sum)
        perc_segdups <- cov_segdups/width(gr[as.numeric(names(cov_segdups))])
    gr$segdup <- 0; gr$segdup[as.numeric(names(perc_segdups))] <- as.numeric(perc_segdups)
    return(gr$segdup)
}

percCov <- function(gr, gr2){
    gr2 <- reduce(gr2)
    ol_segdups <- findOverlaps(gr, gr2)
    pint_segdups <- pintersect(gr[queryHits(ol_segdups)], gr2[subjectHits(ol_segdups)])
        cov_segdups <- by(width(pint_segdups), queryHits(ol_segdups), sum)
        perc_segdups <- cov_segdups/width(gr[as.numeric(names(cov_segdups))])
    gr$segdup <- 0; gr$segdup[as.numeric(names(perc_segdups))] <- as.numeric(perc_segdups)
    return(gr$segdup)
}

extractCallsgCNV <- function(file, chr=FALSE){
    vcf <- readVcf(file)    
    out <- rowRanges(vcf)
        end(out) <- info(vcf)$END
    out$QA <- geno(vcf)$QA
    out$QS <- geno(vcf)$QS
    out$QSS <- geno(vcf)$QSS
    out$QSE <- geno(vcf)$QSE
    out$CN <- as.numeric(unlist(geno(vcf)$CN))
    if(chr){
        calls_x <- out[seqnames(out)=="chrX"]; ploidy_x <- round(sum(width(calls_x)*calls_x$CN)/sum(width(calls_x)))
        calls_y <- out[seqnames(out)=="chrY"]; ploidy_y <- round(sum(width(calls_y)*calls_y$CN)/sum(width(calls_y)))
    }else{
        calls_x <- out[seqnames(out)=="X"]; ploidy_x <- round(sum(width(calls_x)*calls_x$CN)/sum(width(calls_x)))
        calls_y <- out[seqnames(out)=="Y"]; ploidy_y <- round(sum(width(calls_y)*calls_y$CN)/sum(width(calls_y)))
    }
    
    ref_ploidy <- data.frame(chr=c(1:22, "X", "Y"), ploidy=c(rep(2, 22), ploidy_x, ploidy_y))
    if(chr) ref_ploidy[,1] <- paste0("chr", ref_ploidy[,1])
        sex_chr_tot <- sum(ref_ploidy[23:24,2])
        if(sex_chr_tot!=2){
            if(ref_ploidy[24,2]>0){ref_ploidy[23,2] <- 1}
            if(ref_ploidy[24,2]==0){ref_ploidy[23,2] <- 2}
            if(ref_ploidy[24,2]>1){ref_ploidy[24,2] <-1}
        }
    ref_ploidy_matched <- match(as.character(seqnames(out)), ref_ploidy$chr)
    out$ref_ploidy <- ref_ploidy[ref_ploidy_matched, 2]
    out <- out[out$CN != out$ref_ploidy]
    if(length(out)==0){
        return(out)
    }
    out$call <- "DEL"
        out$call[out$CN>out$ref_ploidy] <- "DUP"
    strand(out) <- "-"
        strand(out)[out$call=="DUP"] <- "+"
    names(out) <- NULL
    return(out)
}

extractCallsgCNVnoSex <- function(file, chr=FALSE){
    vcf <- readVcf(file)    
    out <- rowRanges(vcf)
        end(out) <- info(vcf)$END
    out$QA <- geno(vcf)$QA
    out$QS <- geno(vcf)$QS
    out$QSS <- geno(vcf)$QSS
    out$QSE <- geno(vcf)$QSE
    out$CN <- as.numeric(unlist(geno(vcf)$CN))
    
    ref_ploidy <- data.frame(chr=c(1:22), ploidy=rep(2, 22))
        if(chr) ref_ploidy[,1] <- paste0("chr", ref_ploidy[,1])
    ref_ploidy_matched <- match(as.character(seqnames(out)), ref_ploidy$chr)
    out$ref_ploidy <- ref_ploidy[ref_ploidy_matched, 2]
    out <- out[which(out$CN != out$ref_ploidy & !is.na(out$ref_ploidy))]
    out$call <- "DEL"
        out$call[out$CN>out$ref_ploidy] <- "DUP"
    strand(out) <- "-"
        strand(out)[out$call=="DUP"] <- "+"
    names(out) <- NULL
    return(out)
}


processTruth <- function(truth_file, master_bin_file="~/downloads/exons_list.txt", gcnv_samps){
    truth <- read.table(truth_file, stringsAsFactors=FALSE)
        truth <- truth[truth[,5] %in% c("DEL", "DUP", "CN=0"),]
        truth[truth[,5]=="DEL",5] <- "-"; truth[truth[,5]=="DUP",5] <- "+"; truth[truth[,5]=="CN=0",5] <- "*"

    ## Finding only calls that overlap all possible exons 
    master_exons <- read.table(master_bin_file, comment="@")
        gr_master <- GRanges(master_exons[,1], IRanges(master_exons[,2], master_exons[,3]))

    gr_truth <- GRanges(truth[,1], IRanges(as.numeric(truth[,2]), as.numeric(truth[,3])))
    ol_bin_truth <- findOverlaps(gr_truth, gr_master)
        bin_counts <- table(queryHits(ol_bin_truth))
        truth$num_probes <- 0; truth$num_probes[as.numeric(names(bin_counts))] <- as.numeric(bin_counts)
        truth <- truth[truth$num_probes>0,]; tmp <- truth; truth <- NULL
    ## Convert to long-form
    for(i in 1:dim(tmp)[1]){
        message(i)
        base <- tmp[i,-6]
        samps <- tmp[i,6]
        samps <- unlist(str_split(samps, ","))
            samps <- str_replace(samps, "SSC_", ""); samps <- str_replace(samps, "_FA", '.fa')
            samps <- str_replace(samps, "_MO", '.mo'); samps <- str_replace(samps, "_P1", '.p1'); samps <- str_replace(samps, "_S1", '.s1')
            samps <- samps[samps %in% gcnv_samps]
        tmp1 <- matrix(rep(as.vector(base), length(samps)), ncol=6, byrow=TRUE)
        tmp2 <- cbind(tmp1, samps, rep(paste0("var_", i), times=length(samps)), rep(length(samps), length(samps)))
        truth <- rbind(truth, tmp2)
    }
    return(truth)
}

grTruth <- function(truth, chrs=1:22){
    truth <- truth[truth[,1] %in% chrs,]
    gr_truth <- GRanges(paste0(truth[,7], "-", truth[,1]), 
    IRanges(as.numeric(truth[,2]), as.numeric(truth[,3])), strand=unlist(truth[,5]))
    gr_truth$name <- truth[,8]
    gr_truth$sample <- truth[,7]
    ac_truth <- table(unlist(gr_truth$name))
    gr_truth$ac <- ac_truth[match(gr_truth$name, names(ac_truth))]
    gr_truth$af <- gr_truth$ac/length(unique(truth[,7]))
    return(gr_truth)
}

toBinSpace <- function(gr, bins, pattern="-"){
    tmp <- GRanges(str_replace(as.character(seqnames(gr)), paste0(".*", pattern), ""), IRanges(start(gr), end(gr)))
    ol <- findOverlaps(tmp, bins)    
    info <- matrix(unlist(by(subjectHits(ol), queryHits(ol), function(x) c(min(x), max(x)))), ncol=2, byrow=TRUE)
    strds <- strand(gr)[unique(queryHits(ol))]
    out <- GRanges(seqnames(gr)[unique(queryHits(ol))], IRanges(info[,1], info[,2]), strand=strds)
    values(out) <- values(gr)[unique(queryHits(ol)),]
    return(out)
}

toBinSpace2 <- function(gr, bins){
    tmp <- GRanges(str_replace(seqnames(gr), ".*-", ""), IRanges(start(gr), end(gr)))
    ol <- findOverlaps(tmp, bins)    
    
    info <- cbind(unique(queryHits(ol)), matrix(unlist(by(subjectHits(ol), queryHits(ol), function(x) c(min(x), max(x)))), ncol=2, byrow=TRUE))

    out <- GRanges(rep(0, length(gr)), IRanges(0, 0))
        ranges(out[info[,1]]) <- IRanges(info[,2], info[,3])
        strand(out) <- strand(gr)
    return(out)
}

toBinSpace3 <- function(gr, bins){
    tmp <- GRanges(str_replace(seqnames(gr), ".*-", ""), IRanges(start(gr), end(gr)))
    ol <- findOverlaps(tmp, bins)    
    
    info <- cbind(unique(queryHits(ol)), matrix(unlist(by(subjectHits(ol), queryHits(ol), function(x) c(min(x), max(x)))), ncol=2, byrow=TRUE))

    out <- GRanges(seqnames(tmp), IRanges(0, 0))
        ranges(out[info[,1]]) <- IRanges(info[,2], info[,3])
        strand(out) <- strand(gr)
    return(out)
}

getPerf <- function(xx, yy, cov_threshold=0.5, af_threshold=0.01, sample_threshold=0.5, sizes=1:25){
    cov <- suppressWarnings(getRols(xx, yy))
    hits <- xx[cov>=cov_threshold]
    misses <- xx[cov<cov_threshold]

    names_hit <- unlist(hits$name)
    names_miss <- unlist(misses$name)

    cts_hit <- table(names_hit)
    cts_miss <- table(names_miss)
    cts_full <- data.frame(name=sort(unique(unlist(xx$name))))
        cts_full$hit <- 0
        cts_full$miss <- 0
        mat_hit <- match(cts_full$name, names(cts_hit))
        mat_miss <- match(cts_full$name, names(cts_miss))
        cts_full$hit[!is.na(mat_hit)] <- cts_hit[mat_hit[!is.na(mat_hit)]]
        cts_full$miss[!is.na(mat_miss)] <- cts_miss[mat_miss[!is.na(mat_miss)]]
        cts_full$wid <- width(xx[match(cts_full$name, xx$name)])

    cts_full$total <- cts_full$hit+cts_full$miss
    # cts_full$af <- cts_full$total/length(samples) # Internal AF
    cts_full$af <- xx$af[match(cts_full$name, xx$name)] # Preexisting AF
    cts_full$hit_rate <- cts_full$hit/cts_full$total
    cts_full$pass <- cts_full$hit_rate>=sample_threshold
    cts_full$type <- as.character(strand(xx)[match(cts_full$name, xx$name)])
    ## Rare
    rare <- cts_full[cts_full$af <= af_threshold,]
    rare_rate <- sapply(sizes, function(x) c(mean(rare$pass[rare$wid>=x]), sum(rare$wid>=x)))
    rare_del <- cts_full[cts_full$af <= af_threshold & cts_full$type=="-",]
    rare_del_rate <- sapply(sizes, function(x) c(mean(rare_del$pass[rare_del$wid>=x]), sum(rare_del$wid>=x)))
    rare_dup <- cts_full[cts_full$af <= af_threshold & cts_full$type=="+",]
    rare_dup_rate <- sapply(sizes, function(x) c(mean(rare_dup$pass[rare_dup$wid>=x]), sum(rare_dup$wid>=x)))

    ## Common
    common <- cts_full[cts_full$af > af_threshold,]
    common_rate <- sapply(sizes, function(x) c(mean(common$pass[common$wid>=x]), sum(common$wid>=x)))

    return(list(df=data.frame(wid=sizes, rare=rare_rate[1,], rare_ct=rare_rate[2,],
        rare_del=rare_del_rate[1,], rare_del_ct=rare_del_rate[2,],
        rare_dup=rare_dup_rate[1,], rare_dup_ct=rare_dup_rate[2,], 
        common=common_rate[1,], common_ct=common_rate[2,]), raw=cts_full))
}

getRols <- function(x, y){
    ol <- suppressWarnings(findOverlaps(x, y))
    int <- pintersect(x[queryHits(ol)], y[subjectHits(ol)])
    cov <- width(int)/width(x[queryHits(ol)])
    cov <- by(cov, queryHits(ol), sum)
    cov_out <- rep(0, length(x))
        cov_out[as.numeric(names(cov))] <- as.numeric(cov)
    return(cov_out)
}

###
plotPerf <- function(sens, ppv, lab1="sens", lab2="ppv"){
    par(mfrow=c(1, 2))
    par(mar=rep(4, 4))
    par(cex.lab=1, cex.axis=1)
    plot(sens$rare~sens$wid, ylim=c(0, 1.1), ty="b", ylab=lab1, xlab="# exons")
        lines(sens$common~sens$wid, ty="b", col=2)
            text(x=1:25-0.3, y=1.1, pos=1, labels=sens$common_ct, col=2, srt=90, cex=0.6)
            axis(side=3, at=sens$wid, labels=sens$rare_ct, las=2, cex=1)
            abline(h=0.95, lty=2)
            abline(h=1)

    plot(ppv$rare~ppv$wid, ylim=c(0, 1.1), ty="b", ylab=lab2, xlab="# exons")
        lines(ppv$common~ppv$wid, ty="b", col=2)
            text(x=1:25-0.3, y=1.1, pos=1, labels=ppv$common_ct, col=2, srt=90, cex=0.6)
            axis(side=3, at=ppv$wid, labels=ppv$rare_ct, las=2, cex=1)
            abline(h=0.95, lty=2)
            abline(h=1)
}

compileResultsArray <- function(samples, gcnv, gr_array_sample, meta_child, QS_del=100, QS_dup=40){
    cov_sens <- NULL; cov_sens_str <- NULL
    clusters <- as.character(unique(meta_child$batch[meta_child$SFARI_ID %in% samples]))
    ### Divide by batch
    for(j in 1:length(clusters)){
        message(j)
        cluster <- clusters[j]
        samples_sub <- intersect(samples, meta_child$SFARI_ID[meta_child$batch==cluster])
        if(length(samples_sub)>0){
            gcnv_use <- gcnv[gcnv$sample %in% samples_sub,]
            gr_gcnv_sample <- GRanges(paste0(gcnv_use$sample, "-", gcnv_use$chr), IRanges(as.numeric(as.character(gcnv_use$start)), as.numeric(as.character(gcnv_use$end))))
                strand(gr_gcnv_sample) <- "-"; strand(gr_gcnv_sample[gcnv_use$call=="DUP"]) <- "+"
            gr_gcnv_sample$sample <- gcnv_use$sample; gr_gcnv_sample$name <- gcnv_use$name
            gr_gcnv_sample$ac <- gcnv_use$vac; gr_gcnv_sample$af <- gcnv_use$vaf
            gr_gcnv_sample$QS <- gcnv_use$QS
            # gr_gcnv_sample_strigent <- gr_gcnv_sample[which((gr_gcnv_sample$QS>=QS_dup & strand(gr_gcnv_sample)=="+") | (gr_gcnv_sample$QS>=QS_del) & strand(gr_gcnv_sample)=="-")]
            gr_array_sub <- gr_array_sample[gr_array_sample$sample %in% meta_child$SFARI_ID[meta_child$batch==cluster & meta_child$filtered_rawCount==FALSE & meta_child$filt_filtCount==F]]

            cov_sens <- c(cov_sens, findCovGr(gr_array_sub, gr_gcnv_sample, list_bins[[cluster]]))
        }
    }
    return(cov_sens)
}

compileResultsSens <- function(samples, gcnv, truth_set, QS_del=100, QS_dup=50, QS_0=400){
    cov_sens <- NULL; cov_ppv <- NULL; cov_sens_str <- NULL; cov_ppv_str <- NULL
    ### Divide by batch
    for(batch in unique(gcnv$batch)){
        message(batch)
        samples_sub <- intersect(samples, gcnv$sample[gcnv$batch==batch])
        gcnv_use <- gcnv[gcnv$sample %in% samples_sub & gcnv$batch==batch,]
        gr_gcnv_sample <- GRanges(paste0(gcnv_use$sample, "-", gcnv_use$chr), IRanges(as.numeric(as.character(gcnv_use$start)), as.numeric(as.character(gcnv_use$end))))
            strand(gr_gcnv_sample) <- "-"; strand(gr_gcnv_sample[gcnv_use$call=="DUP"]) <- "+"
        gr_gcnv_sample$sample <- gcnv_use$sample; gr_gcnv_sample$name <- gcnv_use$name
        gr_gcnv_sample$ac <- gcnv_use$vac; gr_gcnv_sample$af <- gcnv_use$vaf
        gr_gcnv_sample$QS <- as.numeric(as.character(gcnv_use$QS))
        gr_gcnv_sample$QA <- as.numeric(as.character(gcnv_use$QA))
        gr_gcnv_sample$CN <- gcnv_use$CN
        gr_gcnv_sample$var_name <- gcnv_use$var_name

        gr_gcnv_sample_strigent <- gr_gcnv_sample[which(( as.numeric(as.character(gr_gcnv_sample$QS))>QS_dup & strand(gr_gcnv_sample)=="+") | 
            (as.numeric(as.character(gr_gcnv_sample$QS))>QS_del & strand(gr_gcnv_sample)=="-" & gr_gcnv_sample$CN>0) | gr_gcnv_sample$QS>QS_0)]

        truth_set_sub <- truth_set[truth_set$sample %in% samples_sub]

        cov_sens <- c(cov_sens, findCovGr(truth_set_sub, gr_gcnv_sample, list_bins[[batch]]))
        cov_sens_str <- c(cov_sens_str, findCovGr(truth_set_sub, gr_gcnv_sample_strigent, list_bins[[batch]]))
    }
    return(list(cov_sens, cov_sens_str))
}

compileResultsPPV <- function(samples, gcnv, truth_set, QS_del=100, QS_dup=50, QS_0=400){
    cov_sens <- NULL; cov_ppv <- NULL; cov_sens_str <- NULL; cov_ppv_str <- NULL
    ### Divide by batch
    for(batch in unique(gcnv$batch)){
        message(batch)
        samples_sub <- intersect(samples, gcnv$sample[gcnv$batch==batch])
        gcnv_use <- gcnv[gcnv$sample %in% samples_sub & gcnv$batch==batch,]
        gr_gcnv_sample <- GRanges(paste0(gcnv_use$sample, "-", gcnv_use$chr), IRanges(as.numeric(as.character(gcnv_use$start)), as.numeric(as.character(gcnv_use$end))))
            strand(gr_gcnv_sample) <- "-"; strand(gr_gcnv_sample[gcnv_use$call=="DUP"]) <- "+"
        values(gr_gcnv_sample) <- gcnv_use[,-(1:3)]

        gr_gcnv_sample_strigent <- gr_gcnv_sample[which(( as.numeric(as.character(gr_gcnv_sample$QS))>QS_dup & strand(gr_gcnv_sample)=="+") | 
            (as.numeric(as.character(gr_gcnv_sample$QS))>QS_del & strand(gr_gcnv_sample)=="-" & gr_gcnv_sample$CN>0) | gr_gcnv_sample$QS>QS_0)]

        truth_set_sub <- truth_set[truth_set$sample %in% samples_sub]

        cov_ppv <- c(cov_ppv, suppressWarnings(findCovGr(gr_gcnv_sample, truth_set_sub, list_bins[[batch]])))
        cov_ppv_str <- c(cov_ppv_str, suppressWarnings(findCovGr(gr_gcnv_sample_strigent, truth_set_sub, list_bins[[batch]])))
    }
    return(list(cov_ppv, cov_ppv_str))
}

compileResults3 <- function(samples, gcnv, truth_set, QS_del=100, QS_dup=50){
    cov_sens <- NULL; cov_ppv <- NULL; cov_sens_str <- NULL; cov_ppv_str <- NULL
    clusters <- as.character(unique(gcnv$cluster[gcnv$sample %in% samples]))
    ### Divide by batch
    for(j in 1:length(clusters)){
        message(j)
        cluster <- clusters[j]
        samples_sub <- intersect(samples, gcnv$sample[gcnv$cluster==cluster])
        gcnv_use <- gcnv[gcnv$sample %in% samples_sub & gcnv$cluster==cluster,]
        gr_gcnv_sample <- GRanges(paste0(gcnv_use$sample, "-", gcnv_use$chr), IRanges(as.numeric(as.character(gcnv_use$start)), as.numeric(as.character(gcnv_use$end))))
            strand(gr_gcnv_sample) <- "-"; strand(gr_gcnv_sample[gcnv_use$call=="DUP"]) <- "+"
        gr_gcnv_sample$sample <- gcnv_use$sample; gr_gcnv_sample$name <- gcnv_use$name
        gr_gcnv_sample$ac <- gcnv_use$vac; gr_gcnv_sample$af <- gcnv_use$vaf
        gr_gcnv_sample$QS <- gcnv_use$QS
        gr_gcnv_sample$segdup <- gcnv_use$segdup
        gr_gcnv_sample_strigent <- gr_gcnv_sample[which((gr_gcnv_sample$QS>=QS_dup & strand(gr_gcnv_sample)=="+") | (gr_gcnv_sample$QS>=QS_del) & strand(gr_gcnv_sample)=="-")]

        truth_set_sub <- truth_set[truth_set$sample %in% samples_sub]
        cov_sens <- c(cov_sens, findCov(truth_set_sub[truth_set_sub$segdup<0.3], gr_gcnv_sample, list_bins[[j]]))
        cov_ppv <- c(cov_ppv, suppressWarnings(findCov(gr_gcnv_sample[gr_gcnv_sample$segdup<0.3], truth_set_sub, list_bins[[j]])))
        cov_sens_str <- c(cov_sens_str, findCov(truth_set_sub[truth_set_sub$segdup<0.3], gr_gcnv_sample_strigent, list_bins[[j]]))
        cov_ppv_str <- c(cov_ppv_str, suppressWarnings(findCov(gr_gcnv_sample_strigent[gr_gcnv_sample_strigent$segdup<0.3], truth_set_sub, list_bins[[j]])))
    }
    return(list(cov_sens, cov_ppv, cov_sens_str, cov_ppv_str))
}


findCovGr <- function(x, y, bins){
    x$id <- 1:length(x)
    xx <- toBinSpace(x, bins)
    x <- x[x$id %in% xx$id]
    yy <- toBinSpace(y, bins)
    cov <- rep(0, length(x))
    wid <- width(xx)
    yy <- reduce(yy)
    ol <- suppressWarnings(findOverlaps(xx, yy))
    int <- suppressWarnings(pintersect(xx[queryHits(ol)], yy[subjectHits(ol)]))
    cov_hit <- by(width(int)/width(xx[queryHits(ol)]), queryHits(ol), sum)
    cov[as.numeric(names(cov_hit))] <- cov_hit
    x$cov <- cov
    x$wid <- wid
    return(x)
}

findCov <- function(x, y, bins=NULL){
    if(length(bins)==0){
        ol <- suppressWarnings(findOverlaps(x, y))
        int <- suppressWarnings(pintersect(x[queryHits(ol)], y[subjectHits(ol)]))
        cov <- rep(0, length(x))
        cov_hit <- by(width(int)/width(x[queryHits(ol)]), queryHits(ol), sum)
        cov[as.numeric(names(cov_hit))] <- cov_hit
        return(cov)
    }else{
        x$id <- 1:length(x)
        xx <- toBinSpace(x, bins)
        x <- x[x$id %in% xx$id]
        yy <- toBinSpace(y, bins)
        cov <- rep(0, length(x))
        wid <- rep(0, length(x))
        yy <- reduce(yy)
        ol <- suppressWarnings(findOverlaps(xx, yy))
        int <- suppressWarnings(pintersect(xx[queryHits(ol)], yy[subjectHits(ol)]))
        cov_hit <- by(width(int)/width(xx[queryHits(ol)]), queryHits(ol), sum)
        cov[as.numeric(names(cov_hit))] <- cov_hit
        return(cov)
    }
}


printStats <- function(gr, cov_thresh, af_thresh){
    gr <- suppressWarnings(do.call(c, gr))
    out <- cbind(t(sapply(1:100, function(x) cbind(mean(gr$cov[as.numeric(gr$wid)>=x & gr$af<=af_thresh]>=cov_thresh), sum(as.numeric(gr$wid)>=x & gr$af<=af_thresh)))),
        t(sapply(1:100, function(x) cbind(mean(gr$cov[which(as.numeric(gr$wid)>=x & gr$af<=af_thresh & strand(gr)=="-")]>=cov_thresh), sum(as.numeric(gr$wid)>=x & gr$af<=af_thresh & strand(gr)=="-")))),
        t(sapply(1:100, function(x) cbind(mean(gr$cov[which(as.numeric(gr$wid)>=x & gr$af<=af_thresh & strand(gr)=="+")]>=cov_thresh), sum(as.numeric(gr$wid)>=x & gr$af<=af_thresh & strand(gr)=="+")))))
    return(out)
}


printStatsMar <- function(gr, cov_thresh, af_thresh){
    gr <- suppressWarnings(do.call(c, gr))
    out <- cbind(t(sapply(1:100, function(x) cbind(mean(gr$cov[as.numeric(gr$wid)==x & gr$af<=af_thresh]>=cov_thresh), sum(as.numeric(gr$wid)==x & gr$af<=af_thresh)))),
        t(sapply(1:100, function(x) cbind(mean(gr$cov[which(as.numeric(gr$wid)==x & gr$af<=af_thresh & strand(gr)=="-")]>=cov_thresh), sum(as.numeric(gr$wid)==x & gr$af<=af_thresh & strand(gr)=="-")))),
        t(sapply(1:100, function(x) cbind(mean(gr$cov[which(as.numeric(gr$wid)==x & gr$af<=af_thresh & strand(gr)=="+")]>=cov_thresh), sum(as.numeric(gr$wid)==x & gr$af<=af_thresh & strand(gr)=="+")))))
    return(out)
}


printStatsCommon <- function(gr, cov_thresh, af_thresh){
    gr <- suppressWarnings(do.call(c, gr))
    out <- cbind(t(sapply(1:100, function(x) cbind(mean(gr$cov[as.numeric(gr$wid)>=x & gr$af>af_thresh]>=cov_thresh), sum(as.numeric(gr$wid)>=x & gr$af>af_thresh)))),
        t(sapply(1:100, function(x) cbind(mean(gr$cov[which(as.numeric(gr$wid)>=x & gr$af>af_thresh & strand(gr)=="-")]>=cov_thresh), sum(as.numeric(gr$wid)>=x & gr$af>af_thresh & strand(gr)=="-")))),
        t(sapply(1:100, function(x) cbind(mean(gr$cov[which(as.numeric(gr$wid)>=x & gr$af>af_thresh & strand(gr)=="+")]>=cov_thresh), sum(as.numeric(gr$wid)>=x & gr$af>af_thresh & strand(gr)=="+")))))
    commons <- gr[as.vector(gr$af)>af_thresh]
    common_stat <- by(unlist(commons$cov), unlist(commons$name), function(x) mean(x>=cov_thresh))
    common_name <- by(unlist(commons$wid), unlist(commons$name), function(x) x[1])

    out2 <- cbind(common_name, common_stat)

    return(list(out, out2))
}

returnVars <- function(x, j, i){
    return(output[[x]][[i]][[j]])
}

divideMultiple <- function(svtk){
    if(sum(str_detect(svtk$call_name, ","))){
        svtk_multi <- svtk[str_detect(svtk$call_name, ","),]
        splits <- str_split(svtk_multi$call_name, ",")
        lens <- sapply(splits, length)
        expanded <- data.frame(name=rep(svtk_multi[,1], times=lens), call_name=unlist(splits))
        svtk <- svtk[-which(str_detect(svtk$call_name, ",")),]
        svtk <- rbind(svtk, expanded)
    }
    return(svtk)
}

performancePlots <- function(ind_stats, sens_cum, sens_mar, ppv_cum, ppv_mar){
    par(mfrow=c(2,2))
    par(mar=c(4, 4, 12, 4))
    plot(sens_cum[ind_stats,1], xlab="# of exons", ylab="% with support", main="", ty="b", pch=19, ylim=c(0, 1))
    abline(h=0.95, lty=2)
    lines(sens_cum[,3], col=2, ty="b", pch=19)
    lines(sens_cum[,5], col=4, ty="b", pch=19)
    par(xpd=NA)
    text(x=ind_stats, y=1.2, pos=1, labels=sens_cum[ind_stats,4], col=2, srt=90, cex=1.2)
    text(x=ind_stats, y=1.4, pos=1, labels=sens_cum[ind_stats,6], col=4, srt=90, cex=1.2)
    legend("bottomright", legend=c("overall", "del", "dup"), col=c(1, 2, 4), pch=19)
    axis(side=4, at=0.5, label='cumulative sensitivity. gCNV vs sfari')
    par(xpd=F)

    par(mar=c(4, 4, 12, 4))
    plot(ppv_cum[ind_stats,1], xlab="# of exons", ylab="% with support", main="", ty="b", pch=19, ylim=c(0, 1))
    abline(h=0.95, lty=2)
    lines(ppv_cum[,3], col=2, ty="b", pch=19)
    lines(ppv_cum[,5], col=4, ty="b", pch=19)
    par(xpd=NA)
    text(x=ind_stats-0.3, y=1.2, pos=1, labels=ppv_cum[ind_stats,4], col=2, srt=90, cex=1.2)
    text(x=ind_stats-0.3, y=1.4, pos=1, labels=ppv_cum[ind_stats,6], col=4, srt=90, cex=1.2)
    legend("bottomright", legend=c("overall", "del", "dup"), col=c(1, 2, 4), pch=19)
    axis(side=4, at=0.5, label='cumulative ppv. gCNV vs sfari')
    par(xpd=F)

    plot(sens_mar[ind_stats,1], xlab="# of exons", ylab="% with support", main="", ty="b", pch=19, ylim=c(0, 1))
    abline(h=0.95, lty=2)
    lines(sens_mar[,3], col=2, ty="b", pch=19)
    lines(sens_mar[,5], col=4, ty="b", pch=19)
    par(xpd=NA)
    text(x=ind_stats, y=1.2, pos=1, labels=sens_mar[ind_stats,4], col=2, srt=90, cex=1.2)
    text(x=ind_stats, y=1.4, pos=1, labels=sens_mar[ind_stats,6], col=4, srt=90, cex=1.2)
    legend("bottomright", legend=c("overall", "del", "dup"), col=c(1, 2, 4), pch=19)
    axis(side=4, at=0.5, label='marginal sensitivity. gcnv vs sfari')
    par(xpd=F)

    par(mar=c(4, 4, 12, 4))
    plot(ppv_mar[ind_stats,1], xlab="# of exons", ylab="% with support", main="", ty="b", pch=19, ylim=c(0, 1))
    abline(h=0.95, lty=2)
    lines(ppv_mar[,3], col=2, ty="b", pch=19)
    lines(ppv_mar[,5], col=4, ty="b", pch=19)
    par(xpd=NA)
    text(x=ind_stats-0.3, y=1.2, pos=1, labels=ppv_mar[ind_stats,4], col=2, srt=90, cex=1.2)
    text(x=ind_stats-0.3, y=1.4, pos=1, labels=ppv_mar[ind_stats,6], col=4, srt=90, cex=1.2)
    legend("bottomright", legend=c("overall", "del", "dup"), col=c(1, 2, 4), pch=19)
    axis(side=4, at=0.5, label='marginal ppv. gcnv vs sfari')
    par(xpd=F)
dev.off()
}

# deFragment <- function(tab, bins, extension=0.3){
#     if(dim(tab)[1]>0){
#         gr <- GRanges(paste0(tab$sample, "::", tab$chr), IRanges(as.numeric(as.character(tab$start)), as.numeric(as.character(tab$end))), strand="-")
#             strand(gr[tab$call=="DUP"]) <- "+"
#             gr$QS <- as.numeric(as.character(tab$QS))
#             gr$QA <- as.numeric(as.character(tab$QA))
#             gr$sample <- tab$sample
#             gr$batch <- as.character(tab$batch)
#             gr$call <- tab$call
#             gr$num_exon <- as.numeric(as.character(tab$num_exon))
#             gr$CN <- as.numeric(as.character(tab$CN))
#         bs <- toBinSpace(gr, bins, pattern="::")
#         bse <- bs; start(bse) <- start(bs) - width(bs)*extension; end(bse) <- end(bs) + width(bs)*extension
#         bse <- reduce(bse)
#         ol <- findOverlaps(bs, bse)
#         gr$defragged <- FALSE

#         tbl <- table(subjectHits(ol))
#             ol_change <- ol[subjectHits(ol) %in% names(tbl[tbl>1])]
#             if(length(ol_change)==0){
#                 out <- GRanges(paste0(tab$sample, "::", tab[,1]), IRanges(as.numeric(as.character(tab$start)), as.numeric(as.character(tab$end))), sample=tab$sample, QS=as.numeric(as.character(tab$QS)), QA=as.numeric(as.character(tab$QA)), batch=tab$batch, CN=as.numeric(as.character(tab$CN)), strand="-", num_exon=as.numeric(as.character(tab$num_exon)), call=tab$call)
#                     strand(out)[out$call=="DUP"] <- "+"
#                 out$defragged=FALSE
#                 return(out)
#             }

#         gre <- by(queryHits(ol_change), subjectHits(ol_change), function(x){
#             sub <- gr[x];
#             out <- paste0(c(as.character(seqnames(sub)[1]), min(start(sub)), max(end(sub)), sub$sample[1], max(sub$QS), sum(sub$QA*sub$num_exon)/sum(sub$num_exon), sub$batch[1], sum(sub$CN*sub$num_exon)/sum(sub$num_exon), sum(sub$num_exon), as.character(strand(sub[1])), as.character(sub$call[1])), collapse="::::")
#             return(out)
#         })
#         gr_out <- as.character(gre)
#         gr_out <- str_split_fixed(gr_out, "::::", n=11)
#         out <- GRanges(gr_out[,1], IRanges(as.numeric(gr_out[,2]), as.numeric(gr_out[,3])), sample=gr_out[,4], QS=round(as.numeric(gr_out[,5]), 1), QA=round(as.numeric(gr_out[,6]),1), batch=gr_out[,7], CN=as.numeric(gr_out[,8]), strand=gr_out[,10], num_exon=as.numeric(gr_out[,9]), call=gr_out[,11])
#         gr$defragged <- FALSE
#         out$defragged <- TRUE
#         out <- c(out, gr[-queryHits(ol_change)])
#         return(out)
#     }else{
#         return(GRanges())
#     }
# }

annotate <- function(x, y, truncate=TRUE){
    if(truncate==TRUE) x <- GRanges(str_replace(as.character(x), ".{9}", ""))
    ol <- findOverlaps(x, y)
    pint <- pintersect(x[queryHits(ol)], y[subjectHits(ol)])
    cov <- by(width(pint), queryHits(ol), sum)
    out <- rep(0, length(x))
    out[as.numeric(names(cov))] <- as.numeric(cov)
    out <- round(out/width(x), 3)
    return(out)    
}


exportCallset <- function(callset, path){
    output <- callset[,c("chr", "start", "end", "QS", "call", "CN", "sample", "batch", "var_name", "num_exon")]
    write.table(output, quote=F, row.names=F, col.names=T, sep="\t", file=path)
}



### Visulization code
visualize <- function(var, path='gs://fc-03f45e25-708d-49e8-a139-bd942a3f7f89/dCR/', workDir="~/downloads/", DN=FALSE, meta, SSC=FALSE, force=FALSE){
    matrix_dir <- paste0(workDir, "matrices/")
    dir.create(matrix_dir)
    target <- GRanges(str_replace(seqnames(var), ".*-", ""), ranges(var))
    queryName <- str_replace_all(paste0(var$sample, "_", as.character(var)), ":", "_")

    if(SSC==TRUE){
        # meta_use <- meta[which(meta$SFARI_ID==var$sample),]
        meta_use <- meta[which(meta$Person==var$sample),]
    }else{
        meta_use <- meta[meta$Person==var$sample,]
    }
    

    # meta_off <- meta_child[meta_child$match_id == var$match_id & meta_child$called==1 & meta_child$duplicate==0,]
    if(DN==TRUE){
        meta_father <- meta[meta$Person == meta_use$Father[1] & meta$called==TRUE,]
        meta_mother <- meta[meta$Person == meta_use$Mother[1] & meta$called==TRUE,]
        # call_info_sub <- rbind(call_info_sub, call_info[call_info$sample %in% c(meta_father$Person, meta_mother$Person),])
    }else{
        meta_father <- NULL
        meta_mother <- NULL
    }

    # batches <- unique(c(as.character(meta_off$gcnv_cohort), as.character(meta_father$gcnv_cohort), as.character(meta_mother$gcnv_cohort)))
    batches <- unique(c(meta_use$batch, meta_mother$batch, meta_father$batch))
    tabs <- list.files(matrix_dir, pattern=queryName)
    tabs <- tabs[str_detect(tabs, ".txt")]
    if(force==TRUE){
        if(length(tabs)>=0){
            file.remove(paste0(matrix_dir, tabs))
            tabs <- list.files(matrix_dir, pattern=queryName)
            tabs <- tabs[str_detect(tabs, ".txt")]
        }
    }
    if(length(tabs)!=length(batches)){
        for(batch in batches){
            script <- paste0("GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` tabix -h gs://fc-03f45e25-708d-49e8-a139-bd942a3f7f89/dCR/", batch, ".bed.gz ", as.character(target), " > ", matrix_dir, queryName, "_", batch, ".txt")
            write.table(script, file=paste0(matrix_dir, "vis_script.sh"), quote=F, row.names=F, col.names=F)
            system2("bash", paste0(matrix_dir, "vis_script.sh"))
            file.remove(paste0(matrix_dir, "vis_script.sh"))
        }
        tabs <- list.files(matrix_dir, pattern=queryName)
        tabs <- tabs[str_detect(tabs, ".txt")]
    }
    list_dCR <- lapply(paste0(matrix_dir, tabs), function(x){tab <- read.table(x, header=TRUE, comment="@", sep="\t")})
    if(length(list_dCR)==1){
        merged <- list_dCR[[1]]
        merged[,1:3] <- merged[,c(2, 1, 3)]
    }else{
        merged <- list_dCR[[1]]; list_dCR <- list_dCR[-1]
        while(length(list_dCR)>0){
            merged <- merge(merged, list_dCR[[1]], by="chromStart", all=TRUE, sort=TRUE)
            list_dCR <- list_dCR[-1]
        }
    }
    # rownames(merged) <- round((merged[,1]+merged[,3])/2)
    rownames(merged) <- merged[,1]
    merged <- merged[,-which(str_detect(colnames(merged), "chrom"))]
    colnames(merged) <- str_replace(colnames(merged), "^X", "")

    dCR <- merged
    dCR[dCR<0] <- 0
    dCR[dCR>5] <- 5
    dCR_collapse <- NULL
    indexes <- round(seq(1, min(10, dim(dCR)[1]), length.out=dim(dCR)[1]))
    rn <- NULL
    for(y in unique(indexes)){
        dCR_collapse <- rbind(dCR_collapse, apply(dCR[which(indexes==y),,drop=FALSE], 2, mean, na.rm=TRUE))
        rn <- c(rn, rownames(dCR)[which(indexes==y)][1])
    }
    rownames(dCR_collapse) <- prettyNum(rn, big.mark=",")

    meta_use$Person <- str_replace_all(meta_use$Person, "-", ".")
    meta_mother$Person <- str_replace_all(meta_mother$Person, "-", ".")
    meta_father$Person <- str_replace_all(meta_father$Person, "-", ".")

    title <- paste0(as.character(var), "\nQS:", var$QS, " QA:", var$QA, " CN:", var$CN, "\nnum_exon:", dim(dCR)[1])
    pdf(file=paste0(workDir, queryName, ".pdf"), height=5, width=5)
    matplot(dCR_collapse, ty="l", col=rgb(0, 0, 0, 0.3), ylim=c(0, 5), xaxt="n", xlab="", main=title)
        text(y=-.6, x=1:dim(dCR_collapse)[1], labels=prettyNum(rownames(dCR_collapse), big.mark=","), srt=45, xpd=TRUE, cex=.7)
        lines(dCR_collapse[,match(meta_use$Person[1], colnames(dCR_collapse))], lwd=3, col=(2+(as.character(strand(var))=="+")*2))
        if(DN==TRUE){
            lines(dCR_collapse[,match(meta_father$Person[1], colnames(dCR_collapse))], col=5, lwd=2)
            lines(dCR_collapse[,match(meta_mother$Person[1], colnames(dCR_collapse))], col=6, lwd=2)
            legend("topleft", legend=c("father", "mother"), pch=19, col=c(5, 6), bg='grey')
        }
    dev.off()
}

visualize38 <- function(var, path='gs://fc-8534373c-3ca7-43ae-8f40-7c757b809514/dCR/', workDir="~/downloads/", DN=FALSE, meta, force=FALSE){
    matrix_dir <- paste0(workDir, "matrices/")
    dir.create(matrix_dir)
    target <- GRanges(str_replace(seqnames(var), ".*-", ""), ranges(var))
    queryName <- str_replace_all(paste0(var$sample, "_", as.character(var)), ":", "_")

    meta_off <- meta[which(meta$spid == var$sample & meta$called==1),]
    if(DN==TRUE){
        meta_father <- meta[which(meta$spid %in% meta_off$father & meta$called==TRUE),]
        meta_mother <- meta[which(meta$spid %in% meta_off$mother & meta$called==TRUE),]
    }else{
        meta_father <- NULL
        meta_mother <- NULL
    }

    batches <- unique(c(as.character(meta_off$gcnv_batch), as.character(meta_father$gcnv_batch), as.character(meta_mother$gcnv_batch)))
    tabs <- list.files(matrix_dir, pattern=queryName)
    tabs <- tabs[str_detect(tabs, ".txt")]
    if(force==TRUE){
        if(length(tabs)>=0){
            file.remove(paste0(matrix_dir, tabs))
            tabs <- list.files(matrix_dir, pattern=queryName)
            tabs <- tabs[str_detect(tabs, ".txt")]
        }
    }
    if(length(tabs)!=length(batches)){
        for(batch in batches){
            script <- paste0("GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` tabix -h gs://fc-8534373c-3ca7-43ae-8f40-7c757b809514/dCR/", batch, ".bed.gz ", as.character(target), " > ", matrix_dir, queryName, "_", batch, ".txt")
            write.table(script, file=paste0(matrix_dir, "vis_script.sh"), quote=F, row.names=F, col.names=F)
            system2("bash", paste0(matrix_dir, "vis_script.sh"))
            file.remove(paste0(matrix_dir, "vis_script.sh"))
        }
        tabs <- list.files(matrix_dir, pattern=queryName)
        tabs <- tabs[str_detect(tabs, ".txt")]
    }
    list_dCR <- lapply(paste0(matrix_dir, tabs), function(x){tab <- read.table(x, header=TRUE, comment="@", sep="\t")})
    if(length(list_dCR)==1){
        merged <- list_dCR[[1]]
        merged[,1:3] <- merged[,c(2, 1, 3)]
    }else{
        merged <- list_dCR[[1]]; list_dCR <- list_dCR[-1]
        while(length(list_dCR)>0){
            merged <- merge(merged, list_dCR[[1]], by="chromStart", all=TRUE, sort=TRUE)
            list_dCR <- list_dCR[-1]
        }
    }
    # rownames(merged) <- round((merged[,1]+merged[,3])/2)
    rownames(merged) <- merged[,1]
    merged <- merged[,-which(str_detect(colnames(merged), "chrom"))]
    colnames(merged) <- str_replace(colnames(merged), "^X", "")

    dCR <- merged
    dCR[dCR<0] <- 0
    dCR[dCR>5] <- 5
    dCR_collapse <- NULL
    indexes <- round(seq(1, min(10, dim(dCR)[1]), length.out=dim(dCR)[1]))
    rn <- NULL
    for(y in unique(indexes)){
        dCR_collapse <- rbind(dCR_collapse, apply(dCR[which(indexes==y),,drop=FALSE], 2, mean, na.rm=TRUE))
        rn <- c(rn, rownames(dCR)[which(indexes==y)][1])
    }
    rownames(dCR_collapse) <- prettyNum(rn, big.mark=",")

    title <- paste0(str_replace(as.character(var), "SF[0-9]*", var$sample), "\nQS:", var$QS, " QA:", var$QA, " CN:", var$CN, "\nnum_exon:", dim(dCR)[1])
    pdf(file=paste0(workDir, queryName, ".pdf"), height=5, width=5)
    matplot(dCR_collapse, ty="l", col=rgb(0, 0, 0, 0.3), ylim=c(0, 5), xaxt="n", xlab="", main=title)
        text(y=-.6, x=1:dim(dCR_collapse)[1], labels=prettyNum(rownames(dCR_collapse), big.mark=","), srt=45, xpd=TRUE, cex=.7)
        lines(dCR_collapse[,match(str_replace_all(meta_off$spid, "-", "."), colnames(dCR_collapse))], lwd=3, col=(2+(as.character(strand(var))=="+")*2))
        if(DN==TRUE){
            lines(dCR_collapse[,match(str_replace_all(meta_father$spid, "-", "."), colnames(dCR_collapse))], col=5, lwd=2)
            lines(dCR_collapse[,match(str_replace_all(meta_mother$spid, "-", "."), colnames(dCR_collapse))], col=6, lwd=2)
            legend("topleft", legend=c("father", "mother"), pch=19, col=c(5, 6), bg='grey')
        }
    dev.off()
}

### Remote tabix firecloud
getMatrix <- function(x, hdr){
    tab <- read.table(x, stringsAsF=F)
    rownames(tab) <- paste0(tab[,1], ":", tab[,2], "-", tab[,3])
    tab <- tab[,-(1:3)]
    colnames(tab) <- hdr[-(1:3)]
    return(tab)
}


evalPlots <- function(inspect, QS=FALSE){
    svtype <- rep("DEL", length(inspect)); svtype[which(strand(inspect)=="+")] <- "DUP"
    bed_geno <- cbind(str_replace(as.character(seqnames(inspect)), ".*-", ""), start(inspect), end(inspect), paste0("CNV_", svtype), as.character(inspect$sample), svtype)
    bc <- read.table('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/binCov_loc.tsv', sep="\t", header=TRUE, stringsAsF=F)
    bc_batch <- read.table('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/binCov_loc_batch.tsv', sep="\t", header=TRUE, stringsAsF=F)
    bc_map <- read.table('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/ref_10065.batches', stringsAsF=F)

    for(i in 1:length(inspect)){
        message(i)
        var <- bed_geno[i,]
        samp <- var[5]

        # bc_sub <- bc_map[bc_map[,2]==samp,]
        #     bc_sub <- bc_sub[str_detect(bc_sub[,1], "phase1"),]
        #     phase <- bc_sub[1,1]

        bc_map_sub <- bc_map[bc_map[,1]==samp,2]
            bc_batch_sub <- bc_batch[bc_batch[,1]==bc_map_sub,]

            val_dir <- "/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/plots/mat/"
            locus_info <- var[1:3]; target_locus <- paste0(locus_info[1], ":", locus_info[2], "-", locus_info[3])
            window=100000
            locus_info[2:3] <- as.numeric(locus_info[2:3])

            file_bc <- bc_batch_sub$bincov_matrix
            cmd_mid <- paste0("GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` tabix ", file_bc, " ", target_locus, ' > ', val_dir, target_locus, "_mid.txt")
                write.table(cmd_mid, file="hdr_script.sh", quote=F, row.names=F, col.names=F)
                system2("bash", "hdr_script.sh")
                file.remove("hdr_script.sh")
            cmd_left <- paste0("GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` tabix ", file_bc, " ", paste0(locus_info[1], ":", max(0, as.numeric(locus_info[2])-window), "-", locus_info[2]), ' > ', val_dir, target_locus, "_left.txt")
                write.table(cmd_left, file="hdr_script.sh", quote=F, row.names=F, col.names=F)
                system2("bash", "hdr_script.sh")
                file.remove("hdr_script.sh")
            cmd_right <- paste0("GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` tabix ", file_bc, " ", paste0(locus_info[1], ":", locus_info[3], "-", as.numeric(locus_info[3])+window), ' > ', val_dir, target_locus, "_right.txt")
                write.table(cmd_right, file="hdr_script.sh", quote=F, row.names=F, col.names=F)
                system2("bash", "hdr_script.sh")
                file.remove("hdr_script.sh")
            script_hdr <- paste0("GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` tabix -H ", file_bc, " > ", val_dir, bc_map_sub, ".hdr.txt")
                write.table(script_hdr, file="hdr_script.sh", quote=F, row.names=F, col.names=F)
                system2("bash", "hdr_script.sh")
                file.remove("hdr_script.sh")

            hdr <- as.character(read.table(paste0(val_dir, bc_map_sub, ".hdr.txt"), stringsAsF=F, comment="@"))

            mid <- getMatrix(paste0(val_dir, target_locus, "_mid.txt"), hdr)
            left <- getMatrix(paste0(val_dir, target_locus, "_left.txt"), hdr)
            right <- getMatrix(paste0(val_dir, target_locus, "_right.txt"), hdr)

            resolution <- 10

            norm_fact <- apply(rbind(left, right), 2, median)
            m <- mid
            red <- NULL
            indexes <- round(seq(1, min(resolution, dim(m)[1]), length.out=dim(m)[1]))
            rn <- NULL
            for(y in unique(indexes)){
                red <- rbind(red, apply(m[which(indexes==y),,drop=FALSE], 2, mean, na.rm=TRUE))
                tmp <- rownames(m)[which(indexes==y)]
                    tmp <- str_split_fixed(tmp, ":|-", n=3)
                rn <- c(rn, paste0(tmp[1,1], ":", tmp[1,2], "-", tmp[dim(tmp)[1],3]))
            }
            rownames(red) <- rn
            m_sub <- t(t(red)/norm_fact)

            global <- list.files('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/binCovMed', pattern=paste0(bc_map_sub, "_median"), full=TRUE)
                global <- read.table(global, header=TRUE, stringsAsF=F, sep="\t")
                global <- unlist(global[1,match(colnames(mid), colnames(global)), drop=TRUE])
            m2 <- t(t(mid)/global)

            red2 <- NULL
            indexes <- round(seq(1, min(resolution, dim(m2)[1]), length.out=dim(m2)[1]))
            for(y in unique(indexes)){
                red2 <- rbind(red2, apply(m2[which(indexes==y),,drop=FALSE], 2, mean, na.rm=TRUE))
            }

                bins <- list_bins[[as.character(gcnv$batch[which(gcnv$sample==samp)[1]])]]
                ol <- findOverlaps(bins, GRanges(var[1], IRanges(as.numeric(var[2]), as.numeric(var[3]))))

            pdf(paste0("/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/plots/", target_locus, "-", samp, ".pdf"), width=10, height=5)

            if(QS==FALSE){
                title <- paste0(target_locus, "\n", samp, "\n", length(ol), " exons")
            }else{
                title <- paste0(target_locus, "\n", samp, "\n", length(ol), " exons\nQS: ", inspect$QS[i], " medQS: ", inspect$medQS[i])
            }

            if(as.character(strand(inspect[i]))=="-"){col="red"}else{col="blue"}

                par(mfrow=c(1,2))
                par(las=3)
                par(mar=c(6, 4, 4, 4))
                coords <- apply(str_split_fixed(rn, ":|-", n=3)[,2:3], 1, function(x) mean(as.numeric(x)))
                
                matplot(x=coords,y=2*m_sub, ty="l", col=1, ylim=c(0, 5), main = title, ylab="bincov 1mb neighborhood normalized CN", xaxt="n", xlab="")
                    lines(2*m_sub[,samp]~coords, col=col)
                rect(xleft=start(bins[queryHits(ol)]), xright=end(bins[queryHits(ol)]), ytop=6, ybot=-2, col=rgb(0, 0, 0, 0.2), border=1, lwd=0.1)
                axis(side=1, coords, labels = prettyNum(str_replace(rownames(m_sub), ".*-", ""), big.mark=","), cex=0.5)

                par(las=3)
                par(mar=c(6, 4, 4, 4))
                coords <- apply(str_split_fixed(rn, ":|-", n=3)[,2:3], 1, function(x) mean(as.numeric(x)))
                matplot(x=coords,y=2*(red2+(1-apply(red2, 1, median))), ty="l", col=1, ylim=c(0, 5), main = title, ylab="bincov 1mb neighborhood normalized CN", xaxt="n", xlab="")
                    lines(2*(red2+(1-apply(red2, 1, median)))[,samp]~coords, col=col)
                rect(xleft=start(bins[queryHits(ol)]), xright=end(bins[queryHits(ol)]), ytop=6, ybot=-2, col=rgb(0, 0, 0, 0.2), border=1, lwd=0.1)
                axis(side=1, coords, labels = prettyNum(str_replace(rownames(m_sub), ".*-", ""), big.mark=","), cex=0.5)
            dev.off()
        }
}

xy2ind <- function(x, y, xdim){
    out <- ((y-1)*xdim)+x
    return(out)
}


deFragment <- function(tab, bins, extension=0.3){
    tab$defragged <- FALSE
    tab$cprotect <- FALSE

    gr_tab <- GRanges(paste0(tab$sample, "::::", tab$chr), IRanges(tab$start, tab$end))
    bs_tab <- toBinSpace(gr_tab, bins, pattern="::::")
        values(bs_tab) <- tab

    wid_bs <- width(bs_tab)
        bs_tab_ext <- GRanges(seqnames(bs_tab), 
            IRanges(start(bs_tab) - round(wid_bs*extension), end(bs_tab) + round(wid_bs*extension)))
    values(bs_tab_ext) <- values(bs_tab)

    ### Find overlaps now produced by extension
    ol <- findOverlaps(bs_tab_ext, bs_tab_ext)
        ol <- ol[-which(queryHits(ol)==subjectHits(ol))]
    
    ### Extract CN of overlaps [If there is overlap introduced by extension]
    if(length(ol)>0){
        pairing <- unique(t(apply(cbind(queryHits(ol), subjectHits(ol)), 1, sort)))
            pairing <- cbind(pairing, bs_tab_ext$CN[pairing[,1]], bs_tab_ext$CN[pairing[,2]])
        ### Subset to pairs of compatible CN when extended
        paired <- pairing[pairing[,3]==pairing[,4],, drop=FALSE]

        if(dim(paired)[1]>0){
            ### Create merged-pair with un-extended boundaries
            bs_paired <- GRanges(seqnames(bs_tab)[paired[,1]], 
                IRanges(apply(cbind(start(bs_tab[paired[,1]]), start(bs_tab[paired[,2]])), 1, min),
                apply(cbind(end(bs_tab[paired[,1]]), end(bs_tab[paired[,2]])), 1, max)))
            bs_paired$CN <- bs_tab_ext$CN[paired[,1]]

            ### Check that merged-pair does not bulldoze inconsistent CN calls in the middle
            ol_bull <- findOverlaps(bs_paired, bs_tab)
                clobber <- sapply(1:length(bs_paired), function(x){
                    sum(bs_tab$CN[subjectHits(ol_bull)[queryHits(ol_bull)==x]] != bs_paired$CN[x])
                    })
                protect <- sapply(unique(subjectHits(ol_bull)), function(x){
                    sum(bs_tab$CN[x] != bs_paired$CN[unique(queryHits(ol_bull)[subjectHits(ol_bull)==x])])
                    })
                ### Annotate for when there was clobber protection
                tab$cprotect[unique(subjectHits(ol_bull))] <- protect>0

            ### Reduce non-bulldozing pairs
            bs_paired_noclobber <- bs_paired[clobber==0]
            bs_paired_red <- reduce(bs_paired_noclobber)

            ### Re-insert call-level information
            ol_reinsert <- findOverlaps(bs_paired_red, bs_tab)
            replacements <- do.call(rbind, sapply(1:length(bs_paired_red), function(x){
                vals <- values(bs_tab)[subjectHits(ol_reinsert)[queryHits(ol_reinsert)==x],]
                vals <- vals[order(vals$start),]
                vals_new <- vals[1,]
                vals_new$start <- min(vals$start)
                vals_new$end <- max(vals$end)
                vals_new$QA <- round(mean(vals$QA))
                vals_new$QS <- max(vals$QS)
                vals_new$QSS <- vals$QSS[1]
                vals_new$QSE <- tail(vals$QSE, 1)
                vals_new$num_exon <- sum(vals$num_exon)
                vals_new$defragged <- TRUE
                return(vals_new)
                }))
            ### Now return to genomic-coordinate space, and annotate for whether defragged
            if(length(ol_reinsert)>0){
                tab_out <- tab[-unique(subjectHits(ol_reinsert)),]
                tab_out <- rbind(tab_out, replacements)
                tab_out <- tab_out[order(tab_out$sample, tab_out$chr, tab_out$start),]
            }
            return(tab_out)
        }
    }
    return(tab)
}

inspectPlot <- function(x, bed_geno, dir, sens=FALSE, ppv=FALSE){
    message(i)
    var <- bed_geno[i,]
    samp <- var[5]

    bc_map_sub <- bc_map[bc_map[,1]==samp,2]
        bc_batch_sub <- bc_batch[bc_batch[,1]==bc_map_sub,]

        val_dir <- "/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/plots/mat/"
        locus_info <- var[1:3]; target_locus <- paste0(locus_info[1], ":", locus_info[2], "-", locus_info[3])
        window=100000
        locus_info[2:3] <- as.numeric(locus_info[2:3])
    
    if(sens){out_file <- paste0(dir, "sens/", target_locus, "-", samp, ".pdf")}
    if(ppv){out_file <- paste0(dir, "ppv/", target_locus, "-", samp, ".pdf")}
    if(sens==F & ppv==F){out_file <- paste0(dir, target_locus, "-", samp, ".pdf")}
    if(file.exists(out_file)==F){

        file_bc <- bc_batch_sub$bincov_matrix
        cmd_mid <- paste0("GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` tabix ", file_bc, " ", target_locus, ' > ', val_dir, target_locus, "_mid.txt")
            write.table(cmd_mid, file=paste0(target_locus, "hdr_script.sh"), quote=F, row.names=F, col.names=F)
            system2("bash", paste0(target_locus, "hdr_script.sh"))
            file.remove(paste0(target_locus, "hdr_script.sh"))
        # cmd_left <- paste0("GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` tabix ", file_bc, " ", paste0(locus_info[1], ":", max(0, as.numeric(locus_info[2])-window), "-", locus_info[2]), ' > ', val_dir, target_locus, "_left.txt")
        #     write.table(cmd_left, file="hdr_script.sh", quote=F, row.names=F, col.names=F)
        #     system2("bash", "hdr_script.sh")
        #     file.remove("hdr_script.sh")
        # cmd_right <- paste0("GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` tabix ", file_bc, " ", paste0(locus_info[1], ":", locus_info[3], "-", as.numeric(locus_info[3])+window), ' > ', val_dir, target_locus, "_right.txt")
        #     write.table(cmd_right, file="hdr_script.sh", quote=F, row.names=F, col.names=F)
        #     system2("bash", "hdr_script.sh")
        #     file.remove("hdr_script.sh")
        script_hdr <- paste0("GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` tabix -H ", file_bc, " > ", val_dir, target_locus, ".hdr.txt")
            write.table(script_hdr, file=paste0(target_locus, "hdr_script.sh"), quote=F, row.names=F, col.names=F)
            system2("bash", paste0(target_locus, "hdr_script.sh"))
            file.remove(paste0(target_locus, "hdr_script.sh"))

        hdr <- as.character(read.table(paste0(val_dir, target_locus, ".hdr.txt"), stringsAsF=F, comment="@"))
        mid <- getMatrix(paste0(val_dir, target_locus, "_mid.txt"), hdr)
        # left <- getMatrix(paste0(val_dir, target_locus, "_left.txt"), hdr)
        # right <- getMatrix(paste0(val_dir, target_locus, "_right.txt"), hdr)

        resolution <- 10

        # norm_fact <- apply(rbind(left, right), 2, median)
        m <- mid
        red <- NULL
        indexes <- round(seq(1, min(resolution, dim(m)[1]), length.out=dim(m)[1]))
        rn <- NULL
        for(y in unique(indexes)){
            red <- rbind(red, apply(m[which(indexes==y),,drop=FALSE], 2, mean, na.rm=TRUE))
            tmp <- rownames(m)[which(indexes==y)]
                tmp <- str_split_fixed(tmp, ":|-", n=3)
            rn <- c(rn, paste0(tmp[1,1], ":", tmp[1,2], "-", tmp[dim(tmp)[1],3]))
        }
        rownames(red) <- rn
        m_sub <- t(t(red)/1)

        global <- list.files('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/binCovMed', pattern=paste0(bc_map_sub, "_median"), full=TRUE)
            global <- read.table(global, header=TRUE, stringsAsF=F, sep="\t")
            global <- unlist(global[1,match(colnames(mid), colnames(global)), drop=TRUE])
        m2 <- t(t(mid)/global)

        red2 <- NULL
        indexes <- round(seq(1, min(resolution, dim(m2)[1]), length.out=dim(m2)[1]))
        for(y in unique(indexes)){
            red2 <- rbind(red2, apply(m2[which(indexes==y),,drop=FALSE], 2, mean, na.rm=TRUE))
        }

        bins <- list_bins[[as.character(gcnv_ppv$batch[which(gcnv_ppv$sample==samp)[1]])]]
        ol <- findOverlaps(bins, GRanges(var[1], IRanges(as.numeric(var[2]), as.numeric(var[3]))))

        if(sens){
            pdf(paste0(dir, "sens/", target_locus, "-", samp, ".pdf"), width=10, height=5)
        }
        if(ppv){
            pdf(paste0(dir, "ppv/", target_locus, "-", samp, ".pdf"), width=10, height=5)
        }
        if(sens==F & ppv==F){
            pdf(paste0(dir, target_locus, "-", samp, ".pdf"), width=10, height=5)
        }
            par(mfrow=c(1,2))
            par(las=3)
            par(mar=c(6, 4, 4, 4))
            coords <- apply(str_split_fixed(rn, ":|-", n=3)[,2:3], 1, function(x) mean(as.numeric(x)))
            title <- paste0(target_locus, "\n", samp, "\n", length(ol), " exons")
            matplot(x=coords,y= log(m_sub, 2), ty="l", col=1, main = title, ylab="bincov 1mb neighborhood normalized CN", xaxt="n", xlab="")
                lines(log(m_sub[,samp], 2)~coords, col=2)
            rect(xleft=start(bins[queryHits(ol)]), xright=end(bins[queryHits(ol)]), ytop=6, ybot=-2, col=rgb(0, 0, 0, 0.2), border=1, lwd=0.1)
            axis(side=1, coords, labels = prettyNum(str_replace(rownames(m_sub), ".*-", ""), big.mark=","), cex=0.5)

            par(las=3)
            par(mar=c(6, 4, 4, 4))
            coords <- apply(str_split_fixed(rn, ":|-", n=3)[,2:3], 1, function(x) mean(as.numeric(x)))
            title <- paste0(target_locus, "\n", samp, "\n", length(ol), " exons")
            matplot(x=coords,y=2*(red2+(1-apply(red2, 1, median))), ty="l", col=1, ylim=c(0, 5), main = title, ylab="bincov 1mb neighborhood normalized CN", xaxt="n", xlab="")
                lines(2*(red2+(1-apply(red2, 1, median)))[,samp]~coords, col=2)
            rect(xleft=start(bins[queryHits(ol)]), xright=end(bins[queryHits(ol)]), ytop=6, ybot=-2, col=rgb(0, 0, 0, 0.2), border=1, lwd=0.1)
            axis(side=1, coords, labels = prettyNum(str_replace(rownames(m_sub), ".*-", ""), big.mark=","), cex=0.5)
        dev.off()

        if(sens){
            coords <- apply(str_split_fixed(rn, ":|-", n=3)[,2:3], 1, function(x) mean(as.numeric(x)))
            info <- 2*(red2+(1-apply(red2, 1, median)))
                out <- coords[which(abs(info[,samp]-2)>.8)]
                if(length(out)>0){
                    out <- c(min(out), max(out))
                    out <- c(samp, str_split_fixed(rn, ":|-", n=3)[1,1], out)
                    write.table("resize", file=paste0(dir, "sens_resized/", out[2], ":", out[3], "-", out[4], "-", out[1], ":", var[6], ".pdf"))    
                }       
        }
        if(ppv){
            coords <- apply(str_split_fixed(rn, ":|-", n=3)[,2:3], 1, function(x) mean(as.numeric(x)))
            info <- 2*(red2+(1-apply(red2, 1, median)))
                out <- coords[which(abs(info[,samp]-2)>.4)]
                if(length(out)>0){
                    out <- c(min(out), max(out))
                    out <- c(samp, str_split_fixed(rn, ":|-", n=3)[1,1], out)
                    write.table("resize", file=paste0(dir, "ppv_resized/", out[2], ":", out[3], "-", out[4], "-", out[1], ":", var[6], ".pdf"))    
                }       
        }
    }
}

inspectPlot2 <- function(var, dir, resize=FALSE){
    samp <- as.character(var$sample)

    val_dir <- "/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/plots/mat/"
    locus_info <-cbind(str_replace(seqnames(var), ".*-", ""), start(var), end(var)); target_locus <- paste0(locus_info[1,1], ":", locus_info[1,2], "-", locus_info[1,3])

    bc_map_sub <- bc_map[bc_map[,1]==samp,2]
    bc_batch_sub <- bc_batch[bc_batch[,1]==bc_map_sub,]
    file_bc <- bc_batch_sub$bincov_matrix
    cmd_mid <- paste0("GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` tabix ", file_bc, " ", target_locus, ' > ', val_dir, target_locus, "_mid.txt")
        write.table(cmd_mid, file="hdr_script.sh", quote=F, row.names=F, col.names=F)
        system2("bash", "hdr_script.sh")
        file.remove("hdr_script.sh")
    script_hdr <- paste0("GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` tabix -H ", file_bc, " > ", val_dir, bc_map_sub, ".hdr.txt")
            write.table(script_hdr, file="hdr_script.sh", quote=F, row.names=F, col.names=F)
            system2("bash", "hdr_script.sh")
            file.remove("hdr_script.sh")

    hdr <- as.character(read.table(paste0(val_dir, bc_map_sub, ".hdr.txt"), stringsAsF=F, comment="@"))
    mid <- getMatrix(paste0(val_dir, target_locus, "_mid.txt"), hdr)

    resolution <- 10
    m_tmp <- mid
    red <- NULL
    indexes <- round(seq(1, min(resolution, dim(m_tmp)[1]), length.out=dim(m_tmp)[1]))
    rn <- NULL
    for(y in unique(indexes)){
        red <- rbind(red, apply(m_tmp[which(indexes==y),,drop=FALSE], 2, mean, na.rm=TRUE))
        tmp <- rownames(m_tmp)[which(indexes==y)]
            tmp <- str_split_fixed(tmp, ":|-", n=3)
        rn <- c(rn, paste0(tmp[1,1], ":", tmp[1,2], "-", tmp[dim(tmp)[1],3]))
    }
    rownames(red) <- rn
    m_sub <- t(t(red)/1)

    global <- list.files('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/binCovMed', pattern=paste0(bc_map_sub, "_median"), full=TRUE)
        global <- read.table(global, header=TRUE, stringsAsF=F, sep="\t")
        global <- unlist(global[1,match(colnames(mid), colnames(global)), drop=TRUE])
    m2 <- t(t(mid)/global)

    red2 <- NULL
    indexes <- round(seq(1, min(resolution, dim(m2)[1]), length.out=dim(m2)[1]))
    for(y in unique(indexes)){
        red2 <- rbind(red2, apply(m2[which(indexes==y),,drop=FALSE], 2, mean, na.rm=TRUE))
    }

    bins <- list_bins[[as.character(var$batch)]]
    ol <- findOverlaps(bins, GRanges(target_locus))

    if(var$call[1]=="DEL"){col <- "red"}else{col <- "blue"}

    pdf(paste0(dir, target_locus, "_", samp, ".pdf"), width=10, height=10)
        par(mfrow=c(2,2))
        par(las=3)
        par(mar=c(6, 4, 4, 4))
        coords <- apply(str_split_fixed(rn, ":|-", n=3)[,2:3], 1, function(x) mean(as.numeric(x)))
        title <- paste0(target_locus, "\n", samp, "\n", var$num_exon, " exons")
        matplot(x=coords,y= log(m_sub, 2), ty="l", col=1, main = title, ylab="bincov 1mb neighborhood normalized CN", xaxt="n", xlab="")
            lines(log(m_sub[,samp], 2)~coords, col=col)
        rect(xleft=start(bins[queryHits(ol)]), xright=end(bins[queryHits(ol)]), ytop=100, ybot=-2, col=rgb(0, 0, 0, 0.1), border=0, lwd=0.1)
        axis(side=1, coords, labels = prettyNum(str_replace(rownames(m_sub), ".*-", ""), big.mark=","), cex=0.5)

        par(las=3)
        par(mar=c(6, 4, 4, 4))
        # if(var[1,1] %in% c("chrX", "chrY")==TRUE){adj <- 0}else{adj <- (1-apply(red2, 1, median))}
        coords <- apply(str_split_fixed(rn, ":|-", n=3)[,2:3], 1, function(x) mean(as.numeric(x)))
        title <- paste0(target_locus, "\n", samp, "\nQS: ", var$QS)
        matplot(x=coords,y=2*(red2), ty="l", col=1, ylim=c(0, 5), main = title, ylab="bincov 1mb neighborhood normalized CN", xaxt="n", xlab="")
            lines(2*(red2)[,samp]~coords, col=col, lwd=2)
        rect(xleft=start(bins[queryHits(ol)]), xright=end(bins[queryHits(ol)]), ytop=6, ybot=-2, col=rgb(0, 0, 0, 0.1), border=0, lwd=0.1)
        axis(side=1, coords, labels = prettyNum(str_replace(rownames(m_sub), ".*-", ""), big.mark=","), cex=0.5)

        xy <- par("usr")
        par(xpd=NA)
        points(x=xy[1], y=mean(2*(red2)[,samp]), pch=19, cex=2, col=col)
        par(xpd=F)
        
    dcr_child <- as.character(dcr_index[dcr_index==var$CNV_ID,2])
    dcr_bg <- str_extract(dcr_child, ".*/call-PostprocessGermlineCNVCalls/")
    dcr_bg <- system2('gsutil', paste0('-m ls ', dcr_bg, '**.tsv'), stdout=TRUE)
        set.seed(1337)
        dcr_sel <- dcr_bg[sample(1:length(dcr_bg),min(25, length(dcr_bg)))]

    system2('gsutil', paste0('-m cp ', paste0(c(dcr_child, dcr_sel), collapse=" "), ' /Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/dcr/'))
    bm_sel <- c(basename(dcr_child), basename(dcr_sel))
    tab_sel <- mclapply(paste0('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/dcr/', bm_sel), function(x) read.table(x, sep="\t", comment="@", header=TRUE, stringsAsF=F)[,4], mc.cores=4)
        tab_sel <- do.call(cbind, tab_sel)
        coords_sel <- read.table(paste0('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/dcr/', bm_sel)[1], comment="@", stringsAsF=F, header=TRUE)[,1:3]
        rownames(tab_sel) <- paste0(coords_sel[,1], ":", coords_sel[,2], "-", coords_sel[,3])
        colnames(tab_sel) <- str_replace(str_replace(bm_sel, "denoised_copy_ratios-", ""), ".tsv", "")

    gr_target <- GRanges(target_locus)
    gr_sel <- GRanges(rownames(tab_sel))
    
    ol <- findOverlaps(gr_target, gr_sel)
        coord <- ((end(gr_sel)+start(gr_sel))/2)[subjectHits(ol)]

    matplot(x=coord, y=tab_sel[subjectHits(ol),,drop=FALSE], ty="l", col=1, ylim=c(0, 5), pch=19, ylab="dCR", main="Exome", xlab="", xaxt="n", xlim=c(min(coords), max(coords)))
        lines(x=coord, y=tab_sel[subjectHits(ol),var$CNV_ID], col=col, lwd=2)
        abline(v=coord, lwd=2, col=rgb(0, 0, 0, 0.1))
    axis(side=1, at=coord, labels=prettyNum(coord, big.mark=","))

    grmast <- gr_master; start(grmast) <- start(grmast)+1
        xids <- match(gr_sel[subjectHits(ol)], grmast)
        
    matplot(y=tab_sel[subjectHits(ol),,drop=FALSE], x=xids, ty="l", col=1, ylim=c(0, 5), pch=19, ylab="dCR", xaxt="n", main="Exome", xlab="", xlim=c(min(xids)-.5, max(xids)+.5))
        lines(x=xids, y=tab_sel[subjectHits(ol),var$CNV_ID], col=col, lwd=2)
        abline(v=xids, col=rgb(0, 0, 0, 0.1), lwd=2)
    axis(side=1, at=xids, labels=1:length(xids), las=0)

    xy <- par("usr")
    par(xpd=NA)
    points(x=xy[1], y=mean(tab_sel[subjectHits(ol),var$CNV_ID]), pch=19, cex=2, col=col)
    par(xpd=FALSE)

    file.remove(list.files('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/dcr/', full=TRUE))
    dev.off()
}

inspectPlotDNWGS <- function(var, dir, resize=FALSE){
    samp <- as.character(var$sample)
    sampm <- as.character(var$idm)
    sampp <- as.character(var$idp)

    val_dir <- "/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/plots/mat/"
    locus_info <- var[1:3]; target_locus <- paste0(locus_info[1,1], ":", locus_info[1,2], "-", locus_info[1,3])

    bc_map_sub <- bc_map[bc_map[,1]==samp,2]
    bc_batch_sub <- bc_batch[bc_batch[,1]==bc_map_sub,]
    file_bc <- bc_batch_sub$bincov_matrix
    cmd_mid <- paste0("GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` tabix ", file_bc, " ", target_locus, ' > ', val_dir, target_locus, "_mid.txt")
        write.table(cmd_mid, file="hdr_script.sh", quote=F, row.names=F, col.names=F)
        system2("bash", "hdr_script.sh")
        file.remove("hdr_script.sh")
    script_hdr <- paste0("GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` tabix -H ", file_bc, " > ", val_dir, bc_map_sub, ".hdr.txt")
            write.table(script_hdr, file="hdr_script.sh", quote=F, row.names=F, col.names=F)
            system2("bash", "hdr_script.sh")
            file.remove("hdr_script.sh")

    hdr <- as.character(read.table(paste0(val_dir, bc_map_sub, ".hdr.txt"), stringsAsF=F, comment="@"))
    mid <- getMatrix(paste0(val_dir, target_locus, "_mid.txt"), hdr)

    resolution <- 10
    m_tmp <- mid
    red <- NULL
    indexes <- round(seq(1, min(resolution, dim(m_tmp)[1]), length.out=dim(m_tmp)[1]))
    rn <- NULL
    for(y in unique(indexes)){
        red <- rbind(red, apply(m_tmp[which(indexes==y),,drop=FALSE], 2, mean, na.rm=TRUE))
        tmp <- rownames(m_tmp)[which(indexes==y)]
            tmp <- str_split_fixed(tmp, ":|-", n=3)
        rn <- c(rn, paste0(tmp[1,1], ":", tmp[1,2], "-", tmp[dim(tmp)[1],3]))
    }
    rownames(red) <- rn
    m_sub <- t(t(red)/1)

    global <- list.files('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/binCovMed', pattern=paste0(bc_map_sub, "_median"), full=TRUE)
        global <- read.table(global, header=TRUE, stringsAsF=F, sep="\t")
        global <- unlist(global[1,match(colnames(mid), colnames(global)), drop=TRUE])
    m2 <- t(t(mid)/global)

    red2 <- NULL
    indexes <- round(seq(1, min(resolution, dim(m2)[1]), length.out=dim(m2)[1]))
    for(y in unique(indexes)){
        red2 <- rbind(red2, apply(m2[which(indexes==y),,drop=FALSE], 2, mean, na.rm=TRUE))
    }

    bins <- list_bins[[as.character(var$batch)]]
    ol <- findOverlaps(bins, GRanges(target_locus))

    if(var$call[1]=="DEL"){col <- "red"}else{col <- "blue"}

    pdf(paste0(dir, samp, "_dn_", target_locus, ".pdf"), width=10, height=10)
        par(mfrow=c(2,2))
        par(las=3)
        par(mar=c(6, 4, 4, 4))
        coords <- apply(str_split_fixed(rn, ":|-", n=3)[,2:3], 1, function(x) mean(as.numeric(x)))
        title <- paste0(target_locus, "\n", samp, "\n", var$num_exon, " exons")
        matplot(x=coords,y= log(m_sub, 2), ty="l", col=1, main = title, ylab="bincov 1mb neighborhood normalized CN", xaxt="n", xlab="")
            lines(log(m_sub[,samp], 2)~coords, col=col)
            tryCatch(lines(log(m_sub[,sampm], 2)~coords, col="orange"), error=function(e){})
            tryCatch(lines(log(m_sub[,sampp], 2)~coords, col="purple"), error=function(e){})
        rect(xleft=start(bins[queryHits(ol)]), xright=end(bins[queryHits(ol)]), ytop=100, ybot=-2, col=rgb(0, 0, 0, 0.1), border=0, lwd=0.1)
        axis(side=1, coords, labels = prettyNum(str_replace(rownames(m_sub), ".*-", ""), big.mark=","), cex=0.5)

        par(las=3)
        par(mar=c(6, 4, 4, 4))
        if(var[1,1] %in% c("chrX", "chrY")==TRUE){adj <- 0}else{adj <- (1-apply(red2, 1, median))}
        coords <- apply(str_split_fixed(rn, ":|-", n=3)[,2:3], 1, function(x) mean(as.numeric(x)))
        title <- paste0(target_locus, "\n", samp, "\nQS: ", var$QS)
        matplot(x=coords,y=2*(red2+adj), ty="l", col=1, ylim=c(0, 5), main = title, ylab="bincov 1mb neighborhood normalized CN", xaxt="n", xlab="")
            lines(2*(red2+adj)[,samp]~coords, col=col, lwd=2)
            tryCatch(lines(2*(red2+adj)[,sampm]~coords, col="orange", lwd=2, lty=2), error=function(e){})
            tryCatch(lines(2*(red2+adj)[,sampp]~coords, col="purple", lwd=2, lty=2), error=function(e){})
        rect(xleft=start(bins[queryHits(ol)]), xright=end(bins[queryHits(ol)]), ytop=6, ybot=-2, col=rgb(0, 0, 0, 0.1), border=0, lwd=0.1)
        legend("topright", legend=c("child", "mother", "father"), col=c(col, "orange", "purple"), pch=19, bg="white")
        axis(side=1, coords, labels = prettyNum(str_replace(rownames(m_sub), ".*-", ""), big.mark=","), cex=0.5)

        xy <- par("usr")
        par(xpd=NA)
        points(x=xy[1], y=mean(2*(red2+adj)[,samp]), pch=19, cex=2, col=col)
        tryCatch(points(x=xy[1], y=mean(2*(red2+adj)[,sampm]), pch=19, cex=2, col="orange"), error=function(e){})
        tryCatch(points(x=xy[1], y=mean(2*(red2+adj)[,sampp]), pch=19, cex=2, col="purple"), error=function(e){})
        par(xpd=F)
        
    dcr_child <- as.character(dcr_index[dcr_index==var$CNV_ID,2])
    dad <- map[match(var$idp, map[,2]),1]
        dcr_dad <- as.character(dcr_index[dcr_index[,1] %in% unique(gcnv_nodup$CNV_ID[gcnv_nodup$sample==dad]),2])
        if(length(dcr_dad)==0){
            dad <- m$CNV_ID[which(m$SSC_ID==var$idp)]
            dcr_dad <- as.character(dcr_index[which(dcr_index[,1]==dad),2])
        }
    mom <- map[match(var$idm, map[,2]),1]
        dcr_mom <- as.character(dcr_index[dcr_index[,1] %in% unique(gcnv_nodup$CNV_ID[gcnv_nodup$sample==mom]),2])
        if(length(dcr_mom)==0){
            mom <- m$CNV_ID[which(m$SSC_ID==var$idm)]
            dcr_mom <- as.character(dcr_index[which(dcr_index[,1]==mom),2])
        }
        
    dcr_bg <- str_extract(dcr_child, ".*/call-PostprocessGermlineCNVCalls/")
    dcr_bg <- system2('gsutil', paste0('-m ls ', dcr_bg, '**.tsv'), stdout=TRUE)
        set.seed(1337)
        dcr_sel <- dcr_bg[sample(1:length(dcr_bg),min(25, length(dcr_bg)))]

    system2('gsutil', paste0('-m cp ', paste0(c(dcr_child, dcr_dad, dcr_mom, dcr_sel), collapse=" "), ' /Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/dcr/'))
    bm_sel <- c(basename(dcr_child), basename(dcr_sel))
    tab_sel <- mclapply(paste0('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/dcr/', bm_sel), function(x) read.table(x, sep="\t", comment="@", header=TRUE, stringsAsF=F)[,4], mc.cores=4)
        tab_sel <- do.call(cbind, tab_sel)
        coords_sel <- read.table(paste0('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/dcr/', bm_sel)[1], comment="@", stringsAsF=F, header=TRUE)[,1:3]
        rownames(tab_sel) <- paste0(coords_sel[,1], ":", coords_sel[,2], "-", coords_sel[,3])
        colnames(tab_sel) <- str_replace(str_replace(bm_sel, "denoised_copy_ratios-", ""), ".tsv", "")

        tab_dad <- lapply(basename(dcr_dad), function(x) read.table(paste0('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/dcr/', x), header=TRUE, stringsAsF=F, comment="@"))
        tab_mom <- lapply(basename(dcr_mom), function(x) read.table(paste0('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/dcr/', x), header=TRUE, stringsAsF=F, comment="@"))

    gr_target <- GRanges(target_locus)
    gr_sel <- GRanges(rownames(tab_sel))
    gr_dad <- lapply(tab_dad, function(x) GRanges(x[,1], IRanges(as.numeric(x[,2]), as.numeric(x[,3]))))
    gr_mom <- lapply(tab_mom, function(x) GRanges(x[,1], IRanges(as.numeric(x[,2]), as.numeric(x[,3]))))

    ol <- findOverlaps(gr_target, gr_sel)
        coord <- ((end(gr_sel)+start(gr_sel))/2)[subjectHits(ol)]

    old <- lapply(gr_dad, function(x) findOverlaps(gr_target, x))
        coordd <- lapply(1:length(old), function(x) ((end(gr_dad[[x]])+start(gr_dad[[x]]))/2)[subjectHits(old[[x]])])
    olm <- lapply(gr_mom, function(x) findOverlaps(gr_target, x))
        coordm <- lapply(1:length(olm), function(x) ((end(gr_mom[[x]])+start(gr_mom[[x]]))/2)[subjectHits(olm[[x]])])

    matplot(x=coord, y=tab_sel[subjectHits(ol),,drop=FALSE], ty="l", col=1, ylim=c(0, 5), pch=19, ylab="dCR", main="Exome", xlab="", xaxt="n", xlim=c(min(coords), max(coords)))
        lines(x=coord, y=tab_sel[subjectHits(ol),var$CNV_ID], col=col, lwd=2)
        sapply(1:length(old), function(x) lines(x=coordd[[x]], y=tab_dad[[x]][subjectHits(old[[x]]),4], col="purple", lwd=2))
        sapply(1:length(olm), function(x) lines(x=coordm[[x]], y=tab_mom[[x]][subjectHits(olm[[x]]),4], col="orange", lwd=2))
        abline(v=coord, lwd=2, col=rgb(0, 0, 0, 0.1))
    axis(side=1, at=coord, labels=prettyNum(coord, big.mark=","))

    grmast <- gr_master; start(grmast) <- start(grmast)+1
        xids <- match(gr_sel[subjectHits(ol)], grmast)
        xids_m <- lapply(1:length(gr_mom), function(x) match(grmast[xids], gr_mom[[x]]))
        xids_p <- lapply(1:length(gr_dad), function(x) match(grmast[xids], gr_dad[[x]]))

    matplot(y=tab_sel[subjectHits(ol),,drop=FALSE], x=xids, ty="l", col=1, ylim=c(0, 5), pch=19, ylab="dCR", xaxt="n", main="Exome", xlab="", xlim=c(min(xids)-.5, max(xids)+.5))
        lines(x=xids, y=tab_sel[subjectHits(ol),var$CNV_ID], col=col, lwd=2)
        sapply(1:length(olm), function(x) lines(x=xids[!is.na(xids_m[[x]])], y=tab_mom[[x]][xids_m[[x]][!is.na(xids_m[[x]])],4], col="orange", lwd=2))
        sapply(1:length(old), function(x) lines(x=xids[!is.na(xids_p[[x]])], y=tab_dad[[x]][xids_p[[x]][!is.na(xids_p[[x]])],4], col="purple", lwd=2))
        abline(v=xids, col=rgb(0, 0, 0, 0.1), lwd=2)
    axis(side=1, at=xids, labels=1:length(xids), las=0)

    xy <- par("usr")
    par(xpd=NA)
    points(x=xy[1], y=mean(tab_sel[subjectHits(ol),var$CNV_ID]), pch=19, cex=2, col=col)
    sapply(1:length(olm), function(x) points(x=xy[1], y=mean(tab_mom[[x]][xids_m[[x]][!is.na(xids_m[[x]])],4]), col="orange", pch=19, cex=2))
    sapply(1:length(old), function(x) points(x=xy[1], y=mean(tab_dad[[x]][xids_p[[x]][!is.na(xids_p[[x]])],4]), col="purple", pch=19, cex=2))
    par(xpd=FALSE)

    st_dad <- lapply(1:length(old), function(x) tab_dad[[x]][xids_p[[x]],4,drop=FALSE])
        for(i in 1:length(old)){colnames(st_dad[[i]]) <- paste0("DAD_", str_replace(str_replace(basename(dcr_dad[i]), "denoised_copy_ratios-", ""), ".tsv", ""))}
        st_dad <- do.call(cbind, st_dad)
    st_mom <- lapply(1:length(olm), function(x) tab_mom[[x]][xids_m[[x]],4,drop=FALSE])
        for(i in 1:length(olm)){colnames(st_mom[[i]]) <- paste0("MOM_", str_replace(str_replace(basename(dcr_mom[i]), "denoised_copy_ratios-", ""), ".tsv", ""))}
        st_mom <- do.call(cbind, st_mom)

    info <- cbind(tab_sel[subjectHits(ol), 1], st_mom, st_dad, tab_sel[subjectHits(ol), -1,drop=FALSE])
        colnames(info)[1] <- var$sample
        write.table(info, file=paste0('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/dn_mat/', samp, "_dn_", target_locus, ".txt"), quote=F, sep="\t")

    file.remove(list.files('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/dcr/', full=TRUE))
    dev.off()
}

inspectPlotDN <- function(var, dir, resize=FALSE, XY=FALSE){
    samp <- var$sample
    sampm <- var$idm
    sampp <- var$idp

    val_dir <- "/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/plots/mat/"
    locus_info <- var[1:3]; target_locus <- paste0(locus_info[1,1], ":", locus_info[1,2], "-", locus_info[1,3])
        
    dcr_child <- as.character(dcr_index[dcr_index==var$CNV_ID,2])
    dad <- gcnv_nodup$CNV_ID[match(sampp, gcnv_nodup$sample)]
        dcr_dad <- as.character(dcr_index[dcr_index[,1] %in% dad,2])
    mom <- gcnv_nodup$CNV_ID[match(sampm, gcnv_nodup$sample)]
        dcr_mom <- as.character(dcr_index[dcr_index[,1] %in% mom,2])

    dcr_bg <- str_extract(dcr_child, ".*/call-PostprocessGermlineCNVCalls/")
    dcr_bg <- system2('gsutil', paste0('-m ls ', dcr_bg, '**.tsv'), stdout=TRUE)
        set.seed(1337)
        dcr_sel <- dcr_bg[sample(1:length(dcr_bg),min(25, length(dcr_bg)))]

    system2('gsutil', paste0('-m cp ', paste0(c(dcr_child, dcr_dad, dcr_mom, dcr_sel), collapse=" "), ' /Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/dcr/'))
    bm_sel <- c(basename(dcr_child), basename(dcr_sel))
    tab_sel <- mclapply(paste0('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/dcr/', bm_sel), function(x) read.table(x, sep="\t", comment="@", header=TRUE, stringsAsF=F)[,4], mc.cores=4)
        tab_sel <- do.call(cbind, tab_sel)
        coords_sel <- read.table(paste0('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/dcr/', bm_sel)[1], comment="@", stringsAsF=F, header=TRUE)[,1:3]
        rownames(tab_sel) <- paste0('chr', str_replace(coords_sel[,1], 'chr', ''), ":", coords_sel[,2], "-", coords_sel[,3])
        colnames(tab_sel) <- str_replace(str_replace(bm_sel, "denoised_copy_ratios-", ""), ".tsv", "")

        tab_dad <- lapply(basename(dcr_dad), function(x) read.table(paste0('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/dcr/', x), header=TRUE, stringsAsF=F, comment="@"))
        tab_mom <- lapply(basename(dcr_mom), function(x) read.table(paste0('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/dcr/', x), header=TRUE, stringsAsF=F, comment="@"))

    gr_target <- GRanges(target_locus)
    gr_sel <- GRanges(rownames(tab_sel))
    gr_dad <- lapply(tab_dad, function(x) GRanges(paste0('chr', str_replace(x[,1], 'chr', '')), IRanges(as.numeric(x[,2]), as.numeric(x[,3]))))
    gr_mom <- lapply(tab_mom, function(x) GRanges(paste0('chr', str_replace(x[,1], 'chr', '')), IRanges(as.numeric(x[,2]), as.numeric(x[,3]))))

    ol <- findOverlaps(gr_target, gr_sel)
        coord <- ((end(gr_sel)+start(gr_sel))/2)[subjectHits(ol)]

    old <- lapply(gr_dad, function(x) findOverlaps(gr_target, x))
        coordd <- lapply(1:length(old), function(x) ((end(gr_dad[[x]])+start(gr_dad[[x]]))/2)[subjectHits(old[[x]])])
    olm <- lapply(gr_mom, function(x) findOverlaps(gr_target, x))
        coordm <- lapply(1:length(olm), function(x) ((end(gr_mom[[x]])+start(gr_mom[[x]]))/2)[subjectHits(olm[[x]])])

    if(var$call[1]=="DEL"){col <- "red"}else{col <- "blue"}
    if(XY){pdf(paste0(dir, samp, "_dn_", target_locus, ".pdf"), width=15, height=5)}else{pdf(paste0(dir, samp, "_dn_", target_locus, ".pdf"), width=10, height=5)}
    if(XY){par(mfrow=c(1,3))}else{par(mfrow=c(1,2))}
    title <- paste0(target_locus, "\n", samp, "\n", var$num_exon, " exons")
    matplot(x=coord, y=tab_sel[subjectHits(ol),], ty="l", col=1, ylim=c(0, 5), pch=19, ylab="dCR", main=title, xlab="", xaxt="n", xlim=c(min(coord), max(coord)))
        lines(x=coord, y=tab_sel[subjectHits(ol),var$CNV_ID], col=col, lwd=2)
        sapply(1:length(old), function(x) lines(x=coordd[[x]], y=tab_dad[[x]][subjectHits(old[[x]]),4], col="purple", lwd=2))
        sapply(1:length(olm), function(x) lines(x=coordm[[x]], y=tab_mom[[x]][subjectHits(olm[[x]]),4], col="orange", lwd=2))
        abline(v=coord, lwd=2, col=rgb(0, 0, 0, 0.1))
    axis(side=1, at=coord, labels=prettyNum(coord, big.mark=","))

    grmast <- gr_master; start(grmast) <- start(grmast)+1
        xids <- match(gr_sel[subjectHits(ol)], grmast)
        xids_m <- lapply(1:length(gr_mom), function(x) match(grmast[xids], gr_mom[[x]]))
        xids_p <- lapply(1:length(gr_dad), function(x) match(grmast[xids], gr_dad[[x]]))

    title <- paste0(target_locus, "\n", samp, "\nQS: ", var$QS)
    matplot(y=tab_sel[subjectHits(ol),], x=xids, ty="l", col=1, ylim=c(0, 5), pch=19, ylab="dCR", xaxt="n", main=title, xlab="", xlim=c(min(xids)-.5, max(xids)+.5))
        lines(x=xids, y=tab_sel[subjectHits(ol),var$CNV_ID], col=col, lwd=2)
        sapply(1:length(olm), function(x) lines(x=xids[!is.na(xids_m[[x]])], y=tab_mom[[x]][xids_m[[x]][!is.na(xids_m[[x]])],4], col="orange", lwd=2))
        sapply(1:length(old), function(x) lines(x=xids[!is.na(xids_p[[x]])], y=tab_dad[[x]][xids_p[[x]][!is.na(xids_p[[x]])],4], col="purple", lwd=2))
        abline(v=xids, col=rgb(0, 0, 0, 0.1), lwd=2)
    axis(side=1, at=xids, labels=1:length(xids), las=0)

    xy <- par("usr")
    par(xpd=NA)
    points(x=xy[1], y=mean(tab_sel[subjectHits(ol),var$CNV_ID]), pch=19, cex=2, col=col)
    sapply(1:length(olm), function(x) points(x=xy[1], y=mean(tab_mom[[x]][xids_m[[x]][!is.na(xids_m[[x]])],4]), col="orange", pch=19, cex=2))
    sapply(1:length(old), function(x) points(x=xy[1], y=mean(tab_dad[[x]][xids_p[[x]][!is.na(xids_p[[x]])],4]), col="purple", pch=19, cex=2))
    par(xpd=FALSE)

    st_dad <- sapply(1:length(old), function(x) tab_dad[[x]][xids_p[[x]],4])
        colnames(st_dad) <- paste0("DAD_", str_replace(str_replace(basename(dcr_dad), "denoised_copy_ratios-", ""), ".tsv", ""))
    st_mom <- sapply(1:length(olm), function(x) tab_mom[[x]][xids_m[[x]],4])
        colnames(st_mom) <- paste0("MOM_", str_replace(str_replace(basename(dcr_mom), "denoised_copy_ratios-", ""), ".tsv", ""))

    x <- apply(tab_sel[str_detect(rownames(tab_sel), "X"),], 2, median, na.rm=TRUE)
    y <- apply(tab_sel[str_detect(rownames(tab_sel), "Y"),], 2, median, na.rm=TRUE)

    plot(y~x, main="Sex Chr CN", xlab="chrX CN", ylab="chrY CN", pch=19, col=rgb(0, 0, 0, 0.1), xaxt="n", yaxt="n", xlim=c(0, 4), ylim=c(0,4))
        axis(side=1, at=0:4)
        axis(side=2, at=0:4)

        points(y[var$CNV_ID]~x[var$CNV_ID], pch=19, col=col, cex=4)

        sapply(tab_dad, function(z){y=median(z[str_detect(z[,1], "Y"),4], na.rm=TRUE); x=median(z[str_detect(z[,1], "X"),4], na.rm=TRUE); points(y~x, pch=19, col="purple")})
        sapply(tab_mom, function(z){y=median(z[str_detect(z[,1], "Y"),4], na.rm=TRUE); x=median(z[str_detect(z[,1], "X"),4], na.rm=TRUE); points(y~x, pch=19, col="orange")})

    info <- cbind(tab_sel[subjectHits(ol), 1], st_mom, st_dad, tab_sel[subjectHits(ol), -1])
        colnames(info)[1] <- var$sample
        write.table(info, file=paste0('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/dn_mat/', samp, "_dn_", target_locus, ".txt"), quote=F, sep="\t")

    legend("topright", legend=c("child", "mother", "father"), col=c(col, "orange", "purple"), pch=19)

    file.remove(list.files('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/dcr/', full=TRUE))
    dev.off()
}

regenotypeDN <- function(x){
        if(str_detect(x, "(chrX)|(chrY)")){ cutoff <- 1 }else{ cutoff<-.5 }
        tab <- read.table(x, header=TRUE, stringsAsF=F, sep="\t")
        mads <- apply(tab, 1, mad)
        mads_fail <- sum(mads>=cutoff, na.rm=TRUE)
        means <- apply(tab[mads<cutoff,], 2, mean, na.rm=TRUE)
        md <- min(c(abs(means[1] - means[str_detect(names(means), "MOM_")]), abs(means[1] - means[str_detect(names(means), "DAD_")])))
        return(c(means[1], mean(means[str_detect(names(means), "MOM_")]), mean(means[str_detect(names(means), "DAD_")]), md, mads_fail, length(mads)))
    }

map <- read.csv('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/data/manifest/nygc_sfari_id_sample_map.csv', header=TRUE, stringsAsF=F)
bc <- read.table('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/binCov_loc.tsv', sep="\t", header=TRUE, stringsAsF=F)
    bc_batch <- read.table('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/binCov_loc_batch.tsv', sep="\t", header=TRUE, stringsAsF=F)
    bc_map <- read.table('/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/validation/ref_10065.batches', stringsAsF=F)
    
    meta_dcr <- data.table::fread("/Users/jackfu/Documents/00_Postdoc/19-10-asc_cnv/data/manifest/reformatted_CN.tsv", sep="\t", header=TRUE, stringsAsFactor=FALSE)
        meta_dcr <- data.frame(meta_dcr)
        dcr_index <- lapply(meta_dcr$denoised_copy_ratios, function(x) unlist(str_split(str_replace_all(x, '(\\[)|]|\\"', ""), ",")))
        dcr_index <- do.call(c, dcr_index)
        dcr_index <- data.frame(sample=as.character(str_replace(str_replace(dcr_index, ".*denoised_copy_ratios-", ""), ".tsv", "")), path=dcr_index)









