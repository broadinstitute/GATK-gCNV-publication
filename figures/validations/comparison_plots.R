### gCNV vs XHMM vs CoNIFER
load('comparison_outcome.rda')
par(mfrow=c(1,3))
window <- 1:50
cols <- c("aquamarine4", "darkorchid1", "darkorange1")
plot(stats_sens_gcnv[window,2], ty="l", ylim=c(0, 1), main="Sensitivity", xlab="Min # exons", ylab="% against WGS", col=cols[1], lwd=2)
    lines(stats_sens_xhmm[window,2], col=cols[2], lwd=2)
    lines(stats_sens_conifer[window,2], col=cols[3], lwd=2)
    grid()
plot(stats_sens_gcnv_j[window,2], ty="l", ylim=c(0, 1), main="Sensitivity \n[Strict Intersect]", xlab="Min # exons", ylab="% against WGS", col=cols[1], lwd=2)
    lines(stats_sens_xhmm_j[window,2], col=cols[2], lwd=2)
    lines(stats_sens_conifer_j[window,2], col=cols[3], lwd=2)
    grid()
plot(stats_ppv_gcnv[window,2], ty="l", ylim=c(0, 1), main="PPV", xlab="Min # exons", ylab="% against WGS", col=cols[1], lwd=2)
    lines(stats_ppv_xhmm[window,2], col=cols[2], lwd=2)
    lines(stats_ppv_conifer[window,2], col=cols[3], lwd=2)
    grid()
legend("bottomright", legend=c("gCNV", "XHMM", "CoNIFER"), pch=19, col=cols, cex=0.8)

### gCNV vs WGS
load('comparison.rda')
window <- 50
par(mfrow=c(2,2))
plot(as.numeric(sen_mar_del)~as.numeric(names(sen_mar_del)), ty="b", pch=19, col=2, xlim=c(0, window), ylim=c(0, 1), main="Marginal Sensitivity vs WGS", ylab="Sensitivity", xlab="# exons")
    lines(as.numeric(sen_mar_dup)~as.numeric(names(sen_mar_dup)), ty="b", pch=19, col=4); grid()
plot(as.numeric(ppv_mar_del)~as.numeric(names(ppv_mar_del)), ty="b", pch=19, col=2, xlim=c(0, window), ylim=c(0, 1), main="Marginal PPV vs WGS", ylab="PPV", xlab="# exons")
    lines(as.numeric(ppv_mar_dup)~as.numeric(names(ppv_mar_dup)), ty="b", pch=19, col=4); grid()

plot(y=sen_cum_del, x=1:length(sen_cum_del), ty="b", pch=19, col=2, xlim=c(0, window), ylim=c(0, 1), main="Cumulative Sensitivity vs WGS", ylab="Sensitivity", xlab=">= # exons")
    lines(y=sen_cum_dup, x=1:length(sen_cum_dup), ty="b", pch=19, col=4); grid()
plot(y=ppv_cum_del, x=1:length(ppv_cum_del), ty="b", pch=19, col=2, xlim=c(0, window), ylim=c(0, 1), main="Cumulative PPV vs WGS", ylab="PPV", xlab=">= # exons")
    lines(y=ppv_cum_dup, x=1:length(ppv_cum_dup), ty="b", pch=19, col=4); grid()
