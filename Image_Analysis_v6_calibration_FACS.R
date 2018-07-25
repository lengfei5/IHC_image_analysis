###################
################### Image Analysis starts
###################
source('functions_images.R')
nucleus.all.fitting = TRUE
cell.fitting = FALSE
calibration.facs = TRUE
Method = 'Normal'
Filter.cells = TRUE
nucleus.size.mean = FALSE

WT = TRUE
plot.version = '0.25_WT_refitting_calibration_FACS_correction_1'

if(WT)
{
    dna = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/FACS_DATA.csv', sep=';', header=TRUE)
    N = 32;
    keep = matrix(NA, nrow=N, ncol=18)
    colnames(keep) = c('image.width', 'image.length', 'nb.cell.detected', 'percent.binucleated',
    'cell.size.mono.mean', 'cell.size.binu.mean', 'cell.size.all.mean',
    'nucleus.size.all.2c', 'nucleus.size.all.4c', 'nucleus.size.all.8c',
    'nucleus.percent.all.2c', 'nucleus.percent.all.4c', 'nucleus.percent.all.8c',
    'cell.percent.mono.2c', 'cell.percent.mono.4c','cell.percent.mono.8c',
    'cell.percent.binu.4c', 'cell.percent.binu.8c')

}else{
    dna = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/CRY_KO_FACS_DATA.csv')
    dna = dna[, c(1, 2, 4, 3, 5)]
    N = 4;
    keep = matrix(NA, nrow=N, ncol=18)
    colnames(keep) = c('image.width', 'image.length', 'nb.cell.detected', 'percent.binucleated',
    'cell.size.mono.mean', 'cell.size.binu.mean', 'cell.size.all.mean',
    'nucleus.size.all.2c', 'nucleus.size.all.4c', 'nucleus.size.all.8c',
    'nucleus.percent.all.2c', 'nucleus.percent.all.4c', 'nucleus.percent.all.8c',
    'cell.percent.mono.2c', 'cell.percent.mono.4c','cell.percent.mono.8c',
    'cell.percent.binu.4c', 'cell.percent.binu.8c')
    #keep = c();
}
time = c()

lambda1 = c()
lambda2 = c()
lambda3 = c()

cutoff = rep(450, 64)
ratio.min = 1.2
ratio.max = 2.0

if(WT)
{
    pdfname = paste('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/results/Image_analysis_preliminary_results_', plot.version, '.pdf', sep='')
}else{
    pdfname = paste('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/results/Image_analysis_preliminary_results_', plot.version, '.pdf', sep='')
}


pdf(pdfname, width=16, height=10)

if(WT){
    NN = c(1:N);
}else{
    NN = c(1:N)
    nn = c(9, 33, 21, 45)
}
for(n in NN)
{
    
    tt = (n-1)*3;
    time = c(time, rep(tt,1));
    
    cat('ZT ', tt, '\n');
    par(mfrow=c(2,3))
    
    #n = 7;
    #m = 1;
    test = c()
    
    m = 1;
    if(WT)
    {
        if(n<10) filename = paste('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/results/5th_try_0.25_all/Image_0', n, '_0', m, '.txt', sep='')
        if(n>=10) filename = paste('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/results/5th_try_0.25_all/Image_', n, '_0', m, '.txt', sep='')
    }else{
        filename = paste('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/results/CRYDKO_try_2nd/Image_', nn[n], '_01', '_KO.txt', sep='')
        
    }
    
    aa = read.table(filename, sep='\t', header=FALSE)
    test = c(test, aa[1,2], aa[1, 3])  ### surface of image
    
    aa = aa[-1, ]
    colnames(aa) = c('cell.index', 'cell.size', 'nucleus.cell')
    
    ### remove cells with 3 nuclei
    bb = cbind(aa, find.repeating(aa[,1]))
    bb = bb[which(bb[,4]<3),]
    
    #### filter outliers of cells
    if(Filter.cells)
    {
        source('functions_images.R')
        cat('filter the cells as outliers \n')
        cc = cell.nucleus.filtering(bb)
        cc = data.frame(cc)
        colnames(cc) = c('cell.index', 'cell.size', 'nuclear.size.s', 'nb.nuclei', 'outlier')
        
        index.outlier = cc[which(cc[,5]==1) ,1]
        mm = match(bb[,1], index.outlier)
        mm = which(is.na(mm)==TRUE)
        bb = bb[mm, ]
        cc = cc[-which(cc[,5]==1), ]
    }
    
    #nucleus.size.mean = TRUE
    if(nucleus.size.mean)
    {
        bb = cc
    }
  
    cat(length(unique(bb[,1])), ' cells detected\n')
    
    test = c(test, length(unique(bb[,1])))
    
    #pecent = c(pecent, ((length(which(bb[,2]==2)))/length(which(bb[,2]<=2))))
    #print((nrow(aa)-length(unique(aa[,1])))/length(unique(aa[,1])))
    #pecent = c(pecent, (nrow(bb)-length(unique(bb[,1])))/length(unique(bb[,1])))
    percent = length(which(bb[,4]==2))/2/(length(which(bb[,4]==2))/2+length(which(bb[,4]==1)))
    test = c(test, length(which(bb[,4]==2))/2/(length(which(bb[,4]==2))/2+length(which(bb[,4]==1))))
    #pecent = c(pecent, length(which(bb[,4]==2))/2/(length(which(bb[,4]==1))))
    
    ### cell size analysis: DOES NOT WORK!!! TOO ARBITRARY!!!
    
    jj = which(bb[,4]==1)
    kk = which(bb[,4]==2)
    #par(mfrow=c(1,3))
    #plot(bb[, c(2,3)], cex=0.25, main=paste('ZT_', tt, '_rep_', m,' Cell - Nuclei (All)', sep=''))
    #plot(bb[jj, c(2,3)], cex=0.25, main=paste('ZT_', tt, '_rep_', m,' Cell - Nuclei (Mono)', sep=''))
    #plot(bb[kk, c(2,3)], cex=0.25, main=paste('ZT_', tt, '_rep_', m,' Cell - Nuclei (Binu)', sep=''))
    test = c(test, median(bb[jj,2]))
    index = match(unique(bb[kk,1]), bb[kk, 1])
    test = c(test, median(bb[kk[index],2]))
    test = c(test, median((bb[c(jj, kk[index]),2])))
    
    xlim = c(1000, 10000)
    hist(unique(bb[,2]), breaks=100, xlim=xlim, col='darkorange', freq=FALSE, xlab='piexls', main=paste('ZT_', tt, '_rep_', m,' Cell Size (All)', sep=''))
    abline(v=median(unique(bb[,2])), col='red', lwd=2.0)
    
    hist(bb[jj,2], breaks=100, xlim=xlim,col='gray80',freq=FALSE, xlab='piexls', main=paste('ZT_', tt, '_rep_', m, ' Celll Size (Mono)', sep=''))
    abline(v=median(bb[jj,2]), col='red')
    
    hist(unique(bb[kk,2]), breaks=50, xlim=xlim, freq=FALSE, col='gray50', xlab='piexls', main=paste('ZT_', tt, '_rep_', m,' Cell Size (Bi)', sep=''))
    abline(v=median(unique(bb[kk,2])), col='red')
    
    #########################
    ### nuclear size analysis
    #########################
    ## all nucleui
    x = bb[,3];
    mu.init = c(sample(seq(300, 400, by=10), 1), sample(seq(500, 600, by=10), 1), sample(seq(700, 800, by=10), 1))
    
    yy = x[which(x>250)]
    res = size.seperation(yy, method=Method, k=3, mu.init=mu.init)
    lambda = res[1:3]
    par1 = res[4:6]
    par2 = res[7:9]
    #ratio = c(par1[2]/par1[1], par1[3]/par1[2])
    refs = dna[c(1, 2, 4), (n+1)]
    alpha = refs/sum(refs)
    res = normalmixEM.test(x, k=3, mu.init=par1, sigma.init=par2, lambda.constr=alpha)
    par1 = res$mu
    
    if(n==31){
        cutoff.size=700
    }else{
        cutoff.size=600
    }
    nitr = 0
    while(par1[3]<cutoff.size)
    {
        mu.init = c(sample(seq(300, 400, by=10), 1), sample(seq(500, 600, by=10), 1), sample(seq(700, 800, by=10), 1));
        res = normalmixEM.test(x, k=3, mu.init=mu.init, lambda.constr=alpha)
        par1 = res$mu
        cat(par1, '\n');
        nitr = nitr+1;
        if(nitr>=20) break;
    }
    
    res1 = c(res$lambda, res$mu, res$sigma)
    lambda = res1[1:3]
    par1 = res1[4:6]
    par2 = res1[7:9]
    par10 = par1
    par20 = par2
    
    lambda.a = lambda
    
    test = c(test, par1)
    #test = c(test, lambda)
    
    xlim = c(100, 1000)
    cut1 = 10;
    cut2 = 20;
    
    xfit<-seq(min(x),max(x),length=100)
    if(Method=='Normal')
    {
        yfit1<-dnorm(xfit, par1[1],par2[1])*lambda[1]
        yfit2<-dnorm(xfit, par1[2],par2[2])*lambda[2]
        yfit3<-dnorm(xfit, par1[3],par2[3])*lambda[3]
    }
    
    hist(bb[,3], breaks=100, xlim = xlim, xlab='piexls', freq=FALSE, col='yellow', main=paste('ZT_', tt, '_rep_', m,' Nucleus Size (All) ', sep=''))
    lines(xfit, yfit1, col='green', lwd=2.)
    lines(xfit, yfit2, col='green', lwd=2.)
    lines(xfit, yfit3, col='green', lwd=2.)
    
    jj = which(bb[,4]==1)
    kk = which(bb[,4]==2)
    
    ## nuclei within mononucleated cells
    x = bb[jj,3]
    res2 = size.seperation(x, method=Method, k=3, m.constr=par10, sigma.constr=par20)
    #res = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20)
    lambda = res2[1:3]
    par1 = res2[4:6]
    par2 = res2[7:9]
    
    lambda.m = lambda
    #test = c(test, (1-percent)*lambda)
    
    xfit<-seq(min(x),max(x),length=100)
    yfit1<-dnorm(xfit, par1[1],par2[1])*lambda[1]
    yfit2<-dnorm(xfit, par1[2],par2[2])*lambda[2]
    yfit3<-dnorm(xfit, par1[3],par2[3])*lambda[3]
    
    hist(x, breaks=100, xlim = xlim, xlab='piexls', freq=FALSE, col='blue', main=paste('ZT_', tt, '_rep_', m, '  Nucleus Size (mononuclate) ', sep=''))
    lines(xfit, yfit1, col='green', lwd=2.0)
    lines(xfit, yfit2, col='green', lwd=2.)
    lines(xfit, yfit3, col='green', lwd=2.)
    
    ## nuclei within binucleated cells
    x = bb[kk,3]
    #res3 = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20, m.constr=par10, sigma.constr=par20)
    #res3 = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20, m.constr=par10)
    res3 = size.seperation(x, method=Method, k=3, m.constr=par10, sigma.constr=par20)
    
    #res5 = res1
    lambda = res3[1:3]
    par1 = res3[4:6]
    par2 = res3[7:9]
    lambda.b = lambda
    
    xfit<-seq(min(x),max(x),length=100)
    yfit1<-dnorm(xfit, par1[1],par2[1])*lambda[1]
    yfit2<-dnorm(xfit, par1[2],par2[2])*lambda[2]
    yfit3<-dnorm(xfit, par1[3],par2[3])*lambda[3]
    
    #test = c(test, percent*lambda[c(1:2)])
    xlim = c(100, 1000)
    hist(x, breaks=50, xlim = xlim, xlab='piexls', freq=FALSE, col='darkblue', main=paste('ZT_', tt, '_rep_', m,' Nucleus Size (binucleate) ', sep=''))
    lines(xfit, yfit1, col='green', lwd=2.0)
    lines(xfit, yfit2, col='green', lwd=2.)
    lines(xfit, yfit3, col='green', lwd=2.)
    
    lambda.mm = lambda.m
    lambda.bb = lambda.b
    
    lambda2 = rbind(lambda2, lambda.mm)
    lambda3 = rbind(lambda3, lambda.bb)
    
    #cat(dna[c(1,2,4), (n+1)]/100, '\n')
    cat(lambda.a, '\n')
    lambda.m = lambda.mm*(1-percent);
    lambda.b = lambda.bb*(percent);
    cat(c((lambda.m[1]+2*lambda.b[1])/(1+percent), (lambda.m[2]+2*lambda.b[2])/(1+percent), (lambda.m[3]+2*lambda.b[3])/(1+percent)), '\n')

    ### save table
    test = c(test, c((lambda.m[1]+2*lambda.b[1])/(1+percent), (lambda.m[2]+2*lambda.b[2])/(1+percent), (lambda.m[3]+2*lambda.b[3])/(1+percent)))
    test = c(test, lambda.m)
    test = c(test, lambda.b[1:2])
    
    keep[n, ] = test
    
}

dev.off()

keep = data.frame(keep, stringsAsFactors=FALSE)
#keep = matrix(NA, nrow=N, ncol=18)
colnames(keep) = c('image.width', 'image.length', 'nb.cell.detected', 'percent.binucleated',
'cell.size.mono.mean', 'cell.size.binu.mean', 'cell.size.all.mean',
'nucleus.size.all.2c', 'nucleus.size.all.4c', 'nucleus.size.all.8c',
'nucleus.percent.all.2c', 'nucleus.percent.all.4c', 'nucleus.percent.all.8c',
'cell.percent.mono.2c', 'cell.percent.mono.4c','cell.percent.mono.8c',
'cell.percent.binu.4c', 'cell.percent.binu.8c')

#kkeep = keep
#load(file='Rdata/Results_images_WT_vf_calibrated_FACS_outliers_1_22_29_32.Rdata')
#keep$cell.size.binu.mean = kkeep$cell.size.binu.mean
#keep$cell.size.all.mean = kkeep$cell.size.all.mean
#save(keep, file='Rdata/Results_images_WT_vf_calibrated_FACS_outliers_1_22_29_32.Rdata')
#save(keep, file='Rdata/Results_images_WT_vf_calibrated_FACS_outliers_1_22_29_32.Rdata')
#keep.ko = keep
#filename = paste('Rdata/Results_images_KO_calibrated_FACS.Rdata', sep='')
#save(keep.ko, file='Rdata/Results_images_KO_calibrated_FACS.Rdata')

Correction = FALSE
if(Correction)
{
    load('Rdata/Results_images_WT_v4_remove_outliers_0.25_3_1.Rdata')
    keep0 = keep
    
    
    WT = TRUE
    dna = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/FACS_DATA.csv', sep=';', header=TRUE)
    nucleus.all.fitting = TRUE
    cell.fitting = FALSE
    calibration.facs = TRUE
    WT = TRUE
    Method = 'Normal'
    
    Filter.cells = TRUE
    nucleus.size.mean = FALSE
   
    
    correlation = FALSE
    
    if(correlation)
    {
        n = 7
        nb.components.mono = 2
        nb.components.bino = 2
        Sigma.Constr = FALSE
        source('functions_images.R')
        reduce.small.size.population = FALSE
        
        fit = nuclear.size.fitting.manual(n, nb.components.mono, nb.components.bino, Sigma.Constr, reduce.small.size.population, keep)
        
        ######
        ### test any improvement
        ######
        PLOT_test = FALSE
        if(PLOT_test)
        {
            keep0[n, ] = fit
            
            par(mfrow=c(1,1))
            dna = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/FACS_DATA.csv', sep=';', header=TRUE)
            ref2 = dna[1, c(2:33)]
            #ref2 = cbind(t(ref2[1:8]), t(ref2[9:16]), t(ref2[17:24]), t(ref2[25:32]))
            #ref2 = cbind(t(ref2[1:8]), t(ref2[9:16]))
            #ref2.mean = apply(ref2, 1, mean)
            #ref2.err = apply(ref2, 1, sd)/sqrt(4)
            ref4 = dna[2, c(2:33)]
            ref8 = dna[4, c(2:33)]
            #ref4 = cbind(t(ref4[1:8]), t(ref4[9:16]), t(ref4[17:24]), t(ref4[25:32]))
            #ref4 = cbind(t(ref4[1:8]), t(ref4[9:16]))
            #ref4.mean = apply(ref4, 1, mean)
            #ref4.err = apply(ref4, 1, sd)/sqrt(4)
            
            nuclear.f11 = keep0[,11]*100
            nuclear.f22 = keep0[,12]*100
            nuclear.f33 = keep0[,13]*100
            #nuclear.f11 = (nuclear.f11[1:8]+nuclear.f11[9:16]+nuclear.f11[17:24]+nuclear.f11[25:32])/4
            #nuclear.f22 = (nuclear.f22[1:8]+nuclear.f22[9:16]+nuclear.f22[17:24]+nuclear.f22[25:32])/4
            time = c(0:31)*3
            plot(time, ref2, type='b', col='green', lwd=2.0, ylim=range(c(ref2, ref4, nuclear.f11, nuclear.f22, ref8, nuclear.f33)), ylab='Percentages of DNA contents', main='Percentages of DNA contents', ,xlab='time (ZT)')
            points(time, ref4, type='b', col='orange', lwd=2.0)
            points(time, ref8, type='b', col='violet', lwd=2.0)
            points(time, nuclear.f11, type='p', col='black', cex=1.2, pch=16)
            text(time, (nuclear.f11+3), paste(c(1:32)), col='green')
            points(time, nuclear.f22, type='p', col='black', cex=1.2, pch=14)
            text(time, (nuclear.f22+3), paste(c(1:32)), col='orange')
            points(time, (nuclear.f33), type='p', col='black', cex=1.2, pch=13)
            text(time, (nuclear.f33+3), paste(c(1:32)), col='violet')
            #points(time, nuclear.f3, type='b', col='gray', pch=1)
            abline(v=24, col='gray', lwd=1.5)
            abline(v=48, col='gray', lwd=1.5)
            abline(v=72, col='gray', lwd=1.5)
            abline(v=96, col='gray', lwd=1.5)
        }
        
    }
    
    2correct = FALSE
    if(2correct)
    {
        keep = keep0
        save(keep, dna, lambda2, lambda3, file='Rdata/Results_images_WT_v4_remove_outliers_correction_n_7_13_15_16_17_19_22_29_31_32_save.Rdata')
    }
}

time = c(0:31)*3
ref = dna[c(1,2,4), c(2:33)]/100
ref = ref/apply(ref, 2, sum)
matplot(time, t(ref), type='b', col='blue')
#matplot(time, t(lambda1), type='p', add=TRUE, pch=c(2,3,4), col='black')
#matplot(time, t(lambda2[, (c(1:32)*2-1)]), type='p', add=TRUE, pch=c(2,3,4), col='black')
#matplot(time, t(lambda2[, (c(1:32)*2)]), type='p', add=TRUE, pch=c(2,3,4), col='green')
matplot(time, keep[,c(11:13)], type='p', add=TRUE, pch=c(2,3,4), col='red')


Test.EM = FALSE
if(Test.EM)
{
    data(faithful)
    attach(faithful)
    set.seed(100)
    out<-normalmixEM(waiting, arbvar = FALSE, epsilon = 1e-03)
    #plot(out)
    
    N <- 10000
    alpha = c(0.6,0.3,0.1)
    components <- sample(1:3,prob=alpha,size=N,replace=TRUE)
    mus <- c(54,80,120)
    sds <- c(6,10,4)
    
    data <- rnorm(n=N,mean=mus[components],sd=sds[components])
    
    source('functions_images.R')
    fit = normalmixEM.test(data, k=3, lambda.constr=alpha)
}



pdfname = paste('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/results/Image_analysis_preliminary_results_summary_', plot.version, '.pdf', sep='')
pdf(pdfname, width=18, height=8)

par(mfrow=c(1,2))
## mean cell size

plot(time, keep$cell.size.all.mean,  xlab='time (ZT)', ylab='cell size (pixel)', lwd=1.5, ylim=range(c(keep$cell.size.all.mean, keep$cell.size.mono.mean, keep$cell.size.binu.mean)), type='b', col='darkblue', main='Cell Size (median))')
points(time, keep$cell.size.mono.mean, type='b', lwd=1.5, col='darkgreen',pch=1)
points(time, keep$cell.size.binu.mean, type='b', lwd=1.5, col='darkorange', pch=1)
legend('topright', legend = c('Mono','Binu','All'), pch=c(1),lty=c(1,1,1,1), cex=1.2,col = c('darkblue','darkgreen', 'darkorange'), border = NA, bty = 'n')
abline(v=24, col='gray', lwd=1.5)
abline(v=48, col='gray', lwd=1.5)
abline(v=72, col='gray', lwd=1.5)
abline(v=96, col='gray', lwd=1.5)

## detected cell numbers
ss = keep$image.width*keep$image.length
ss = ss/ss[1]
plot(time, keep$nb.cell.detected/ss, type='b', xlab='time (ZT)', lwd=2.0, col='black',ylab='Nb of cells detected', main='Nb of cells in a fixed surface')
#points(time, nb.cell, type='b', col='black')
abline(v=24, col='gray', lwd=1.5)
abline(v=48, col='gray', lwd=1.5)
abline(v=72, col='gray', lwd=1.5)
abline(v=96, col='gray', lwd=1.5)

## percentages of mono- and binucleated cells
plot(time, keep$percent.binucleated, type='b',lwd=2.0, col='darkred', xlab='ZT', main='percentage of binucleate cells')
abline(v=24, col='gray', lwd=1.5)
abline(v=48, col='gray', lwd=1.5)
abline(v=72, col='gray', lwd=1.5)
abline(v=96, col='gray', lwd=1.5)

## control parts: compare percentages of 2C and 4C nuclei with FACS
dna = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/FACS_DATA.csv', sep=';', header=TRUE)
ref2 = dna[1, c(2:33)]
#ref2 = cbind(t(ref2[1:8]), t(ref2[9:16]), t(ref2[17:24]), t(ref2[25:32]))
#ref2 = cbind(t(ref2[1:8]), t(ref2[9:16]))
#ref2.mean = apply(ref2, 1, mean)
#ref2.err = apply(ref2, 1, sd)/sqrt(4)
ref4 = dna[2, c(2:33)]
ref8 = dna[4, c(2:33)]
#ref4 = cbind(t(ref4[1:8]), t(ref4[9:16]), t(ref4[17:24]), t(ref4[25:32]))
#ref4 = cbind(t(ref4[1:8]), t(ref4[9:16]))
#ref4.mean = apply(ref4, 1, mean)
#ref4.err = apply(ref4, 1, sd)/sqrt(4)

nuclear.f11 = keep[,11]*100
nuclear.f22 = keep[,12]*100
nuclear.f33 = keep[,13]*100
#nuclear.f11 = (nuclear.f11[1:8]+nuclear.f11[9:16]+nuclear.f11[17:24]+nuclear.f11[25:32])/4
#nuclear.f22 = (nuclear.f22[1:8]+nuclear.f22[9:16]+nuclear.f22[17:24]+nuclear.f22[25:32])/4
time = c(0:31)*3
plot(time, ref2, type='b', col='green', lwd=2.0, ylim=range(c(ref2, ref4, nuclear.f11, nuclear.f22, ref8, nuclear.f33)), ylab='Percentages of DNA contents', main='Percentages of DNA contents', ,xlab='time (ZT)')
points(time, ref4, type='b', col='orange', lwd=2.0)
points(time, ref8, type='b', col='violet', lwd=2.0)
points(time, nuclear.f11, type='p', col='black', cex=1.2, pch=16)
text(time, (nuclear.f11+3), paste(c(1:32)), col='green')
points(time, nuclear.f22, type='p', col='black', cex=1.2, pch=14)
text(time, (nuclear.f22+3), paste(c(1:32)), col='orange')
points(time, (nuclear.f33), type='p', col='black', cex=1.2, pch=13)
text(time, (nuclear.f33+3), paste(c(1:32)), col='violet')
#points(time, nuclear.f3, type='b', col='gray', pch=1)
abline(v=24, col='gray', lwd=1.5)
abline(v=48, col='gray', lwd=1.5)
abline(v=72, col='gray', lwd=1.5)
abline(v=96, col='gray', lwd=1.5)

one.day.dna.content= TRUE
if(one.day.dna.content)
{
    source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')
    #remv = c(1, 7, 15,16, 17, 22,26, 29, 31, 32)
    remv = c(1, 26)
    kkeep = keep
    kkeep[remv, ] = NA
    source('functions_nuclear.R')
    ref2 = dna[1, c(2:33)]/100
    ref4 = dna[2, c(2:33)]/100
    ref8 = dna[4, c(2:33)]/100
    averg0 = mean.err(ref2)[1,]
    err0 = mean.err(ref2)[2,]
    averg1 = mean.err(ref4)[1,]
    err1 = mean.err(ref4)[2,]
    averg2 = mean.err(ref8)[1,]
    err2 = mean.err(ref8)[2,]
    lims = range(c(averg0+err0, averg0-err0, averg1+err1, averg1-err1, averg2+err2, averg2-err2));
    pch=1
    ttime = c(0:7)*3
    tt = seq(0, 24, length.out=100)
    lims = c(0., 0.7)
    col = 'darkblue'
    plot(ttime, averg0, type='p', ylim=lims, xlim=c(0,21), col=col, lwd=1.0,  lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA, cex=0.5)
    arrows(ttime, averg0-err0, ttime, averg0+err0, length=0.05, angle=90, code=3, col=col, lty=1, pch=pch,lwd=1.0)
    res = f24_R2_alt2(as.numeric(ref2), t=c(0:31)*3)
    yy = res[2]*(1+res[4]*cos(2*pi/24*(tt-res[5])))
    points(tt, yy, type='l', col=col, lwd=1.)
    col='orange'
    points(ttime, averg1, type='p', col=col, lwd=1., lty=1, pch=pch, cex=0.5)
    arrows(ttime, averg1-err1, ttime, averg1+err1, length=0.05, angle=90, code=3,  lty=1, col=col, pch=pch,lwd=1.)
    res = f24_R2_alt2(as.numeric(ref4), t=c(0:31)*3)
    yy = res[2]*(1+res[4]*cos(2*pi/24*(tt-res[5])))
    points(tt, yy, type='l', col=col, lwd=1.)
    
    
    col='darkgray'
    points(ttime, averg2, type='p', col=col, lwd=1.,  lty=1, pch=pch, cex=0.5)
    arrows(ttime, averg2-err2, ttime, averg2+err2, length=0.05, angle=90, code=3,  lty=1,col=col, pch=pch,lwd=1.)
    res = f24_R2_alt2(as.numeric(ref8), t=c(0:31)*3)
    yy = res[2]*(1+res[4]*cos(2*pi/24*(tt-res[5])))
    points(tt, yy, type='l', col=col, lwd=1.0)
    
    axis(1,at=3*c(0:7),cex.axis =0.8)
    #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
    lims = c(signif(lims[1], d=1), signif(lims[2], d=2))
    axis(2, at= seq(0.3,0.7, 0.1), las=1,cex.axis = 0.8)
    box()
    
    points(ttime, mean.err(kkeep[,11])[1,], type='p', pch=1, col='red')
    points(ttime, mean.err(kkeep[,12])[1,], type='p', pch=2, col='red')
    points(ttime, mean.err(kkeep[,13])[1,], type='p', pch=3, col='red')
    
    
    ### different cell population
    averg0 = mean.err(kkeep[,14])[1,]
    err0 = mean.err(kkeep[,14])[2,]
    averg1 = mean.err(kkeep[,15])[1,]
    err1 = mean.err(kkeep[,15])[2,]
    averg2 = mean.err(kkeep[,17])[1,]
    err2 = mean.err(kkeep[,17])[2,]
    
    #lims = range(c(averg0+err0, averg0-err0, averg1+err1, averg1-err1, averg2+err2, averg2-err2));
    #lims = range(c(averg0+err0, averg0-err0, averg1+err1, averg1-err1));
    pch=1
    ttime = c(0:7)*3
    tt = seq(0, 24, length.out=100)
    #lims = c(0., 0.7)
    col = 'darkblue'
    plot(ttime, averg0, type='p', ylim=lims, xlim=c(0,21), col=col, lwd=1.0,  lty=1, pch=pch, main=NA, axes=FALSE, xlab=NA, ylab=NA, cex=0.5)
    arrows(ttime, averg0-err0, ttime, averg0+err0, length=0.05, angle=90, code=3, col=col, lty=1, pch=pch,lwd=1.0)
    #res = f24_R2_alt2(as.numeric(ref2), t=c(0:31)*3)
    #yy = res[2]*(1+res[4]*cos(2*pi/24*(tt-res[5])))
    #points(tt, yy, type='l', col=col, lwd=1.)
    col='orange'
    points(ttime, averg1, type='p', col=col, lwd=1., lty=1, pch=pch, cex=0.5)
    arrows(ttime, averg1-err1, ttime, averg1+err1, length=0.05, angle=90, code=3,  lty=1, col=col, pch=pch,lwd=1.)
    #res = f24_R2_alt2(as.numeric(ref4), t=c(0:31)*3)
    #yy = res[2]*(1+res[4]*cos(2*pi/24*(tt-res[5])))
    #points(tt, yy, type='l', col=col, lwd=1.)
    
    col='darkgray'
    points(ttime, averg2, type='p', col=col, lwd=1.,  lty=1, pch=pch, cex=0.5)
    arrows(ttime, averg2-err2, ttime, averg2+err2, length=0.05, angle=90, code=3,  lty=1,col=col, pch=pch,lwd=1.)

    axis(1,at=3*c(0:7),cex.axis =0.8)
    #axis(1,at=24,'ZT[hr]',tick=FALSE,cex.axis =1.0)
    lims = c(signif(lims[1], d=1), signif(lims[2], d=2))
    axis(2, at= seq(0.3,0.7, 0.1), las=1,cex.axis = 0.8)
    box()

}
#require('plotrix')
#plotCI(time, averg, err, scol="darkblue",lwd=1.5, pch=16
#plotCI(c(0:7)*3, ref2.mean, ref2.err, lwd=2.0, scol='green', col='green', xlab='ZT', ylim=range(c(nuclear.f11, ref2)), main='percentages of 2C nuclei', ylab='DNA contents per nucleus')
#points(c(0:7)*3, ref2.mean, lwd=2.0, col='green', type='l')
#points(c(0:7)*3, nuclear.f11, type='p', col='black', cex=1.2, pch=16)

#### percentages of nucleus in binu and mono cells
matplot(time, lambda2, type='b', main='Percentages of nucleus in mono cells', log='')
abline(v=24, col='gray', lwd=1.5)
abline(v=48, col='gray', lwd=1.5)
abline(v=72, col='gray', lwd=1.5)
abline(v=96, col='gray', lwd=1.5)

matplot(time, lambda3, type='b', main='Percentages of nucleus in bino cells', log='')
abline(v=24, col='gray', lwd=1.5)
abline(v=48, col='gray', lwd=1.5)
abline(v=72, col='gray', lwd=1.5)
abline(v=96, col='gray', lwd=1.5)

### nucleus size for different DNA centents
mm = grep('nucleus.size', colnames(keep))
matplot(time, keep[,mm], type='b',  main='Nuclear sizes')
abline(v=24, col='gray', lwd=1.5)
abline(v=48, col='gray', lwd=1.5)
abline(v=72, col='gray', lwd=1.5)
abline(v=96, col='gray', lwd=1.5)
abline(h=650, col='black', lwd=2.0)
abline(h=400, col='black', lwd=2.0)

### percentages of cells with polyploidy
xx = rbind( c(1, 1, 1, 1, 1),
            c(1, 1, 1, 0, 0),
            c(1, 0, 0, 2, 0),
            c(0, 1, 0, 0, 2),
            c(0, 0, 1, 0, 0)
            )
n = 1
s = keep$percent.binucleated[n]
x = (1+s)*as.numeric(dna[1, (n+1)])/100
y = (1+s)*as.numeric(dna[2, (n+1)])/100
z = (1+s)*(1-as.numeric(dna[1, (n+1)])/100-as.numeric(dna[2, (n+1)])/100)
yy = c(1, 1-s, x, y, z)
#solve(xx, yy)

mm = grep('cell.percent', colnames(keep))
matplot(time, keep[, mm], type='b', main='Percentages of cells in liver', log='')
abline(v=24, col='gray', lwd=1.5)
abline(v=48, col='gray', lwd=1.5)
abline(v=72, col='gray', lwd=1.5)
abline(v=96, col='gray', lwd=1.5)


source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')

f24_R2_alt2(as.numeric(dna[c(1), c(2:33)]/100), t=c(0:31)*3)
f24_R2_alt2(as.numeric(dna[c(2), c(2:33)]/100), t=c(0:31)*3)
f24_R2_alt2(as.numeric(dna[c(4), c(2:33)]/100), t=c(0:31)*3)

f24_R2_alt2(as.numeric(keep[,11]), t=c(0:31)*3)
f24_R2_alt2(as.numeric(keep[,12]), t=c(0:31)*3)
f24_R2_alt2(as.numeric(keep[,13]), t=c(0:31)*3)

f24_R2_alt2(keep[,14], t=c(0:31)*3)
f24_R2_alt2(keep[,15], t=c(0:31)*3)
f24_R2_alt2(keep[,17], t=c(0:31)*3)
f24_R2_alt2(keep[,16], t=c(0:31)*3)
f24_R2_alt2(keep[,18], t=c(0:31)*3)
#f24_R2_alt2(keep[,23], t=c(0:31)*3)

dev.off()
