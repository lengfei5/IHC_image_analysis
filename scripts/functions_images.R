## Analyze and plot the results of image processsing
library(mixtools)
############ functions used in this file
find.repeating = function(x=c(1:10))
{
    bb = x;
    rept  = c();
    
    for(b in bb)
    {
        rept = c(rept, length(which(x==b)))
    }
    return(rept)
}

cell.nucleus.processing = function(bb)
{
    index = unique(bb[,1])
    cc = matrix(NA, nrow = length(index), ncol = 4)
    for(n in 1:length(index))
    {
        kk = which(bb[, 1]==index[n])
        cc[n, ] = c(n, mean(bb[kk,2]), mean(bb[kk, 3]), length(kk))
    }
    
    return(cc)
}

cell.nucleus.filtering = function(bb)
{
    index = unique(bb[,1])
    cc = matrix(NA, nrow = length(index), ncol = 4)
    for(n in 1:length(index))
    {
        kk = which(bb[, 1]==index[n])
        cc[n, ] = c(index[n], mean(bb[kk,2]), sum(bb[kk, 3]), length(kk))
    }
    
    nb.bins = 20
    bins = seq(range(cc[,3])[1], range(cc[,3])[2], length.out=(nb.bins+1))
    outliers = rep(NA, nrow(cc))
    
    for(n in 1:nb.bins)
    {
        if(n ==1)
        {
            kk = which(cc[,3]>=bins[n] & cc[,3]<=bins[n+1])
        }else{
            kk = which(cc[,3]>bins[n] & cc[,3]<=bins[n+1])
        }
        
        if(length(kk)>50)
        {
            #cat(length(kk), '\n')
            test = cc[kk, 2]
            index = index.outliers(test)
            outliers[kk[index]] = 1
        }
        
    }
    
    cc = cbind(cc, outliers=outliers)
    #colnames(cc) = c('cell.index', 'cell.size', 'nuclear.size.s', 'outlier')
    
    plot.test =FALSE
    
    if(plot.test)
    {
        ii = which(cc[,4]==1)
        jj = which(cc[,4]==2)
        plot(cc[ii, 2:3], cex=0.2, ylim=range(cc[,3]))
        points(cc[jj, 2:3], cex=0.2, col='blue')
        kk = which(cc[,5]==1)
        points(cc[kk, 2:3], cex=0.2, col='red')
    }
    
    return(cc)
}

find.cutoff.size = function(x)
{
    # x = bb[,3]
    test = hist(x, plot=FALSE)
    ll = length(test$counts)
    derive.forward = test$counts[2:(ll-1)]- test$counts[1:(ll-2)]
    derive.backward = test$counts[2:(ll-1)]- test$counts[3:(ll)]
    peaks.index = which(derive.forward>=0 & derive.backward>=0)+1
    peaks = test$counts[peaks.index]
    kk = order(-peaks)
    peaks = peaks[kk]
    peaks.index = peaks.index[kk]
    peaks = peaks[c(1:2)]
    peaks.index = peaks.index[1:2]
    
    index = which(derive.forward<=0 & derive.backward<=0)+1
    index = index[which(index>min(peaks.index) & index<max(peaks.index))]
    if(length(index)==1)
    cut = 0.5*(test$breaks[index] + test$breaks[(index+1)])
    
    cut
    
}

size.seperation = function(x, method='Normal',k=3, mu.init=NULL, sigma.init=NULL, lambda.init=c(0.595,0.385,0.02), m.constr=NULL, sigma.constr=NULL)
{
    if(method=='Normal')
    {
        if(k==3)
        {
            fit <- normalmixEM(x, lambda = lambda.init, mu = mu.init, sigma = sigma.init, mean.constr=m.constr, sd.constr=sigma.constr, k=3, maxrestarts=20, maxit = 1500)
            #plot(fit, density=TRUE)
            lambda = fit$lambda;
            par1 = fit$mu
            par2 = fit$sigma
            #o1 = order(par1)
            #lambda = lambda[o1]
            #par1 = par1[o1]
            #par2 = par2[o1]
            loglike = fit$loglik
            return(c(lambda, par1, par2))
            
        }
        
        if(k==4)
        {
            
            fit <- normalmixEM(x, lambda = lambda.init, mu = mu.init, sigma = sigma.init, k=4, maxrestarts=20, maxit = 1500)
            #plot(fit, density=TRUE)
            lambda = fit$lambda;
            par1 = fit$mu
            par2 = fit$sigma
            o1 = order(par1)
            lambda = lambda[o1]
            par1 = par1[o1]
            par2 = par2[o1]
            loglike = fit$loglik
            
            #prob = data.frame(lambda[1]*dnorm(x, par1[1], par2[1], log=FALSE), lambda[2]*dnorm(x, par1[2], par2[2], log=FALSE), lambda[3]*dnorm(x, par1[3], par2[3], log=FALSE))
            #fam = apply(prob, 1, whichmax)
            #return(c(lambda1=lambda[1], lambda2=lambda[2], lambda3=lambda[3], par1, par2, family1=length(which(fam==1)), family2=length(which(fam==2)), family3=length(which(fam==3)),  cutoff1 = 0.5*(max(x[which(fam==1)])+min(x[which(fam==2)])), cutoff2 = 0.5*(max(x[which(fam==2)])+min(x[which(fam==3)]))))
            #return(c(lambda, par1, par2, cutoff1 = 0.5*(max(x[which(fam==1)])+min(x[which(fam==2)])), cutoff2 = 0.5*(max(x[which(fam==2)])+min(x[which(fam==3)]))) )
            
            return(c(lambda, par1, par2))
        }

        if(k==2)
        {
            fit <- normalmixEM(x, lambda = c(0.6, 0.4), mu = mu.init, sigma = sigma.init, mean.constr=m.constr, sd.constr=sigma.constr, k=2, maxrestarts=20, maxit = 2000)
            #plot(fit, density=TRUE)
            lambda = fit$lambda;
            par1 = fit$mu
            par2 = fit$sigma
            o1 = order(par1)
            lambda = lambda[o1]
            par1 = par1[o1]
            par2 = par2[o1]
            loglike = fit$loglik
            
            #prob = data.frame(lambda[1]*dnorm(x, par1[1], par2[1], log=FALSE), lambda[2]*dnorm(x, par1[2], par2[2], log=FALSE), lambda[3]*dnorm(x, par1[3], par2[3], log=FALSE))
            #fam = apply(prob, 1, whichmax)
            #return(c(lambda1=lambda[1], lambda2=lambda[2], lambda3=lambda[3], par1, par2, family1=length(which(fam==1)), family2=length(which(fam==2)), family3=length(which(fam==3)),  cutoff1 = 0.5*(max(x[which(fam==1)])+min(x[which(fam==2)])), cutoff2 = 0.5*(max(x[which(fam==2)])+min(x[which(fam==3)]))))
            return(c(lambda, par1, par2))
            
        }
    }
    
    if(method=='Gamma')
    {
        fit <- gammamixEM(x, lambda = c(0.5, 0.4, 0.1), alpha = NULL, beta = NULL, k=3)
        #plot(fit, density=TRUE)
        lambda = fit$lambda;
        par1 = fit$gamma.pars[1,]
        par2 = fit$gamma.pars[2,]
        o1 = order(par1*par2)
        lambda = lambda[o1]
        par1 = par1[o1]
        par2 = par2[o1]
        
        PLOT = FALSE
        if(PLOT)
        {
            xfit<-seq(min(x),max(x),length=100)
            yfit1<-dgamma(xfit, shape=par1[1],scale=par2[1])*lambda[1]
            yfit2<-dgamma(xfit, shape=par1[2],scale=par2[2])*lambda[2]
            yfit3<-dgamma(xfit,shape=par1[3],scale=par2[3])*lambda[3]
            hist(x, breaks=40, xlim = xlim, xlab='piexls', freq=FALSE, col='blue', main=paste('ZT_', tt, '_rep_', m, '  Nucleus Size (mononuclate) ', sep=''))
            #abline(v=cut1, col='red',lwd=1.5)
            #abline(v=cut2, col='red',lwd=1.5)
            lines(xfit, yfit1, col='green', lwd=2.0)
            lines(xfit, yfit2, col='green', lwd=2.)
            lines(xfit, yfit3, col='green', lwd=2.)
            
        }
        
        #prob = data.frame(lambda[1]*dgamma(x, shape=par1[1],scale=par2[1], log=FALSE), lambda[2]*dgamma(x, shape=par1[2],scale=par2[2], log=FALSE), lambda[3]*dgamma(x, shape=par1[3],scale=par2[3], log=FALSE))
        #fam = apply(prob, 1, whichmax)
        #return(c(lambda1=lambda[1], lambda2=lambda[2], lambda3=lambda[3], par1, par2, family1=length(which(fam==1)), family2=length(which(fam==2)), family3=length(which(fam==3)),  cutoff1 = 0.5*(max(x[which(fam==1)])+min(x[which(fam==2)])), cutoff2 = 0.5*(max(x[which(fam==2)])+min(x[which(fam==3)]))))
        return(c(lambda, par1, par2))
    }
    
    if(method=='NP')
    {
        if(k==3)
        {
            fit <- npEM(x, mu0=3, verb=FALSE)
            #plot(fit, density=TRUE)
            lambda = fit$lambda;
            par1 = fit$mu
            par2 = fit$sigma
            o1 = order(par1)
            lambda = lambda[o1]
            par1 = par1[o1]
            par2 = par2[o1]
            loglike = fit$loglik
            
            #prob = data.frame(lambda[1]*dnorm(x, par1[1], par2[1], log=FALSE), lambda[2]*dnorm(x, par1[2], par2[2], log=FALSE), lambda[3]*dnorm(x, par1[3], par2[3], log=FALSE))
            
            #fam = apply(prob, 1, whichmax)
            #return(c(lambda1=lambda[1], lambda2=lambda[2], lambda3=lambda[3], par1, par2, family1=length(which(fam==1)), family2=length(which(fam==2)), family3=length(which(fam==3)),  cutoff1 = 0.5*(max(x[which(fam==1)])+min(x[which(fam==2)])), cutoff2 = 0.5*(max(x[which(fam==2)])+min(x[which(fam==3)]))))
            return(c(lambda, par1, par2, loglike))
        }
    }
    
}

whichmax = function(x)
{
    return(which(x==max(x)))
}

index.outliers = function(data.xx)
{
    c = 1.5
    #data.xx = c(2, 3, 6, 9, 13, 18, 21, 106)
    Q1 = quantile(data.xx, 0.25,type=5)
    Q3 = quantile(data.xx, 0.75, type=5)
    IQD = Q3 - Q1
    lower = Q1 - c*IQD
    upper = Q3 + c*IQD
    index = which(data.xx<lower|data.xx>upper)
    #boxplot(data.xx);abline(h=Q1);abline(h=Q3);
}

correction.size.fitting = function(index=1, nb.components=3, WT=TRUE)
{
    n = index;
    
    #n =22;
    #nb.components=4;
    test = c()
    
    if(WT)
    {
        if(n<10) filename = paste('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/results/5th_try_0.25_all/Image_0', n, '_0', m, '.txt', sep='')
        if(n>=10) filename = paste('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/results/5th_try_0.25_all/Image_', n, '_0', m, '.txt', sep='')
    }else{
        filename = paste('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/results/E2f_KO_1st/Image_', (n+17), '_01', '_E2f1.txt', sep='')
        
    }
    aa = read.table(filename, sep='\t', header=FALSE)
    
    #images.size = rbind(images.size, aa[1,])
    test = c(test, aa[1,2], aa[1, 3])
    aa = aa[-1, ]
    
    colnames(aa) = c('cell.index', 'cell.size', 'nucleus.cell')
    
    bb = cbind(aa, find.repeating(aa[,1]))
    
    kk = which(bb[,4]<3 & bb[,2]<8000)
    bb = bb[kk,]
    
    
    #cat(length(unique(bb[,1])), ' cells detected\n')
    test = c(test, length(unique(bb[,1])))
    
    percent = length(which(bb[,4]==2))/2/(length(which(bb[,4]==2))/2+length(which(bb[,4]==1)))
    
    test = c(test, length(which(bb[,4]==2))/2/(length(which(bb[,4]==2))/2+length(which(bb[,4]==1))))
    
    ### cell size analysis: DOES NOT WORK!!! TOO ARBITRARY!!!
    jj = which(bb[,4]==1)
    kk = which(bb[,4]==2)
    
    test = c(test, median(bb[jj,2]))
    test = c(test, median(unique(bb[kk,2])))
    test = c(test, median(unique(bb[,2])))
    
    ## all nucleui
    x = bb[,3];
    mu.init = c(sample(seq(300, 400, by=20), 1), sample(seq(500, 600, by=20), 1), sample(seq(700, 800, by=20), 1))
    
    yy = x[which(x>250)]
    res = size.seperation(yy, method=Method, k=3, mu.init=mu.init)
    lambda = res[1:3]
    par1 = res[4:6]
    
    ### all nuclei
    res = size.seperation(x, method=Method, k=3, mu.init=par1, m.constr=par1)
    lambda = res[1:3]
    par1 = res[4:6]
    par2 = res[7:9]
    nitr = 0
    while(par1[1]<260 | par1[3]<600 |lambda[2]>lambda[1]|lambda[3]>lambda[2])
    {
        mu.init = c(sample(seq(300, 400, by=10), 1), sample(seq(500, 600, by=10), 1), sample(seq(700, 800, by=10), 1));
        res = size.seperation(x, method=Method, k=3, mu.init=mu.init)
        par1 = res[4:6]
        #cat(par1, '\n');
        nitr = nitr+1;
        if(nitr>=5) break;
    }
    
    res1 = res
    lambda = res1[1:3]
    par1 = res1[4:6]
    par2 = res1[7:9]
    
    par10 = par1
    par20 = par2
    
    test = c(test, par1)
    
    jj = which(bb[,4]==1)
    kk = which(bb[,4]==2)

    if(nb.components==3)
    {
        
        ## nuclei within mononucleated cells
        x = bb[jj,3]
        #res2 = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20, m.constr=par10, sigma.constr=par20)
        res = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20)
        lambda = res[1:3]
        par1 = res[4:6]
        par2 = res[7:9]
        nitr = 0
        while(par1[1]<260 | par1[3]<600 |lambda[2]>lambda[1]|lambda[3]>lambda[2])
        {
            #mu.init = c(sample(seq(300, 400, by=10), 1), sample(seq(500, 600, by=10), 1), sample(seq(700, 800, by=10), 1));
            mu.init = c(sample(seq(ceiling(par10[1]*0.9), ceiling(par10[1]*1.1), by=10), 1), sample(seq(ceiling(par10[2]*0.9), ceiling(par10[2]*1.1), by=10), 1),
            sample(seq(ceiling(par10[3]*0.9), ceiling(par10[3]*1.1), by=10), 1));
            res = size.seperation(x, method=Method, k=3, mu.init=mu.init, sigma.init=par20)
            par1 = res[4:6]
            #cat(par1, '\n');
            nitr = nitr+1;
            if(nitr>=5) break;
        }
        res2 = res
        
        lambda = res2[1:3]
        par1 = res2[4:6]
        par2 = res2[7:9]
        
        lambda.m = lambda
        
        xfit<-seq(min(x),max(x),length=100)
        yfit1<-dnorm(xfit, par1[1],par2[1])*lambda[1]
        yfit2<-dnorm(xfit, par1[2],par2[2])*lambda[2]
        yfit3<-dnorm(xfit, par1[3],par2[3])*lambda[3]
        #test = c(test, percent*lambda[c(1:2)])
        hist(x, breaks=50, xlim = xlim, xlab='piexls', freq=FALSE, col='darkblue', main=paste(' Nucleus Size (binucleate) ', sep=''))
        lines(xfit, yfit1, col='green', lwd=2.0)
        lines(xfit, yfit2, col='green', lwd=2.)
        lines(xfit, yfit3, col='green', lwd=2.)

        #test = c(test, (1-percent)*lambda)
        
        ## nuclei within binucleated cells
        x = bb[kk,3]
        
        res = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20, m.constr=par10, sigma.constr=par20)
        #res3 = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20, m.constr=par10)
        #res = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20)
        
        lambda = res[1:3]
        par1 = res[4:6]
        par2 = res[7:9]
        nitr = 0
        while(par1[1]<260 | par1[3]<600 |lambda[3]>lambda[2])
        {
            mu.init = c(sample(seq(ceiling(par10[1]*0.9), ceiling(par10[1]*1.1), by=10), 1), sample(seq(ceiling(par10[2]*0.9), ceiling(par10[2]*1.1), by=10), 1),
            sample(seq(ceiling(par10[3]*0.9), ceiling(par10[3]*1.1), by=10), 1));
            #res = size.seperation(x, method=Method, k=3, mu.init=mu.init, sigma.init=par20)
            res = size.seperation(x, method=Method, k=3, mu.init=mu.init, sigma.init=par20, m.constr=par10, sigma.constr=par20)
            par1 = res[4:6]
            #cat(par1, '\n');
            nitr = nitr+1;
            if(nitr>=5) break;
        }
        res3 = res
        
        lambda = res3[1:3]
        par1 = res3[4:6]
        par2 = res3[7:9]
        lambda.b = lambda
        
        xfit<-seq(min(x),max(x),length=100)
        yfit1<-dnorm(xfit, par1[1],par2[1])*lambda[1]
        yfit2<-dnorm(xfit, par1[2],par2[2])*lambda[2]
        yfit3<-dnorm(xfit, par1[3],par2[3])*lambda[3]
        #test = c(test, percent*lambda[c(1:2)])
        hist(x, breaks=50, xlim = xlim, xlab='piexls', freq=FALSE, col='darkblue', main=paste(' Nucleus Size (binucleate) ', sep=''))
        lines(xfit, yfit1, col='green', lwd=2.0)
        lines(xfit, yfit2, col='green', lwd=2.)
        lines(xfit, yfit3, col='green', lwd=2.)

        
        lambda.mm = lambda.m
        lambda.bb = lambda.b
        
        lambda2 = rbind(lambda2, lambda.mm)
        lambda3 = rbind(lambda3, lambda.bb)
        
        lambda.m = lambda.mm*(1-percent);
        lambda.b = lambda.bb*(percent);
        
        test = c(test, c((lambda.m[1]+2*lambda.b[1])/(1+percent), (lambda.m[2]+2*lambda.b[2])/(1+percent), (lambda.m[3]+2*lambda.b[3])/(1+percent)))
        test = c(test, lambda.m)
        test = c(test, lambda.b[1:2])

        return(list(test=test, lambda.m=lambda.mm, lambda.b=lambda.bb))

    }else{
        
        print('model selection')
        ## nuclei within mononucleated cells
        x = bb[jj,3]
        #res2 = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20, m.constr=par10, sigma.constr=par20)
        par10.s = c(sample(seq(150, 250, by=20), 1), par10)
        par20.s = c(sample(seq(50, 100, by=10), 1), par20)
        
        res = size.seperation(x, method=Method, k=4, mu.init=par10.s, sigma.init=par20.s)
        lambda = res[1:4]
        par1 = res[5:8]
        par2 = res[9:12]
        
        nitr = 0
        while(par1[1]>250 |par1[2]<260 | par1[3]<600 |lambda[3]>(lambda[1]+lambda[2])|lambda[4]>lambda[3])
        {
            #mu.init = c(sample(seq(300, 400, by=10), 1), sample(seq(500, 600, by=10), 1), sample(seq(700, 800, by=10), 1));
            mu.init = c(sample(seq(ceiling(par1[1]*0.9), ceiling(par1[1]*1.1), by=10), 1), sample(seq(ceiling(par1[2]*0.9), ceiling(par1[2]*1.1), by=10), 1),
            sample(seq(ceiling(par1[3]*0.9), ceiling(par1[3]*1.1), by=10), 1), sample(seq(ceiling(par1[4]*0.9), ceiling(par1[4]*1.1), by=10), 1));
            res = size.seperation(x, method=Method, k=4, mu.init=mu.init, sigma.init=par2)
            lambda = res[1:4]
            par1 = res[5:8]
            
            #cat(par1, '\n');
            nitr = nitr+1;
            if(nitr>=5) break;
        }
        res2 = res
        
        lambda = res2[1:4]
        par1 = res2[5:8]
        par2 = res2[9:12]
        
        lambda.m2 = c((lambda[1]+lambda[2]), lambda[3], lambda[4])
        #test = c(test, (1-percent)*lambda)
        
        ## nuclei within binucleated cells
        x = bb[kk,3]
        #res3 = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20, m.constr=par10, sigma.constr=par20)
        #res3 = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20, m.constr=par10)
        par10.s = c(sample(seq(150, 250, by=20), 1), par10)
        par20.s = c(sample(seq(50, 100, by=10), 1), par20)
        
        res = size.seperation(x, method=Method, k=4, mu.init=par10.s, sigma.init=par20.s)
        lambda = res[1:4]
        par1 = res[5:8]
        par2 = res[9:12]
        
        nitr = 0
        while(par1[1]>250 |par1[2]<260 | par1[3]<600 |lambda[3]>(lambda[1]+lambda[2])|lambda[4]>lambda[3])
        {
            #mu.init = c(sample(seq(300, 400, by=10), 1), sample(seq(500, 600, by=10), 1), sample(seq(700, 800, by=10), 1));
            mu.init = c(sample(seq(ceiling(par1[1]*0.9), ceiling(par1[1]*1.1), by=10), 1), sample(seq(ceiling(par1[2]*0.9), ceiling(par1[2]*1.1), by=10), 1),
            sample(seq(ceiling(par1[3]*0.9), ceiling(par1[3]*1.1), by=10), 1), sample(seq(ceiling(par1[4]*0.9), ceiling(par1[4]*1.1), by=10), 1));
            res = size.seperation(x, method=Method, k=4, mu.init=mu.init, sigma.init=par2)
            lambda = res[1:4]
            par1 = res[5:8]
            
            #cat(par1, '\n');
            nitr = nitr+1;
            if(nitr>=5) break;
        }
        res3 = res
        lambda = res3[1:4]
        par1 = res3[5:8]
        par2 = res3[9:12]
        lambda.b2 = c((lambda[1]+lambda[2]), lambda[3], lambda[4])

        lambda.mm2 = lambda.m2
        lambda.bb2 = lambda.b2
        
        lambda2 = rbind(lambda2, lambda.mm2)
        lambda3 = rbind(lambda3, lambda.bb2)
        
        lambda.m2 = lambda.mm2*(1-percent);
        lambda.b2 = lambda.bb2*(percent);
        
        test = c(test, c((lambda.m2[1]+2*lambda.b2[1])/(1+percent), (lambda.m2[2]+2*lambda.b2[2])/(1+percent), (lambda.m2[3]+2*lambda.b2[3])/(1+percent)))
        test = c(test, lambda.m)
        test = c(test, lambda.b2[1:2])
        
        return(list(test=test, lambda.m=lambda.mm2, lambda.b=lambda.bb2))

    }

}

nuclear.size.fitting.manual = function(n=1, nb.components.mono=3, nb.components.bino=3, Sigma.Constr=FALSE, reduce.small.size.population=FALSE, keep)
{
    WT = TRUE
    dna = read.csv('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/FACS_DATA.csv', sep=';', header=TRUE)
    nucleus.all.fitting = TRUE
    cell.fitting = FALSE
    calibration.facs = TRUE
    Method = 'Normal'
    
    Filter.cells = TRUE
    nucleus.size.mean = FALSE
    
    test = c()
    
    cat('n = ', n, '\n');
    tt = (n-1)*3;
    #time = c(time, rep(tt,1));
    
    par(mfrow=c(1,3))
    m = 1;
    if(WT)
    {
        if(n<10) filename = paste('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/results/5th_try_0.25_all/Image_0', n, '_0', m, '.txt', sep='')
        if(n>=10) filename = paste('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/results/5th_try_0.25_all/Image_', n, '_0', m, '.txt', sep='')
    }else{
        filename = paste('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/results/E2f_KO_1st/Image_', (n+17), '_01', '_E2f1.txt', sep='')
        
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
        cat('filter outliers of cells \n')
        cc = cell.nucleus.filtering(bb)
        cc = data.frame(cc)
        colnames(cc) = c('cell.index', 'cell.size', 'nuclear.size.s', 'nb.nuclei', 'outlier')
        
        index.outlier = cc[which(cc[,5]==1) ,1]
        mm = match(bb[,1], index.outlier)
        mm = which(is.na(mm)==TRUE)
        bb = bb[mm, ]
        cc = cc[-which(cc[,5]==1), ]
    }
    
    jj = which(bb[,4]==1)
    kk = which(bb[,4]==2)
    
    test = c(test, length(unique(bb[,1])))
    percent = length(which(bb[,4]==2))/2/(length(which(bb[,4]==2))/2+length(which(bb[,4]==1)))
    test = c(test, length(which(bb[,4]==2))/2/(length(which(bb[,4]==2))/2+length(which(bb[,4]==1))))
    test = c(test, median(bb[jj,2]))
    test = c(test, median(unique(bb[kk,2])))
    test = c(test, median(unique(bb[,2])))
    
    nucleus.size.mean = FALSE
    
    if(!nucleus.size.mean)
    {
        ## all nucleui
        x = bb[,3];
        mu.init = c(sample(seq(300, 400, by=10), 1), sample(seq(500, 600, by=10), 1), sample(seq(700, 800, by=10), 1))
        
        #res = size.seperation(x, method=Method, k=3, mu.init=mu.init)
        #lambda1 = cbind(lambda1, res[1:3])
        #res = size.seperation(x, method=Method, k=3, mu.init=mu.init)
        
        yy = x[which(x>250)]
        res = size.seperation(yy, method=Method, k=3, mu.init=mu.init)
        lambda = res[1:3]
        par1 = res[4:6]
        par2 = res[7:9]
        if(reduce.small.size.population)
        {
            #par20 = par2
            #Sigma.Constr = TRUE
            #x[which(x>250)]
            x = yy
        }
        #ratio = c(par1[2]/par1[1], par1[3]/par1[2])
        
        #cat(res[1:3], '\n')
        #res11 = size.seperation(x, method=Method, k=3, mu.init=mu.init)
        #cat(res11[1:3], '\n')
        
        ### all nuclei
        res = size.seperation(x, method=Method, k=3, mu.init=par1, m.constr=par1)
        lambda = res[1:3]
        par1 = res[4:6]
        par2 = res[7:9]
        
        nitr = 0
        while(par1[1]<260 | par1[3]<600 |lambda[2]>lambda[1]|lambda[3]>lambda[2])
        {
            mu.init = c(sample(seq(300, 400, by=10), 1), sample(seq(500, 600, by=10), 1), sample(seq(700, 800, by=10), 1));
            res = size.seperation(x, method=Method, k=3, mu.init=mu.init)
            par1 = res[4:6]
            lambda = res[1:3]
            #cat(par1, '\n');
            nitr = nitr+1;
            if(nitr>=10) break;
        }
        
        res1 = res
        lambda = res1[1:3]
        par1 = res1[4:6]
        par2 = res1[7:9]
        par10 = par1
        par20 = par2
        
        lambda.a = lambda
        
        test = c(test, par1)
        #test = c(test, lambda)
        
        xlim = c(100, 1000)
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
        
        if(reduce.small.size.population)
        {
            #par20 = par2
            #Sigma.Constr = TRUE
            x = x[which(x>250)]
            #x = yy
        }

        if(nb.components.mono==3)
        {
            if(!Sigma.Constr)
            {
                res = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20)
            }else{
                res = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20, sigma.constr=par20)
            }
            lambda = res[1:3]
            par1 = res[4:6]
            par2 = res[7:9]
            nitr = 0
            while(par1[1]<260 | par1[3]<600 |lambda[2]>lambda[1]|lambda[3]>lambda[2])
            {
                #mu.init = c(sample(seq(300, 400, by=10), 1), sample(seq(500, 600, by=10), 1), sample(seq(700, 800, by=10), 1));
                mu.init = c(sample(seq(ceiling(par10[1]*0.9), ceiling(par10[1]*1.1), by=10), 1), sample(seq(ceiling(par10[2]*0.9), ceiling(par10[2]*1.1), by=10), 1),
                sample(seq(ceiling(par10[3]*0.9), ceiling(par10[3]*1.1), by=10), 1));
                if(!Sigma.Contr)
                {
                    res = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20)
                }else{
                    res = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20, sigma.constr=par20)
                }
                
                par1 = res[4:6]
                lambda = res[1:3]
                #cat(par1, '\n');
                nitr = nitr+1;
                if(nitr>=6) break;
            }
            
            res2 = res
            lambda = res2[1:3]
            par1 = res2[4:6]
            par2 = res2[7:9]
            lambda.m = lambda
            
        }else{
            if(!Sigma.Constr)
            {
                #res = size.seperation(x, method=Method, k=2, mu.init=par10[1:2], sigma.init=par20[1:2])
                res = size.seperation(x, method=Method, k=2, mu.init=par10[1:2], sigma.init=par20[1:2])
            }else{
                res = size.seperation(x, method=Method, k=2, mu.init=par10[1:2], sigma.init=par20[1:2], sigma.constr=par20[1:2])
                #res = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20, sigma.constr=par20)
            }
            
            res2 = res
            #res5 = res1
            lambda = c(res2[1:2], 0)
            par1 = c(res2[3:4], 0)
            par2 = c(res2[5:6], 0)
            lambda.m = lambda
        }
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
        
        if(reduce.small.size.population)
        {
            #par20 = par2
            #Sigma.Constr = TRUE
            x = x[which(x>=250)]
            #x = yy
        }
        
        if(nb.components.bino==3)
        {
            if(!Sigma.Contr)
            {
                res = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20)
            }else{
                res = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20, sigma.constr=par20)
            }
            lambda = res[1:3]
            par1 = res[4:6]
            par2 = res[7:9]
            
            nitr = 0
            while(par1[1]<260 | par1[3]<600 |lambda[3]>lambda[2])
            {
                mu.init = c(sample(seq(ceiling(par10[1]*0.9), ceiling(par10[1]*1.1), by=10), 1), sample(seq(ceiling(par10[2]*0.9), ceiling(par10[2]*1.1), by=10), 1),
                sample(seq(ceiling(par10[3]*0.9), ceiling(par10[3]*1.1), by=10), 1));
                if(!Sigma.Constr)
                {
                    res = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20)
                }else{
                    res = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20, sigma.constr=par20)
                }
                
                par1 = res[4:6]
                lambda = res[1:3]
                cat(c(lambda, par1), '\n');
                nitr = nitr+1;
                if(nitr>=10) break;
            }
            res3 = res
            #res5 = res1
            lambda = res3[1:3]
            par1 = res3[4:6]
            par2 = res3[7:9]
            lambda.b = lambda
            
        }else{
            if(!Sigma.Constr)
            {
                #res = size.seperation(x, method=Method, k=2, mu.init=par10[1:2], sigma.init=par20[1:2])
                res = size.seperation(x, method=Method, k=2, mu.init=par10[1:2], sigma.init=par20[1:2])
            }else{
                res = size.seperation(x, method=Method, k=2, mu.init=par10[1:2], sigma.init=par20[1:2], sigma.constr=par20[1:2])
                #res = size.seperation(x, method=Method, k=3, mu.init=par10, sigma.init=par20, sigma.constr=par20)
            }
            #res = size.seperation(x, method=Method, k=2, mu.init=par10[1:2], sigma.init=par20[1:2])
            res3 = res
            #res5 = res1
            lambda = c(res3[1:2], 0)
            par1 = c(res3[3:4], 0)
            par2 = c(res3[5:6], 0)
            lambda.b = lambda
        }
        
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
        
        cat(dna[c(1,2,4), (n+1)]/100, '\n') ## FACS
        cat(unlist(keep[n, c(11:13)]), '\n') ## previous fitting
        cat(lambda.a, '\n') ## new fitting with all
        lambda.m = lambda.mm*(1-percent); ## new percentages of mononucleated cell populations
        lambda.b = lambda.bb*(percent); ## new percentages of binucleated cell populations
        cat(c((lambda.m[1]+2*lambda.b[1])/(1+percent), (lambda.m[2]+2*lambda.b[2])/(1+percent), (lambda.m[3]+2*lambda.b[3])/(1+percent)), '\n')
        
    }
    
    test = c(test, c((lambda.m[1]+2*lambda.b[1])/(1+percent), (lambda.m[2]+2*lambda.b[2])/(1+percent), (lambda.m[3]+2*lambda.b[3])/(1+percent)))
    test = c(test, lambda.m)
    test = c(test, lambda.b[1:2])
    
    return(test)
}

####################
########## fitting percentages of cell populaitons with different polyploidy
####################
library(emdbook)
library(deSolve)
library(fdrtool)
library(circular)
library(preprocessCore)
library(gtools)
library(biomaRt)

f2min = function(pars, data, lambda=1/c(0.0019486686, 0.0010855907, 0.0004689430, 0.0010686874, 0.0001258006), m.index = rep(1, 8))
{
    #source('functions_images.R')
    #ptm <- proc.time()
    
    #sim.percents = compute.percentages(t=c(0:31)*3, pars, rep(1,8), Decay=TRUE)
    #proc.time() - ptm
    #sim.percents = sim.percents[, -c(1,26)]
    #cat(sim.percents, '\n')
    #ptm <- proc.time()
    #nn = 2
    #m.index = as.numeric(model.index[nn, c(1:8)])
    #sim.percents =compute.percentages.approx(t=c(0:31)*3, pars, m.index=m.index, Decay=TRUE)
    #proc.time() - ptm
    #cat(sim.percents, '\n')
    #if(any(Re(sim.percents)<0) |any(Im(sim.percents)!=0))
    Decay = FALSE
    simulation = FALSE
    
    if(simulation & Decay)
    {
        t=c(0:31)*3;
        #sim.percents = compute.percentages(t=c(0:31)*3, pars=pars, m.index=m.index, Decay=TRUE)
        if(((max(t)-min(t))>24) & (max(t)+(t[2]-t[1]))%%24 == 0)
        {
            zt = t[1:(length(t)/((max(t)+(t[2]-t[1]))/24))]
        }
        
        nbs = simulate.nbs(t = zt, pars=pars, m.index=m.index, Decay=Decay);
        nb.total = apply(nbs, 1, sum)
        nb.fc = max(nb.total)/min(nb.total)
        
        x = nbs/apply(nbs, 1, sum);
        sim.percents = t(rbind(x, x, x, x))
        
        if(any(sim.percents<0) | any(nbs<10^(-8)) | any(nb.fc<0) |any(nb.fc>5))
        {
            err = Inf
        }else{
            err = sum((percents-Re(sim.percents))^2*lambda)
        }

    }else{
        sim.percents =compute.percentages.approx(t=c(0:31)*3, pars, m.index=m.index, Decay=Decay)
        if(any(Re(sim.percents)<0))
        {
            err = Inf
        }else{
            err = sum((percents-Re(sim.percents))^2*lambda)
        }

    }
    #cat('here\n')
    #cat(err, '\n' )
    return(err)
}

compute.percentages.approx = function(t=c(0:31)*3, pars, m.index=rep(1, 8), w=2*pi/24, Decay=FALSE)
{
    par.s = pars[1:4]
    par.n = pars[5:8]
    par.m = pars[9:12]
    par.d = pars[13:19]
    
    #par.d = pars[13:19]
    t = t[1:(length(t)/((max(t)+(t[2]-t[1]))/24))]
    keep = matrix(NA, nrow=5, ncol=8)
    
    for(kk in 1:length(t))
    {
        tt = t[kk];
        
        if(m.index[1]==0){
            s0 = 1*(1+m.index[2]*par.s[3]*cos(w*(tt-par.s[4])));
            s1 = s0;
            s2 = s0;
            
        }else{
            s0 = 1*(1+m.index[2]*par.s[3]*cos(w*(tt-par.s[4])));
            s1 = par.s[1]*(1+m.index[2]*par.s[3]*cos(w*(tt-par.s[4])));
            s2 = par.s[2]*(1+m.index[2]*par.s[3]*cos(w*(tt-par.s[4])));
        }
        if(m.index[3]==0){
            n0 = par.n[1]*(1+m.index[4]*par.n[3]*cos(w*(tt-par.n[4])));
            n1 = n0;
        }else{
            n0 = par.n[1]*(1+m.index[4]*par.n[3]*cos(w*(tt-par.n[4])));
            n1 = par.n[2]*(1+m.index[4]*par.n[3]*cos(w*(tt-par.n[4])));
        }
        if(m.index[5]==0){
            m0 = par.m[1]*(1+m.index[6]*par.m[3]*cos(w*(tt-par.m[4])));
            m1 = m0;
        }else{
            m0 = par.m[1]*(1+m.index[6]*par.m[3]*cos(w*(tt-par.m[4])));
            m1 = par.m[2]*(1+m.index[6]*par.m[3]*cos(w*(tt-par.m[4])));
        }
        if(Decay)
        {
            if(m.index[7]==0){
                d1 = par.d[1]*(1+m.index[8]*par.d[6]*cos(w*(tt-par.d[7])));
                d2 = d1;
                d3 = d1;
                d4 = d1;
                d5 = d1;
            }else{
                d1 = par.d[1]*(1+m.index[8]*par.d[6]*cos(w*(tt-par.d[7])));
                d2 = par.d[2]*(1+m.index[8]*par.d[6]*cos(w*(tt-par.d[7])));
                d3 = par.d[3]*(1+m.index[8]*par.d[6]*cos(w*(tt-par.d[7])));
                d4 = par.d[4]*(1+m.index[8]*par.d[6]*cos(w*(tt-par.d[7])));
                d5 = par.d[5]*(1+m.index[8]*par.d[6]*cos(w*(tt-par.d[7])));
            }
            
            M=rbind(c(-s0-d1, 0, 2*m0, 0 , 0),
            c(s0, -n0-s1-d2, 0, 0, 2*m1),
            c(0, n0, -m0-s2-d3, 0, 0),
            c(0, s1, 0, -n1-d4, 0),
            c(0, 0, s2, n1, -m1-d5))
            #print(M)
        }else{
            M=rbind(c(-s0, 0, 2*m0, 0 , 0),
            c(s0, -n0-s1, 0, 0, 2*m1),
            c(0, n0, -m0-s2, 0, 0),
            c(0, s1, 0, -n1, 0),
            c(0, 0, s2, n1, -m1))
        }
        e=eigen(M)
        #cat(e$values, '\n');
        #print(e$vectors)
        
        o=order(-Re(e$val))
        #cat(e$val[o], '\n');
        v=e$vec[,o[1]]
        
        keep[, kk] = v/sum(v)
    }
    
    keep = cbind(keep, keep, keep, keep)
    return(keep)
}

dpdt = function(t, y, pars, m.index=rep(1, 8), Decay=FALSE, w = 2*pi/24)
{
    par.s0 = 1;
    par.s = pars[1:4]
    par.n = pars[5:8]
    par.m = pars[9:12]
    par.d = pars[13:19]
    
    #m.index = pars[20:27]
    
    tt = t;
    if(m.index[1]==0){
        s0 = par.s0*(1+m.index[2]*par.s[3]*cos(w*(tt-par.s[4])));
        s1 = s0;
        s2 = s0;
        
    }else{
        s0 = par.s0*(1+m.index[2]*par.s[3]*cos(w*(tt-par.s[4])));
        s1 = par.s[1]*(1+m.index[2]*par.s[3]*cos(w*(tt-par.s[4])));
        s2 = par.s[2]*(1+m.index[2]*par.s[3]*cos(w*(tt-par.s[4])));
    }
    if(m.index[3]==0){
        n0 = par.n[1]*(1+m.index[4]*par.n[3]*cos(w*(tt-par.n[4])));
        n1 = n0;
    }else{
        n0 = par.n[1]*(1+m.index[4]*par.n[3]*cos(w*(tt-par.n[4])));
        n1 = par.n[2]*(1+m.index[4]*par.n[3]*cos(w*(tt-par.n[4])));
    }
    if(m.index[5]==0){
        m0 = par.m[1]*(1+m.index[6]*par.m[3]*cos(w*(tt-par.m[4])));
        m1 = m0;
    }else{
        m0 = par.m[1]*(1+m.index[6]*par.m[3]*cos(w*(tt-par.m[4])));
        m1 = par.m[2]*(1+m.index[6]*par.m[3]*cos(w*(tt-par.m[4])));
    }

    #Decay = TRUE
    if(Decay)
    {
        #print('Here !')
        if(m.index[7]==0){
            d1 = par.d[1]*(1+m.index[8]*par.d[6]*cos(w*(t-par.d[7])));
            d2 = d1;
            d3 = d1;
            d4 = d1;
            d5 = d1;
        }else{
            d1 = par.d[1]*(1+m.index[8]*par.d[6]*cos(w*(t-par.d[7])));
            d2 = par.d[2]*(1+m.index[8]*par.d[6]*cos(w*(t-par.d[7])));
            d3 = par.d[3]*(1+m.index[8]*par.d[6]*cos(w*(t-par.d[7])));
            d4 = par.d[4]*(1+m.index[8]*par.d[6]*cos(w*(t-par.d[7])));
            d5 = par.d[5]*(1+m.index[8]*par.d[6]*cos(w*(t-par.d[7])));
        }
        
        dxdt = 2*m0*y[3] - s0*y[1] - d1*y[1];
        dydt = s0*y[1] + 2*m1*y[5] - n0*y[2] - s1*y[2] - d2*y[2]
        dzdt = n0*y[2] - m0*y[3] - s2*y[3] - d3*y[3]
        dudt = s1*y[2] - n1*y[4] - d4*y[4]
        dvdt = n1*y[4] + s2*y[3] - m1*y[5] - d5*y[5]

    }else{
        dxdt = 2*m0*y[3] - s0*y[1];
        dydt = s0*y[1] + 2*m1*y[5]- n0*y[2] - s1*y[2]
        dzdt = n0*y[2] - m0*y[3] - s2*y[3]
        dudt = s1*y[2] - n1*y[4]
        dvdt = n1*y[4] + s2*y[3] - m1*y[5]
    }
    #cat('here \n')
    list(c(dxdt, dydt, dzdt, dudt, dvdt),NULL)
}

simulate.nbs = function(t, pars, m.index=rep(1, 8), Decay=FALSE)
{
    #t=c(0:7)*3; par.s=c(0.5, 1, 0., 22); par.n=c(0.1, 0.1, 0, 22); par.m=c(1, 1, 0., 9); par.d =c(rep(0.2, 5), 0, 10)
    offset.period = 20
    #pars = pars*100;
    
    Tstable =  24*offset.period ## burning time is two period
    t.res = 3;
    if(length(t)!=1){t.res = (t[2]-t[1])}
    t.sup = seq(0, Tstable+max(t),by= t.res)
    
    y.init = c(0.6, 0.332, 0.106, 0.034, 0.048)
    y.init = y.init/sum(y.init)*1000
    
    #pars = c(pars, m.index)
    
    ptm <- proc.time()
    simu = lsoda(y = y.init, #init.conditions
    times = t.sup, ## times
    dpdt,
    pars, m.index=m.index, Decay=Decay, rtol = 1e-10, atol = 1e-10) ## parameter values
    proc.time() - ptm
    
    #matplot(simu1[,1], simu[, c(2:6)], type = 'l',col=c(1:5), log='y')
    
    #soln[match(t+48, soln[,1]),2]
    i.last = nrow(simu);
    i.keep = seq(i.last - length(t)+1,i.last,by = 1)
    
    nbs = simu[i.keep, c(2:6)];
    
    #cat(percents[1,])
    #matplot(c(0:7)*3,percents, type = 'b',col=c(1:5))
    #mean = mean(m);m = m/mean
    #cat(soln[,2],'\n')
    #cat(parametrization,'\n')
    return(nbs)
}

compute.percentages = function(t=c(0:31)*3, pars = c(c(1, 1, 0.5, 10), c(1, 1, 0.5, 20), c(1, 1, 0.5, 2), c(1, 1, 1, 1, 1, 0.5, 10)), w=2*pi/24, m.index=rep(1, 8), Decay=FALSE, simulation=TRUE)
{
    #t=c(0:31)*3; par.s=c(1, 1, 0.5, 10); par.n=c(1, 1, 0.5, 20); par.m=c(1, 1, 0.5, 2); par.d =c(0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 10)
    #### because m is periodic, so just compute the first period and then repeat it.
    #cat(m.index, '\n');
    
    if(((max(t)-min(t))>24) & (max(t)+(t[2]-t[1]))%%24 == 0)
    {
        zt = t[1:(length(t)/((max(t)+(t[2]-t[1]))/24))]
    }
    
    nbs = simulate.nbs(t = zt, pars=pars, m.index=m.index, Decay=Decay);
    percents = nbs/apply(nbs, 1, sum);
    x = percents;
    
    return(t(rbind(x, x, x, x)));
}

compute.nbs = function(t=c(0:31)*3, pars = c(c(1, 1, 0.5, 10), c(1, 1, 0.5, 20), c(1, 1, 0.5, 2), c(1, 1, 1, 1, 1, 0.5, 10)), w=2*pi/24, m.index=rep(1, 8), Decay=FALSE, simulation=TRUE)
{
    #t=c(0:31)*3; par.s=c(1, 1, 0.5, 10); par.n=c(1, 1, 0.5, 20); par.m=c(1, 1, 0.5, 2); par.d =c(0.1, 0.1, 0.1, 0.1, 0.1, 0.5, 10)
    #### because m is periodic, so just compute the first period and then repeat it.
    if(((max(t)-min(t))>24) & (max(t)+(t[2]-t[1]))%%24 == 0)
    {
        zt = t[1:(length(t)/((max(t)+(t[2]-t[1]))/24))]
    }
    
    nbs = simulate.nbs(t = zt, pars=pars, m.index=m.index, Decay=Decay);
    #percents = nbs/apply(nbs, 1, sum);
    #x = percents;
    
    x = nbs;
    return(t(rbind(x, x, x, x)));
}

find.variance = function(x, t, period=24)
{
    #x=percents[1,];
    #t=c(2:25, 27:32)*3;
    #period = 24
    
    n=length(x)
    
    rss0 = sum((x-mean(x))^2)
    sig2.m0 = rss0/n
    bic0 = n*log(rss0/n) + 1*log(n);
    w0 = exp(-0.5*bic0)

    c=cos(2*pi*t/period)
    s=sin(2*pi*t/period)
    fit = lm(x~c+s)
    
    rss1 = sum(fit$residuals^2)
    sig2.m1 = rss1/n
    bic1 = n*log(rss1/n) + 3*log(n);
    w1 = exp(-0.5*bic1)
    
    return((sig2.m0*w0/(w0+w1)+sig2.m1*w1/(w0+w1)))
}

neighbor.detection = function(aa)
{
    test = c()
    
    x = as.matrix(aa[, c(6, 7)])
    dist <- as.matrix(dist(x, method = "euclidean", diag = FALSE, upper = FALSE))
    
    for(n in 1:nrow(aa))
    {
        #cat(n, '\n')
        #x0 = aa[n, 6];
        #y0 = aa[n, 7];
        r0 = aa[n, 3];
        test = c(test, (length(which(dist[n,]<=r0))-1))
        #if(length(which(dist[n,]<=r0))>1) {
        #    test = c(test, 1)
        #}else{
        #    test = c(test, 0)
        #}
    }
    return(test)
}

weights = function(data, index, lambda, mu, sigma)
{
    # mu = mu.0; lambda = lambda.0; sigma = sigma.0; x = data;
    prob = lambda[1]*dnorm(data, mean=mu[1], sd=sigma[1]) + lambda[2]*dnorm(data, mean=mu[2], sd=sigma[2]) + lambda[3]*dnorm(data, mean=mu[3], sd=sigma[3])
    #index =1
    part = lambda[index]*dnorm(data, mean=mu[index], sd=sigma[index])
    weight = part/prob
    return(weight)
}

loglike = function(data, lambda, mu, sigma)
{
    # mu = mu.0; lambda = lambda.0; sigma = sigma.0
    prob = lambda[1]*dnorm(data, mean=mu[1], sd=sigma[1]) + lambda[2]*dnorm(data, mean=mu[2], sd=sigma[2]) +lambda[3]*dnorm(data, mean=mu[3], sd=sigma[3])
    ll = sum(log(prob))
    return(ll)
}

normalmixEM.test = function(data, k=3, lambda.init=NULL,mu.init=NULL, sigma.init=NULL, lambda.constr=NULL, maxit = 1000, maxrestarts=10)
{
    if(!is.null(lambda.constr))
    {
        lambda.0 = lambda.constr;
    }else{
        lambda.0 = c(0.5, 0.3, 0.2)
    }
    if(!is.null(mu.init))
    {
        mu.0 = mu.init
    }else{
        mu.0=c(50, 100, 150);
    }
    if(!is.null(sigma.init))
    {
        sigma.0 = sigma.init
    }else{
        sigma.0 = c(10, 20, 30);
    }
    loglikelihood = c()
    loglike.0 = loglike(data, lambda.0, mu.0, sigma.0)
    loglikelihood = c(loglikelihood, loglike.0)
    
    epsilon = 100
    iteration = 1
    while(epsilon>10^(-10) & iteration<maxit)
    {
        #cat(iteration, '\n');
        
        lambda.1 = lambda.0
        mu.1 = rep(0, k)
        sigma.1 = rep(0, k)
        
        for(index in c(1:k))
        {
            pl = weights(data, index, lambda.0, mu.0, sigma.0)
            #lambda.1[index] = sum(pl)/length(data)
            
            pl.x = data*pl;
            mu.1[index] = sum(pl.x)/sum(pl)
            
            pl.x2 = (data-mu.1[index])^2*pl
            sigma.1[index] = sqrt(sum(pl.x2)/sum(pl))
        }
    
        #epsilon = max(c(sum((lambda.1-lambda.0)^2), sum((mu.1-mu.0)^2), sum((sigma.1-sigma.0)^2)))
        
        loglike.1 = loglike(data, lambda.1, mu.1, sigma.1)
        epsilon = abs(loglike.1 - loglike.0)
        loglikelihood = c(loglikelihood, loglike.1)
        loglike.0 = loglike.1;
        
        #cat(c(loglike.0, loglike.1, epsilon), '\n')
        #pars0 = c(lambda.0, mu.0, sigma.0)
        #pars1 = c(lambda.1, mu.1, sigma.1)
        #cat(pars0, '\n')
        #cat(pars1, '\n')
        lambda.0 = lambda.1;
        mu.0 = mu.1;
        sigma.0 = sigma.1
        
        iteration = iteration +1;
    }
    
    return(list(lambda=lambda.0, mu=mu.0, sigma=sigma.0))
}

mean.err = function(x, period=24, interval=4)
{
  x = as.numeric(x)
  kk = which(!is.na(x)==TRUE)
  index = period/interval;
  nn = length(x)/(period/interval)
  test = c()
  for(mm in 1:nn) test = cbind(test, x[c(1:index)+index*(mm-1)])
  test = data.frame(test)
  
  mean = apply(test, 1, mean.nona)
  err = apply(test, 1, sme.nona)
  
  return(rbind(mean=mean, err=err))
}
mean.nona = function(x)
{	
  xx = mean(x[which(!is.na(x)==TRUE)])
  return(xx)
}
sme.nona = function(x)
{
  kk = which(!is.na(x))
  if(length(kk)>1)
  {
    return(sd(x[kk])/sqrt(length(kk)))
  }else{
    return(0.0)
  }
}

