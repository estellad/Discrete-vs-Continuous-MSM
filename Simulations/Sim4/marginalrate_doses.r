
# Proportional hazard functions:

# eta=-1.5
eta=1.5
b1=-1.0
b2=2*b1

lambda <- function(x, zt=Inf, xt=Inf, d1, d2) {
    return(0.02 * exp(0.0 * x + b1 * (x < zt) * (d1==1) + b2 * (x < zt) * (d1==2) + 
                                b1 * (x >= zt) * (d2==1) + b2 * (x >= zt) * (d2==2) + 1.5 * I(x >= xt)))
}
mu <- function(x, zt=Inf, d1, d2) {
    return(0.05 * exp(eta * (x < zt) * (d1==1) + eta * (x < zt) * (d1==2) + 
                      eta * (x >= zt) * (d2==1) + eta * (x >= zt) * (d2==2)))
}

# Cumulative hazards:

clambda <- function(x, zt=Inf, xt=Inf, d1, d2) {
    int <- rep(0, length(x))
    for (k in 1:length(x)) {
        times <- sort(c(0, x[k], zt, xt))
        for (i in 2:(2 + I(zt < x[k]) + I(xt < x[k]))) {
            int[k] <- int[k] + integrate(lambda, lower=times[i-1], upper=times[i], zt=zt, xt=xt, d1=d1, d2=d2)$value
        }
    }
    return(int)
}

cmu <- function(x, zt=Inf, d1, d2) {
    int <- rep(0, length(x))
    for (k in 1:length(x)) {
        times <- sort(c(0, x[k], zt))
        for (i in 2:(2 + I(zt < x[k]))) {
            int[k] <- int[k] + integrate(mu, lower=times[i-1], upper=times[i], zt=zt, d1=d1, d2=d2)$value
        }
    }
    return(int)
}

# Plot the functions:

tgrid <- seq(0,5,by=0.01)
zt=Inf; xt=Inf
d1=1;d2=1

cl <- clambda(tgrid, zt=zt, xt=xt, d1=d1, d2=d2)
plot(tgrid, lambda(tgrid, zt=zt, xt=xt, d1=d1, d2=d2), type='l', ylim=c(0,0.2))
lines(tgrid, cl, col='red')
lines(tgrid, 1-exp(-cl), col='blue')
abline(h=0, lty='dotted')

cm <- cmu(tgrid, zt=zt, d1=d1, d2=d2)
plot(tgrid, mu(tgrid, zt=zt, d1=d1, d2=d2), type='l', ylim=c(0,1))
lines(tgrid, cm, col='red')
lines(tgrid, 1-exp(-cm), col='blue')
abline(h=0, lty='dotted')

xprob <- function(x, t, zt=Inf, d1, d2) {
    int <- rep(0, length(x))
    for (k in 1:length(x)) {
        int[k] <- exp(-clambda(t, zt=zt, xt=x[k], d1=d1, d2=d2)) * mu(x[k], zt=zt, d1=d1, d2=d2) * exp(-cmu(x[k], zt=zt, d1=d1, d2=d2))
    }    
    return(int)
}

xp <- xprob(tgrid, t=5, zt=zt, d1=d1, d2=d2)
plot(tgrid, xp, type='l', ylim=c(0,1))

xint <- function(x, zt=zt, d1, d2) {
    int <- 0
    times <- sort(c(0, x, zt))
    for (i in 2:(2 + I(zt < x))) {
        int <- int + integrate(xprob, lower=times[i-1], upper=times[i], t=x, zt=zt, d1=d1, d2=d2)$value
    }
    return(int)
}

eps <- 0.001
getlm <- function(tgrid, zt, d1, d2) {
    lm <- rep(NA, length(tgrid))
    for (i in 1:length(lm)) {
        c0 <- exp(-clambda(tgrid[i], zt=zt, xt=(tgrid[i] + eps), d1=d1, d2=d2)) * exp(-cmu(tgrid[i], zt=zt, d1=d1, d2=d2))
        c1 <- xint(tgrid[i], zt=zt, d1=d1, d2=d2)
        lm[i] <- (lambda(tgrid[i], zt=zt, xt=(tgrid[i] + eps), d1=d1, d2=d2) * c0 + lambda(tgrid[i], zt=zt, xt=(tgrid[i] - eps), d1=d1, d2=d2) * c1)/(c0+c1)
        if (i %% 10 == 0)
            print(i)
    }
    return(lm)
}
lm0A0 <- getlm(tgrid, Inf, d1=1, d2=1)
lm1A0 <- getlm(tgrid, 2, d1=1, d2=0)
# lm2 <- getlm(tgrid, 1)

lm0A2 <- getlm(tgrid, Inf, d1=1, d2=1)
lm1A2 <- getlm(tgrid, 2, d1=1, d2=2)

if (1) {
    outpath <- 'C:/Users/dongy30/Desktop/Results/Scenario4'
    pdf(file.path(outpath, paste('rates.pdf', sep='')), height=6, width=6, paper='special')
    op <- par(mar=c(4.5,4.5,1,1), oma=c(0,0,0,0))
    layout(matrix(c(1,1,2), 3, 1, byrow = TRUE))
    # layout.show(n = 2)

    plot(tgrid, lambda(tgrid, zt=2, xt=Inf), type='l', ylim=c(0,0.1), lty='dashed', xlab='Follow-up years', ylab='Hazard')
    lines(tgrid, lambda(tgrid, zt=2, xt=0), lty='dotdash')
    # lines(tgrid, lm2, col='red', lwd=2, lty='dashed')
    lines(tgrid, lm1, col='red', lwd=2, lty='solid')
    lines(tgrid, lm0, col='black', lwd=2, lty='solid')
    legend('topright', legend=c('Marginal hazard (always treat)',
    'Marginal hazard (stop at t=2)','Conditional hazard (stop at t=2, X(0)=1)', 'Conditional hazard (stop at t=2, X(5)=0)'), 
    col=c('black','red','black','black'), lwd=c(2,2,1,1), lty=c('solid','solid','dotdash','dashed'), bg='white')

    plot(tgrid, lm1/lm0, type='l', ylim=c(0.5,3), ylab='Hazard ratio', col='red', xlab='Follow-up years', lwd=2, lty='solid')
    # lines(tgrid, lm1/lm0, col='red', lwd=2, lty='solid')
    abline(h=seq(0.5,3,by=0.5), lty='dotted')
    legend('topright', legend=c('Hazard ratio (stop at t=2 vs always treat)'), 
           col=c('red'), lwd=2, bg='white', lty=c('solid'))

    par(op)
    dev.off()

    # op <- par(mar=c(4.5,4.5,1,1), oma=c(0,0,0,0))
    # layout(matrix(c(1,1,2), 3, 1, byrow = TRUE))
    # # layout.show(n = 2)
    # 
    # pdf(file.path(outpath, paste('logrates.pdf', sep='')), height=6, width=6, paper='special')
    # op <- par(mar=c(4.5,4.5,1,1), oma=c(0,0,0,0))
    # layout(matrix(c(1,1,2), 3, 1, byrow = TRUE))
    # layout.show(n = 2)

    # plot(tgrid, log(lambda(tgrid, zt=2, xt=Inf)), type='l', ylim=c(-5,-2), lty='dashed', xlab='Follow-up years', ylab='Log hazard')
    # lines(tgrid, log(lambda(tgrid, zt=2, xt=0)), lty='dotdash')
    # # lines(tgrid, log(lm2), col='red', lwd=2, lty='dashed')
    # lines(tgrid, log(lm1), col='red', lwd=2, lty='solid')
    # lines(tgrid, log(lm0), col='black', lwd=2, lty='solid')
    # legend('topright', legend=c('Marginal log hazard (always treat)',
    # 'Marginal log hazard (stop at t=2)','Conditional log hazard (stop at t=2, X(0)=1)', 'Conditional log hazard (stop at t=2, X(5)=0)'), 
    # col=c('black','red','black','black'), lwd=c(2,2,1,1), lty=c('solid','solid','dotdash','dashed'), bg='white')

    ################### Put A=0 A=2 together ####################
    pdf(file.path(outpath, paste('logrates_sim4.pdf', sep='')), height=6, width=6, paper='special')
    op <- par(mar=c(4.5,4.5,1,1), oma=c(0,0,0,0))
    layout(matrix(c(1,2), 2, 1, byrow = TRUE))
    
    plot(tgrid, log(lm1A0/lm0A0), type='l', ylim=c(-0.5,2.5), ylab='Log hazard ratio', col='red', xlab='Follow-up years', lwd=2, lty='solid', yaxt="n", xaxt="n", cex.lab=0.7)
    abline(h=seq(-0.5, 2.5,by=0.5), lty='dotted')
    axis(2,cex.axis=0.7)
    axis(1,cex.axis=0.7)
    legend('topright', legend=c('Log hazard ratio (stop treatment A=0 at t=2 vs always treat with A=1)'), 
           col=c('red'), lwd=2, bg='white', lty=c('solid'), cex = 0.70)
    
    plot(tgrid, log(lm1A2/lm0A2), type='l', ylim=c(-1.5,1.5), ylab='Log hazard ratio', col='red', xlab='Follow-up years', lwd=2, lty='solid', yaxt="n", xaxt="n", cex.lab=0.7)
    abline(h=seq(-1.5, 1.5,by=0.5), lty='dotted')
    axis(2,cex.axis=0.7)
    axis(1,cex.axis=0.7)
    legend('topright', legend=c('Log hazard ratio (change dosage A=2 at t=2 vs always treat with A=1)'), 
           col=c('red'), lwd=2, bg='white', lty=c('solid'), cex = 0.70)

    par(op)
    dev.off()
}




