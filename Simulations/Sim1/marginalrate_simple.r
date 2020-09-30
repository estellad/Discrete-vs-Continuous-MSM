
# Additive hazard functions:

if (0) {
    lambda <- function(x, zt=Inf, xt=Inf) {
        return(0.02 - 0.01 * I(x >= zt) + 0.01 * I(x >= xt))
    }
    mu <- function(x, zt=Inf) {
        return(0.3 - 0.15 * (x >= zt))
    }
}

# Proportional hazard functions:

eta=-1.5
#eta=0.0

lambda <- function(x, zt=Inf, xt=Inf) {
    return(0.02 * exp(0.0 * x - 1.0 * I(x >= zt) + 1.5 * I(x >= xt)))
}
mu <- function(x, zt=Inf) {
    return(0.05 * exp(eta * (x >= zt)))
}

# Cumulative hazards:

clambda <- function(x, zt=Inf, xt=Inf) {
    int <- rep(0, length(x))
    for (k in 1:length(x)) {
        times <- sort(c(0, x[k], zt, xt))
        for (i in 2:(2 + I(zt < x[k]) + I(xt < x[k]))) {
            int[k] <- int[k] + integrate(lambda, lower=times[i-1], upper=times[i], zt=zt, xt=xt)$value
        }
    }
    return(int)
}

cmu <- function(x, zt=Inf) {
    int <- rep(0, length(x))
    for (k in 1:length(x)) {
        times <- sort(c(0, x[k], zt))
        for (i in 2:(2 + I(zt < x[k]))) {
            int[k] <- int[k] + integrate(mu, lower=times[i-1], upper=times[i], zt=zt)$value
        }
    }
    return(int)
}

# Plot the functions:

tgrid <- seq(0,5,by=0.01)
zt=Inf; xt=Inf

cl <- clambda(tgrid, zt=zt, xt=xt)
plot(tgrid, lambda(tgrid, zt=zt, xt=xt), type='l', ylim=c(0,0.2))
lines(tgrid, cl, col='red')
lines(tgrid, 1-exp(-cl), col='blue')
abline(h=0, lty='dotted')

cm <- cmu(tgrid, zt=zt)
plot(tgrid, mu(tgrid, zt=zt), type='l', ylim=c(0,1))
lines(tgrid, cm, col='red')
lines(tgrid, 1-exp(-cm), col='blue')
abline(h=0, lty='dotted')

xprob <- function(x, t, zt=Inf) {
    int <- rep(0, length(x))
    for (k in 1:length(x)) {
        int[k] <- exp(-clambda(t, zt=zt, xt=x[k])) * mu(x[k], zt=zt) * exp(-cmu(x[k], zt=zt))
    }    
    return(int)
}

xp <- xprob(tgrid, t=5, zt=zt)
plot(tgrid, xp, type='l', ylim=c(0,1))

xint <- function(x, zt=zt) {
    int <- 0
    times <- sort(c(0, x, zt))
    for (i in 2:(2 + I(zt < x))) {
        int <- int + integrate(xprob, lower=times[i-1], upper=times[i], t=x, zt=zt)$value
    }
    return(int)
}

eps <- 0.001
getlm <- function(tgrid, zt) {
    lm <- rep(NA, length(tgrid))
    for (i in 1:length(lm)) {
        c0 <- exp(-clambda(tgrid[i], zt=zt, xt=(tgrid[i] + eps))) * exp(-cmu(tgrid[i], zt=zt))
        c1 <- xint(tgrid[i], zt=zt)
        lm[i] <- (lambda(tgrid[i], zt=zt, xt=(tgrid[i] + eps)) * c0 + lambda(tgrid[i], zt=zt, xt=(tgrid[i] - eps)) * c1)/(c0+c1)
        if (i %% 10 == 0)
            print(i)
    }
    return(lm)
}
lm0 <- getlm(tgrid, Inf)
lm1 <- getlm(tgrid, 2)
# lm2 <- getlm(tgrid, 1)

if (0) {
    outpath <- 'C:/Users/dongy30/Desktop/Olli_Paper_Try_Sim'
    pdf(file.path(outpath, paste('rates.pdf', sep='')), height=6, width=6, paper='special')
    op <- par(mar=c(4.5,4.5,1,1), oma=c(0,0,0,0))
    layout(matrix(c(1,1,2), 3, 1, byrow = TRUE))
    # layout.show(n = 2)

    plot(tgrid, lambda(tgrid, zt=2, xt=Inf), type='l', ylim=c(0,0.1), lty='dashed', xlab='Follow-up years', ylab='Hazard')
    lines(tgrid, lambda(tgrid, zt=2, xt=0), lty='dotdash')
    # lines(tgrid, lm2, col='red', lwd=2, lty='dashed')
    lines(tgrid, lm1, col='red', lwd=2, lty='solid')
    lines(tgrid, lm0, col='black', lwd=2, lty='solid')
    legend('topright', legend=c('Marginal hazard (never treat)',
    'Marginal hazard (treat at t=2)','Conditional hazard (treat at t=2, Z(0)=1)', 'Conditional hazard (treat at t=2, Z(5)=0)'), 
    col=c('black','red','black','black'), lwd=c(2,2,1,1), lty=c('solid','solid','dotdash','dashed'))

    plot(tgrid, lm1/lm0, type='l', ylim=c(0,1.5), ylab='Hazard ratio', col='red', xlab='Follow-up years', lwd=2, lty='solid')
    # lines(tgrid, lm1/lm0, col='red', lwd=2, lty='solid')
    abline(h=c(0,0.5,1,1.5), lty='dotted')
    legend('topright', legend=c('Hazard ratio (treat at t=2 vs never treat)'), 
           col=c('red'), lwd=2, bg='white', lty=c('solid'))

    par(op)
    dev.off()

    op <- par(mar=c(4.5,4.5,1,1), oma=c(0,0,0,0))
    layout(matrix(c(1,1,2), 3, 1, byrow = TRUE))
    # layout.show(n = 2)

    pdf(file.path(outpath, paste('logrates.pdf', sep='')), height=6, width=6, paper='special')
    op <- par(mar=c(4.5,4.5,1,1), oma=c(0,0,0,0))
    layout(matrix(c(1,1,2), 3, 1, byrow = TRUE))
    # layout.show(n = 2)

    plot(tgrid, log(lambda(tgrid, zt=2, xt=Inf)), type='l', ylim=c(-5,-2), lty='dashed', xlab='Follow-up years', ylab='Log hazard')
    lines(tgrid, log(lambda(tgrid, zt=2, xt=0)), lty='dotdash')
    # lines(tgrid, log(lm2), col='red', lwd=2, lty='dashed')
    lines(tgrid, log(lm1), col='red', lwd=2, lty='solid')
    lines(tgrid, log(lm0), col='black', lwd=2, lty='solid')
    legend('topright', legend=c('Marginal log hazard (never treat)',
    'Marginal log hazard (treat at t=2)','Conditional log hazard (treat at t=2, X(0)=1)', 'Conditional log hazard (treat at t=2, X(5)=0)'), 
    col=c('black','red','black','black'), lwd=c(2,2,1,1), lty=c('solid','solid','dotdash','dashed'))

    plot(tgrid, log(lm1/lm0), type='l', ylim=c(-2,0), ylab='Log hazard ratio', col='red', xlab='Follow-up years', lwd=2, lty='solid')
    # lines(tgrid, log(lm1/lm0), col='red', lwd=2, lty='solid')

    # lines(tgrid-1, log(lm1/lm0), col='blue', lwd=1)
    abline(h=seq(-2,0,by=0.5), lty='dotted')
    legend('topright', legend=c('Log hazard ratio (treat at t=2 vs never treat)'), 
           col=c('red'), lwd=2, bg='white', lty=c('solid'))

    par(op)
    dev.off()
}

# plot(tgrid, log(lm2/lm0), type='l', ylim=c(-1.5,0), ylab='Log hazard ratio', col='green', xlab='Follow-up years', lwd=2)
# lines(tgrid, log(lm1/lm0), col='blue', lwd=2)
# zt <- 1
# lines(tgrid, (tgrid >= zt) * (-1) - 0.11 * (tgrid >= zt) * (tgrid-zt)^(0.65) + 0.0 * (tgrid >= zt) * (tgrid-zt)^(0.65) * zt, col='black', lwd=1)
# zt <- 2
# lines(tgrid, (tgrid >= zt) * (-1) - 0.11 * (tgrid >= zt) * (tgrid-zt)^(0.65) + 0.0 * (tgrid >= zt) * (tgrid-zt)^(0.65) * zt, col='black', lwd=1)
# abline(h=c(-2,-1,0), lty='dotted')



