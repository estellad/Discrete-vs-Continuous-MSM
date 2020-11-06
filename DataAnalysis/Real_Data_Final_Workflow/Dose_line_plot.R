plot_dose_line_all <- function(Dinit, dat_expo, string){
  
  d0tab<- table(Dinit$D0)
  d0tabdf <- c(0, d0tab[1], d0tab[2], d0tab[3], d0tab[4], d0tab[5])
  
  tab <- table(dat_expo$Dtrt, dat_expo$cycle)
  data <- rbind(tab[1, ], tab[2, ], tab[3, ], tab[4, ], tab[5, ], tab[6, ])
  data <- cbind(d0tabdf, data)
  
  colnames(data) <- c("0", as.character(as.factor(sort(unique(dat_expo$cycle)))))
  rownames(data) <- c(as.character(levels(as.factor(dat_expo$Dtrt))))
  
  data_percentage <- apply(data, 2, function(x){x/sum(x, na.rm = T)})
  
  #write.csv(data_percentage, '/users/edong/Output/doselineall.csv')
  
  ## Line plot
  cycle <- c(0,dtimes_V)
  ndose <- nrow(data_percentage)
  xrange <- c(-0.5,365)
  yrange <- c(0,1)
  colors <- c(7,1,2,3,4,5)
  linetype <- c(1:ndose)
  plotchar <- seq(18, 18+ndose, 1)
  
  ## Percentage plot
  pdf(paste0('/users/edong/Output/Descriptive/doseline', string, '.pdf'), width = 9, height = 6)
  plot(xrange, yrange, type = 'n', xlab="Interval (every 5 days)", ylab = "Proportion", xaxt = "n", yaxt = "n")
  for(i in 1: nrow(data_percentage)){
    lines(cycle, data_percentage[i, ], type = "b", lwd=1.5, lty=linetype[i], col = colors[i], pch=plotchar[i])
  }
  
  title("Patient Percentage in Each Dose Category Over Time")
  legend(270, 0.7, legend=c('Unexposed', 
                            "0 < daily dose <= 5mg", 
                            "5 < daily dose <= 10mg", 
                            "10 < daily dose < 30mg",
                            "30 <= daily dose < 50mg", 
                            "50mg <= daily dose"), 
         cex = 0.8, col=colors, pch=plotchar, lty=linetype, title="Dose")
  axis(side = 1, at = cycle)
  axis(side = 2, at = seq(0, 1, by=0.1))
  
  dev.off()
}

# All Dose Line
plot_dose_line_all(Dinit, dat_expo, string= 'all')


if(nonadmcensplt){
  
  plot_dose_line_nonadm <- function(Dinit, dat_expo, string, sumevents){
    
    d0tab<- table(Dinit$D0)
    d0tabdf <- c(0, d0tab[1], d0tab[2], d0tab[3], d0tab[4], d0tab[5])
    
    tab <- table(dat_expo$Dtrt, dat_expo$cycle)
    data <- rbind(tab[1, ], tab[2, ], tab[3, ], tab[4, ], tab[5, ], tab[6, ])
    data <- cbind(d0tabdf, data)
    data <- rbind(data, sumevents - colSums(data))
    
    subleg <- ifelse(string == 'dead', 'Death', 'Fracture')
    subcol <- ifelse(string == 'dead', 6, 'purple')
    colnames(data) <- c("0", as.character(as.factor(sort(unique(dat_expo$cycle)))))
    rownames(data) <- c(as.character(levels(as.factor(dat_expo$Dtrt))), subleg)
    
    data_percentage <- apply(data, 2, function(x){x/sum(x, na.rm = T)})
    data_percentage <- data_percentage[-1, ] # do not plot unexposed
    data_percentage_last_row <- data_percentage[nrow(data_percentage), ] # lag one to get difference
    lead_data_percentage_last_row <- ifelse(is.na(lead(data_percentage_last_row)), 1, lead(data_percentage_last_row))
    data_percentage_last_row <- lead_data_percentage_last_row - data_percentage_last_row
    data_percentage[nrow(data_percentage), ] <- data_percentage_last_row
    
    #write.csv(data_percentage, '/users/edong/Output/doselinedead.csv')
    #write.csv(data_percentage, '/users/edong/Output/doselinefx.csv')
    
    ## Line plot
    cycle <- c(0,dtimes_V)
    ndose <- nrow(data_percentage)
    xrange <- c(-0.5,365)
    yrange <- c(0, max((max(data_percentage)+0.1), 0.4))
    colors <- c(1,2,3,4,5,subcol)
    linetype <- c(1:ndose)
    plotchar <- seq(18, 18+ndose, 1)
    
    ## Percentage plot
    pdf(paste0('/users/edong/Output/Descriptive/doseline', string, '.pdf'), width = 9, height = 6)
    plot(xrange, yrange, type = 'n', xlab="Interval (every 5 days)", ylab = "Percentage (%)", xaxt = "n", yaxt = "n")
    for(i in 1: nrow(data_percentage)){
      lines(cycle, data_percentage[i, ], type = "b", lwd=1.5, lty=linetype[i], col = colors[i], pch=plotchar[i])
    }
    
    tlt <- ifelse(string == 'dead', 'Dead Patients', 'Fracture Patients')
    title(paste0("Patient Percentage in Each Dose Category Over Time Among ", tlt))
    legend(270, 0.3, legend=c("0<daily dose <= 5mg", 
                              "5 < daily dose <= 10mg", 
                              "10 < daily dose < 30mg",
                              "30 <= daily dose < 50mg", 
                              "50mg <= daily dose",
                              subleg), 
           cex = 0.8, col=colors, pch=plotchar, lty=linetype, title="Dose")
    axis(side = 1, at = cycle)
    axis(side = 2, at = seq(0, yrange[2], by=0.1))
    
    dev.off()
  }
  
  
  
  
  
  ############ Check among those who died, what's the percentage of high dose in initial dose ###########
  # Among ___ Initial Dose: gc_init
  dead_pt <- Dinit %>%
    filter(yt<365) %>%
    filter(status == 0) 
  # nrow(dead_pt) 
  #table(dead_pt$sex, dead_pt$D0)
  
  dat_expo_dead <- dat_expo %>%
    right_join(dead_pt[, 'study_id'])
  
  
  fx_pt <- Dinit %>%
    filter(status == 1) 
  # nrow(fx_pt) 
  #table(fx_pt$sex, fx_pt$D0)
  
  dat_expo_fx <- dat_expo %>%
    right_join(fx_pt[, 'study_id'])
  
  
  # Death Pt Line
  plot_dose_line_nonadm(dead_pt, dat_expo_dead, 'dead', nrow(dead_pt))
  
  # Fx Pt Line
  plot_dose_line_nonadm(fx_pt, dat_expo_fx, 'fx', nrow(fx_pt))
}




## Count plot
#   plot(xrange, c(0,nobs), type = 'n', xlab="Interval (every 5 days)", ylab = "Number Pt (n)", xaxt = "n", yaxt = "n")
#   for(i in 1: nrow(data)){
#     lines(cycle, data[i, ], type = "b", lwd=1.5, lty=linetype[i], col = colors[i], pch=plotchar[i])
#   }
#   
#   title("Patient Count in Each Dose Category Over Time")
#   legend(270, 58000, rownames(data), cex = 1, col=colors, pch=plotchar, lty=linetype, title="Dose")
#   axis(side = 1, at = cycle)
#   axis(side = 2, at = seq(0, 80000, by=10000))
