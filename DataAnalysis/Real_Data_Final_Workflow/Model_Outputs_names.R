###############################################
####             Exposure Models           ####
###############################################
d0baseorder <- c('agescale', 'falldr1', 'sexhor1', 'totprevBPdayscless1', 'totprevBPdaysc1to3y', 
                 'totprevBPdaysc3to5y', 'totprevBPdayscmore5', 'inhalegc1', 'thiazide1',
                 'arthritis1', 'fallcomorb1', 'emerhosp1', 'pfx1')
expobaseorder <- c('agescale', 'falldr', 'sexhor', 'totprevBPdayscless1', 'totprevBPdaysc1to3y', 
                   'totprevBPdaysc3to5y', 'totprevBPdayscmore5', 'inhalegc', 'thiazide',
                   'arthritis', 'fallcomorb', 'emerhosp', 'pfx')
expobaselabel <- c('Age', 'Fall related drugs use', 'Sex hormone', 'Prev. BP duration 0-1 year (vs. 0 days)', 
                   'Prev. BP duration 1-3 years (vs. 0 days)', 'Prev. BP duration 3-5 years (vs. 0 days)', 
                   'Prev. BP duration $\\geq$ 5 year (vs. 0 days)', 'Inhale GC', 'Thiazide',
                   'Arthritis', 'Fall related conditions', 'Emergency or hospitalization', 'Prev. fracture')
expobasenota <- c('$Z_{i1}$', '$Z_{i2}$', '$Z_{i3}$', '$Z_{i4}$', '$Z_{i5}$', 
                  '$Z_{i6}$', '$Z_{i7}$', '$Z_{i8}$', '$Z_{i9}$', '$Z_{i10}$', 
                  '$Z_{i11}$', '$Z_{i12}$', '$Z_{i13}$')


##############################################
####             Outcome Models           ####
##############################################
outcomebaseorder <- c('(Intercept)', 'timebasis_epi1', 'timebasis_epi2', d0baseorder)
outcomebaselabel <- c('Intercept', 'Baseline hazard basis 1', 'Baseline hazard basis 2', expobaselabel)
outcomebasenota <- c('$\\eta_0$', '\\theta_{01}', '\\theta_{02}', expobasenota)


############################### Cumdose in model #################################
if(!(female)){
  # D0 Model - Male
  d0order <- c(d0baseorder, '1|2', '2|3', '3|4', '4|5')
  d0label <- c(expobaselabel, 'Intercept $d_0 \\leq 1 \\mid d_0 \\geq 2$', 'Intercept $d_0 \\leq 2 \\mid d_0 \\geq 3$',
               'Intercept $d_0 \\leq 3 \\mid d_0 \\geq 4$', 'Intercept $d_0 \\leq 4 \\mid d_0 \\geq 5$')
  d0nota <- c(expobasenota, '$\\xi_{01}$', '$\\xi_{02}$', '$\\xi_{03}$', '$\\xi_{04}$')
  
  # V Model - Male
  vorder <- c('(Intercept)', expobaseorder, 'Ltc', 'cumBPdayspscale', 'Dtrtp0', 'Dtrtp2', 'Dtrtp3', 'Dtrtp4', 'Dtrtp5',
              'gccumdosescale', 'timebasis_epi1', 'timebasis_epi2')
  vlabel <- c('Intercept', expobaselabel, 'LTC', 'Cumulative BP days', 'Prev. dose $d = 0$ (vs. $d = 1$)', 'Prev. dose $d = 2$ (vs. $d = 1$)', 'Prev. dose $d = 3$ (vs. $d = 1$)', 
              'Prev. dose $d = 4$ (vs. $d = 1$)', 'Prev. dose $d = 5$ (vs. $d = 1$)', 'Cumulative GC dose', 'Time basis 1', 'Time basis 2')
  vnota <- c('$\\eta_0$', expobasenota, '$X_{i1,k}$', '$X_{i2,k}$', '$\\mathbf{1}\\{D_{i,k-1} = 0\\}$', '$\\mathbf{1}\\{D_{i,k-1} = 2\\}$', '$\\mathbf{1}\\{D_{i,k-1} = 3\\}$',
             '$\\mathbf{1}\\{D_{i,k-1} = 4\\}$', '$\\mathbf{1}\\{D_{i,k-1} = 5\\}$', '$\\textrm{cum}(\\bar{a})$', 'b_{1k}', 'b_{2k}')
  
  # D Model - Male Dtrtp=0
  dp0order <- c(expobaseorder, 'Ltc', 'cumBPdayspscale', 'gccumdosescale', 'timebasis_epi1', 'timebasis_epi2', '1|2', '2|3', '3|4', '4|5')
  dp0label <- c(expobaselabel, 'LTC', 'Cumulative BP days', 'Cumulative GC dose', 'Time basis 1', 'Time basis 2', 'Intercept $d \\leq 1 \\mid d \\geq 2$', 
                'Intercept $d \\leq 2 \\mid d \\geq 3$', 'Intercept $d \\leq 3 \\mid d \\geq 4$', 'Intercept $d \\leq 4 \\mid d \\geq 5$')
  dp0nota <- c(expobasenota, '$X_{i1,k}$', '$X_{i2,k}$', '$\\textrm{cum}(\\bar{a})$', 'b_{1k}', 'b_{2k}', '\\phi_{01}', '\\phi_{02}', '\\phi_{03}', '\\phi_{04}')
  
  
  # D Model - Male Dtrtp>0
  dpnot0order <- c(expobaseorder, 'Ltc', 'cumBPdayspscale', 'Dtrtp2', 'Dtrtp3', 'Dtrtp4', 'Dtrtp5',
                   'gccumdosescale', 'timebasis_epi1', 'timebasis_epi2', '0|1', '1|2', '2|3', '3|4', '4|5')
  dpnot0label <- c(expobaselabel, 'LTC', 'Cumulative BP days', 'Prev. dose $d = 2$ (vs. $d = 1$)', 'Prev. dose $d = 3$ (vs. $d = 1$)', 
                   'Prev. dose $d = 4$ (vs. $d = 1$)', 'Prev. dose $d = 5$ (vs. $d = 1$)', 'Cumulative GC dose', 'Time basis 1', 'Time basis 2',
                   'Intercept $d \\leq 0 \\mid d \\geq 1$', 'Intercept $d \\leq 1 \\mid d \\geq 2$', 'Intercept $d \\leq 2 \\mid d \\geq 3$',
                   'Intercept $d \\leq 3 \\mid d \\geq 4$', 'Intercept $d \\leq 4 \\mid d \\geq 5$')
  dpnot0nota <- c(expobasenota, '$X_{i1,k}$', '$X_{i2,k}$', '$\\mathbf{1}\\{D_{i,k-1} = 2\\}$', '$\\mathbf{1}\\{D_{i,k-1} = 3\\}$',
                  '$\\mathbf{1}\\{D_{i,k-1} = 4\\}$', '$\\mathbf{1}\\{D_{i,k-1} = 5\\}$', '$\\textrm{cum}(\\bar{a})$',
                  'b_{1k}', 'b_{2k}', '\\phi_{00}', '\\phi_{01}', '\\phi_{02}', '\\phi_{03}', '\\phi_{04}')  
  
  
  ### Outcome A model - Male
  aorder <- c(outcomebaseorder, 'A0', 'A2', 'A3', 'A4', 'A5')
  alabel <- c(outcomebaselabel, 'Treatment $A = 0$ (vs. $A = 1$)', 'Treatment $A = 2$ (vs. $A = 1$)', 
              'Treatment $A = 3$ (vs. $A = 1$)', 'Treatment $A = 4$ (vs. $A = 1$)', 'Treatment $A = 5$ (vs. $A = 1$)')
  anota <- c(outcomebasenota, '$\\mathbf{1}\\{A_i(t) = 0\\}$', '$\\mathbf{1}\\{A_i(t) = 2\\}$', 
             '$\\mathbf{1}\\{A_i(t) = 3\\}$', '$\\mathbf{1}\\{A_i(t) = 4\\}$', '$\\mathbf{1}\\{A_i(t) = 5\\}$')
  
  ### Outcome cumd model - Male
  cumdorder <- c(outcomebaseorder, 'gccumdosescale')
  cumdlabel <- c(outcomebaselabel, 'Cumulative GC dose')
  cumdnota <- c(outcomebasenota, '$\\textrm{cum}(\\bar{a})}$')
  
  ### Outcome spcumd model - Male
  spcumdorder <- c(outcomebaseorder, 'timebasis_cumdosescale1', 'timebasis_cumdosescale2')
  spcumdlabel <- c(outcomebaselabel, 'Cumulative GC dose spline basis 1', 'Cumulative GC dose spline basis 2')
  spcumdnota <- c(outcomebasenota, '$\\theta_{21}\\textrm{Basis}(\\textrm{cum}(\\bar{a}))}$', '$\\theta_{22}\\textrm{Basis}(\\textrm{cum}(\\bar{a}))}$')
  
}else{
  # D0 Model - Female
  d0order <- c(d0baseorder, 'pralcal1', 'diabetes1', '1|2', '2|3', '3|4', '4|5')
  d0label <- c(expobaselabel, 'Prev. raloxifene or calcitonin', 'Diabetes', 'Intercept $d_0 \\leq 1 \\mid d_0 \\geq 2$', 'Intercept $d_0 \\leq 2 \\mid d_0 \\geq 3$',
               'Intercept $d_0 \\leq 3 \\mid d_0 \\geq 4$', 'Intercept $d_0 \\leq 4 \\mid d_0 \\geq 5$')
  d0nota <- c(expobasenota, '$Z_{i14}$', '$Z_{i15}$', '$\\xi_{01}$', '$\\xi_{02}$', '$\\xi_{03}$', '$\\xi_{04}$')
  
  # V Model - Female
  vorder <- c('(Intercept)', expobaseorder, 'pralcal', 'diabetes', 'Ltc', 'cumBPdayspscale', 'Dtrtp0', 'Dtrtp2', 'Dtrtp3', 'Dtrtp4', 'Dtrtp5',
              'gccumdosescale', 'timebasis_epi1', 'timebasis_epi2')
  vlabel <- c('Intercept', expobaselabel, 'Prev. raloxifene or calcitonin', 'Diabetes', 'LTC', 'Cumulative BP days', 
              'Prev. dose $d = 0$ (vs. $d = 1$)', 'Prev. dose $d = 2$ (vs. $d = 1$)', 'Prev. dose $d = 3$ (vs. $d = 1$)', 
              'Prev. dose $d = 4$ (vs. $d = 1$)', 'Prev. dose $d = 5$ (vs. $d = 1$)', 'Cumulative GC dose', 'Time basis 1', 'Time basis 2')
  vnota <- c('$\\eta_0$', expobasenota, '$Z_{i14}$', '$Z_{i15}$', '$X_{i1,k}$', '$X_{i2,k}$',
             '$\\mathbf{1}\\{D_{i,k-1} = 0\\}$', '$\\mathbf{1}\\{D_{i,k-1} = 2\\}$', '$\\mathbf{1}\\{D_{i,k-1} = 3\\}$',
             '$\\mathbf{1}\\{D_{i,k-1} = 4\\}$', '$\\mathbf{1}\\{D_{i,k-1} = 5\\}$', '$\\textrm{cum}(\\bar{a})$', 'b_{1k}', 'b_{2k}')
  
  # D Model - Female p=0
  dp0order <- c(expobaseorder, 'pralcal', 'diabetes', 'Ltc', 'cumBPdayspscale', 'gccumdosescale', 'timebasis_epi1', 'timebasis_epi2', '1|2', '2|3', '3|4', '4|5')
  dp0label <- c(expobaselabel, 'Prev. raloxifene or calcitonin', 'Diabetes', 'LTC', 'Cumulative BP days', 'Cumulative GC dose', 'Time basis 1', 'Time basis 2',
                'Intercept $d \\leq 1 \\mid d \\geq 2$', 'Intercept $d \\leq 2 \\mid d \\geq 3$', 'Intercept $d \\leq 3 \\mid d \\geq 4$', 'Intercept $d \\leq 4 \\mid d \\geq 5$')
  dp0nota <- c(expobasenota, '$Z_{i14}$', '$Z_{i15}$', '$X_{i1,k}$', '$X_{i2,k}$', '$\\textrm{cum}(\\bar{a})$', 'b_{1k}', 'b_{2k}', '\\phi_{01}', '\\phi_{02}', '\\phi_{03}', '\\phi_{04}')
  
  # D Model - Female p>0
  dpnot0order <- c(expobaseorder, 'pralcal', 'diabetes', 'Ltc', 'cumBPdayspscale', 'Dtrtp2', 'Dtrtp3', 'Dtrtp4', 'Dtrtp5',
                   'gccumdosescale', 'timebasis_epi1', 'timebasis_epi2', '0|1', '1|2', '2|3', '3|4', '4|5')
  dpnot0label <- c(expobaselabel, 'Prev. raloxifene or calcitonin', 'Diabetes', 'LTC', 'Cumulative BP days', 'Prev. dose $d = 2$ (vs. $d = 1$)', 'Prev. dose $d = 3$ (vs. $d = 1$)', 
                   'Prev. dose $d = 4$ (vs. $d = 1$)', 'Prev. dose $d = 5$ (vs. $d = 1$)', 'Cumulative GC dose', 'Time basis 1', 'Time basis 2',
                   'Intercept $d \\leq 0 \\mid d \\geq 1$', 'Intercept $d \\leq 1 \\mid d \\geq 2$', 'Intercept $d \\leq 2 \\mid d \\geq 3$',
                   'Intercept $d \\leq 3 \\mid d \\geq 4$', 'Intercept $d \\leq 4 \\mid d \\geq 5$')
  dpnot0nota <- c(expobasenota, '$Z_{i14}$', '$Z_{i15}$', '$X_{i1,k}$', '$X_{i2,k}$', '$\\mathbf{1}\\{D_{i,k-1} = 2\\}$', 
                  '$\\mathbf{1}\\{D_{i,k-1} = 3\\}$', '$\\mathbf{1}\\{D_{i,k-1} = 4\\}$', '$\\mathbf{1}\\{D_{i,k-1} = 5\\}$', 
                  '$\\textrm{cum}(\\bar{a})$', 'b_{1k}', 'b_{2k}', '\\phi_{00}', '\\phi_{01}', '\\phi_{02}', '\\phi_{03}', '\\phi_{04}')
  
  
  ### Outcome A model - Female
  aorder <- c(outcomebaseorder, 'pralcal1', 'diabetes1', 'A0', 'A2', 'A3', 'A4', 'A5')
  alabel <- c(outcomebaselabel, 'Prev. raloxifene or calcitonin', 'Diabetes', 'Treatment $A = 0$ (vs. $A = 1$)', 'Treatment $A = 2$ (vs. $A = 1$)', 
              'Treatment $A = 3$ (vs. $A = 1$)', 'Treatment $A = 4$ (vs. $A = 1$)', 'Treatment $A = 5$ (vs. $A = 1$)')
  anota <- c(outcomebasenota, '$Z_{i14}$', '$Z_{i15}$', '$\\mathbf{1}\\{A_i(t) = 0\\}$', '$\\mathbf{1}\\{A_i(t) = 2\\}$', 
             '$\\mathbf{1}\\{A_i(t) = 3\\}$', '$\\mathbf{1}\\{A_i(t) = 4\\}$', '$\\mathbf{1}\\{A_i(t) = 5\\}$')
  
  ### Outcome cumd model - Female
  cumdorder <- c(outcomebaseorder, 'pralcal1', 'diabetes1', 'gccumdosescale')
  cumdlabel <- c(outcomebaselabel, 'Prev. raloxifene or calcitonin', 'Diabetes', 'Cumulative GC dose')
  cumdnota <- c(outcomebasenota, '$Z_{i14}$', '$Z_{i15}$', '$\\textrm{cum}(\\bar{a})}$')
  
  ### Outcome spcumd model - Female
  spcumdorder <- c(outcomebaseorder, 'pralcal1', 'diabetes1', 'timebasis_cumdosescale1', 'timebasis_cumdosescale2')
  spcumdlabel <- c(outcomebaselabel, 'Prev. raloxifene or calcitonin', 'Diabetes', 'Cumulative GC dose spline basis 1', 'Cumulative GC dose spline basis 2')
  spcumdnota <- c(outcomebasenota, '$Z_{i14}$', '$Z_{i15}$', '$\\theta_{21}\\textrm{Basis}(\\textrm{cum}(\\bar{a}))}$', 
                  '$\\theta_{22}\\textrm{Basis}(\\textrm{cum}(\\bar{a}))}$')
}
