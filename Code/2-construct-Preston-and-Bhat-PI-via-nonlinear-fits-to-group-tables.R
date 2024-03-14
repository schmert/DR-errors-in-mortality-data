#..................................................................
# Carl Schmertmann
# 21 Nov 2023
# 
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ACHTUNG!
#   The parameter fits are VERY sensitive to the 
# choice of fitting algorithm, and to exact way that Pmat()
# is calculated. That's ok, because we only want
# one example that fits reasonably well. But be aware
# that there are many other plausible patterns of
# net imports, mx errors, etc. that could also be
# "Preston" or "Bhat" patterns.
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
# 
# finds a set of single-year age misstatement probabilities
# for ages {0,1,...,110} x {0,1,...,110} that closely
# (but not exactly) match a given set of 5-year group
# misstatement probabilities. 
# 
# In this implementation we match either the estimated
# Male misstatement probabilities for India from Table 3 of 
# 
#  Bhat, Estimating Transition Probabilities of Age Misstatement
#  Demography, Vol. 27, No. 1 (Feb., 1990), pp. 149-163 (15 pages)
#  https://doi.org/10.2307/2061559
#
# or those for African-American females in Table 2 of
# 
# PRESTON, SAMUEL & Elo, Irma. (1999). 
# Effects of Age Misreporting on Mortality Estimates at Older Ages. 
# Population Studies. 53. 165-177. 10.1080/00324720308075. 
# 
# For the Preston table, the program uses the single-year distribution of *deaths* 
# in a stationary population with São Paulo male 2010 mortality rates to weight 
# within age groups. For the Bhat table it uses the single-year distribution
# of *exposure* as the within-group weights.
#..................................................................

library(tidyverse)
library(here)
library(kableExtra)

rm(list=ls())
graphics.off()

# SP data ----

age  = 0:110
nage = length(age)

# approximate distribution of male exposure (eta) and deaths (delta) by 
# single-year age 0,1,...110 for a stationary
# population with Sao Paulo state 2010 mortality

# "true" non-stationary population

eta = c(269285, 268487, 270377, 275873, 281490, 283945, 284819, 287323, 
        298616, 316883, 330783, 336906, 335200, 336711, 338696, 337500, 
        334031, 330425, 332426, 340674, 351293, 361943, 367144, 372828, 
        373930, 372627, 373693, 378392, 382437, 382962, 370419, 359287, 
        341154, 333843, 329087, 322005, 314128, 304822, 303283, 306637, 
        302076, 296625, 282824, 279392, 278976, 277463, 270325, 257455, 
        250511, 249664, 244512, 239161, 226744, 219626, 212211, 205289, 
        195219, 183815, 173928, 168460, 158396, 149950, 138939, 131246, 
        122948, 114511, 106753, 98477, 92243, 88735, 83637, 78994, 72250, 
        68461, 63655, 58560, 53068, 48750, 45170, 42277, 38003, 33905, 
        29273, 25512, 21740, 18248, 15038, 12212, 9952, 8300, 6689, 5151, 
        3764, 2759, 2047, 1563, 1153, 835, 591, 406, 275, 183, 119, 75, 
        47, 30, 20, 12, 7, 4, 2)

# "true" non-stationary deaths

delta = c(4081, 272.3, 148.5, 112.6, 90.2, 80.4, 76.2, 72.7, 71.2, 72.3, 
          74.6, 78.5, 87.3, 108.6, 150.8, 226.6, 328.5, 415.6, 483, 545.8, 
          593.2, 637.8, 654.4, 667.4, 676.4, 690.5, 716.4, 752.2, 784.2, 
          803.5, 801.9, 820.4, 825.4, 844.8, 868.9, 896.8, 928, 953.2, 
          1002.2, 1074.3, 1139.4, 1209.1, 1250, 1317.6, 1393.1, 1475.5, 
          1545.8, 1596.1, 1672.3, 1764.7, 1845, 1984.6, 2075.3, 2161.8, 
          2220.2, 2299.3, 2355.3, 2388.4, 2415, 2473.5, 2476.8, 2550.4, 
          2570.7, 2589.2, 2582.1, 2610.1, 2668.1, 2690.1, 2739.7, 2857.2, 
          2934.1, 3051.7, 3069.9, 3160.3, 3173.2, 3162, 3119.2, 3130.8, 
          3162.8, 3208.5, 3141.4, 3093.2, 2938.5, 2768.5, 2543.2, 2332.2, 
          2107.7, 1858, 1617.4, 1420.8, 1221, 1035.8, 838.5, 664, 520.6, 
          414.9, 316.4, 235, 171.4, 126.4, 97.5, 68.5, 47, 31.4, 20.9, 
          13.7, 9.6, 6.1, 3.6, 2.1, 1.3)

names(eta) = names(delta) = age

# convert to stationary exposure and events

mx = delta/eta
Hx = cumsum(c(0,mx))
lx = exp(-Hx)

exposure = eta['0'] * (head(lx,-1)+tail(lx,-1))/2 
events   = exposure * mx

names(exposure) = names(events) = age

# create age group misstatment tables ----
source( here('Code','make-age-group-misstatement-tables.R') )

# add some table-specific information ---- 
table_info = tribble(
  ~author, ~wt_var,
  'Preston', 'events',
  'Bhat'   , 'exposure'
) %>% 
  mutate( 
   data = map(author, function(a) get(paste0(a,'_table'))),
   gnames = map(data, function(d) sort(unique(d$Y)))
  )      
      
# construct w vectors for age groups ----

# W[i,] has 1s and 0s such that
# W[i,] %*% delta  = total deaths in i-th age group
# W[i,] %*% Q %*% diag(delta) %*% t(W[j,]) = deaths in true group j rep. as group i

make_W = function(gnames,Ymin=40) {

  rr = cut(age, c(gnames,Inf), right=FALSE) %>% 
         as.numeric()
  cc = seq(age)
  
  ix = cbind(rr,cc)

  W     = matrix(0, nrow= length(gnames), ncol=nage,
               dimnames=list(gnames, age))
  W[ix] = 1

  return( W[as.numeric(rownames(W)) >= Ymin,])
}

Y0 = 0  # we'll only consider groups with starting age >= Y0

table_info = table_info %>% 
  mutate(W  = map(gnames, make_W, Ymin= 0),
         W0 = map(gnames, make_W, Ymin=Y0)
         )

# table_info is a list data frame
table_info 


# probability functions ---- 

# calculate the probs of correct, understated, overstated ages for 
# true ages y in {age}
# 
# The model for ages y=40+ is 
#   log( P[understate]/P[correct]) at age y = alpha_u + beta_u * (y-40)
#   log( P[overstate] /P[correct]) at age y = alpha_v + beta_v * (y-40)
#   
# for y= 0-39 we just assume correct age reporting (this won't matter because
# we focus on results for ages 60+)  

# returns an nage x 3 table of probabilties for (correct,under,over) reports
# of age, 1 triple per age

calc_p = function(alpha_u, beta_u, alpha_v, beta_v) {
  
  v = cbind(1, 
          exp(alpha_u + beta_u * (age-Y0)),
          exp(alpha_v + beta_v * (age-Y0)))
  v[age < Y0, 1] = 1
  v[age < Y0, 2] = 0
  v[age < Y0, 3] = 0
  dimnames(v) = list(age, c('diag','under','over'))
  return( prop.table(v,1) )
}


# construct an nage x nage matrix of single-year misstatement probabilities
# using the model above for (correct,under,over) misstatement, plus
#  P(error of k=1,2,... years | UNDERstatement) prop to (rho_u)^k
#  P(error of k=1,2,... years | OVERstatement)  prop to (rho_v)^k
 
Pmat = function(alpha_u, beta_u, 
                alpha_v, beta_v, 
                rho_u, rho_v, 
                debug=FALSE) {
  
  mat = matrix(0, nrow=nage, ncol=nage, dimnames=list(age,age))
  
  if (debug) print( sapply(age, function(a) p(alpha_u,beta_u,alpha_v,beta_v,a)) )
  
  pr = calc_p(alpha_u,beta_u,alpha_v,beta_v)
  
  diag(mat) = pr[,'diag']
  
  for (y in age) {
    
    xu = age[age<y]    
    #    xu = age[age<y & abs(age-y) <= 15]    
    
    if (length(xu) > 0) {
      multip = prop.table( dgeom(y-xu-1, 1-rho_u))
      mat[paste(xu), paste(y)] = multip * pr[paste(y),'under']
    }
    
    xv = age[age>y]        
    #    xv = age[age>y & abs(age-y) <= 15]    
    
    if (length(xv) > 0) {
      multip = prop.table( dgeom(xv-y-1, 1-rho_v))
      mat[paste(xv), paste(y)] = multip * pr[paste(y),'over']
    }
    
  } # for y
  
  # adjust corner cells 
  mat['0','0']     = mat['0','0']     + (1-colSums(mat)['0']) 
  mat['110','110'] = mat['110','110'] + (1-colSums(mat)['110']) 
  
  return(mat)
}  

# create list to save estimated parameters 
par        = vector('list', nrow(table_info))
names(par) = table_info$author

# main loop over tables/cases ----

for (case in 1:nrow(table_info)) {

  target_author = table_info$author[case]
  target_table  = table_info$data[[case]]
  W             = table_info$W[[case]]
  W0            = table_info$W0[[case]]
  gnames        = table_info$gnames[[case]]
  
  print("=========================================")
  print(target_author)
  print("=========================================")
  
  # select either exposure or events as single-year weighting variable ----
  wt_var  = table_info$wt_var[[case]]
  this_wt = get(wt_var)

  
  # construct a target df with reported vs. true age groups
  # for the selected measure (either exposure or events)
  # in the SP male stationary population

  target_df =   expand_grid(Y    = gnames, 
                            X    = gnames, 
                            frac = 0)   %>% 
    left_join(target_table, by=c('Y','X')) %>% 
    mutate( frac   = map2_dbl(frac.x,frac.y, function(x,y) {ifelse(is.na(y),x,y)})) %>% 
    select(-frac.x, -frac.y)

  count = as.vector(W %*% this_wt)
  names(count) = gnames
  
  count_df = tibble(Y = as.numeric(names(count)), count=count)
  
  target_df = target_df %>% 
               left_join(count_df) %>% 
               mutate(target = frac * count) %>% 
               filter(X >= Y0, Y >= Y0)
  
  # predicted deaths by GROUP at ages 40+
  model_estimate = function(alpha_u,beta_u, alpha_v, beta_v, rho_u, rho_v) {
  
  # print(c(   alpha_u = alpha_u,
  #            beta_u  = beta_u,
  #            alpha_v = alpha_v,
  #            beta_v  = beta_v,
  #            rho_u   = rho_u,
  #            rho_v   = rho_v))
    
    PI      = Pmat(alpha_u,beta_u, alpha_v, beta_v, rho_u, rho_v) 
    est1    = PI %*% diag(this_wt)
    estG    = W0 %*% est1 %*% t(W0)
    return(as.vector(estG))
  }

  obs = target_df %>% pull(target)
  
  
  err = function(theta) {
    pred = model_estimate(alpha_u = theta[1],
                          beta_u  = theta[2],
                          alpha_v = theta[3],
                          beta_v  = theta[4],
                          rho_u   = theta[5],
                          rho_v   = theta[6])
    
    return( sum( (obs - pred)^2)  )
  }
  
  fit = optim(par=c(0,0,0,0,.50,.50), err,
              method='BFGS',
              control=list(trace  = 4,
                           maxit  = 5000,
                           fnscale= 1e5,
                           reltol=1e-12))
  
  summary(fit)

  theta = fit$par
  names(theta) = c('alpha_u','beta_u',
                   'alpha_v','beta_v',
                   'rho_u', 'rho_v')
  
  for (k in names(theta)) assign(k, theta[k])
  
  par[[target_author]] = theta
  
  # add model fits to target df
  
  target_df = target_df  %>% 
    add_column(estimate = model_estimate(alpha_u,beta_u,
                                         alpha_v,beta_v,
                                         rho_u,rho_v)) %>% 
    group_by(Y) %>% 
    mutate(frac = prop.table(target),
           est_frac = prop.table(estimate)) %>%
    mutate(across(ends_with('frac'),round,4))
  
  target_df %>% 
    select(Y,X,target,estimate,frac,est_frac) %>% 
    print(n=999)

  # lots of plots ....
  
  graphics.off()
  
  pdf(here('Plots', paste0('final-',target_author, '-nonlinear-fit.pdf')))
  
  G = ggplot(data=target_df) +
    aes(x=target, y=estimate, color=Y, label=Y) +
    geom_text() +
    geom_abline(intercept = 0, slope=1) +
    theme_bw() +
    labs(title=target_author, y=wt_var) +
    guides(color=FALSE)
  
  G = G +
    geom_point(data=filter(target_df,X==Y),
               shape=1,size=10)
    
  print(G)
  
  G = ggplot(data=target_df) + 
    aes(x=Y,y=X, fill=(estimate-target)) + 
    geom_tile() + 
    scale_fill_distiller(type='seq', palette='Spectral') + 
    labs(title=target_author,
         x='True Age Group',y='Reported Age Group')
  
  G = G + geom_text(data=filter(target_df,X==Y),
                aes(label=X))

  print(G)
  
  PI = Pmat(alpha_u, beta_u, alpha_v, beta_v, rho_u, rho_v)
  
  U               = 0*PI
  U[upper.tri(U)] = PI[upper.tri(PI)]
  diag(U)         = -colSums(U) 
  
  V               = 0*PI
  V[lower.tri(V)] = PI[lower.tri(PI)]
  diag(V)         = -colSums(V) 
  
  dd = diag(PI)
  uu = -diag(U)
  vv = -diag(V)
  
  z = tibble(age,dd,uu,vv) %>% 
    pivot_longer(cols=dd:vv) %>% 
    filter(age %in% 60:105)
  
  G = ggplot(data=z) +
    aes(x=age,y=value,color=name) +
    geom_line(size=1) +
    theme_bw() +
    labs(title=target_author)
  
  print(G)
  
  # net imports of exposure
  xx = 60:100
  nn           = PI %*% diag(exposure)
  dimnames(nn) = list(age,age)
  nimport      = rowSums(nn[paste(xx),])/colSums(nn[,paste(xx)]) - 1
  
  # net imports of exposure
  dd           = PI %*% diag(events)
  dimnames(dd) = list(age,age)
  dimport      = rowSums(dd[paste(xx),])/colSums(dd[,paste(xx)]) - 1
  
  matplot(xx, cbind(nimport, dimport), pch=c('N','D'),
          main=target_author, type='o', col=c('green3','blue'))
  
  abline(h=0)
  
  import_df = tibble(age=xx,
               nimport,
               dimport) %>% 
        pivot_longer(cols=ends_with('import'),
                     names_to='cat',
                     values_to='net_import') 
  
G=  ggplot(data=import_df) +
    aes(x=age,y=net_import,color=cat) +
    geom_point() +
    geom_smooth(se=FALSE,span=.60, alpha=.80) +
    theme_bw() +
    geom_hline(yintercept = 0) +
    scale_color_manual(values=c('blue','green3')) +
    labs(title=target_author)
  
print(G)

  plot(xx, dimport-nimport, type='o',
       main=paste(target_author, 'mx bias'))
  abline(h=0,lty='longdash')
  abline(v=seq(60,100,5),lty='dotted')
             
  mx_true = delta[paste(xx)]/eta[paste(xx)]  
  mx_est  = rowSums(dd[paste(xx),]) / rowSums(nn[paste(xx),])
  ratio   = mx_est/mx_true
  
  plot(xx,ratio, type='o')
  abline(h=1)
  
 G= ggplot() +
    aes(x=xx,y=ratio) +
    geom_point() +
    geom_smooth(se=FALSE) +
    theme_bw() +
    scale_y_continuous(limits=c(0.5,1.2)) + 
    geom_hline(yintercept = 1) +
    labs(title=target_author,
         subtitle='Ratio of Estimated/True mx') 
  
  
 print(G)
 
  for (y in age[age>Y0]) {
    yy = paste(y)
    plot(age, PI[,yy], type='h', main=y, ylim=c(0,1))
    segments(y,0,y,PI[yy,yy], col=2, lwd=2)
    abline(v=seq(Y0,105,5)-.50, lty='dotted')
  }
  
  dev.off()

  I = diag(nrow(PI))
  
  U               = 0*PI
  U[upper.tri(U)] = PI[upper.tri(PI)]
  diag(U)         = -colSums(U) 
  
  V               = 0*PI
  V[lower.tri(V)] = PI[lower.tri(PI)]
  diag(V)         = -colSums(V) 
  
  save(PI,I,U,V, file=here('Data',paste0(target_author,'-PI-IUV.Rdata')))
  

} # for case

tab = sapply(par, round, 3) %>% 
         kable(format='latex')

write_lines(tab, file=here('Tables','nonlinear-parameter-fit-Preston-and-Bhat.tex'))
