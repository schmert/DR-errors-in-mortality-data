#...............................................
# Carl Schmertmann
# created 17 Apr 2023
# updated 13 Mar 2024
# 
# main script for producing plots and tables
# for paper on analytics if age misreporting
# and under-registration
#...............................................

#.......................................
# housekeeping ----
# ......................................

graphics.off()
rm(list=ls())

library(tidyverse)
library(here)
library(kableExtra)
library(cowplot)

# subdirectory for plot output

output_dir = 'Plots'

# custom theme
theme_carl <- function () { 
  theme_bw(base_size=13) %+replace% 
    theme(
      axis.text        = element_text(size=15, face='bold'),
      axis.title       = element_text(size=15, face='bold'),
      legend.title     = element_text(size=13,face='bold'),
      legend.text      = element_text(size=12),
      panel.grid.major = element_line(color='grey75'),
      panel.grid.minor = element_line(color='grey90'),
      strip.text       = element_text(size=14, face='bold'),
      strip.background = element_rect(fill='grey95')
     )
}

DDhue = 'royalblue'    # color for death reg errors
CChue = 'green3'       # color for census enumeration errors
BBhue = 'red'          # color for combinations of errors
AAhue = 'orangered'    # color for Afr-Am age misreporting
CRhue = 'purple'       # color for Costa Rica age misreporting
INhue = 'turquoise'    # color for Indian age misreporting

# function to make labels for age misreporting patterns from authors' names
make_pattern= function(author) {
 factor(author,levels=c('Palloni','Preston','Bhat'),
               labels=paste(c('Costa Rica','African-American','India'),
                            'Age Misreporting'))
}  
  

#...............................................
# plot list and on/off flags ----
#...............................................

UNDERREG            = TRUE
MISREPORT           = TRUE
DERIVATIVES         = TRUE
NET_ERRORS          = TRUE
PXY_PROBS           = TRUE
COMPARE_MX_BIAS     = TRUE
COMPARE_EX_BIAS     = TRUE
COMPARATIVE_EFFECTS = TRUE
CROSSOVER           = TRUE
CROSSOVER_ERRORS    = TRUE

#...............................................
# get SP data, smooth and extrapolate mx ----
#...............................................

# read Nx_SP and Dx_SP, ages0...99 for SP males


source( here('Code', 'SP-male-2010-Nx-Dx.R'))

age = 0:119

source( here('Code','construct-smooth-logmx-schedule-for-SP.R'))

#..........................................................
# convert to a STATIONARY population at ages 0..119 ----
# .........................................................

mx = exp(logmx_smooth)
Hx = cumsum( c(0,mx))
lx = exp(-Hx)

Nx = Nx_SP[1] * head(lx,-1)
Dx = Nx * mx

names(Dx) = names(Nx) = age

# define the mortality-info function ----
source(here('Code', 'mortality-info-function.R'))

M = mortality_info(Dx, Nx)

#...................................................
# under-registration calculations and plot ----
# ..................................................

if (UNDERREG) {
  
  source(here('Code','under-registration-calculations.R'))
  
  # plots use points for the theoretical derivative calculations
  # and lines for the empirical calculations.
  # If they align then all is well 
  
  hue1 = CChue
  hue2 = DDhue

  G = ggplot(data=df) +
    geom_point(aes(x=y, y=Ceff),     color=hue1, size=2) +
    geom_line(aes(x=y,  y=Ceff_emp), color=hue1) +
    geom_point(aes(x=y, y=Deff),     color=hue2, size=2) +
    geom_line(aes(x=y,  y=Deff_emp), color=hue2) +
    labs(title='',
         x = 'Age at which omitted',
         y=  'Bias in e0 (years)') +
    geom_segment(x=0,y=0,xend=99, yend=0, color='grey40', lwd=1) +
    scale_y_continuous(limits=c(-.003, +.003)) +
    scale_x_continuous(breaks=seq(0,100,10), minor_breaks = NULL) +
    geom_point(x=0, y=.0022, color=hue2) +
    geom_segment(x=0, y=.0022, xend=0, yend=.0032, color=hue2,
                 arrow=arrow(length=unit(0.2,'cm'))) +
    geom_text (x=3, y=.003,size=6, color=hue2, 
               label=paste(round(df$Deff[1],3),'at age 0'),
               hjust=0) +
    geom_point(x=0, y=-.0022, color=hue1) +
    geom_segment(x=0, y=-.0022, xend=0, yend=-.0032, color=hue1, 
                 arrow=arrow(length=unit(0.2,'cm'))) +
    geom_text (x=3, y=-.003, size=6,color=hue1,  
               label=paste(round(df$Ceff[1],3),'at age 0'),
               hjust=0) +
    geom_text(x=60, y=+.003, size=6, label='Death Under-registration', color=hue2) +
    geom_text(x=60, y=-.003, size=6, label='Census Under-enumeration', color=hue1) +
    geom_text(x=60, y=+.0002,size=6, label='Equal Census & Death Undercounting', color='grey40') +
    theme_carl()
  
  print(G)
  
  ggsave(here(output_dir,'Fig-underregistration-effects-on-e0.pdf'), plot=G, 
         height=8, width=8, units='in')

} # if UNDERREG
#...................................................
# condensed plot for y -> x age misreporting  ----
# ..................................................

if (MISREPORT) {

  # Theoretical calculation with derivs, assuming
  # P=\diag{c}=Q=\diag{v}=I, N=n, D=d at starting point
  
  effect = function(x,y,kp=1,kq=1) {
    Tx = mean( M$Tx[paste(x+0:1)] ) 
    Ty = mean( M$Tx[paste(y+0:1)] )
  
    ix = paste(x)
    iy = paste(y)
    
    eff1 = -Ty * (kp-kq) * M$mx[iy]
    eff2 = -Tx *  (M$Nx[iy]/M$Nx[ix]) * 
                   (-kp * M$mx[ix] + kq * M$mx[iy])
      
    return(eff1 + eff2)
  }
  
  # empirical calculations of the same quantities
  
  empirical_effect = function(x,y,kp=1, kq=1) {
    ix = paste(x)
    iy = paste(y)
    
    eps = .01
    
    P = C = Q = V = diag(length(age))
    dimnames(P) = list(age,age)
    dimnames(C) = list(age,age)
    dimnames(Q) = list(age,age)
    dimnames(V) = list(age,age)
    
    P[ix,iy] = P[ix,iy] + kp*eps
    Q[ix,iy] = Q[ix,iy] + kq*eps
    
    P[iy,iy] = P[iy,iy] - kp*eps
    Q[iy,iy] = Q[iy,iy] - kq*eps
  
    nn = P %*% C %*% Nx
    dd = Q %*% V %*% Dx
    
    M2 = mortality_info(dd,nn)
    
    tmp = (M2$e0 - M$e0)/eps
    
    return(tmp)
  } 
  
  ## derivative calculations for 
  ## census and death undercounting
  ## (both theoretical and empirical [emp_])
  
  D_colors = c( '55'='lightsteelblue',
                '65'='cornflowerblue',
                '75'='royalblue',
                '85'='navy')
  
  CD_colors = c('55'='grey80',
                '65'='grey60',
                '75'='grey40',
                '85'='grey20')
  
  df = expand_grid(y = seq(55,85,10),
                   k = seq(-9,+9)) %>% 
           mutate(true_age = factor(y),
                  Dcolor   = D_colors[paste(y)],
                  CDcolor  = CD_colors[paste(y)],
                  x        = y+k,
                  Ceff     = map2_dbl(x,y, effect, kp=.01, kq=  0),
                  Deff     = map2_dbl(x,y, effect, kp=  0, kq=.01),
                  CDeff    = map2_dbl(x,y, effect, kp=.01, kq=.01),
                  emp_Ceff = map2_dbl(x,y,empirical_effect,kp=.01, kq=  0),
                  emp_Deff = map2_dbl(x,y,empirical_effect,kp=  0, kq=.01),
                  emp_CDeff = map2_dbl(x,y,empirical_effect,kp=.01, kq=.01)
                  )
  
  
  
  YL = range(df$Ceff, df$Deff, df$CDeff, na.rm=TRUE)
  
  
  G = ggplot(data=df) +
    aes(x=x, y=CDeff, color=true_age) +
    geom_line(aes(y=emp_CDeff, color=CDcolor),
              lwd=3,alpha=.80) +
    labs(x='Reported Age (x)',
         y = 'Bias in e0 (years)',
         color='True Age (y)') +
    theme_carl() +
    scale_y_continuous(limits=YL) +
    scale_x_continuous(limits=c(40,95),breaks=seq(45,95,10)) +
    geom_hline(yintercept = 0) +
    geom_line(aes(y=Deff,color=Dcolor), 
              lwd=1, lty='longdash') +
    scale_color_identity()
  
  txt_df = df %>% 
    filter(y==55) %>% 
    slice_min(x,n=1) %>% 
    select(x,y,Deff,CDeff, true_age) %>% 
    pivot_longer(cols=c(Deff,CDeff)) %>% 
    mutate(id=c('Deff'='Age Misreporting\nin Deaths Only', 
                'CDeff'='in Both\nDeaths & Census')[name])
  
  G = G +
    geom_text(data=txt_df, aes(y=value, label=id),
              color='black', hjust=0.5, nudge_y = -.00015,
              size=4.5)
  
  # add y-points
  
  pt_df = tibble(x=seq(55,85,10),
                 true_age = x,
                 CDeff=0,
                 id = c('True Age\ny=55',
                        '65',
                        '75',
                        '85'))
  
  G = G + geom_point(data=pt_df,
                 color='black',shape=1,size=6) +
      geom_text(data=pt_df, aes(label=id), 
                color='black',size=5,
                nudge_x = -1.5,
                nudge_y = .0002,
                hjust=c(1,0,0,0))
      
  
  print(G)
  
  ggsave(here(output_dir,'Fig-misreporting-effects-on-e0.pdf'), plot=G, 
         height=8, width=8, units='in')

} # if MISREPORT

#......................................
# read PI matrices ----
# .....................................
author_list = c('Palloni','Preston','Bhat')

P = vector('list', length(author_list))
names(P) = author_list

for (author in author_list) {
  load(here('Data', paste0(author,'-PI-IUV.Rdata')))
  P[[author]] = PI
}  

if (PXY_PROBS) {
  
  sel_ages = c(55,65,75,85,95)
  
  p = lapply(P, function(M) M[,paste(sel_ages)])
  
  tmp = expand_grid(author = author_list,
                    y      = sel_ages,
                    x      = 0:110) %>% 
         add_column(pxy    = unlist(p)) %>% 
         mutate(pattern=make_pattern(author)) %>% 
         group_by(author,y) %>% 
         mutate( cumul_pxy = cumsum(pxy)) %>% 
         filter( cumul_pxy > .00001, cumul_pxy < .99999)
  
  G = ggplot(data=tmp) +
    aes(x=x,y=cumul_pxy,color=pattern,shape=pattern,group=factor(y)) +
    geom_step(linewidth=1, alpha=.60) +
#    geom_point(size=2) +
    facet_wrap(~pattern) +
    geom_vline(aes(xintercept=y), lty='dotted',linewidth=.50) +
    theme_carl() +
    scale_x_continuous(limits=c(45,110)) +
#    scale_y_continuous(breaks=seq(0,1,.50), minor_breaks = NULL) +
    scale_color_manual(values=c(CRhue,AAhue,INhue)) +
    theme(legend.position = 'bottom',
          panel.grid = element_line(linewidth = 0)) +
    labs(x='Reported Age', y='Cumulative Prob')
          
  
  print(G)
  
  ggsave(filename=here(output_dir,paste0('Fig-pxy-examples.pdf')),
         height=6, width=12, units='in')
  
  
}

#...................................................
# derivative calculations for P,Q matrices  ----
# ..................................................

if (DERIVATIVES) {

  #......................................
  # plot derivatives ----
  # .....................................
  
  nn        = head(Nx,111)
  names(nn) = 0:110
  
  dd        = head(Dx,111)
  names(dd) = 0:110
  
  mm = as.vector( head(M$mx,111) )
  
  # calculate scaling factor for 
  # Delta such that
  # P = I + (0.01) * Delta
  # will corresp. to 1% of
  # ages misreported among 60+ yr olds
  
  
  # mx_deriv = function(author,kP,kQ) {
  #   
  #   PI = P[[author]]
  #   old = (0:110) > 59
  #   frac_wrong = weighted.mean((1-diag(PI))[old], nn[old])
  #   
  #   # derivs of log mx wrt Costa Rica errors
  #   
  #   Delta = (PI - diag(nrow(PI))) / frac_wrong
  #   
  #   # deriv formula
  #   part1 = as.vector( diag(mm/nn) %*% Delta %*% nn)
  #   part2 = as.vector( diag( 1/nn) %*% Delta %*% dd)
  #   
  #   
  #   diag(1/mm) %*% (-kP * part1 + kQ * part2)  
  # }

  mx_deriv = function(author,kP,kQ) {
    
    PI = P[[author]]
  
    # derivs of log mx wrt Costa Rica errors
    
    Delta = (PI - diag(nrow(PI))) 
    
    # deriv formula
    part1 = as.vector( diag(mm/nn) %*% Delta %*% nn)
    part2 = as.vector( diag( 1/nn) %*% Delta %*% dd)
    
    
    diag(1/mm) %*% (-kP * part1 + kQ * part2)  
  }
  
  
  df = tibble( kP = c(.01 ,  0 ,.01),
               kQ = c(  0 ,.01 ,.01)) %>% 
       expand_grid(author=author_list) %>% 
         mutate( eff = pmap(list(author,kP,kQ),mx_deriv),
               age = list(0:110),
               expmt = case_when(
                 kP >0 & kQ==0 ~ 'C',
                 kP==0 & kQ >0 ~ 'D',
                 kP>0  & kQ >0 ~ 'B'),
               errors = factor(expmt, 
                               levels=c('D','C','B'),
                               labels=c('Death Records',
                                        'Census Records',
                                        'Both')),
               pattern=make_pattern(author)) %>% 
       unnest(cols=c(eff,age)) %>% 
       filter(age %in% 60:100)
  
  G = ggplot(data=df) +
    aes(x=age, y=eff, color=errors, linetype=errors, size=errors, alpha=errors) +
    geom_hline(yintercept = 0, lwd=0.3) +
    geom_line( ) +
    theme_carl() +
    theme(legend.key.width = unit(0.8,'in'),
          legend.position='bottom') +
    theme(strip.text = element_text(size=12)) +
    labs(x='Age',
         y=expression(paste(Delta,'ln ',m[x])),
         color   ='Type of Age Error', 
         linetype='Type of Age Error',
         size    ='Type of Age Error',
         alpha   ='Type of Age Error') +
    scale_color_manual(values=c(DDhue,CChue,BBhue)) +
    scale_alpha_manual(values=c(1,1,.50)) +
    scale_linetype_manual(values=c('solid','dashed','solid')) +
    scale_size_manual(values=c(1, 1.5, 3)) +
    scale_x_continuous(limits=c(60,102)) +
    scale_y_continuous(labels=function(x) sprintf('%5.3f',x),
                       limits=c(-.0125,+.006)) +
    facet_wrap(~pattern)
  
  print(G)
  
  ggsave(filename=here(output_dir,
            paste0('Fig-deriv-effects-on-mx.pdf')),
         plot=G,height=8, width=15, units='in')

} # if DERIVATIVES

#...................................................
# net errors in exposure and events by age ----
# ..................................................

if (NET_ERRORS) {

  #...............................................
  # percent exp & death errors by rep. age ----
  #...............................................
  #
  nn        = head(Nx,111)
  names(nn) = 0:110
  
  dd        = head(Dx,111)
  names(dd) = 0:110
  
  mm = as.vector( head(M$mx,111) )

  a   = 0:110
  old = (a %in% 60:100)
  
  calc_net_n = function(author) {
  
    # proportional import of census reports at each age
    PI = P[[author]]
    
    net_n = ((PI - diag(nrow(PI))) %*% nn)[old]
    pct_n = 100 * net_n/nn[old]
    
    return(pct_n)
  }

  calc_net_d = function(author) {  
    PI = P[[author]]  
    # proportional import of death reports at each age
    net_d = ((PI - diag(nrow(PI))) %*% dd)[old]
    pct_d = 100 * net_d/dd[old]
    
    return(pct_d)
  }
  
  df = tibble(age=a[old]) %>% 
        expand_grid(author=author_list) %>% 
        mutate( age = list(a[old]),
                pct_n = map(author, calc_net_n),
                pct_d = map(author, calc_net_d),
                pattern=make_pattern(author)) %>% 
        unnest(cols=c(age,pct_n,pct_d)) %>% 
    pivot_longer(cols=starts_with('pct'))  

  txt_df = expand_grid(author=author_list,
                       name = unique(df$name)) %>% 
    mutate(name = factor(name, levels=c('pct_n','pct_d'),
                         labels=c('Exposure','Deaths')),
           pattern= make_pattern(author)) %>% 
    add_column(age = c(  81,87, 67, 73, 73, 83),
               value = c(+8,-5,-10,+10,+21,-11))
  
  
  
  G = ggplot(data=df) +
    aes(x=age, y=value, color=name) +
    geom_line(lwd=1.5) +
    geom_hline(yintercept = 0) +
    theme_carl() +
    labs(x='Age', 
         y = 'Net Percent Error') +
    scale_color_manual(values=c('royalblue',CChue)) +
    guides(color='none') +
    scale_y_continuous(limits=c(-25,75),breaks=seq(-20,80,20),
                       minor_breaks = NULL) +
    facet_wrap(~pattern) +
    geom_text(data=txt_df, aes(label=name),
              size=5,fontface='bold',color='black',
              hjust=0)
  

  print(G)
  
  
  ggsave(filename=here(output_dir,paste0('Fig-net-D-and-N-pct-errors-by-age.pdf')),
         height=8, width=12, units='in')
  

} # if NET_ERRORS

#...................................................
# calculate P,Q effects on estimated rates ----
# ..................................................

age  = 0:110
nage = length(age)

exposure = Nx[paste(age)]
events   = Dx[paste(age)]

names(exposure) = names(events) = age

df = tibble( author=author_list,
             pattern=make_pattern(author)) 


calc_Nxy = function(author) {
  tmp = P[[author]] %*% diag(exposure)
  dimnames(tmp) = list(age,age)
  return( tmp)
}

calc_Dxy = function(author) {
  tmp = P[[author]] %*% diag(events)
  dimnames(tmp) = list(age,age)
  return( tmp)
}

calc_ex = function(mx) {
  Hx = cumsum(c(0,mx))
  lx = exp(-Hx)
  Lx = (head(lx,-1) + tail(lx,-1))/2
  names(Lx) = age
  Tx = rev( cumsum( rev(Lx)))
  ex = Tx/head(lx,-1)
}


df = df %>% 
  mutate(age     = list(age),
         Nxy     = map(author, calc_Nxy),
         Dxy     = map(author, calc_Dxy),
         Nx      = map(Nxy, rowSums),
         Ny      = map(Nxy, colSums),
         Dx      = map(Dxy, rowSums),
         Dy      = map(Dxy, colSums),
         mx      = map2(Dx,Nx,function(d,n) d/n),
         true_mx = map2(Dy,Ny,function(d,n) d/n),
         ex      = map(mx,calc_ex),
         true_ex = map(true_mx, calc_ex)
  )

#...................................................
# compare bias in mx across ages and misreporting patterns ----
# ..................................................


if (COMPARE_MX_BIAS) {
  
  bias_df = df %>% 
    select(author,pattern,age,mx,true_mx) %>% 
    unnest(cols=age:true_mx) %>% 
    mutate(ratio = mx/true_mx) %>% 
    filter(age %in% 60:100)
  
  
  G=  ggplot(data=bias_df) +
    aes(x=age,y=ratio) +
    geom_hline(yintercept = 1, lwd=1) +
#    geom_point(size=2, color='black') +
    geom_smooth(aes(color=pattern),se=FALSE,span=.50, alpha=.80, lwd=1.5) +
    theme_carl() +
    facet_wrap(~pattern) +
    labs(x='Reported Age',y='Estimated/True Mortality Rate') +
    scale_color_manual(values=c(CRhue,AAhue,INhue)) +
    guides(color='none')
  
  
  print(G)
  
  
  ggsave(filename=here(output_dir,'Fig-compare-mx-bias.pdf'),
         width=14, height=8, units='in')
  
  
} # if COMPARE_MX_BIAS

#...................................................
# compare bias in ex across ages and misreporting patterns ----
# ..................................................

if (COMPARE_EX_BIAS) {
  
  
  ex_df = df %>% 
    select(author,pattern,mx,true_mx) %>% 
    mutate(true_ex = map(true_mx, calc_ex),
           ex      = map(mx, calc_ex),
           age     = list(age)) %>% 
    select(author,pattern,age,ex,true_ex) %>% 
    unnest(cols=c(age,ex,true_ex)) %>% 
    filter(age %in% 60:100)
  
  true_df = ex_df %>% 
    filter(author=='Palloni') %>% 
    select(age,pattern,author,ex=true_ex) %>% 
    mutate(pattern='TRUE')
  
  G = ggplot(data=ex_df) +
    aes(x=age, y=ex, color=pattern) +
    geom_point(data=true_df, color='black',size=4, 
               shape=0, stroke=1) +
    geom_line(size=2, alpha=.80) +
    theme_carl() +
    scale_y_continuous(limits=c(0,20.5), expand=c(.02,.02)) +
    scale_color_manual(values=c(CRhue,AAhue,INhue)) +
    labs(color='Pattern', x='Age',y=expression(paste(e[x],' (years)')))
  
  G = G +
    geom_text(aes(x=63.5,y=17.0,label='Costa Rica'),hjust=0,
              size=4.5,fontface='bold',color='black') +
    geom_text(aes(x=65,y=12.0,label='African-\nAmerican'),hjust=0,
              size=4.5,fontface='bold',color='black') +
    geom_text(aes(x=61.0,y=19.8,label='India Age Misreporting'),hjust=0,
              size=4.5,fontface='bold',color='black') +
    guides(color='none')
  
  print(G)
  
  ggsave(filename=here(output_dir,'Fig-compare-ex-bias.pdf'),
         width=9, height=8, units='in')
  
  
} # if COMPARE_EX_BIAS

#...................................................
# compare effects of undercounts vs. age misreports ----
# ..................................................

if (COMPARATIVE_EFFECTS) {
  
  # figure out the fraction misreported at 60+ under the
  # difft misreporting patterns
  
  old = (age > 60)
  frac_wrong = sapply(P, function(p) weighted.mean((1-diag(p))[old], exposure[old] ))

  # Compare 5 types of errors: 
  #    census underenumeration
  #    death  underregistration
  #    age misstatement (Palloni/Costa Rica)
  #    age misstatement (Preston/Afr-Am)
  #    age misstatement (Bhat/India)
  #    
  # for each error type we'll calculate ex and mx values for
  # cases with 0,1,...,25% of records with errors

  etype = c('Uncounted Deaths',
            'Unenumerated Population',
            'Costa Rica Misreporting',
            'African-American Misreporting',
            'India Misreporting')
    
  underreg_df = expand_grid(etype=etype[1], pct=0:25) %>% 
     mutate( age  = list(0:110),
             frac = 1-pct/100,
             mx = map(frac, function(f) { (f * events)/exposure }),
             ex = map(mx,calc_ex))

  underenum_df = expand_grid(etype=etype[2], pct=0:25) %>% 
    mutate( age  = list(0:110),
            frac = 1-pct/100,
            mx = map(frac, function(f) { events/(f * exposure) }),
            ex = map(mx,calc_ex))
  
  CR_df = expand_grid(etype=etype[3], pct=0:25) %>% 
    mutate( age = list(0:110),
            ev  = map(pct, function(p) {
                    k = (p/100)/frac_wrong['Palloni']
                    this_P = I + k*(P[['Palloni']] - I)
                    return(as.vector(this_P %*% events))
                  }),
            ex  = map(pct, function(p) {
                    k = (p/100)/frac_wrong['Palloni']
                    this_P = I + k*(P[['Palloni']] - I)
                    return(as.vector(this_P %*% exposure))
            }),
            mx = map2(ev,ex, function(n,d) n/d),
            ex = map(mx,calc_ex))

  AA_df = expand_grid(etype=etype[4], pct=0:25) %>% 
    mutate( age = list(0:110),
            ev  = map(pct, function(p) {
              k = (p/100)/frac_wrong['Preston']
              this_P = I + k*(P[['Preston']] - I)
              return(as.vector(this_P %*% events))
            }),
            ex  = map(pct, function(p) {
              k = (p/100)/frac_wrong['Preston']
              this_P = I + k*(P[['Preston']] - I)
              return(as.vector(this_P %*% exposure))
            }),
            mx = map2(ev,ex, function(n,d) n/d),
            ex = map(mx,calc_ex))

  IN_df = expand_grid(etype=etype[5], pct=0:25) %>% 
    mutate( age = list(0:110),
            ev  = map(pct, function(p) {
              k = (p/100)/frac_wrong['Bhat']
              this_P = I + k*(P[['Bhat']] - I)
              return(as.vector(this_P %*% events))
            }),
            ex  = map(pct, function(p) {
              k = (p/100)/frac_wrong['Bhat']
              this_P = I + k*(P[['Bhat']] - I)
              return(as.vector(this_P %*% exposure))
            }),
            mx = map2(ev,ex, function(n,d) n/d),
            ex = map(mx,calc_ex))
  

  # pick an age (eg 80) and compare biases in mx and ex across error types
  sel_age = 80
  
  big = underenum_df %>% 
         bind_rows(underreg_df) %>% 
         bind_rows(CR_df) %>% 
         bind_rows(AA_df) %>% 
         bind_rows(IN_df) %>% 
         select(etype,pct,age,mx,ex) %>% 
         unnest(cols=c(age,mx,ex)) %>% 
         group_by(etype,age) %>% 
         mutate(mx_bias = 100*(mx/mx[pct==0] - 1),
                ex_bias = 100*(ex/ex[pct==0] - 1)) %>% 
         filter(age == sel_age) %>% 
         select(age,etype,pct,ends_with('bias')) %>% 
         pivot_longer(cols=ends_with('bias'), names_to = 'var', values_to = 'bias') %>% 
         mutate(index = factor(var,
                               levels=c('ex_bias','mx_bias'),
                               labels=paste0(c('e','m'), sel_age)))
  
  
  G = ggplot(data=big) +
    aes(pct, y=bias, color=etype) +
    geom_line(lwd=1.5) +
    geom_hline(yintercept = 0) +
    facet_wrap(~index) +
    theme_carl() +
    labs(x='Fraction of Records with Errors',
         y='Percent Error') +
    scale_color_manual(values=c(AAhue,CRhue,INhue,DDhue,CChue)) +
    scale_x_continuous(minor_breaks = NULL) +
    guides(color='none') 

  

  txt_df = tribble(
    ~index, ~etype, ~pct, ~bias, ~txt,
    'e80', 'Unenumerated Population',        18, -20,  'Uncounted\nPopulation',
    'e80', 'Uncounted Deaths',               18, +21,  'Uncounted\nDeaths',
    'e80', 'Costa Rica Misreporting',        18, +2  , 'CR Age Errors',
    'e80', 'India Misreporting',             18, +5.5, 'IN Age Errors',
    'e80', 'African-American Misreporting',  18, -3  , 'AA Age Errors',
    'm80', 'Unenumerated Population',        13, +22,  'Uncounted\nPopulation',
    'm80', 'Uncounted Deaths',               13, -19,  'Uncounted\nDeaths',
    'm80', 'Costa Rica Misreporting',        18, -1.5, 'CR Age Errors',
    'm80', 'India Misreporting',             18, -4,   'IN Age Errors',
    'm80', 'African-American Misreporting',  18, +8  , 'AA Age Errors'
  )
  
  G = G + geom_text(data=txt_df, aes(label=txt), hjust=0, color='black',fontface='bold')

  
  print(G)
  
  ggsave(filename=here(output_dir,'Fig-comparative-effect-sizes.pdf'),
         width=14, height=8, units='in')
  
  
  
  
    
} # if COMPARATIVE_EFFECTS


# process Rio Grande do Norte data for crossover plots ----

# raw RN death and exposure counts for ages 0:99

Nx_RN = c(71180.4, 71149.2, 71788.7, 73387.2, 75269.7, 76233.2, 77081.4, 
          78049.3, 80318.8, 83910.6, 86584.6, 88883.3, 89075.5, 90808, 
          93001.6, 93757.3, 92742, 89670.3, 86802.5, 86327.5, 88445.9, 
          92617.3, 94278.4, 94222.9, 91841.7, 88514.9, 86260.4, 85467.3, 
          85471.1, 84315.7, 80530.5, 77215.5, 73107.5, 70973.3, 69789.1, 
          68702.8, 67269.6, 64765.6, 62798.8, 63555.9, 63405.1, 63716.7, 
          61150.5, 60256.9, 59496.4, 58898.8, 56979.2, 54037.9, 51601.5, 
          49532.8, 45704, 43034.6, 40232.4, 39330.3, 37263.4, 35130.5, 
          32778, 30864.8, 30149.5, 31208.9, 30654, 29934.5, 27559.3, 26387.8, 
          24942.2, 23016.4, 21173.3, 19860.7, 19555, 19943.5, 18778.5, 
          17531, 15533.6, 14845, 13941.1, 12819.4, 11064.8, 9831.8, 9232.2, 
          9420.5, 8943.8, 8367.2, 7687.7, 7345.9, 6959, 6307.5, 5441.8, 
          4523.7, 3797, 3295, 2805.8, 2287.4, 1780.3, 1356.5, 1100.4, 911.2, 
          730.8, 545.5, 416.1, 300.9)

Dx_RN = c(1089L, 73L, 36L, 35L, 18L, 17L, 21L, 31L, 28L, 13L, 29L, 25L, 
          33L, 44L, 51L, 85L, 121L, 160L, 179L, 196L, 244L, 244L, 220L, 
          230L, 219L, 219L, 212L, 215L, 229L, 243L, 202L, 200L, 206L, 205L, 
          220L, 218L, 209L, 195L, 224L, 216L, 215L, 232L, 254L, 249L, 278L, 
          256L, 250L, 281L, 319L, 299L, 287L, 310L, 289L, 351L, 308L, 317L, 
          310L, 301L, 331L, 353L, 363L, 413L, 349L, 424L, 416L, 425L, 447L, 
          397L, 366L, 440L, 469L, 483L, 464L, 474L, 451L, 495L, 491L, 416L, 
          447L, 503L, 506L, 520L, 589L, 561L, 572L, 612L, 533L, 539L, 483L, 
          427L, 404L, 352L, 292L, 252L, 214L, 185L, 148L, 134L, 107L, 74L)

names(Nx_RN) = names(Dx_RN) = seq(Nx_RN)-1

# extrapolate from (0:99) to (0:110), using the 
# reported N[99] and the reported [99,100) mortality rate
mm = Dx_RN/Nx_RN

extended_Nx_RN = Nx_RN['99'] * exp(-mm['99'] * (100:110 - 99))
extended_Dx_RN = mm['99']    * extended_Nx_RN

Nx_RN = c(Nx_RN, extended_Nx_RN)
Dx_RN = c(Dx_RN, extended_Dx_RN)

names(Nx_RN) = names(Dx_RN) = seq(Nx_RN)-1

mx_RN = Dx_RN/Nx_RN


if (CROSSOVER) {
  
  # illustrate crossover ----
  
  mx_SP = exp(logmx_smooth)
  
  df = expand_grid( UF=c('SP','RN'),
                    age = 0:99) %>% 
    add_column( mx = c(head(mx_SP,100), head(mx_RN,100))) %>% 
    filter(age %in% 30:100)
  
  hue1 = 'orange'
  hue2 = 'purple'
  
  G = ggplot(data=df) +
    aes(x=age, y=mx, color=UF) +
    geom_line(lwd=1.5) +
    scale_y_log10(breaks=c(.001,.010,.100),
                  minor_breaks=NULL,
                  limits=c(.001,.350)) +
    scale_x_continuous(breaks=seq(30,100,10),
                       minor_breaks = NULL) +
    theme_carl() +
    labs(color='State', y='mx') +
    scale_color_manual(values=c(hue1,hue2)) +
    guides(color='none') +
    geom_text(aes(x=76,y=.025,label='Rio Grande\ndo Norte (RN)'),
              inherit.aes = FALSE, hjust=0, color=hue1, size=6) +
    geom_text(aes(x=57,y=.060,label='São Paulo (SP)'),
              inherit.aes = FALSE, hjust=0, color=hue2,size=6)
  
  print(G)
  
  ggsave(filename=here(output_dir,'Fig-SP-RN-crossover.pdf'),
         height=8, width=8, units='in')
  
} # if CROSSOVER

# errors that could generate a RN-SP crossover ----
 

if (CROSSOVER_ERRORS) {

  D  = diff(diag(nage), diff=1)
  DD = crossprod(D)
  
  # approx_solve finds an approximate solution to Ax=b,
  # adding the condition that x should have small-ish
  # first differences. It finds good solutions for
  # true exposure and deaths even when the misreporting
  # matrix (as in Bhat) is nearly singular
  
  approx_solve = function(A,b,lambda=1) {
    x =solve(lambda*DD + crossprod(A), crossprod(A,b) )
    return(as.vector(x))
  }
  
  target_mx = exp(logmx_smooth[paste(0:110)])
  
  case_df = expand_grid(author=author_list,
                        census_complete = c(FALSE,TRUE),
                        age_errors      = c(FALSE,TRUE))
  
  calc_cv = function(author,census_complete, age_errors, var='c') {
    
    if (age_errors) {
      this_P = this_Q = P[[author]]
    } else {
      this_P = this_Q = I
    }  
    
    Qinv_d    = approx_solve(this_Q, Dx_RN)
    Pinv_n    = approx_solve(this_P, Nx_RN)

    R = diag( Qinv_d / (target_mx * Pinv_n))
    
    oneA = rep(1, 111)  
    if (census_complete) {
      cc = oneA
    } else {
      cc = solve( I + crossprod(R), (I + t(R)) %*% oneA)
    }
    
    vv = R %*% cc
    
    if (var == 'c') return(as.vector(cc))
    if (var == 'v') return(as.vector(vv))
  }
  
  calc_eta = function(author,c,age_errors) {
    
    if (age_errors) {
      this_P = this_Q = P[[author]]
    } else {
      this_P = this_Q = I
    }  
    
    eta  = approx_solve(this_P, Nx_RN) / c

    return(as.vector(eta))
  }
  
  calc_delta = function(author,v,age_errors) {
    
    if (age_errors) {
      this_P = this_Q = P[[author]]
    } else {
      this_P = this_Q = I
    }  
    
    delta  = approx_solve(this_Q, Dx_RN) / v
  
    return(as.vector(delta))
  }
  
  case_df = case_df %>% 
    mutate( age=list(0:110),
            c  = pmap(list(author,census_complete,age_errors),
                      calc_cv, var='c'),
            v  = pmap(list(author,census_complete,age_errors),
                      calc_cv, var='v'),
            eta   = pmap(list(author,c,age_errors),calc_eta),
            delta = pmap(list(author,v,age_errors),calc_delta))
  
  
  long_df = case_df %>% 
    unnest(cols=c(age,c,v)) %>% 
    pivot_longer(cols=c('c','v'), names_to = 'var') %>% 
    mutate(pattern=make_pattern(author),
           census = factor(census_complete,
                           levels=c(TRUE, FALSE),
                           labels=c('No Enumeration Errors','With Census Enumeration Errors')),
           categ = factor(var,
                          levels=c('v','c'),
                          labels=c('Death Reg','Census Enum')))
  
  # smoothed only, without points
  
  G = ggplot(data=long_df) +
    aes(x=age,y=value*100, color=categ,
        lty=age_errors, size=categ) +
    geom_smooth(se=FALSE, span=.5) +
    facet_grid(census~pattern, scales='fixed') +
    theme_carl() +
    theme(legend.position = 'bottom') +
    geom_hline(yintercept = 100) +
    geom_hline(yintercept = 0, lty='dotted') +
    scale_x_continuous(limits=c(40,100)) +
    scale_y_continuous(limits=c(60,120),breaks=seq(60,120,20)) +
    scale_color_manual(values=c(DDhue,CChue)) +
    scale_size_discrete(range=c(0.6,1.2)) +
    labs(y='Coverage Level (%)', color='Coverage',
         lty='Age Errors', size='Coverage')
  
  print(G)
  
  ggsave(filename=here(output_dir,'Fig-errors-for-RN-SP-crossover.pdf'),
         height=8, width=14, units='in')
  
  
  # check whether we have really reproduced the RN
  # rates
  
  tmp = case_df %>% 
    select(author:age, eta,delta) %>% 
    unnest(cols=age:delta) %>% 
    mutate(mx = delta/eta)
  
  target_df = tibble(age=0:110, mx=target_mx, author='TARGET')
  
  raw_df = tibble(age=0:110, mx = Dx_RN/Nx_RN, author='RAW')
  
  
  ggplot(data=tmp) +
    aes(x=age,y=mx, color=interaction(author,census_complete,age_errors),
        shape=author) +
    geom_point(size=2,alpha=.60) +
    scale_y_log10() +
    scale_x_continuous(limits=c(60,110)) +
    theme_bw() +
    scale_color_manual(values=rainbow(12)) +
    geom_line(data=target_df,lwd=4,color='steelblue2',
              alpha=.30) +
    geom_line(data=raw_df,color='red',lwd=2,
              alpha=.30) +
    geom_line( data=filter(tmp,author=='Bhat',census_complete,age_errors)) 
  
  
} # if CROSSOVER_ERRORS
