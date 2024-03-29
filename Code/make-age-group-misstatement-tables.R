#.....................................................................
# Male misstatement probabilities for India from Table 3 of 
# 
#  Bhat, Estimating Transition Probabilities of Age Misstatement
#  Demography, Vol. 27, No. 1 (Feb., 1990), pp. 149-163 (15 pages)
#  https://doi.org/10.2307/2061559
#
# African-American females in Table 2 of
# 
# PRESTON, SAMUEL & Elo, Irma. (1999). 
# Effects of Age Misreporting on Mortality Estimates at Older Ages. 
# Population Studies. 53. 165-177. 10.1080/00324720308075. 
# 
# Y is the TRUE age group
# X is the REPORTED age group
# frac is P(X|Y)
#.....................................................................

library(tidyverse)

Bhat_table = tribble(
  ~Y,~X,~frac,
  0,0,.898,
  0,5,.102,
  5,0,.053,
  5,5,.812,
  5,10,.135,
  10,5,.138,
  10,10,.721,
  10,15,.140,
  15,10,.213,
  15,15,.578,
  15,20,.209,
  20,15,.187,
  20,20,.520,
  20,25,.293,
  25,15,.001,
  25,20,.176,
  25,25,.477,
  25,30,.344,
  25,35,.002,
  30,20,.001,
  30,25,.174,
  30,30,.418,
  30,35,.399,
  30,40,.007,
  35,25,.006,
  35,30,.163,
  35,35,.392,
  35,40,.408,
  35,45,.031,
  40,30,.012,
  40,35,.171,
  40,40,.374,
  40,45,.379,
  40,50,.063,
  45,35,.021,
  45,40,.183,
  45,45,.341,
  45,50,.386,
  45,55,.069,
  50,40,.035,
  50,45,.192,
  50,50,.360,
  50,55,.283,
  50,60,.131,
  55,45,.048,
  55,50,.212,
  55,55,.256,
  55,60,.377,
  55,65,.107,
  60,50,.069,
  60,55,.164,
  60,60,.352,
  60,65,.251,
  60,70,.165,
  65,50,.001,
  65,55,.067,
  65,60,.242,
  65,65,.241,
  65,70,.323,
  65,75,.122,
  65,80,.004,
  70,55,.004,
  70,60,.107,
  70,65,.165,
  70,70,.297,
  70,75,.196,
  70,80,.218,
  70,85,.014,
  75,60,.011,
  75,65,.079,
  75,70,.204,
  75,75,.174,
  75,80,.311,
  75,85,.222,
  80,65,.013,
  80,70,.107,
  80,75,.126,
  80,80,.283,
  80,85,.472,
  85,70,.026,
  85,75,.074,
  85,80,.218,
  85,85,.682
) 

# target table from Preston 1999, Table 2, females

Preston_table = tribble(
  ~Y,~X,~frac,
  0,0,1,
  5,5,1,
  10,10,1,
  15,15,1,
  20,20,1,
  25,25,1,
  30,30,1,
  35,35,1,
  40,40,.9987,
  40,45,.0003,
  45,40,.0127,
  45,45,.9867,
  45,50,.0005,
  50,40,.0003,
  50,45,.0211,
  50,50,.9776,
  50,55,.0010,
  55,45,.0216,
  55,50,.0450,
  55,55,.9298,
  55,60,.0035,
  60,50,.0429,
  60,55,.0690,
  60,60,.8820,
  60,65,.0061,
  65,55,.0642,
  65,60,.0929,
  65,65,.8342,
  65,70,.0087,
  70,60,.0642,
  70,65,.1088,
  70,70,.8131,
  70,75,.0102,
  70,80,.0037,
  75,65,.0463,
  75,70,.1331,
  75,75,.8060,
  75,80,.0093,
  75,85,.0053,
  80,70,.0672,
  80,75,.1350,
  80,80,.7757,
  80,85,.0176,
  80,90,.0045,
  85,75,.0923,
  85,80,.1429,
  85,85,.6891,
  85,90,.0593,
  85,95,.0164,
  90,80,.0923,
  90,85,.1429,
  90,90,.6891,
  90,95,.0593,
  90,100,.0164,
  95,85,.1179,
  95,90,.0727,
  95,95,.7507,
  95,100,.0587,
  100,100,1
) 

# tweak so that columns sum to 1 (keep proportions)
Preston_table = Preston_table %>% 
  group_by(Y) %>% 
  mutate(frac=prop.table(frac)) %>% 
  ungroup()

