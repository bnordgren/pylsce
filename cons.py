#!/usr/bin/python

import pb
import numpy as np


#PFTs in ORCHIDEE
#    1 - Bare soil
#    2 - tropical  broad-leaved evergreen
#    3 - tropical  broad-leaved raingreen
#    4 - temperate needleleaf   evergreen
#    5 - temperate broad-leaved evergreen
#    6 - temperate broad-leaved summergreen
#    7 - boreal    needleleaf   evergreen
#    8 - boreal    broad-leaved summergreen
#    9 - boreal    needleleaf   summergreen
#   10 -           C3           grass
#   11 -           C4           grass
#   12 -           C3           agriculture
#   13 -           C4           agriculture

pftdic= \
{\
1: 'Bare soil',\
2: 'tropical broad-leaved evergreen',\
3: 'tropical broad-leaved raingreen',\
4: 'temperate needleleaf evergreen',\
5: 'temperate broad-leaved evergreen',\
6: 'temperate broad-leaved summergreen',\
7: 'boreal needleleaf evergreen',\
8: 'boreal broad-leaved summergreen',\
9: 'boreal needleleaf summergreen',\
10: 'C3 grass',\
11: 'C4 grass',\
12: 'C3 agriculture',\
13: 'C4 agriculture',\
}

pftlist = [pftdic[num+1] for num in range(13)]

cal = pb.calendar()
cal.get_month_doy()
index_first_day_month_noleap = cal.index_first_day_month_noleap
index_first_day_month_leap = cal.index_first_day_month_leap

class levels(object):
    biomass_gC_m2 = [0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000, 35000, 40000]
    biomass_kgC_m2 = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 35.0, 40.0]
    biomass_tC_ha = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 350, 400]





