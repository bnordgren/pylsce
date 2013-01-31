#!/usr/bin/python

import pb

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







