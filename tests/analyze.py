import tests.restart as rs

def report_stomate(variable) : 
  print "Stomate Variable: %s" % variable
  for yr in range(2000,2011):
    data = rs.load_stomate(yr, variable)
    print '%d - %s' % (yr, str(rs.sumByPFT(data)))

def report_sechiba(variable) : 
  print "Sechiba Variable: %s" % variable
  for yr in range(2000,2011):
    data = rs.load_sechiba(yr, variable)
    print '%d - %s' % (yr, str(rs.sumByPFT(data)))

report_sechiba('veget_max')
report_sechiba('lai')
report_sechiba('gpp')

report_stomate('PFTpresent')
report_stomate('firelitter')
report_stomate('fuel_1hr_met')
report_stomate('fuel_1hr_str')
report_stomate('litterpart_met')
report_stomate('litterpart_str')

