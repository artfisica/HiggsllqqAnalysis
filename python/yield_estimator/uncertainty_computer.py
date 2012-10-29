# takes care of efficiency computation
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

def sum_sq(the_list):
  result = 0

  for an_element in the_list:
    result += an_element * an_element

  return sqrt(result)


def get_bayes_eff_unc(num, den):
# returns the bayesian uncertainty over the num/den ratio
  result = [0, 0]

  h_num = TH1D('h_num', '', 1, 0, 1)
  h_den = TH1D('h_den', '', 1, 0, 1)

  h_num.SetBinContent(1, num)
  h_den.SetBinContent(1, den)

  h_eff = TGraphAsymmErrors()
  h_eff.BayesDivide(h_num, h_den)

  result[0] = h_eff.GetErrorYlow(0)
  result[1] = h_eff.GetErrorYhigh(0)

  h_num.Delete()
  h_den.Delete()
  h_eff.Delete()

  return result

def get_bayes_eff_unc_relative(num, den):
  result = get_bayes_eff_unc(num, den)

  if (den != 0):
    if (num != 0):
      return (result[0]+result[1])*0.5 * den/num
    else:
      return 0.
  else:
    print 'ERROR: called get_bayes_eff_unc_relative with num=%lf den=%lf'
    return 0.
