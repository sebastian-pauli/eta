# Copyright 2026 by Torre Caparatta and Sebastian Pauli
#
# This program is free software: you can redistribute it and/or modify it under the terms of the 
# GNU General Public License as published by the Free Software Foundation, either version 3 
# of the License, or any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# <https://www.gnu.org/licenses/>.

# Functions for the fractional derivatives of the  Dirichlet/Euler eta function for Sage (www.sagemath.org)

RR = RealField(53)

import mpmath
tolerance = 10**-7
zetaprec = 200 
mpmath.mp.prec = zetaprec
mpnull = mpmath.mpf(0)
mpone = mpmath.mpf(1)
mptwo = mpmath.mpf(2)
mppi = mpmath.mpf(pi())
mpnegone = mpmath.mpf(-1)
mphalf = mpmath.mpf(0.5)
mpf = mpmath.mpf
mpc = mpmath.mpc
ff = mpmath.ff
rf = mpmath.rf
findroot = mpmath.findroot

def euler_eta_sub(s, alpha, m = 200, v = 20, bound = False):
    """
    alpha-th Grunwald-Letnikov fractional derivative of the Dirichlet/Euler eta function minus 1 evaluated at s.
    """
    k = alpha

    @cached_function
    def ll(k):
        lk = log(mpf(k))
        llk = log(lk)
        return lk, llk

    s = s                              #s
    mm = m                             #m
    mm1 = mpf(2*m)                          #2m
    mm2 = mm1+mpone                         #2m+1
    lm1,llm1 = ll(mm1)                      #log(2m) and log(log(2m))
    lm2, llm2 = ll(mm2)                     #log(2m+1) and log(log(2m+1))
    mms1 = mpmath.exp(lm1*(-s))             #(2m)^(-s)
    mms2 = mpmath.exp(lm2*(-s))             #(2m+1)^(-s)
    malpha = alpha                         #alpha
    lmalpha1 = mpmath.exp(llm1*malpha)      #log^(alpha)(2m)
    lmalpha2 = mpmath.exp(llm2*malpha)      #log^(alpha)(2m+1)
    lminv1 = lm1**(-1)                      #log^(-1)(2m)
    lminv2 = lm2**(-1)                      #log^(-1)(2m+1)


    @cached_function
    def stirling_non_central1_sub(n,i,s):
        if n == 1 and i == 0:
            return -s
        if n == 1 and i == 1:
            return 1
        if i == n:
            return stirling_non_central1_sub(n-1,n-1,s)
        if n >= 1 and i == 0:
            return (-s-(n-1))*stirling_non_central1_sub(n-1,0,s)
        else:   ### i>=1 and i<=n-1
            return  (-s-(n-1))*stirling_non_central1_sub(n-1,i,s)+stirling_non_central1_sub(n-1,i-1,s)

    def term_deriv_eval(alpha, k, s, m):
        # kth derivative of log^alpha x/x^s evaluated at x=m
            m = mpf(m)              #m
            malpha = alpha     #alpha
            #malpha = mpc(alpha)     #alpha
            #malpha = mpf(alpha)     #alpha
            tde = mpnull            #zero
            fact = mpone            #one
            lmpow = mpmath.exp(malpha*log(log(m))) #log(m)^(alpha)
            lminv = mpmath.exp(-1*log(log(m)))     #log(m)^(-1)
            # kth derivative of log^alpha(x)/x^s evaluated at x=m
            for i in range(0,k+1):
              tde += stirling_non_central1_sub(k,i,s)*fact*lmpow
              fact *= malpha-i
              lmpow *= lminv
            return tde*(m**(-s-k))

    def stirling_non_central_sum_sub(v):
           stir = 0
           thisff = 1
           for j in range(0,v+1):
             stir += abs(stirling_non_central1_sub(v,j,s))*thisff
             thisff *= malpha-j
           return stir
   
    sum1 = mpnull #zero
    
    ## \sum\limits_{k=2}^{m-1} -\frac{\log^\alpha 2k}{(2k)^s}+\frac{\log^\alpha (2k+1)}{(2k+1)^s}

    for k in range(1,m):   #1 to m-1
        lmk1, llmk1 = ll(2*k)       #log(2k) and log(log(2k))
        lmk2, llmk2 = ll(2*k+mpone) #log(2k+1) and log(log(2k+1))
        sum1 += -mpmath.exp(llmk1*malpha-lmk1*s)+mpmath.exp(llmk2*malpha-lmk2*s)
    

    if s == 1:
        alpha1=alpha+mpone
        sum2=-(mphalf*(1/alpha1))*(lm2**alpha1-lm1**alpha1)
    else:
        ao = malpha+mpone
        ao1 = -1*ao
        sao = s-mpone
        genic1 = mpmath.gammainc((ao),sao*lm2)
        genic2 = mpmath.gammainc((ao),sao*lm1)
        sum2 = mphalf*(genic1-genic2)*(sao**ao1)
    
    sum3 = mphalf*((-lmalpha1/(mm1**s))+(lmalpha2/(mm2**s))) 

    sum4 = mpnull
    for j in range(2,v+1,2):
           t1 = term_deriv_eval(alpha,j-1,s,mm1)
           t2 = term_deriv_eval(alpha,j-1,s,mm2)
           sum4 += mpf(mpmath.bernoulli(j))/(mpf(mpmath.fac(j)))*(mpf(2^(j-1)))*(-t1+t2)
           #print("t1 =", t1, "\n","t2 =", t2, "\n","sum4 =", sum4, "\n")
    res = ((-1)**malpha)*(sum1+sum2+sum3-sum4)
    
    return res


def euler_eta(s, alpha = 0, m = 200, v = 20):
    """
    alpha-th Grunwald-Letnikov fractional derivative of the Dirichlet/Euler eta function evaluated at s.
    m and v are the parameters of the Euler-MacLaurin summation used in its approximation.
    The algorithm used is an adaption of the algorithm for evaluating fractional derivatives of the Riemann zeta function fro,
    R. Farr, S. Pauli, and F. Saidak: Evaluating fractional derivatives of the Riemann zeta function, 
    Mathematical software—ICMS 2020, Lecture Notes in Comput. Sci., vol. 12097, Springer
    """
    z = euler_eta_sub(s,alpha)
    if alpha == 0:
        return z+mpone
    else:
        return z

def euler_eta_number_of_zeros(z1,z2,alpha):
  """
  The number of zeros counting multiplicities of the alpha-th Gruenwald-Letnikov 
  fractional derivative of the Dirichlet/Euler eta function in the rectangle given by z1 and z2.
  """
  userule=1
  useeps_rel=mpf(1e-12)
  k=mpf(alpha)
  i=mpc(I)

  z1=mpc(z1)
  z2=mpc(z2)

  dre = z2.real-z1.real
  dim = z2.imag-z1.imag

  I1im=integral_numerical(lambda x:(euler_eta(z2+x,k+1)/euler_eta(z2+x,k)).imag, 0 , -dre, rule=userule, eps_rel=useeps_rel)

  I2im=integral_numerical(lambda x:(i*euler_eta(z1+i*x,k+1)/euler_eta(z1+i*x,k)).imag, dim, 0, rule=userule, eps_rel=useeps_rel)

  I3im=integral_numerical(lambda x:(euler_eta(z1+x,k+1)/euler_eta(z1+x,k)).imag, 0,dre, rule=userule, eps_rel=useeps_rel)

  I4im=integral_numerical(lambda x:(i*euler_eta(z2+i*x,k+1)/euler_eta(z2+i*x,k)).imag, -dim, 0, rule=userule, eps_rel=useeps_rel)

  Icont = I1im[0]+I2im[0]+I3im[0]+I4im[0]
  ret = Icont/(2*RR(pi()))
  
  return ret

def euler_eta_frac_next_zero(s,alpha,tol=10**-14):
   """
   Find zero of Gruenwald-Letnikov fractional derivative eta^(alpha) starting at s
   """
   tolerance = min(tol,(abs(euler_eta((ceil(s.real)))*tol))) # smaller tolerence further right
   z = findroot(lambda x:euler_eta(x,mpf(alpha)),s,tol=mpf(tolerance))
   return z

def euler_eta_c_next_zero(s,c,tol=10**-14):
   """
   Find zero of eta-c starting at s
   """
   tol = mpf(tol)
   # adding 0.03i keeps search on right track
   z = findroot(lambda x:euler_eta(x)-c,mpc(s+0.03*I),tol=mpf(tol))
   return z

def euler_eta_c_path(start,delta=10^-2,max_real=12):
  """
  Path of zeros of eta(s)-c for c in [0,1].
  start should be a(n approximation to) zero of eta.
  """
  c = delta
  s = start+delta*I
  L = [(0,s)]
  while (c < 1) and (s.real < max_real):
    s = euler_eta_c_next_zero(s,c)
    L.append((c,s))
    c += delta
  if (s.real < 1.95):
    s = euler_eta_c_next_zero(s,1)
    L.append((1,s))
  return L

def euler_eta_frac_path(start,delta=0.05,max_real=12):
  """
  Path of zeros of Gruenwald-Letnikov fractional derivatives eta^(alpha)(s) with Re(s)<max_real.
  start should be a(n approximation to) zero of eta(s)-1.
  """
  s = start
  L = []
  alpha = delta
  while s.real<max_real:
    s = euler_eta_frac_next_zero(s,alpha)
    L.append((alpha,s))
    alpha += delta
  return L

def plot_labels(L,delta=0.01):
  """
  Label the integral derivatives in a plot of path of derivatives.
  """
  G = Graphics()
  for s in L:
    r = RR(s[0]).round()
    if abs(r-s[0]) < delta:
      G += text(str(r),(RR(s[1].real),RR(s[1].imag)),color="black")
  return G

