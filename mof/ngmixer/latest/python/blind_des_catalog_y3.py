import sys
import hashlib
import numpy as np

code_phrase = "DES is modest and humble"


def get_factor(phrase):
    #hex number derived from code phrase
    m = hashlib.md5(phrase).hexdigest()
    #convert to decimal
    s = int(m, 16)
    # last 8 digits
    f = s%100000000
    # turn 8 digit number into value between 0 and 1
    g = f*1e-8
    #get value between 0.9 and 1.1
    return 0.9 + 0.2*g


def eta_to_g(eta1,eta2):
    eta=np.sqrt(eta1*eta1 + eta2*eta2)
    g = np.tanh(0.5*eta)
    fac = g/eta
    g1 = fac*eta1
    g2 = fac*eta2
    #correct for the zero case
    g1[eta==0.0] = 0.0
    g2[eta==0.0] = 0.0
    return g1,g2


def g_to_eta(g1, g2):
    g=np.sqrt(g1*g1 + g2*g2)
    eta = 2*np.arctanh(g)
    fac = eta/g
    eta1 = fac*g1
    eta2 = fac*g2
    
    eta1[g==0.0]=0.0
    eta2[g==0.0]=0.0
    
    return eta1,eta2

def blind_arrays(e1, e2):
    factor = get_factor(code_phrase)
    eta1, eta2 = g_to_eta(e1, e2)
    eta1 *= factor
    eta2 *= factor
    e1, e2 = eta_to_g(eta1, eta2)
    return e1, e2

def unblind_arrays(e1, e2):
    factor = get_factor(code_phrase)
    eta1, eta2 = g_to_eta(e1, e2)
    eta1 /= factor
    eta2 /=factor
    e1, e2 = eta_to_g(eta1,eta2)
    return e1, e2
