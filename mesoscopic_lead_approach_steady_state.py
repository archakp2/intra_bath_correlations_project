import numpy as np
from scipy.linalg import solve_continuous_lyapunov

# Gx + xG_dagger  = Q

def NESS(G, Q, Ns):

    c_tot_NESS=solve_continuous_lyapunov(G,Q)

    # why do we  even have this
    # here we are solving CGamma_dagger + GammaC = Q
    # and solving it for C
    # Gamma = M - iC

    '''Why have we taken M == 0'''
        
    c_s_NESS = c_tot_NESS[:Ns,:Ns]

    return c_s_NESS



