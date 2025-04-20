import numpy as np
from scipy.linalg import solve_continuous_lyapunov

# Gx + xG_dagger  = Q

def NESS(G, Q, Ns):

    c_tot_NESS=solve_continuous_lyapunov(G,Q)

    c_s_NESS = c_tot_NESS[:Ns,:Ns]

    return c_s_NESS



