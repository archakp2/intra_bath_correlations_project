import numpy as np


#von_Neumann_entropy_from_correlation_matrix
def S_from_C(C):

    # eigvals=np.abs(np.linalg.eigvalsh(C))
    eigvals=np.linalg.eigvalsh(C)
    
    # directly using the eigvals to find the entropy
    S=-np.sum( eigvals*np.log(eigvals+1e-40) + (1-eigvals)*np.log(1-eigvals+1e-40) )

    return S