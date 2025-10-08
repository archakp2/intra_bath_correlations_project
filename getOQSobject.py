import numpy as np
import matplotlib.pyplot as plt

# def plotCorrels(SB, BB, g):
    
#     fig, ax = plt.subplots()
    
    
    
#     for gi in g:
#         ax.plot(gi, SB, label='SB', marker='.')
#         ax.plot(gi, BB, label='BB', marker='+')

#         ax.xlabel('g')
#         ax.ylabel('Correlation')
#         ax.title('Correlations vs g')
#         ax.legend()
#         ax.grid(True)
#         ax.show()
    
#     return ax

def plotCorrels(SB, BB, g, ax=None, label_prefix='',logx=False,logy = False):
    if ax is None:
        fig, ax = plt.subplots()

    ax.plot(g, SB, label=f'{label_prefix}SB', marker='.')
    ax.plot(g, BB, label=f'{label_prefix}BB', marker='+')

    ax.set_xlabel('g')
    ax.set_ylabel('Correlation')
    ax.set_title('Correlations vs g')
    ax.legend()
    ax.grid(True)
    
    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')


    return ax


def nf(omega,beta,mu):

	if omega*(beta-mu)<1e2:
		return 1/(np.exp(beta*(omega-mu))+1)
	else:
		return 0
	


def getCoupling(gInfo):
    
    g_series = {}
    gList = []

    i = 0
    while i < len(gInfo):
        g_name = gInfo[i]
        start = float(gInfo[i+1])
        end = float(gInfo[i+2])
        step = float(gInfo[i+3])
        gList.append(np.arange(start, end + step, step))
        i += 4
    
    return gList




def getFixedHamiltonian(Hinfo,n):

    # create zero Hamiltonian
    H = np.zeros((n, n), dtype=np.double)

    # parse entries
    for idx in range(0, len(Hinfo), 2):
        # Hij
        entry = Hinfo[idx]
        # value in Hij
        value = float(Hinfo[idx+1])
        if entry.startswith('h') and len(entry) == 3:
            i = int(entry[1]) - 1
            j = int(entry[2]) - 1
            H[i, j] = value    

    return H


def getGamma(gammaInfo,Ns):

    gamma = np.zeros(Ns-2, dtype=np.complex128)
    for i in range(Ns-2):
        gamma[i] = gammaInfo[i]
    
    return gamma


def getBetas(betaInfo,Ns):

    betas = np.zeros(Ns-2, dtype=np.complex128)

    for i in range(Ns-2):
        betas[i] = betaInfo[i]
    
    return betas



def getMyu(myuInfo,Ns):

    myus = np.zeros(Ns-2, dtype=np.complex128)

    for i in range(Ns-2):
        myus[i] = myuInfo[i]
    
    return myus


def getFermiDis(Hs,betas,myus):
    
    n = len(betas)
    
    nfs = []
    
    for i in range(n):
        nfs.append(nf(Hs[2+i,2+i],betas[i],myus[i]))
    
    return nfs


def getQ(N,gammas,nfs):
    
    # a matrix Q is initialized here
    Q=np.zeros((N,N), dtype=np.complex128)    

    for i in range(N-2):
        Q[2+i,2+i]=gammas[i]*nfs[i]
        
    
    return Q

def modifyHamiltonian(H_s, gammas):
    
    nGammas = len(gammas)
    Hshape = H_s.shape
    
    H_modif = np.zeros(Hshape, dtype=complex)
    
    H_modif = np.array(H_s,dtype=complex)    
    gammas = np.array(gammas, dtype=np.complex128)
    
    for i in range(nGammas):
        
        H_modif[2+i,2+i] = H_s[2+i,2+i] -1j*gammas[i]/2
        
    return H_modif




def plotCijVSg(C_list, g_list, i, j, ax=None, part='real',logx = False, logy = False):
    """
    Plots C[i, j] vs g for a list of matrices and corresponding g values.

    Parameters:
    - C_list: list of N x N complex matrices
    - g_list: list of scalar values (same length)
    - i, j: indices for the matrix element to extract
    - part: 'real', 'imag', or 'abs' to select component to plot
    """
    if ax is None:
        fig, ax = plt.subplots()
    
    C_list = np.array(C_list)

    # Extract the (i, j) element from each matrix
    values = [C[i, j] for C in C_list]
    print(C_list.shape)
    # Select part of the complex number
    if part == 'real':
        values = [v.real for v in values]
    elif part == 'imag':
        values = [v.imag for v in values]
    elif part == 'abs':
        values = [abs(v) for v in values]
    else:
        raise ValueError("part must be 'real', 'imag', or 'abs'")

    print(g_list)
    print(values)

    ax.plot(g_list,values,label = f'C{i}{j} vs g' )
    
    ax.set_xlabel('g')
    ax.set_ylabel('C{i}{j}')
    ax.set_title('C{i}{j} vs g')
    ax.legend()
    ax.grid(True)
    
    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')


    return ax



# H11 can be tweeked
# interested in g1 
# there should be a region where it matches, and where it does not

# plot both on the same plot, look up for overlap

# expectation of c1c2 c2c2 c1c1 variation with g1

# in all of this there should be some regime where they match
# need to see if there MIGHT be same g value where things get intersting for all of them 

# try to understand why is this happening, they deviation? 

'''HW'''
# just try to make these plots, 
# or a program which can plot al of this easily

