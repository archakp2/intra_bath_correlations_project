import argparse
import numpy as np
import matplotlib.pyplot as plt
from two_baths_construction_for_steady_state import buildFullSetUp
from mesoscopic_lead_approach_steady_state import NESS
import spectral_functions
import time
from von_Neumann_entropy import S_from_C
import getOQSobject




def getMesoEntropies(H, Ns, bath_pos, myus, betas , gammas , Ws , gb , N = 200,blim =2):
    
    
    bath_modes_correlation = []
    system_bath_correlation = []
    J = []
    wLst = []

    specFunc = spectral_functions.constant
    nbath = Ns - 2

    system_bath_correlation = []
    bath_modes_correlation= []

    C_info  = []

    for i in range(nbath):
        
        w = np.linspace(H[2+i,2+i]-blim,H[2+i,2+i] + blim,N)
        J.append(specFunc(w,min(Ws),max(Ws)))
        wLst.append(w)
    
    for g in gb:
        
        nbs = []
        Sbs = []
        
        for i in range(Ns-2):
    
            # choosing the 'g' value here
            H[1, i + 2] = g
            H[2 + i, 1] = g
            
        G, Q = buildFullSetUp(H,bath_pos, myus,betas, gammas, J , wLst)
          
        # print(G,Q)
        # print()
        c_s_NESS=NESS(G, Q, Ns)

        for i in range(Ns-2):
            
            nb = c_s_NESS[i+2,i+2].real
            nbs.append(nb)

        # upper 2x2 square
        c_s=c_s_NESS[:2,:2]
        # rest of the square
        c_b=c_s_NESS[2:,2:]
        # print("shape of cs is ", c_b.shape)

        # Create a 4x4 zero matrix
        block_diag = np.zeros((Ns, Ns))

        # Place c_s in the top-left 2x2 block
        block_diag[:2, :2] = c_s

        # Place c_b in the bottom-right 2x2 block
        block_diag[2:, 2:] = c_b

        # Append the block diagonal array to Cs
        C_info.append(block_diag)

        # entropy total from given C matrix
        S_tot = S_from_C(c_s_NESS)

        # entropy of system from the upper square
        S_s = S_from_C(c_s)
        # entropy of the bath
        S_b = S_from_C(c_b)

        
        for i in range(Ns-2):
            Sbs.append(-nbs[i]*np.log(nbs[i]+1e-40) - (1-nbs[i])*np.log(1-nbs[i]+1e-40))

        system_bath_correlation.append(S_s + S_b - S_tot)        
        bath_modes_correlation.append(sum(Sbs)- S_b)
    
    
    
    return system_bath_correlation, bath_modes_correlation , np.array(C_info)


def getBathPos(bathPosInfo,Ns):
    

    bathPos = np.zeros(Ns-2, dtype=np.complex128)

    for i in range(Ns-2):
        bathPos[i] = bathPosInfo[i]
    
    return bathPos


def getWs(Winfo,Ns):

    W = np.zeros(2, dtype=np.complex128)

    for i in range(2):
        W[i] = Winfo[i]

    
    return W

def main():
    
    
    parser = argparse.ArgumentParser()
    

    # dimension
    parser.add_argument('-n', type=int, required=True, help='Dimension of Hamiltonian')
    
    parser.add_argument('-hs', nargs='*', help='Matrix entries, like h11 1 h23 21, where hij = H[i,j]', default=[])
    
    parser.add_argument('-gs', nargs='*', help='Coupling constants, like g1 0.1 1.1 0.1. Starting from 0.1 till 1.1 with a gap of 0.1 between each entry', default=[])

    parser.add_argument('-gamma', nargs='*', help='', default = [0,0])
    
    parser.add_argument('-betas', nargs='*', help='an array of Betas containing information about temperature of bath', default = [0,0])
    
    parser.add_argument('-mus', nargs='*', help='an array of myus containing information about chemical potential of each bath', default = [0,0])
    
    parser.add_argument('-bathPos', nargs='*', help='an array of myus containing information about chemical potential of each bath', default = [0,0])
    # has to be an int array
    
    parser.add_argument('-W', nargs='*', help='input W = [wmin, wmax]', default = [0.0,0.0])
    
    
    args = parser.parse_args()

    Ns = args.n
    N = np.arange(20,201,20,dtype=int)
    
    H = getOQSobject.getFixedHamiltonian(args.hs,Ns)
    gb = getOQSobject.getCoupling(args.gs)
    
    print(len(gb[0]))
    
    gammas = getOQSobject.getGamma(args.gamma,Ns)
    betas = getOQSobject.getBetas(args.betas,Ns)
    myus = getOQSobject.getMyu(args.mus,Ns)
    
    bathPos = getBathPos(args.bathPos,Ns)
    Ws  = getWs(args.W,Ns)
    
        
    nBeta = len(betas)
    nMyus = len(myus)
    nGammas = len(gammas)


    
    if nBeta != nMyus or nMyus !=nGammas or nGammas!=nBeta or nBeta + 2 != Ns:
        print("Warning : dimensions are inconsistent")
        exit(-1)    
    
    H_modified = H
    
    
    # C --> list of different C matrix corresponding to each gb
    SBcorrel,BBcorrel,C = getMesoEntropies(H_modified, Ns ,bathPos , myus, betas , gammas , Ws, gb[0])
    C = np.array(C)
    print(C.shape)
    Cs = np.array([c[0] for c in C])
    Cb = np.array([c[1] for c in C])
    print(Cs.shape,Cb.shape)
    getOQSobject.plotCorrels(SBcorrel,BBcorrel,gb[0])
    
    # getOQSobject.plotCijVSg(Cs,gb,0,1)
    
    plt.show()
    
    return




if __name__ == "__main__":
    
    main()


# python3 mesoscopic_EntropyTrend.py -n 4 -hs h11 0.5 h22 0 h33 0.2 h44 -0.2 h21 1 h12 1 -gs g1 0.1 1.01 0.1 -gamma 0.01 0.01 -betas 30.0 30.0 -mus 0.0 0.0 -W -20 20 -bathPos 2 3
