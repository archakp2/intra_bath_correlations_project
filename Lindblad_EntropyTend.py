import argparse
import numpy as np
import matplotlib.pyplot as plt
from mesoscopic_lead_approach_steady_state import NESS
from von_Neumann_entropy import S_from_C
import getOQSobject


def getEntropies(H,Q,Ns,gb):
    
    bath_modes_correlation = []
    system_bath_correlation = []
    
    for g in gb:
        
        nbs = []
        
        for i in range(Ns-2):
        
            # choosing the 'g' value here
            H[1,i+2] = g
            H[2+i,1] = g
          
        G=-1j*H.conj()
        c_s_NESS=NESS(G, Q, Ns)

        for i in range(Ns-2):
            
            nb = c_s_NESS[i+2,i+2].real
            nbs.append(nb)

        # upper 2x2 square
        c_s=c_s_NESS[:2,:2]
        
        # rest of the square
        c_b=c_s_NESS[2:,2:]

        # entropy total from given C matrix
        S_tot = S_from_C(c_s_NESS)

        # entropy of system from the upper square
        S_s = S_from_C(c_s)
        # entropy of the bath
        S_b = S_from_C(c_b)

        S_bi = []
    
        
        for i in range(Ns-2):
            S_bi.append(-nbs[i]*np.log(nbs[i]+1e-40) - (1-nbs[i])*np.log(1-nbs[i]+1e-40))

        bath_modes_correlation.append(sum(S_bi)- S_b)
        
        system_bath_correlation.append(S_s + S_b - S_tot)


    return  system_bath_correlation,bath_modes_correlation


def main():
    
    
    parser = argparse.ArgumentParser()
    

    # dimension
    parser.add_argument('-n', type=int, required=True, help='Dimension of Hamiltonian')
    
    parser.add_argument('-hs', nargs='*', help='Matrix entries, like h11 1 h23 21, where hij = H[i,j]', default=[])
    
    parser.add_argument('-gs', nargs='*', help='Coupling constants, like g1 0.1 1.1 0.1. Starting from 0.1 till 1.1 with a gap of 0.1 between each entry', default=[])

    parser.add_argument('-gamma', nargs='*', help='', default = [0.0,0.0])
    
    parser.add_argument('-betas', nargs='*', help='an array of Betas containing information about temperature of bath', default = [0.0,0.0])
    
    parser.add_argument('-mus', nargs='*', help='an array of myus containing information about chemical potential of each bath', default = [0.0,0.0])
    
    args = parser.parse_args()

    Ns = args.n
    
    H = getOQSobject.getFixedHamiltonian(args.hs,Ns)
    gb = getOQSobject.getCoupling(args.gs)
    
    
    gammas = getOQSobject.getGamma(args.gamma,Ns)
    betas = getOQSobject.getBetas(args.betas,Ns)
    myus = getOQSobject.getMyu(args.mus,Ns)
    
    
    nBeta = len(betas)
    nMyus = len(myus)
    nGammas = len(gammas)

    
    if nBeta != nMyus or nMyus !=nGammas or nGammas!=nBeta or nBeta + 2 != Ns:
        print("Warning : dimensions are inconsistent")
        exit(-1)    
    
    H_modified = getOQSobject.modifyHamiltonian(H,gammas)

    nfs = getOQSobject.getFermiDis(H_modified, betas, myus)
 
    Q = getOQSobject.getQ(Ns,gammas,nfs)

    SBcorrel,BBcorrel = getEntropies(H_modified,Q,Ns,gb[0])
    
    getOQSobject.plotCorrels(SBcorrel,BBcorrel,gb[0])
    plt.show()
    
    return




if __name__ == "__main__":
    
    main()



# python3 yourscript.py -n 4 -hs h11 1 h23 2 h34 3 -gs g1 0.1 1.1 0.1 -gamma 0.5 0.5 0.5 0.5 -betas 1.0 1.0 1.0 1.0 -mus 0.0 0.0 0.0 0.0
# python3 yourscript.py -n 4 -hs h11 1 h23 2 h34 3 -gs g1 0.1 1.1 0.1 -gamma 0.5 0.5 -betas 1.0 1.0 1.0 -mus 0.0 0.0 


# python3 Lindblad_EntropyTend.py -n 4 -hs h11 0.5 h22 0 h33 0.2 h44 -0.2 h21 1 h12 1 -gs g1 0.1 1.01 0.1 -gamma 0.01 0.01 -betas 30.0 30.0 -mus 0.0 0.0 

