import argparse
import numpy as np
import matplotlib.pyplot as plt
from two_baths_construction_for_steady_state import buildFullSetUp
from mesoscopic_lead_approach_steady_state import NESS
import spectral_functions
import time
from von_Neumann_entropy import S_from_C
import getOQSobject

import Lindblad_EntropyTend
import mesoscopic_EntropyTrend







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
    
    
    gammas = getOQSobject.getGamma(args.gamma,Ns)
    betas = getOQSobject.getBetas(args.betas,Ns)
    myus = getOQSobject.getMyu(args.mus,Ns)
    
    bathPos = mesoscopic_EntropyTrend.getBathPos(args.bathPos,Ns)
    Ws  = mesoscopic_EntropyTrend.getWs(args.W,Ns)
    
        
    nBeta = len(betas)
    nMyus = len(myus)
    nGammas = len(gammas)


    
    if nBeta != nMyus or nMyus !=nGammas or nGammas!=nBeta or nBeta + 2 != Ns:
        print("Warning : dimensions are inconsistent")
        exit(-1)    
    
    H_modifiedLindBlad = getOQSobject.modifyHamiltonian(H,gammas)
    
    H_mesocscopic = H
    
    SBcorrelMeso,BBcorrelMeso = mesoscopic_EntropyTrend.getMesoEntropies(H_mesocscopic, Ns ,bathPos , myus, betas , gammas , Ws, gb[0])


    nfsLindblad = getOQSobject.getFermiDis(H_modifiedLindBlad, betas, myus)
    Qlindblad = getOQSobject.getQ(Ns,gammas,nfsLindblad)
    SBcorrelLindblad, BBcorrelLindblad = Lindblad_EntropyTend.getEntropies(H_modifiedLindBlad,Qlindblad,Ns,gb[0])
        

    fig, ax = plt.subplots()
    
    getOQSobject.plotCorrels(SBcorrelMeso, BBcorrelMeso, gb[0], ax=ax, label_prefix='Meso ')
    getOQSobject.plotCorrels(SBcorrelLindblad, BBcorrelLindblad, gb[0], ax=ax, label_prefix='Lindblad ')
    
    plt.show()
    
    return




if __name__ == "__main__":
    
    main()


# python3 EntropyCompare.py -n 4 -hs h11 0.5 h22 0 h33 0.2 h44 -0.2 h21 1 h12 1 -gs g1 0.1 1.01 0.1 -gamma 0.01 0.01 -betas 30.0 30.0 -mus 0.0 0.0 -W -20 20 -bathPos 2 3
