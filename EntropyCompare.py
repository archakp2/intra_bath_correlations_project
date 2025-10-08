import argparse
import numpy as np
import matplotlib.pyplot as plt
from two_baths_construction_for_steady_state import buildFullSetUp
from mesoscopic_lead_approach_steady_state import NESS
from von_Neumann_entropy import S_from_C
import getOQSobject
import pandas as pd

import Lindblad_EntropyTend
import mesoscopic_EntropyTrend
import csv
import importlib
importlib.reload(Lindblad_EntropyTend)





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
    
    parser.add_argument( '-o1', '--output1',
            default="",
            help='Specify output name. If not given, uses system default of Geometry+DateTime.root' )
    parser.add_argument( '-o2', '--output2',
            default="",
            help='Specify output name. If not given, uses system default of Geometry+DateTime.root' )
    
    
    args = parser.parse_args()
    
    output1 = args.output1 + ".csv"
    output2 = args.output2 + ".csv"
    
    
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
    
    SBcorrelMeso,BBcorrelMeso, Cmeso = mesoscopic_EntropyTrend.getMesoEntropies(H_mesocscopic, Ns ,bathPos , myus, betas , gammas , Ws, gb[0])


    nfsLindblad = getOQSobject.getFermiDis(H_modifiedLindBlad, betas, myus)
    Qlindblad = getOQSobject.getQ(Ns,gammas,nfsLindblad)
    SBcorrelLindblad, BBcorrelLindblad, Clindblad = Lindblad_EntropyTend.getEntropies(H_modifiedLindBlad,Qlindblad,Ns,gb[0])

    data1 = np.array([gb[0], SBcorrelLindblad, BBcorrelLindblad, SBcorrelMeso, BBcorrelMeso]).T
    output_file1 = output1
    
    # Save as space-separated text file with header
    header1 = "g SBcorrelLindblad BBcorrelLindblad SBcorrelMeso BBcorrelMeso"
    np.savetxt(output_file1, data1, header=header1, fmt="%.6e", delimiter=" ")
    
    # --- Second dataset ---
    Cmeso10 = np.array([C[0, 1] for C in Cmeso])
    Clindblad10 = np.array([C[0, 1] for C in Clindblad])
    print(Cmeso10.shape)
    
    data2 = np.array([gb[0], Cmeso10, Clindblad10]).T
    output_file2 = output2
    
    header2 = "g Cmeso Clindblad"
    np.savetxt(output_file2, data2, header=header2, fmt="%.6e", delimiter=" ")
    
    # data1 = np.array([gb[0], SBcorrelLindblad, BBcorrelLindblad, SBcorrelMeso, BBcorrelMeso])
    # data1 = data1.T
    # output_file1 = output1

    # df1 = pd.DataFrame(data1, columns=['g', 'SBcorrelLindblad',"BBcorrelLindblad","SBcorrelMeso", "BBcorrelMeso"])
    # df1.to_csv(output_file1, index=False)
    
    
    # Cmeso10 = np.array([C[0, 1] for C in Cmeso])
    # Clindblad10 = np.array([C[0, 1] for C in Clindblad])
    # print(Cmeso10.shape)
    
    # data2 = np.array([gb[0],Cmeso10,Clindblad10])

    # data2 = data2.T
    # output_file2 = output2
    # df2 = pd.DataFrame(data2, columns=['g', 'Cmeso',"Clindblad"])
    # df2.to_csv(output_file2, index=False)
    
    
    # fig1, ax1 = plt.subplots()
    
    # getOQSobject.plotCorrels(SBcorrelMeso, BBcorrelMeso, gb[0], ax=ax1, label_prefix='Meso ')
    # getOQSobject.plotCorrels(SBcorrelLindblad, BBcorrelLindblad, gb[0], ax=ax1, label_prefix='Lindblad ')
    
    # ax1.set_title("Correlations")
    # plt.legend()
    # plt.grid(True)
    
    # Second canvas: C[i,j] vs g plots
    # fig2, ax2 = plt.subplots()
    
    # getOQSobject.plotCijVSg(Cmeso, gb[0], 0, 1, ax2)
    # getOQSobject.plotCijVSg(Clindblad, gb[0], 0, 1, ax2)
    
    # ax2.set_title("C[0,1] vs g")
    # plt.grid(True)
    
    # # Show both plots
    # plt.show()
       
    return




if __name__ == "__main__":
    
    main()


# python3 EntropyCompare.py -n 4 -hs h11 0.5 h22 0 h33 0.2 ch44 -0.2 h21 1 h12 1 -gs g1 0.1 1.01 0.1 -gamma 0.01 0.01 -betas 30.0 30.0 -mus 0.0 0.0 -W -20 20 -bathPos 2 3
