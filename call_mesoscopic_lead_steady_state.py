import numpy as np
import matplotlib.pyplot as plt
from two_baths_construction_for_steady_state import build_full_set_up
from mesoscopic_lead_approach_steady_state import NESS
import spectral_functions
import time
from von_Neumann_entropy import S_from_C



def oqs_steady_state_convergence_with_no_of_bath_modes():

    #--------- bath parameters -----------------#

    beta1=30
    beta2=30
    mu1=1
    mu2=1
    #--------------------------------------------#
    
    #---------- system Hamiltonian --------------#

    Ns=4
    H_s=np.zeros((Ns, Ns))
    g=1
    e0=0.5
    e1=0
    eb1=0.2
    eb2=-0.2
    gb1=1
    gb2=1
    H_s_diag=np.array([e0,e1,eb1,eb2])
    H_s[np.arange(4),np.arange(4)]=H_s_diag
    H_s[0,1]=g 
    H_s[1,0]=g
    H_s[1,2]=gb1
    H_s[2,1]=gb1
    H_s[1,3]=gb2
    H_s[3,1]=gb2

    print (H_s)
    #---------------------------------------------#

    #---------------- system-bath coupling params----------#

    w_min=-20
    w_max=20
    
    

    bath_spectral1=spectral_functions.constant
    

    Gamma1=0.01
    

    
    bath_spectral2=spectral_functions.constant
    
    Gamma2=0.01

    bath1_pos=2
    bath2_pos=3

    N=np.arange(20,201,20,dtype=int)
    

    
    nb1=np.zeros(len(N))
    nb2=np.zeros(len(N))
    esb=np.zeros(len(N))
    isb=np.zeros(len(N))
    b1b2_real=np.zeros(len(N))
    b1b2_imag=np.zeros(len(N))
    system_bath_correlation=np.zeros(len(N))
    bath_modes_correlation=np.zeros(len(N))

    for i in range(len(N)):

        print ('N=',N[i])

        b1_lim=2
        w=np.linspace(eb1-b1_lim,eb1+b1_lim,N[i])
        J1=bath_spectral1(w,w_min,w_max)
        w1=w


        b2_lim=2
        w=np.linspace(eb2-b2_lim, eb2+b2_lim,N[i])
        J2=bath_spectral2(w,w_min,w_max)
        w2=w

    

        #--------------------------------------------#
    
    

       

        G, Q=build_full_set_up(H_s,bath1_pos,bath2_pos,mu1,mu2,beta1,beta2, Gamma1, Gamma2, J1,w1,J2,w2)


    
        c_s_NESS=NESS(G, Q, Ns)

        nb1[i]=c_s_NESS[2,2].real
        nb2[i]=c_s_NESS[3,3].real
        esb[i]=2*(gb1*c_s_NESS[1,2].real + gb2*c_s_NESS[1,3].real)
        isb[i]=2*(gb1*c_s_NESS[1,2].imag + gb2*c_s_NESS[1,3].imag)
        b1b2_real[i]=c_s_NESS[2,3].real
        b1b2_imag[i]=c_s_NESS[2,3].imag

        c_s=c_s_NESS[:2,:2]
        c_b=c_s_NESS[2:,2:]

        S_tot=S_from_C(c_s_NESS)
        S_s=S_from_C(c_s)
        S_b=S_from_C(c_b)
        S_b1 = -nb1[i]*np.log(nb1[i]+1e-40) - (1-nb1[i])*np.log(1-nb1[i]+1e-40)
        S_b2 = -nb2[i]*np.log(nb2[i]+1e-40) - (1-nb2[i])*np.log(1-nb2[i]+1e-40)


        system_bath_correlation[i]= S_s + S_b - S_tot
        bath_modes_correlation[i]= S_b1 + S_b2 - S_b



    
    
    plt.plot(N, esb, '*')
    plt.plot(N, system_bath_correlation, '.')
    plt.plot(N, bath_modes_correlation, '+')
    plt.plot(N, b1b2_real,'^')


    
#oqs_steady_state_convergence_with_no_of_bath_modes()




def correlations_vs_coupling():

    #--------- bath parameters -----------------#

    beta1=30
    beta2=30
    mu1=0
    mu2=0
    #--------------------------------------------#
    
    #---------- system Hamiltonian --------------#

    Ns=4
    H_s=np.zeros((Ns, Ns))
    g=1
    e0=0.5
    e1=0
    eb1=0.2
    eb2=-0.2
    H_s_diag=np.array([e0,e1,eb1,eb2])
    H_s[np.arange(4),np.arange(4)]=H_s_diag
    H_s[0,1]=g 
    H_s[1,0]=g

    #gb to be specified later
    
    #---------------------------------------------#

    #---------------- system-bath coupling params----------#

    w_min=-20
    w_max=20
    
    

    bath_spectral1=spectral_functions.constant
    

    Gamma1=0.01
    

    
    bath_spectral2=spectral_functions.constant
    
    Gamma2=0.01

    bath1_pos=2
    bath2_pos=3
    
    N=200


    b1_lim=2
    w=np.linspace(eb1-b1_lim,eb1+b1_lim,N)
    J1=bath_spectral1(w,w_min,w_max)
    w1=w


    b2_lim=2
    w=np.linspace(eb2-b2_lim, eb2+b2_lim,N)
    J2=bath_spectral2(w,w_min,w_max)
    w2=w

    #--------------------------------------------------------#



    gb=np.arange(0.1,1.01,0.1)

    system_bath_correlation=np.zeros(len(gb))
    bath_modes_correlation=np.zeros(len(gb))


    for i in range(len(gb)):

        H_s[1,2]=gb[i]
        H_s[2,1]=gb[i]
        H_s[1,3]=gb[i]
        H_s[3,1]=gb[i]


        G, Q=build_full_set_up(H_s,bath1_pos,bath2_pos,mu1,mu2,beta1,beta2, Gamma1, Gamma2, J1,w1,J2,w2)


    
        c_s_NESS=NESS(G, Q, Ns)

        nb1=c_s_NESS[2,2].real
        nb2=c_s_NESS[3,3].real
        
        c_s=c_s_NESS[:2,:2]
        c_b=c_s_NESS[2:,2:]

        S_tot=S_from_C(c_s_NESS)
        S_s=S_from_C(c_s)
        S_b=S_from_C(c_b)
        S_b1 = -nb1*np.log(nb1+1e-40) - (1-nb1)*np.log(1-nb1+1e-40)
        S_b2 = -nb2*np.log(nb2+1e-40) - (1-nb2)*np.log(1-nb2+1e-40)


    
        system_bath_correlation[i]= S_s + S_b - S_tot
        bath_modes_correlation[i]= S_b1 + S_b2 - S_b
        



    
    
    plt.plot(gb, system_bath_correlation, '.')
    plt.plot(gb, bath_modes_correlation, '+')
    


correlations_vs_coupling()


plt.show()