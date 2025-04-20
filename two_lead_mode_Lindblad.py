import numpy as np
import matplotlib.pyplot as plt
from mesoscopic_lead_approach_steady_state import NESS
from von_Neumann_entropy import S_from_C

def nf(omega,beta,mu):

	if omega*(beta-mu)<1e2:
		return 1/(np.exp(beta*(omega-mu))+1)
	else:
		return 0
	
nf = np.vectorize(nf)

def two_lead_mode():

    #--------- bath parameters -----------------#

    beta1=30
    beta2=30
    mu1=0
    mu2=0
    #--------------------------------------------#
    
    #---------- system Hamiltonian --------------#

    # dimension
    Ns=4
    # initialising the Hamiltonian
    H_s=np.zeros((Ns, Ns), dtype=np.complex128)
    g=1
    e0=0.5
    e1=0
    eb1=0.2
    eb2=-0.2
    
    H_s_diag=np.array([e0,e1,eb1,eb2])
    H_s[np.arange(4),np.arange(4)]=H_s_diag
    H_s[0,1]=g 
    H_s[1,0]=g
    

    
    #---------------------------------------------#


    gamma1=0.01
    gamma2=0.01

    H_s[2,2]=H_s[2,2]-1j*gamma1/2
    H_s[3,3]=H_s[3,3]-1j*gamma2/2

    Q=np.zeros((Ns,Ns), dtype=np.complex128)
	
    nf1=nf(eb1,beta1,mu1)
    nf2=nf(eb2,beta2,mu2)

    Q[2,2]=gamma1*nf1
    Q[3,3]=gamma2*nf2
	

    print (Q)
	


     #--------------------------------------------------------#



    gb=np.arange(0.1,1.01,0.1)

    system_bath_correlation=np.zeros(len(gb))
    bath_modes_correlation=np.zeros(len(gb))


    for i in range(len(gb)):

        H_s[1,2]=gb[i]
        H_s[2,1]=gb[i]
        H_s[1,3]=gb[i]
        H_s[3,1]=gb[i]
		

        
        G=-1j*H_s.conj()

    
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



    print("archak  : this is sys bath corr", system_bath_correlation)
    print("archak  : this is bath bath corr", bath_modes_correlation)

    
    plt.plot(gb, system_bath_correlation, '.')
    plt.plot(gb, bath_modes_correlation, '+')


two_lead_mode()
plt.show()