import numpy as np



def n(omega,beta,mu):

    if omega*(beta-mu)<1e2:
        return 1/(np.exp(beta*(omega-mu))+1).real
    else:
    	return 0
	

n = np.vectorize(n)


def construct_bath(J,w, Gamma, beta, mu):

	gamma=np.diff(w)

	kappa=(Gamma*J[:-1]*gamma/(2*np.pi))**0.5

	w=w[:-1]

	nf=n(w,beta,mu)


	return w,nf,kappa.real,gamma.real





def build_full_set_up(H_s,bath1_pos,bath2_pos,mu1,mu2,beta1,beta2, Gamma1, Gamma2, J1,w1,J2,w2):
    
    w1,nf1,kappa1,gamma1=construct_bath(J1,w1, Gamma1, beta1, mu1)
    w2,nf2,kappa2,gamma2=construct_bath(J2,w2, Gamma2, beta2, mu2)
    
    
    Nb1=len(w1)
    Nb2=len(w2)
    Ns=H_s.shape[0]
    N_tot=Ns+Nb1+Nb2

    H_tot = np.zeros((N_tot,N_tot), dtype=np.complex128)
	
	
    Omega1=np.diag(w1-1j*gamma1/2)
    Omega2=np.diag(w2-1j*gamma2/2)
    
    
    system_bath_mat1 = np.zeros((Ns,Nb1))
    system_bath_mat2 = np.zeros((Ns,Nb2))
	
    
 
    system_bath_mat1[bath1_pos]=kappa1
    system_bath_mat2[bath2_pos]=kappa2
    
    
    H_tot[:Ns,:Ns]=H_s
    H_tot[Ns:Ns+Nb1,Ns:Ns+Nb1]=Omega1
    H_tot[Ns+Nb1:,Ns+Nb1:]=Omega2
    H_tot[:Ns,Ns:Ns+Nb1]=system_bath_mat1
    H_tot[:Ns,Ns+Nb1:]=system_bath_mat2
    H_tot[Ns:Ns+Nb1,:Ns]=system_bath_mat1.T.conj()
    H_tot[Ns+Nb1:,:Ns]=system_bath_mat2.T.conj()
 
 
 
    G=-1j*H_tot.conj()

	


    Q=np.zeros((N_tot,N_tot), dtype=np.complex128)

    Q[Ns:Ns+Nb1,Ns:Ns+Nb1]=np.diag(gamma1*nf1)
    Q[Ns+Nb1:,Ns+Nb1:]=np.diag(gamma2*nf2)
	

    return G, Q


def getOmega(ws,gammas):
    
    omega = []
    
    for w,gamma in zip(ws,gammas):
        omega.append(np.diag(w-1j*gamma/2))
    
    
    return omega 



def getSysBathMat(Ns,Nb):
    
    systemBathMats = []   
    for N in Nb:
        systemBathMat = np.zeros((Ns, N), dtype=np.double)
        systemBathMats.append(systemBathMat)
    

    return systemBathMats


def buildFullSetUp(H_s, bathPosns, mus, betas, gammas, J, wlist):

    Nb = []
    ws = []
    gammaLst = []
    kappaLst = []
    nfLst = []
    
    n = len(J)
    
    for i in range(n):
        w, nf, kappa, gamma = construct_bath(J[i], wlist[i], gammas[i], betas[i], mus[i])
        ws.append(w)
        gammaLst.append(gamma)
        kappaLst.append(kappa)  
        nfLst.append(nf)
    
    #  ques : are all w of same dimensions?
    nWs = 0
    for w in ws:
        Nbi = len(w)
        Nb.append(Nbi)


    Ns = H_s.shape[0]
    N_tot = Ns + sum(Nb)

    
    
    H_tot = np.zeros((N_tot, N_tot), dtype=np.complex128)
    
    omega = getOmega(ws, gammaLst)
    

    
    system_bath_mats = getSysBathMat(Ns, Nb)
    
    for (systemBathMat,bathPos,kappa) in zip(system_bath_mats, bathPosns,kappaLst):
        systemBathMat[int(bathPos)] = kappa
    
    # These variables must be prepared properly before:
    # Omega1, Omega2, system_bath_mat1, system_bath_mat2, gamma1, gamma2, nf1, nf2
    
        
    
    H_tot[:Ns, :Ns] = H_s
    Nstart = Ns 

    Q = np.zeros((N_tot, N_tot), dtype=np.complex128)
    
    for i in range(n):
        Ni = Nstart + Nb[i]
        
        H_tot[Nstart:Ni, Nstart:Ni] = omega[i]
        H_tot[:Ns, Nstart:Ni] = system_bath_mats[i]
        H_tot[Nstart:Ni, :Ns] = system_bath_mats[i].T.conj()
        
        Q[Nstart:Ni, Nstart:Ni] = np.diag(gammaLst[i] * nfLst[i])
        
        Nstart = Ni
    
    G = -1j * H_tot.conj()

    return G, Q













