import numpy as np



def n(omega,beta,mu):

	if omega*(beta-mu)<1e2:
		return 1/(np.exp(beta*(omega-mu))+1)
	else:
		return 0
	
n = np.vectorize(n)


def construct_bath(J,w, Gamma, beta, mu):
	

	gamma=np.diff(w)

	kappa=(Gamma*J[:-1]*gamma/(2*np.pi))**0.5

	w=w[:-1]

	nf=n(w,beta,mu)


	return w,nf,kappa,gamma





def build_full_set_up(H_s,bath1_pos,bath2_pos,mu1,mu2,beta1,beta2, Gamma1, Gamma2, J1,w1,J2,w2):
	
	w1,nf1,kappa1,gamma1=construct_bath(J1,w1, Gamma1, beta1, mu1)
	w2,nf2,kappa2,gamma2=construct_bath(J2,w2, Gamma2, beta2, mu2)

	Nb1=len(w1)
	Nb2=len(w2)

	#============== Build full Lyapunov matrices ================================================#	
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













