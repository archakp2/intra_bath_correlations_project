from __future__ import division
import numpy as np


def constant(w,w_min,w_max):
	
	Gamma=1
	
	if w>=w_min and w<=w_max:
		
		return Gamma
	else:
		
		return 0
	
constant=np.vectorize(constant)


def semi_circle(w,w_min,w_max):
	
	
	tB=(w_max-w_min)/4
	w0=(w_max+w_min)/2
	gamma=(tB/2)**0.5
	
	return ((2*gamma**2)/tB)*(1-((w+w0)/(2*tB))**2)**0.5
	
semi_circle=np.vectorize(semi_circle)


def Lorentzian(w,lamda, e, w_min,w_max):
	
	
	
	if w>=w_min and w<=w_max:
		
		return lamda/((w-e)**2+lamda**2)
	else:
		
		return 0

Lorentzian=np.vectorize(Lorentzian)	




def Lorentzian1(w,w_min,w_max):
	
	lamda=0.1
	e3=1
	
	if w>=w_min and w<=w_max:
		
		return lamda/((w-e3)**2+lamda**2)
	else:
		
		return 0

Lorentzian1=np.vectorize(Lorentzian1)	


def Lorentzian2(w,w_min,w_max):
	
	lamda=0.1
	e4=-1
	
	if w>=w_min and w<=w_max:
		
		return lamda/((w-e4)**2+lamda**2)
	else:
		
		return 0

Lorentzian2=np.vectorize(Lorentzian2)


def Ohmic_type_Gaussian_cutoff(w,s,w_0,w_max):
	
	
	return (w/w_0)**s * np.exp(-(w/w_max)**2)
	
