#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This code is under MIT license. See the License.txt file.
This module contains the functions useful to numerically solve the model

Boris Sauterey
boris.sauterey@ens.fr
"""

import numpy as np
from Metabolisme.Energy import *
from Metabolisme.Rates import *
from PBE.Balancefun import *
#from Environment import *
#from Traits import *
from scipy.stats import truncnorm

def Step_Profile(NC,X0,traits,S,gamma,T=TS,dt = 0.01):
	"""
	Computes one timestep of the profile evolution with profile at time t N
	where S is the substrate concentration vector [H,C,N,G]
	returns the profile after t+dt without touching to
	nutrients concentrations

	traits should be a vector with following order:
	[rc,Vc,ks,qmax,mg,kd,thresh,slope]
	for more information on these traits, see module traits
	"""

	## Extraction of substrate concentrations !! Memory efficient?
	for i in range(0,len(S)):
		if S[i] < 1e-100:
			S[i] = 1e-100
	H = S[0]
	C = S[1]
	N = S[2]
	G = S[3]

	## Traits extraction  !! This could be a waste of time and memory as I don't know how memory is managed  when it comes to put variables in a list and then out
	rc       = traits[0]
	Vc       = traits[1]
	Qc       = traits[2]
	ks       = traits[3]
	qmax     = traits[4]
	mg       = traits[5]
	kd       = traits[6]
	mort     = traits[7]
	thresh   = traits[8]
	slope    = traits[9]
	gmax     = traits[10]
    
	## Computing energetical values that are constant across x
	dgcat    = DeltaGcat(T,H,C,G)			# Energy that a run of metabolic reaction yields
	qcat     = QCat(dgcat,H,C,qmax,ks)		# Rate at which catabolic reaction occurs
	mreq     = Mreq(mg,dgcat) 				# Minimum rate of catabolic reaction for cell maintenance
	decay    = Decay(mreq,qcat,dgcat,kd)	# Decay rate of the cells, that only depends on energy available in the environment

	## Adjusting the delta t so that the numerical scheme remains stable
	lim = np.min([H,C]) 
	if qcat > 0:
		dt = lim/(qcat*NC*1000)
	else:
		dt = 0.01
        
	## Cell dynamics
	dgana    = DeltaGana(T,H,C,N,X0)                               # Energy requirements for anabolic reaction
	Lam      = -((dgana+dgdiss)/dgcat)                             # Metabolic coupling
	Y        = Yl(Lam)                                             # Metabolic stochiometry
	slim     = Slim([H,C,N],Y[:-2])                                # Limiting substrate
	QMet_t   = QMet(dgcat,qmax,ks,slim)                            # Metabolic rate
	qana     = QAna(dgcat,dgana,Lam,qcat,QMet_t,mreq,qmax,ks,slim) # Anabolic rate
	qcat     = qcat                                                # Catabolic rates

	new_cell = Gamma(thresh,slope,gmax,X0) 
	nNC      = NC + (new_cell - decay - mort)*NC*dt 				# First part of time derivative addition
	if nNC < 0: nNC = 0

	nX0 = (X0 + qana*dt) / (1+new_cell*dt)   
	return(nNC,nX0,qana,qcat,decay,mort,dt) # It is critical to note that qanamap, Decaymap and qcatmap are extracted from N at t and that nNc is N at t+dt

def Step_Substrates(S,Hinf,Cinf,Ninf,Ginf,QH,QC,QN,QG,NC,qana,qcat,dt,Vc):
	"""
	Computes the new S substrates vector after dt
	if several cell populations are competing, one should put as arguments:
		Nc = sum(Nci)
		qanamap = sum(qanamapi)
		qcatmap = sum(qcatmapi)
	"""
	H  = S[0]
	C  = S[1]
	N  = S[2]
	G  = S[3]

	nH = H + (QH*(Hinf-H)+(qcat*Catabolism[0]+qana*Anabolism[0])*NC)*dt
	nC = C + (QC*(Cinf-C)+(qcat*Catabolism[1]+qana*Anabolism[1])*NC)*dt
	nN = N + (QN + (qcat*Catabolism[2]+qana*Anabolism[2])*NC)*dt
	nG = G + (QG*(Ginf-G)+(qcat*Catabolism[3]+qana*Anabolism[3])*NC)*dt

	nS = np.array([nH,nC,nN,nG])
	nS[np.where(nS <= 1e-100)] = 1e-100

	return(nS)

def Step_DeadBiomass(Xo,Hinf,Cinf,Ninf,Ginf,QH,QC,QN,QG,Nc,decay,mort,Qc,X,dt,Vc):
	"""
	Computes the increase in dead biomass between t and t+dt
	"""
	return(Xo + (-0.1*Xo + (decay+mort)*Nc*(Qc+X))*dt) #Here the term with Q can be replaced with a specific biomass sedimentation flux

def Run_Profile(init,traits,Env,sig = 0.0001,Ntot0 = 10,tmax = 100,T=TS,dt = 0.01,mu=0.005):
	"""
	This function runs the profile evolution with high output volume because it computes
	and save the whole profile evolution across time tmax with initial conditions init 
	and for microbial population with traits traits

	for a single population?

	init should be [H0,C0,N0,G0]
	"""

	## Environmental conditions
	Hinf   = Env[0]
	Cinf   = Env[1]
	Ninf   = Env[2]
	Ginf   = Env[3]
	QH     = Env[4]
	QC     = Env[5]
	QN     = Env[6]
	QG     = Env[7]
    
	## Traits 
	thresh = traits[7]
	slope  = traits[8]
	gmax   = traits[9]
	Vc     = traits[1]
	Qc     = traits[2]

	## Calculation of constants over timescale of interest (here, the temperature is constant)
	DeltaG0catT = DeltaG0(T,deltaG0Cat,deltaH0Cat)
	DeltaG0anaT = DeltaG0(T,deltaG0Ana,deltaH0Ana)
    
	## Initialization
	HT    = []
	CT    = []
	NT    = []
	GT    = []
	XoT   = []
	NCT   = []
	XT    = []
	D     = []
	time  = []
	t=1

	HT.append(init[0])
	CT.append(init[1])
	NT.append(init[2])
	GT.append(init[3])
	XoT.append(init[4])
	NCT.append(init[5])
	XT.append(init[6])
	D.append(0)
	time.append(0)
	t=1
    
	while time[t-1] < tmax:   
		H  = HT[t-1]
		C  = CT[t-1]
		N  = NT[t-1]
		G  = GT[t-1]
		Xo = XoT[t-1]
		NC = NCT[t-1]
		X0 = XT[t-1]

		nNCT,nXT,qana,qcat,decay,mort,dt = Step_Profile(NC,X0,traits,[H,C,N,G],gamma,T,dt)
		NCT.append(nNCT)
		XT.append(nXT)
		D.append(decay+mort)
		nS     = Step_Substrates([H,C,N,G],Hinf,Cinf,Ninf,Ginf,QH,QC,QN,QG,NCT[t-1],qana,qcat,dt,Vc)
		HT.append(nS[0])
		CT.append(nS[1])
		NT.append(nS[2])
		GT.append(nS[3])

		nXo = Step_DeadBiomass(Xo,Hinf,Cinf,Ninf,Ginf,QH,QC,QN,QG,NCT[t-1],decay,mort,Qc,XT[t-1],dt,Vc)
		XoT.append(nXo)
		time.append(time[t-1] + dt)
		t=t+1        
	return(NCT,XT,HT,CT,NT,GT,XoT,D,time)




