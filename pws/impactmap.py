#!/usr/bin/python

#1) find the fixed pt by NRaphson method
#2) jacobian

from math import *
import numpy as np
import numpy.linalg as linalg
import commands
import sys
from random import uniform

sigma=1
gamma=0.08
n=4
w=1
wg=n*w/2.0
k=wg**2+(gamma/2.0)**2
w0=sqrt(k)
F=0.39
T=2*pi/w

Sign=1

Amp=F/sqrt((w0**2-w**2)**2+(w*gamma)**2)
print "initamp: ",Amp

def M(t):
	return exp(-gamma*t/2)/wg*np.matrix([[wg*cos(wg*t)+gamma*sin(wg*t)/2, sin(wg*t)],\
						[-k*sin(wg*t), wg*cos(wg*t)-gamma*sin(wg*t)/2]])

M2T=M(T)

def G1(y):
	x=float(y[0][0])
	v=float(y[1][0])
	return M2T*np.matrix([x,-v-2*w*sqrt(Amp**2-(sigma-x)**2)]).transpose()-y



def Jackmap1(y):
	x=float(y[0])
	if abs(x-sigma)<Amp:
		return M2T*np.matrix([[1,0],\
			[2*w*(x-sigma)/sqrt(Amp**2-(sigma-x)**2),-1]])
	else:
		return None
			

def Jackeqn1(y):
	j1=Jackmap1(y)
	if j1!=None:
		return j1-np.matrix(np.ma.identity(2))
	else:
		return None	


def nraphsonsolve1(x,v,tol,maxsteps):
	y=np.matrix([x,v]).transpose()

	nsteps=0
	absdy=1
	
	while absdy>tol and nsteps<maxsteps:
		Je=Jackeqn1(y)
		if Je==None:
			return
		dy=-linalg.inv(Je)*G1(y)
		absdy=float((dy.transpose()*dy)[0][0])
		y=y+dy
		nsteps+=1
		
	
	if nsteps<maxsteps:
		return y


def fp_exact():
	global Sign

	sys.stderr.write("FP_exact_Amp "+Amp.__str__()+'\n')
	a=float(M2T[0,0])
	b=float(M2T[0,1])
	c=float(M2T[1,0])
	d=float(M2T[1,1])
	
	alpha=(a-d+a*d-b*c-1)/(2*b*w)

	discrim=sigma**2-(alpha**2+1)*(sigma**2-Amp**2)
	if discrim>0:
		fx=(sigma+sqrt(discrim))/(alpha**2+1)
		fv=(d-a*d+b*c)*fx/b
		
		Sign=(fx*(1-a)+b*fv)/(2*b*w*sqrt(Amp**2-(sigma-fx)**2))
		return np.matrix([fx,fv]).transpose()
#		sys.stderr.write("Sign "+str(Sign)+"\n")
	else:
		return None	



def plot_fp_eigvals(fmin,fmax,npts,newn):
	global sigma,gamma,n,w,wg,k,w0,F,T,Sign,Amp,M2T
	sigma=1
	gamma=0.08
	n=float(newn)
	w=1
	wg=n*w/2.0
	k=wg**2+(gamma/2.0)**2
	w0=sqrt(k)
	F=0.39
	T=2*pi/w
	Sign=1
	M2T=M(T)
	
	df=(fmax-fmin)/float(npts)
	f=fmin
	
	Amp=F/sqrt((w0**2-w**2)**2+(w*gamma)**2)
	print "plot_fp_eigvals_Amp: ", Amp

	while f<fmax:
#		sys.stderr.write("inside n:"+str(n)+"\n")
		F=f
		Amp=F/sqrt((w0**2-w**2)**2+(w*gamma)**2)

		yexact=fp_exact()
		if yexact==None:
			f+=df		
			continue
		eigs=linalg.eigvals(Jackmap1(yexact))
		abseigs=[abs(i) for i in eigs]
		abseigs.sort()
#		print F,float(y[0]),float(y[1]),abs(eigs[0]),abs(eigs[1]), 0
		print F,float(yexact[0]),float(yexact[1]),abseigs[0],abseigs[1]
		
		f+=df		

Fgraz=sigma*sqrt((w0**2-w**2)**2+(w*gamma)**2)
print "Fgraz: ",Fgraz
plot_fp_eigvals(Fgraz*1,Fgraz*1.9,100,float(sys.argv[1]))
#	print fp_exact()
#	for i in range(100):
#		try:
#			y=nraphsonsolve(uniform(sigma-Amp,sigma),uniform(0,3),0.001,5000)
#			print "y"
#			print y
#			print "eigenvals"
#			print Jackmap(y)
#		except:
#			pass	
#	for x in np.arange(-2,0,0.01):
#		try:
#			print x, float(G(np.matrix([x,0.02]).transpose())[0][0])
#		except:
#			pass	
