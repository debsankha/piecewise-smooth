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
gamma=0.062
n=2.2
w=1
wg=n*w/2.0
k=wg**2+(gamma/2.0)**2
w0=sqrt(k)
F=0.39
T=4*pi/w

Amp=F/sqrt((w0**2-w**2)**2+(w*gamma)**2)

def M(t):
	return exp(-gamma*t/2)/wg*np.matrix([[wg*cos(wg*t)+gamma*sin(wg*t)/2, sin(wg*t)],\
						[-k*sin(wg*t), wg*cos(wg*t)-gamma*sin(wg*t)/2]])

M2T=M(T)

def G1(y):
	x=float(y[0][0])
	v=float(y[1][0])
	return M2T*np.matrix([x,-v-2*w*sqrt(Amp**2-(sigma-x)**2)]).transpose()-y

def G2(y):
	x=float(y[0][0])
	v=float(y[1][0])
	return M2T*np.matrix([x,-v+2*w*sqrt(Amp**2-(sigma-x)**2)]).transpose()-y



def Jackmap1(y):
	x=float(y[0])
	if abs(x-sigma)<Amp:
		return M2T*np.matrix([[1,0],\
			[2*w*(x-sigma)/sqrt(Amp**2-(sigma-x)**2),-1]])
	else:
		return None
			
def Jackmap2(y):
	x=float(y[0])
	if abs(x-sigma)<Amp:
		return M2T*np.matrix([[1,0],\
			[-2*w*(x-sigma)/sqrt(Amp**2-(sigma-x)**2),-1]])
	else:
		return None

def Jackeqn1(y):
	j1=Jackmap1(y)
	if j1!=None:
		return j1-np.matrix(np.ma.identity(2))
	else:
		return None	

def Jackeqn2(y):
	j1=Jackmap2(y)
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


def nraphsonsolve2(x,v,tol,maxsteps):
	y=np.matrix([x,v]).transpose()

	nsteps=0
	absdy=1
	
	while absdy>tol and nsteps<maxsteps:
		Je=Jackeqn2(y)
		if Je==None:
			return
		dy=-linalg.inv(Je)*G2(y)
		absdy=float((dy.transpose()*dy)[0][0])
		y=y+dy
		nsteps+=1
		
	
	if nsteps<maxsteps:
		return y




def plot_fp_eigvals(fmin,fmax):
	global F,Amp
	df=(fmax-fmin)/50
	f=fmax

	while f>fmin:
		F=f
		Amp=F/sqrt((w0**2-w**2)**2+(w*gamma)**2)
		for i in range(50):
			if uniform(0,1)<0.5:
				xmin=sigma-Amp
				xmax=sigma
			else:
				xmin=sigma
				xmax=sigma+Amp	

			y=nraphsonsolve1(uniform(xmin,xmax),uniform(0,8),0.001,5000)

			if y!=None:
			#	print "%dth attempt"%i
				eigs=linalg.eigvals(Jackmap1(y))
				print F,float(y[0]),float(y[1]),abs(eigs[0]),abs(eigs[1]), 0
				break
		
			y=nraphsonsolve2(uniform(xmin,xmax),uniform(0,8),0.001,5000)

			if y!=None:
			#	print "%dth attempt"%i
				eigs=linalg.eigvals(Jackmap2(y))
				print F,float(y[0]),float(y[1]),abs(eigs[0]),abs(eigs[1]), 1
				break



		f-=df		

if __name__=='__main__':
	Fgraz=sigma*F/Amp
	print Fgraz
	plot_fp_eigvals(Fgraz*1,Fgraz*1.9)
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
