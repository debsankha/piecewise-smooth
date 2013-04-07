#!/usr/bin/python

#1) find the fixed pt by NRaphson method
#2) jacobian

#Problem:
# Must ensure that the fixed pts are physical, i.e. check that the orbit is always <sigma


from math import *
import numpy as np
import numpy.linalg as linalg
import commands
import sys
from random import uniform
import commands

sigma=None
gamma=None
n=None
w=None
wg=None
k=None
w0=None
F=None
T=None

M2T=None
Signs=None

Amp=None


def M(t):
	return exp(-gamma*t/2)/wg*np.matrix([[wg*cos(wg*t)+gamma*sin(wg*t)/2, sin(wg*t)],\
						[-k*sin(wg*t), wg*cos(wg*t)-gamma*sin(wg*t)/2]])


def G(y):
	x=float(y[0][0])
	v=float(y[1][0])
	return M2T*np.matrix([x,-v-2*w*sqrt(Amp**2-(sigma-x)**2)]).transpose()-y



def Jackmap(y,pm):
	x=float(y[0])
	if abs(x-sigma)<Amp:
		return M2T*np.matrix([[1,0],\
			[2*w*pm*(sigma-x)/sqrt(Amp**2-(sigma-x)**2),-1]])
	else:
		return None
			

def Jackeqn(y):
	j1=Jackmap(y)
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
	global Signs

	a=float(M2T[0,0])
	b=float(M2T[0,1])
	c=float(M2T[1,0])
	d=float(M2T[1,1])
	
	alpha=(a-d+a*d-b*c-1)/(2*b*w)

	discrim=sigma**2-(alpha**2+1)*(sigma**2-Amp**2)
	if discrim>0:
		fx1=(sigma+sqrt(discrim))/(alpha**2+1)	#Why not -ve sign?
		fv1=(d-a*d+b*c)*fx1/b

		fx2=(sigma-sqrt(discrim))/(alpha**2+1)	#Why not -ve sign?
		fv2=(d-a*d+b*c)*fx2/b

		Signs=((fx1*(1-a)+b*fv1)/(2*b*w*sqrt(Amp**2-(sigma-fx1)**2)), (fx2*(1-a)+b*fv2)/(2*b*w*sqrt(Amp**2-(sigma-fx2)**2)))
		return (np.matrix([fx1,fv1]).transpose(),np.matrix([fx2,fv2]).transpose())
	else:
		print "Aww"
		return None	


def check_solution(x0,v0,pm):
	phasediff=atan(w*gamma/(w**2-w0**2))
	xp=sigma-x0
	vp=-pm*w*sqrt(Amp**2-xp**2)

	ph=acos(xp/Amp)

	if pm<0:
		ph=-1*ph
	
	tmin=(ph-phasediff)/w

	nexty=commands.getoutput("./hardcol.out traj %f %f %f %f %f %f %f 2>/dev/null | tail -n 1 | awk '{print $1\"\t\"$2\"\t\"$3}' "%(sigma,v0+vp,T,F,gamma,n,tmin))
	
	sys.stderr.write("yexact: %f, %f\n"%(x0,v0))


	tcol,nextx,nextv=(float(i) for i in nexty.split('\t'))
	
	nextx-=xp
	nextv-=vp
		
	sys.stderr.write("tcol: %f, nexty: %f, %f\n"%(tcol,nextx,nextv))

	if sqrt((x0-nextx)**2+(v0-nextv)**2)>0.001:
		return 0
	else:
		return 1	


	

def plot_fp_eigvals(fmin,fmax,npts):
	global Amp,F
	
	df=(fmax-fmin)/float(npts)
	F=fmin

	Amp=F/sqrt((w0**2-w**2)**2+(w*gamma)**2)

	while F<fmax:
		Amp=F/sqrt((w0**2-w**2)**2+(w*gamma)**2)

		yexacts=fp_exact()
		if yexacts==None:
			F+=df		
			continue
		
	#	abseigs=[]
		for num in [0,1]:
			pm=Signs[num]

			yexact=yexacts[num]
			(x0,v0)=(float(i) for i in yexact)
			
#			check_solution(x0,v0,pm)

#	Insert code to ckeck the fixed point does not lead to orbits crossing the wall
		
				
			eigs=linalg.eigvals(Jackmap(yexact,pm))
			abseigs=[abs(i) for i in eigs]
			abseigs.sort()

			print F, abseigs[0],abseigs[1],pm

	#	print F,abseigs[0][0],abseigs[0][1],abseigs[1][0],abseigs[1][1],float(yexacts[0][0]),float(yexacts[1][0]),float(yexacts[0][1]),float(yexacts[1][1]),pm
		F+=df	
			

def probe_stability_vs_F(nval,pval):	
	global sigma,gamma,n,w,wg,k,w0,T,M2T,Signs
	sigma=1
	gamma=0.08
	n=nval
	w=1
	wg=n*w/2.0
	k=wg**2+(gamma/2.0)**2
	w0=sqrt(k)
	T=pval*pi/w
	
	
	M2T=M(T)
	Signs=(1,1)
	
	Fgraz=sigma*sqrt((w0**2-w**2)**2+(w*gamma)**2)
	print "#Fgraz: ",Fgraz
	plot_fp_eigvals(Fgraz*1.001,Fgraz*1.9,100)

def probe_stability_vs_n(nmin,nmax):
	global sigma,gamma,n,w,wg,k,w0,T,M2T,Amp,F
	sigma=1
	gamma=0.08
	n=nmin
	w=1
	n=nmin
	
	pval=4
	Amp=sigma
	while n<nmax:
		wg=n*w/2.0
		k=wg**2+(gamma/2.0)**2
		w0=sqrt(k)
		T=pval*pi/w
		
		M2T=M(T)


		a=float(M2T[0,0])
		b=float(M2T[0,1])
		c=float(M2T[1,0])
		d=float(M2T[1,1])
	
		alpha=(a-d+a*d-b*c-1)/(2*b*w)
		fx=2*sigma/(alpha**2+1)
		fv=(d-a*d+b*c)*fx/b
			
		pm=1 if (fx*(1-a)+b*fv)*(2*b*w*sqrt(sigma**2-(sigma-fx)**2))>0 else -1
	
		Ja=Jackmap([fx,fv],pm)
		
		if Ja ==None:
			n+=0.1
			continue

		eigs=linalg.eigvals(Ja)
		abseigs=[abs(i) for i in eigs]
		abseigs.sort()
#		print F,float(y[0]),float(y[1]),abs(eigs[0]),abs(eigs[1]), 0
		print n,fx,fv,abseigs[0],abseigs[1]

		n+=0.001
	

if __name__=='__main__':
	probe_stability_vs_F(float(sys.argv[1]),int(sys.argv[2]))
#	probe_stability_vs_n(float(sys.argv[1]),int(sys.argv[2]))

