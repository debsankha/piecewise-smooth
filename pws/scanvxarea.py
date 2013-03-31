#!/usr/bin/python
from math import *
import sys
import random
#The constants needed
NPTS=300
mcpts=100000

F=0.29	#only a default value
Gamma=0.062	#only a default value
m=1	#mass
sigma=1	#switching manifold
w=1	#force freq
n_val=2.3556	#2*w_g/w
K1=(w*n_val/2)*(w*n_val/2.0)+Gamma*Gamma/4.0
w_0=sqrt(K1/m)
w_g=n_val*w/2	#sqrt(w_0**2-Gamma**2/4)



#The necessary functions

C=lambda x,v:-atan(((Gamma/2*x+v)/w_g)/x)

B=lambda x,v:sqrt(x**2+((Gamma/2*x+v)/w_g)**2)

def nextcol(x,v):
	c=C(x,v)
	b=B(x,v)

	res=(pi-c)/(w_g-w) if w_g>w else (pi+c)/(w-w_g) 

	nextpeak_height=(Ast+b*exp(-Gamma*res/2))
	heightnow=-Ast+b*cos(c)		#(Ast+b)*sin(c/2)

	if nextpeak_height>heightnow:	#max height reached at the next peak of the envelope
		tc=res
		max_height=nextpeak_height

		return max_height,tc,1
	else:				#max height reached at next peak of the high freq. signal itself, problem with calc
		tc=(2*pi)/(w+w_g)		#very rough estimate
		max_height=Ast+b*exp(-Gamma*tc/2)
		return max_height,tc,0
	

def vx_area(f,g):
	global F,Gamma,K1,w_0,Ast

	F=f
	Gamma=g
	K1=(w*n_val/2)*(w*n_val/2.0)+Gamma*Gamma/4.0
	w_0=sqrt(K1/m)
	Ast=F/(m*sqrt((w**2-w_0**2)**2+(Gamma*w)**2))


	chi=(sigma-Ast)*exp(3*Gamma*pi/(2*(w_g-w)))
	
	sys.stderr.write("chi= "+str(chi)+'\n')
	xmin=-1*chi
	xmax=min(chi,1)
	
	bigbox_area=(xmax-xmin)*2*chi*w_g

	sofar=0
	istrue=0
	
	while sofar<mcpts:
		x=random.uniform(xmin,xmax)
		ymin=-1*chi*w_g-Gamma*x/2
		ymax=chi*w_g-Gamma*x/2
		y=random.uniform(ymin,ymax)

		nexthei,tc,boo=nextcol(x,y)
		if nexthei<sigma:
			istrue+=1
#			if not boo:
#				print "bool: 0"

		sofar+=1
	return (istrue*1.0/mcpts)*bigbox_area, bigbox_area


def scan_f(fmin,fmax,g,npts):
	f=fmin
	dF=(fmax-fmin)/npts

	while f<fmax:
		area=vx_area(f,g)
		print f,area[0],area[1]
		f+=dF
		sys.stderr.write(str(f)+'\n')


#print vx_area(float(sys.argv[1]),float(sys.argv[2]))
scan_f(float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]),int(sys.argv[4]))
