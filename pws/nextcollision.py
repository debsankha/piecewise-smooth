from math import *
import sys

NPTS=300
x=float(sys.argv[1])
xmax=float(sys.argv[2])

F=float(sys.argv[3])
gamma=float(sys.argv[4])
m=1
sigma=1
w=1
n_val=2.3556
K1=(w*n_val/2)*(w*n_val/2.0)+gamma*gamma/4.0
w_0=sqrt(K1/m)
w_g=n_val*w/2


Ast=F/(m*sqrt((w**2-w_0**2)**2+(gamma*w)**2))

C=lambda x,v:-atan(((gamma/2*x+v)/w_g)/x)

B=lambda x,v:sqrt(x**2+((gamma/2*x+v)/w_g)**2)

def nextcol(x,v):
	c=C(x,v)
	b=B(x,v)

	res=(pi-c)/(w_g-w) if w_g>w else (pi+c)/(w-w_g) 

	nextpeak_height=(Ast+b*exp(-gamma*res/2))
	heightnow=-Ast+b*cos(c)#(Ast+b)*sin(c/2)

	if nextpeak_height>heightnow:	#max height reached at the next peak of the envelope
		tc=res
		max_height=nextpeak_height

		return max_height,tc,1
	else:				#max height reached at next peak itself, problem with calc
		tc=(2*pi)/(w+w_g)		#very rough estimate
		max_height=-Ast+b*exp(-gamma*tc/2)
		return max_height,tc,0
	
	


v=3.0
#v=float(sys.argv[3])
#vmax=float(sys.argv[4])

dx=(xmax-x)/NPTS

sys.stderr.write("Ast: "+str(Ast)+" B: "+str(B(x,v))+"\n")

while x<xmax:
	nexthei,col_time,boo=nextcol(x,v)
	print x, nexthei,boo, col_time, B(x,v), C(x,v)
	x+=dx
