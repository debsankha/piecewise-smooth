from math import *
import sys

F=0.393
m=1
gamma=0.062
sigma=1
w=1
n_val=2.3556
K1=(w*n_val/2)*(w*n_val/2.0)+gamma*gamma/4.0
w_0=sqrt(K1/m)
w_g=n_val*w/2


Ast=F/(m*sqrt((w**2-w_0**2)**2+(gamma*w)**2))
C=lambda v:atan(((gamma/2*(Ast+sigma)-v)/w_g)/(Ast+sigma))
B=lambda v:sqrt((Ast+sigma)**2+((gamma/2*(Ast+sigma)-v)/w_g)**2)

def t_c(v):
	res=(pi-C(v))/(w_g-w)
	if w_g<w:
		res=res*-1
	return res

def maxnext(v):
	tc=t_c(v)
#	return -Ast*cos(w*tc)+exp(-gamma*tc/2)*B(v)*cos(w_g*tc+C(v))
	return Ast+exp(-gamma*tc/2)*B(v)


v=float(sys.argv[1])
vmax=float(sys.argv[2])

dv=(vmax-v)/200

sys.stderr.write("Ast: "+str(Ast)+" B: "+str(B(v))+" t_c: "+str(t_c(v))+"\n")

while v<vmax:
	print v, maxnext(v), t_c(v), B(v)
	v+=dv

