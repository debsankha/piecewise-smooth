from math import *

a=-0.3

f=lambda x:(a-x)**(1.5) if (x<a) else 0.7*(x-a)

x=3

dx=1
n=0

while dx>0.001 and n<100:
	xnew=f(x)
	dx=abs(xnew-x)
	print x,xnew
	print xnew, xnew
	x=xnew
	n+=1


