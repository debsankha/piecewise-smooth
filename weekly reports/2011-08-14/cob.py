from math import sqrt
import sys

mu=0.1
nu=0.5

def iter(x):
	global mu,nu
	if x<mu:
		return sqrt(mu-x)+mu*nu
	else:
		return nu*x	

x=float(sys.argv[1])
print x,0,'\n',x,

if (x<mu):
	x=iter(x)
	print x, '\n', x, x, '\n', x,
	i=1
else:
	i=0
while mu<x or x<mu*nu:
	x = iter(x)
	print x, '\n', x, x, '\n', x,
	i+=1

sys.stderr.write("i=%d\n"%i)

while i<20:
	x=iter(x)
	print x, '\n', x, x, '\n', x,
	i+=1
	
print iter(x)
