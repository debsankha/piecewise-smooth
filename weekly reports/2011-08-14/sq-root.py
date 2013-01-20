from math import *
import sys
from random import uniform

def f(x,mu,nu):
	if x<mu:
		return sqrt(mu-x)+mu*nu
	else:
		return nu*x	

nu=float(sys.argv[1])

mu=0.00001
dm=0.00001
while mu<5:
	i=0
	x=uniform(-10,10)
	while i<400:
		x=f(x,mu,nu)
		i=i+1
	i=0
	while i<100:
		x=f(x,mu,nu)
		print mu,x
		i=i+1
	mu+=dm	
	
