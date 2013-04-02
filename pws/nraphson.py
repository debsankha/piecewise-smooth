from math import *
import numpy as np
import numpy.linalg as linalg

sigma=1
gamma=0.062
n=2.3556
w=1.0
wg=n*w/2.0
k=wg**2+(gamma/2.0)**2
w0=sqrt(k)
F=0.54
T=4*pi/w

def M(t):
	return exp(-gamma*t/2)/wg*np.matrix([[wg*cos(wg*t)+gamma*sin(wg*t)/2, sin(wg*t)],\
						[-k*sin(wg*t), wg*cos(wg*t)-gamma*sin(wg*t)/2]])
def Mdot(t):
	return -gamma/2*exp(-gamma*t/2)/wg*np.matrix([[wg*cos(wg*t)+gamma*sin(wg*t)/2, sin(wg*t)],\
						[-k*sin(wg*t), wg*cos(wg*t)-gamma*sin(wg*t)/2]])+\
			exp(-gamma*t/2)*np.matrix([[-wg*sin(wg*t)+gamma*cos(wg*t)/2, cos(wg*t)],\
						[-k*cos(wg*t), -wg*sin(wg*t)-gamma*cos(wg*t)/2]])

def xp(t):
	return F/sqrt((w0**2-w**2)**2+(w*gamma)**2)*cos(w*t+atan(w*gamma/(w**2-w0**2)))

def vp(t):
	return -F*w/sqrt((w0**2-w**2)**2+(w*gamma)**2)*sin(w*t+atan(w*gamma/(w**2-w0**2)))

def ap(t):
	return -w**2*F/sqrt((w0**2-w**2)**2+(w*gamma)**2)*cos(w*t+atan(w*gamma/(w**2-w0**2)))

def phi(tf,ti,x0,v0,Mdt):
	xpti=xp(ti)
	vpti=vp(ti)
	xptf=xp(tf)
	vptf=vp(tf)

	return np.matrix([[xptf],[vptf]])+Mdt*np.matrix([[x0-xpti],[v0-vpti]])

def delphideltf(tf,ti,x0,v0,Mdt,Mdotdt):
	xpti=xp(ti)
	vpti=vp(ti)
	apti=ap(ti)
	vptf=vp(tf)
	aptf=ap(tf)

	parta=np.matrix([vptf,aptf]).transpose()
	partb=Mdotdt*np.matrix([x0-xpti,v0-vpti]).transpose()

	return parta+partb

def delphidelti(tf,ti,x0,v0,Mdt,Mdotdt):
	xpti=xp(ti)
	vpti=vp(ti)
	apti=ap(ti)
	return -Mdotdt*np.matrix([x0-xpti,v0-vpti]).transpose()-Mdt*np.matrix([vpti,apti]).transpose()



def J1(y):
	[x0,v0,x1,v1,t1]=[float(i[0]) for i in y]
	
#	T=2*pi/w
	Mdt=M(t1)
	Mdotdt=Mdot(t1)
	row12=np.ma.hstack((-Mdt,np.ma.identity(2),-delphideltf(t1,0,x0,v0,Mdt,Mdotdt)))

	Mdt=M(T-t1)
	Mdotdt=Mdot(T-t1)
	row34=np.ma.hstack((np.ma.identity(2),-Mdt*np.matrix([[1,0],[0,-1]]),-delphidelti(T,t1,x1,-v1,Mdt,Mdotdt)))

	row5=np.matrix([0,0,1,0,0])

	return np.matrix(np.ma.vstack((row12,row34,row5)))


def G(y):
#	T=2*pi/w
	[x0,v0,x1,v1,t1]=[float(i[0]) for i in y]
	Mdt=M(t1)
	row12=np.matrix([x1,v1]).transpose()-phi(t1,0,x0,v0,Mdt)

	Mdt=M(T-t1)
	row34=np.matrix([x0,v0]).transpose()-phi(T,t1,x1,-v1,Mdt)

	row5=np.matrix([x1-sigma])
	
	return np.matrix(np.ma.vstack((row12,row34,row5)))



def nraphsonsolve(y0,tol,maxstep):
	y=y0
	dmody=1
	
	step=0
	while dmody>tol*tol and step<maxstep:
		yold=y
		dy=np.matrix(linalg.inv(J1(y)))*G(y)

#		dmody=float(dy.transpose()*dy)
		
		y-=dy
#desperate attempt	
		[x0,v0,x1,v1,t1]=[float(i[0]) for i in y]
		Mdt=M(t1)
		nexty=np.matrix([[1,0],[0,-1]])*phi(t1,0,x0,v0,Mdt)
		print "should be sigma: ",nexty[0]
		Mdt=M(T-t1)
		nexty=phi(T,t1,float(nexty[0]),float(nexty[1]),Mdt)
		diff=np.matrix([x0,v0]).transpose()-nexty
		dmody=float(diff.transpose()*diff)
		print "y"
		print y

		print "nexty"
		print nexty

		print "dmody"
		print dmody


#desperate attempt end

		step+=1
	return 	y

def test_funcs():
	print "M(1)"
	print M(1)
	print "Mdot(1)"
	print Mdot(1)

	Mdt=M(2)
	Mdotdt=Mdot(2)

	print "phi(2,0,-1,0,Mdt)"
	print phi(2,0,-1,0,Mdt)

	print "delphideltf(2,0,-1,0,Mdt,Mdotdt)"
	print delphideltf(2,0,-1,0,Mdt,Mdotdt)
	print "delphidelti(2,0,-1,0,Mdotdt)"
	print delphidelti(2,0,-1,0,MdtMdt,Mdotdt)
	
	y=np.matrix([1,2,3,4,5]).transpose()
	print "J1(y)"
	print J1(y)
	
	print "G(y)"
	print G(y)

	t=0
	while t<1000:
		Mdt=M(t)
		ph=phi(t,0,-1,0,Mdt)
		print t, float(ph[0]), float(ph[1])
		t+=0.1



def test_phi():
	t=0
	tmax=10
	dt=0.01
	
	r=np.matrix([-1.23,0.45]).transpose()
	x=-1.23
	v=0.45

	while t<tmax:
		Mdt=M(t)
		ph=phi(t,0,-1.23,0.45,Mdt)

		v+=(-gamma*v-k*x+F*cos(w*t))*dt
		x+=v*dt

		print t,float(ph[0]), float(ph[1]), x,v
		t+=dt

def test_nraph():
	y=np.matrix([-0.699318,0.738832,0.99,0.706029,1.58]).transpose()
	y=nraphsonsolve(y,0.001,2000)
	print "y"
	print y
	import commands
	commands.getoutput("./hardcol.out traj %f %f %f %f %f > traj.dat"%(float(y[0]),float(y[1]),100,F,gamma))

if __name__=='__main__':
#	test_funcs()
	test_nraph()
#	test_phi()
