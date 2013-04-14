import sys
from nraphson import *
from math import *
from random import uniform

def test_nraph(f,g):
	fps=[]
	for i in range(50):
#		y=np.matrix([-0.699318,0.738832,0.99,0.706029,1.58]).transpose()
		y=np.matrix([uniform(-2,1),uniform(-3,3),1,uniform(-1,1),uniform(0,4*pi)]).transpose()
		y=nraphsonsolve(y,0.001,2000,f,g)
		
		isnew=1
				
		if y!=None:	
			[x0,v0,x1,v1,t1]=[float(i[0]) for i in y]
			if x0<sigma and abs(x1-sigma)<0.3 and t1>0 and t1<4*pi:
				for (fx,fv) in fps:
					if sqrt((fx-x0)**2+(fv-v0)**2)<0.01:
						isnew=0
						break
				if isnew:
					fps.append((x0,v0))

					commands.getoutput("./hardcol.out traj %f %f %f %f %f 2>/dev/null | tail -n 2 | grep -v '#' | awk '{print $2\"\t\"$3}' "%(x0,v0,20,f,g))
				#	nextx,nextv=(float(i) for i in nexty.split('\t'))
				#	print "Jacobian eigenvalues: "
					print float(y[0]), float(y[1]),Jfp(float(y[0]),float(y[1]),f,g)
	
	#	commands.getoutput("./hardcol.out traj %f %f %f %f %f > traj.dat"%(float(y[0]),float(y[1]),100,F,gamma))


def plot_fp_eigvals(fmin,fmax):
	g=0.062
	df=(fmax-fmin)/50
	f=fmax

	while f>fmin:
		fps=[]
		for i in range(50):
			y=np.matrix([uniform(-2,1),uniform(-3,3),1,uniform(-1,1),uniform(0,4*pi)]).transpose()
			y=nraphsonsolve(y,0.001,2000,f,g)
			
			isnew=1
					
			if y!=None:	
				[x0,v0,x1,v1,t1]=[float(i[0]) for i in y]

				for (fx,fv) in fps:
					if sqrt((fx-x0)**2+(fv-v0)**2)<0.01:
						isnew=0
						break
				
				if isnew:
					fps.append((x0,v0))
				else:
					sys.stderr.write("old one\n")
					continue

				if x0<sigma and abs(x1-sigma)<0.1 and t1>0 and t1<4*pi:
					nexty=commands.getoutput("./hardcol.out traj %f %f %f %f %f 2>/dev/null | tail -n 2 | grep -v '#' | awk '{print $2\"\t\"$3}' "%(x0,v0,4*pi,f,g))
					nextx,nextv=(float(i) for i in nexty.split('\t'))
					if sqrt((x0-nextx)**2+(nextv-v0)**2)<0.01:
						abseigs=[abs(i) for i in Jfp(x0,v0,f,g)]
						print f,min(abseigs),max(abseigs)
					else:
						sys.stderr.write("false positive\n")
				else:
					sys.stderr.write("unphysical\n")
			else:
				sys.stderr.write("NR failed\n")
		
		f-=df

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


if __name__=='__main__':
#	test_funcs()
#	test_nraph(float(sys.argv[1]),float(sys.argv[2]))
#	test_phi()
	plot_fp_eigvals(0.38,0.54)
