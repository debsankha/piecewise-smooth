import re


xli=[]
yli=[]


numlines=0
while True:
	try:
		line=raw_input()
		x,y=[float(i) for i in re.split("\s+",line)]
		xli.append(x)
		yli.append(y)
		numlines+=1

	except ValueError:
		continue

	except EOFError:
		break


i=2
while i<numlines-2:
	print xli[i],yli[i],(-yli[i+2]+8*yli[i+1]-8*yli[i-1]+yli[i-2])/(12*(xli[i+1]-xli[i]))
	i+=1	
