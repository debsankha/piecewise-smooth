#!/usr/bin/python
import re,sys
freq=0
yval=0
ysqval=0
oldx=0

while True:
	try:
		line=raw_input()
		dats=[float(i) for i in re.split("\s+",line)]
		if len(dats)<2:
			continue
		linex,liney=dats[0],dats[1]

		if freq==0 or abs(linex-oldx)<10e-5:
			yval+=liney
			ysqval+=liney**2
			freq+=1
		else:	
			mean=yval/freq
			sys.stderr.write(str(freq)+"\n")
			try:
				sd=(ysqval/freq-mean**2)**(0.5)	#sometimes ysqval/freq-mean**2 gets <0, ~10e-13. flotaing pt error?
				print oldx,mean,sd
			except:
				print oldx, mean, 0	
			yval=liney
			ysqval=liney**2
			freq=1
		oldx=linex

	except EOFError:
		mean=yval/freq
		try:
			sd=(ysqval/freq-mean**2)**(0.5)
			print oldx,mean,sd
		except:
			print oldx, mean, 0	

		break
