#!/usr/bin/python
import re
freq=0
yval=0
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
			freq+=1
		else:	
			print oldx,yval/freq
			yval=liney
			freq=1
		oldx=linex

	except EOFError:
		print oldx,yval/freq
		break
