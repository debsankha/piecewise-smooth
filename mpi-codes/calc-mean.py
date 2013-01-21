#!/usr/bin/python
freq=0
yval=0
oldx=0

while True:
	try:
		line=raw_input()
		if line[0]=='#':
			continue
		
		linex,liney=[float(i) for i in line.split('\t')]

		if freq==0 or abs(linex-oldx)<10e-5:
			yval+=liney
			freq+=1
		else:	
			print oldx,yval/freq
			yval=liney
			freq=1
		oldx=linex
	except:
		if freq>0:
			print linex,yval/freq
		break
