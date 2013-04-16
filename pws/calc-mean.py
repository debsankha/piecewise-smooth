n=0
y=0

while True:
	try:
		line=raw_input()
		if line[0]=='#':
			continue
		n+=1
		linex,liney=[float(i) for i in line.split('\t')]
		y+=liney
		if n==100:
			n=0
			print linex,y/10.0
			y=0
		else:	
			pass
	except:
		break
