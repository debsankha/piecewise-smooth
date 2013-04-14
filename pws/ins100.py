i=0
while True:
	try:
		line=raw_input()
		i+=1
		print line
		if i==100:
			i=0
			print 
	except:
		break		


