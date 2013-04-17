system sprintf("python nextcollision.py -10 1 %f %g > next.dat",f,g)
p [:][0:] "next.dat" w l title sprintf("f=%f,g=%f",f,g), 1
