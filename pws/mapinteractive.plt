#define g,maxg before loading this file

system sprintf("./hardcol.out plotmap 200 1 0.393094 %f >  map.dat",a)
p "./map.dat" w lp title sprintf("g=%f",a), x
a=a+da
pause mouse 
print a
if (a>mina) reread
