system sprintf("./hardcol.out plotmap 200 %d %f %f 2>err.log | sort -n | python calc_deriv.py >  map%d.dat",n,f,g,n)
f_graz="`head -n 1 err.log | awk '{print $4}'`" 
g_graz="`head -n 1 err.log | awk '{print $2}'`"
amp="`head -n 1 err.log | awk '{print $6}'`"

set arrow 1 from amp,graph 0 to amp,graph 1 nohead

p sprintf("./map%d.dat",n) w l title sprintf("f=%f,f_graz=%s,g=%f,g_graz=%s",f,f_graz,g,g_graz), x

