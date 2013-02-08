system sprintf("./hardcol.out plotmap 200 %d 0.393094 %f | sort -n | python calc_deriv.py >  map.dat",n,a)

p "./map.dat" w p title sprintf("g=%f",a), x

