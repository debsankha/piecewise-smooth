f=f+0.02
load "plotmap.plt"
pause 0.1
f_graz="`head -n 1 err.log | awk '{print $4}'`"
if (f<f_graz) reread
