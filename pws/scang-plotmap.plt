g=g-0.02
load "plotmap.plt"
pause 0.01
g_graz="`head -n 1 err.log | awk '{print $2}'`"
if (g>loat(g_graz)) reread
