f(x,y)=(abs(x)<1) && (abs(y)<1) ? 1 : 0
p [0.7:6.5] "./n.dat" w l title "Eigenvalue 1", '' u 1:3 w l title "Eigenvalue 2", '' u ($1):(f($2,$3)) w boxes fs solid  0.7 noborder lc rgb "#9ad8ec" title "Region where FP is attractor"
