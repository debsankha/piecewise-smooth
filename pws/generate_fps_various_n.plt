set object 9 rect fc rgb "#9ad8ec" fs solid 0.3
set object 9 from 0.08, -1 to 0.26, 1
set grid
p [:0.18][-1.5:3.5] "./fpsnew_n_2.0.dat" w p ps 1 title "n=2.0", '' u 1:3 w p ps 1 title "n=2.0", "./fpsnew_n_2.01.dat" w p  title "n=2.01", '' u 1:3 w p title "n=2.01", "./fpsnew_n_2.1.dat" w p title "n=2.1", '' u 1:3 w p ls 18 title "n=2.1"
