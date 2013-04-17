set terminal epslatex color solid linewidth 3
set output 'thesis-gnuplottex-fig1.tex'
mu=0.4
unset ytics
set ytics add ('$\frac{\mu}{2}$' mu/2)
unset xtics
set arrow 1
set xtics 0, 0.5, 1
f(x)=x<0.5?mu*x:mu*(1-x)
set samples 1000
set arrow from 0.5, graph 0 to 0.5, graph 1 nohead ls 3 lw 0.5
p [0:1][:mu/1.4]  f(x)
