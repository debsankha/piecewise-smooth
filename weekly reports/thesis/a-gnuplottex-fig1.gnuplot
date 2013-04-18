set terminal epslatex color solid linewidth 3
set output 'a-gnuplottex-fig1.tex'
mu=0.4
unset ytics
set ytics add ('$\frac{\mu}{2}$' mu/2)
unset xtics
set xtics 0, 0.5, 1
f(x)=x<0.5?mu*x:mu*(1-x)
set samples 1000
p [0:1][:mu/1.4]  f(x)
