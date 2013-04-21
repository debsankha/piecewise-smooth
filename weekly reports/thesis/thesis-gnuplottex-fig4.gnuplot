set terminal epslatex color solid linewidth 3
set output 'thesis-gnuplottex-fig4.tex'
r=3
phi=1.5
m=0.1
w=5
set samples 2000
g=-0.1

set grid
se xlabel '$t$'
se ylabel '$x$'

f(x)=-cos(w*x)+r*exp(g*x/2)*cos(m*w*x/2+phi)
g(x)=1+r*exp(g*x/2)*cos(m*w*x/2+phi)
h(x)=-1+r*exp(g*x/2)*cos(m*w*x/2+phi)

p [0:50] f(x) title '$-\cos{\omega x}+re^{-\gamma x/2}\cos{(m\omega x/2+\varphi)}$', g(x) lt 3 title "envelope", h(x) lt 3 title ""
