set terminal epslatex color solid linewidth 3
set output 'thesis-gnuplottex-fig3.tex'
r=1.3
phi=1.5
m=2.2
w=3
set samples 2000
g=-0.1

set grid
se xlabel '$t$'
se ylabel '$x$'

f(x)=-cos(w*x)+r*exp(g*x/2)*cos(m*w*x/2+phi)
g(x)=(1+r*exp(g*x/2))*sin(phi/2+(m/2-1)*w*x/2)
h(x)=-(1+r*exp(g*x/2))*sin(phi/2+(m/2-1)*w*x/2)

p [0:90] f(x) title '$-\cos{\omega x}+re^{-\gamma x/2}\cos{(m\omega x/2+\varphi)}$', g(x) lt 3 title "envelope", h(x) lt 3 title ""
