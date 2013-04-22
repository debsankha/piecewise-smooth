set terminal epslatex color solid linewidth 3
set output 'presentation-gnuplottex-fig1.tex'
r=1.3; phi=1.5; m=2.2; w=3; set samples 2000; g=-0.1; set grid; se xlabel '$t$'; se ylabel '$x$'; f(x)=-cos(w*x)+r*exp(g*x/2)*cos(m*w*x/2+phi); g(x)=(1+r*exp(g*x/2))*sin(phi/2+(m/2-1)*w*x/2); h(x)=-(1+r*exp(g*x/2))*sin(phi/2+(m/2-1)*w*x/2); p [0:90] f(x) title '$-\cos {\omega x}+re^{-\gamma x/2}\cos {(m\omega x/2+\varphi )}$', g(x) lt 3 title "envelope", h(x) lt 3 title ""; \end {gnuplot} \end {center} \end {figure} \end {beamer@frameslide}\ifbeamer@twoscreenstext \beamer@dosecondscreennow {\begin {figure}[!htb] \begin {center} \caption {An envelope of \eqref {eq-shm-sol} for $\omega _g\approx \omega $} \label {fig-envelope-beats} \begin {gnuplot}[terminal=epslatex,terminaloptions=color solid linewidth 3,scale=0.7] r=1.3; phi=1.5; m=2.2; w=3; set samples 2000; g=-0.1; set grid; se xlabel '$t$'; se ylabel '$x$'; f(x)=-cos(w*x)+r*exp(g*x/2)*cos(m*w*x/2+phi); g(x)=(1+r*exp(g*x/2))*sin(phi/2+(m/2-1)*w*x/2); h(x)=-(1+r*exp(g*x/2))*sin(phi/2+(m/2-1)*w*x/2); p [0:90] f(x) title '$-\cos {\omega x}+re^{-\gamma x/2}\cos {(m\omega x/2+\varphi )}$', g(x) lt 3 title "envelope", h(x) lt 3 title ""; \end {gnuplot} \end {center} \end {figure} }\fi \ifbeamer@anotherslide \advance \beamer@slideinframe by 1\relax \relax \expandafter \iterate \fi \let \iterate \relax \beamer@writeslidentry \beamer@reseteecodes 

\begin{frame}
\begin{figure}[!htb]
\begin{center}
\caption{An envelope of \eqref{eq-shm-sol} for $\omega_g\gg\omega$}
\begin{gnuplot}[terminal=epslatex,terminaloptions=color solid linewidth 3,scale=0.7]
r=3;
phi=1.5;
m=0.1;
w=5;
set samples 2000;
g=-0.1;
set grid;
se xlabel '$t$';
se ylabel '$x$';
f(x)=-cos(w*x)+r*exp(g*x/2)*cos(m*w*x/2+phi);
g(x)=1+r*exp(g*x/2)*cos(m*w*x/2+phi);
h(x)=-1+r*exp(g*x/2)*cos(m*w*x/2+phi);
p [0:50] f(x) title '$-\cos{\omega x}+re^{-\gamma x/2}\cos{(m\omega x/2+\varphi)}$', g(x) lt 3 title "envelope", h(x) lt 3 title "";
