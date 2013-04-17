thesis: thesis.tex thesis.bib 
	pdflatex --shell-escape thesis.tex

thesis.bib: thesis.aux 
	pdflatex --shell-escape thesis.tex && bibtex thesis

thesis.aux:
	pdflatex --shell-escape thesis.tex	

