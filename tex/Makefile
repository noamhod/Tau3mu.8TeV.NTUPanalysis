project=cutflow.muons.mva

all	: ltx \
	  pdf

allbib	: ltx \
	  bib \
	  pdf

pdf     :
	dvips $(project).dvi
	ps2pdf $(project).ps
	open $(project).pdf

bib     :
	bibtex $(project)
	latex $(project).tex
	latex $(project).tex

ltx     :
	latex $(project).tex

dvi     :
	open $(project).dvi

clean:
	rm -f *.ps *.dvi *.pdf *.snm *.out *.log *.nav *.aux *.dvi *.toc *.blg *.bbl *.gz* tex/*.aux
	find . -maxdepth 1 -type l -exec rm -f {} \;
