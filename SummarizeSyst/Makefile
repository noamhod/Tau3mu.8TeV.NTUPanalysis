project=SummarizeSyst

all	: ltx \
	  pdf

pdf     :
	dvips $(project).dvi
	ps2pdf $(project).ps
	open $(project).pdf

ltx     :
	latex $(project).tex

dvi     :
	open $(project).dvi

clean:
	rm -f *.ps *.dvi *.pdf *.snm *.out *.log *.nav *.aux *.dvi *.toc *.blg *.bbl *.gz*
