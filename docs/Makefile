# Makefile for manual

MAN_EN=manual_en
MAN_RU=manual_ru
TEX=pdflatex
BIB=bibtex

all: MAN_EN MAN_RU

MAN_EN:
	$(TEX) $(MAN_EN).tex
	$(BIB) $(MAN_EN).aux
	$(TEX) $(MAN_EN).tex
	$(TEX) $(MAN_EN).tex

MAN_RU:
	$(TEX) $(MAN_RU).tex
	$(BIB) $(MAN_RU).aux
	$(TEX) $(MAN_RU).tex
	$(TEX) $(MAN_RU).tex

clean:
	rm -rf *.{aux,log,out,gz,toc,blg,bbl}

