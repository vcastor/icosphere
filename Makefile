#Makefile

ALL=$(wildcard *.tex */*.tex bibl/*.bib)
MAIN=homework.tex
LATEX=pdflatex
SHELL=/bin/zsh
PYTHON=python3
PYPLOT=$(wildcard ./img/*.py)

all: plots tex

plots:                            ## Build python plots
	$(foreach var,$(PYPLOT),$(PYTHON) $(var);)

tex:                               ## here we go! 🥳
	$(LATEX) $(MAIN)                # main run
	bibtex $(MAIN:.tex=)            # bibliography
	makeglossaries $(MAIN:.tex=)    # list of abbreviations, nomenclature
	$(LATEX) $(MAIN)
	$(LATEX) $(MAIN)

clean:
	rm -rf *.aux *.bbl *.blg *.glg *.glo *.gls *.ist *.log *.not *.ntt *.out *.sbl *.sym *.tld *.toc */*.aux
