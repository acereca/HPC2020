doc: document.pdf

%.pdf:%.tex preamble.sty
	latexrun $< -Wall --latex-args="-shell-escape -synctex=1"

img/%.pdf: eval/%.py
	python $<
