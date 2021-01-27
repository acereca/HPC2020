doc: document.pdf

%.pdf:%.tex img/comm_overlap.pdf
	latexrun $< -Wall --latex-args="-shell-escape -synctex=1"

img/%.pdf: eval/%.py
	python $<

data/%.tex: eval/%.py
	python $<
