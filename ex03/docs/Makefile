doc: document.pdf

%.pdf:%.tex preamble.sty img/latency.pdf img/bandwidth.pdf
	latexrun $< -Wall --latex-args="-shell-escape -synctex=1"

img/%.pdf: eval/%.py
	python $<
