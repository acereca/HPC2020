doc: document.pdf

%.pdf:%.tex data/heat_s.tex data/heat_t.tex data/heat_e.tex
	latexrun $< -Wall --latex-args="-shell-escape -synctex=1"

img/%.pdf: eval/%.py
	python $<

data/%.tex: eval/%.py
	python $<
