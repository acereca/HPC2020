doc: document.pdf

%.pdf:%.tex preamble.sty data/scaling_proc.tex data/scaling_prob.tex
	latexrun $< -Wall --latex-args="-shell-escape -synctex=1"

img/%.pdf: eval/%.py
	python $<

data/%.tex: eval/%.py
	python $<
