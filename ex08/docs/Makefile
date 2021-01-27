doc: document.pdf

%.pdf:%.tex img/comm_overlap.pdf
	latexrun $< -Wall --latex-args="-shell-escape -synctex=1"

img/%.pdf: eval/%.py
	python $<

data/%.tex: eval/%.py
	python $<

vis:
	python eval/nbody_vis.py && mpv --fs img/nbody_vis.mp4
