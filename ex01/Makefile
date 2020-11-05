all: document.pdf

%.pdf:%.tex img/amdahl.png
	latexrun $< -Wall

img/amdahl.png: scripts/amdahl.py
	python scripts/amdahl.py
