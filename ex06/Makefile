EX = $$(basename $$(pwd))

########
# EXPERIMENT

run/heat: push
	ssh hpc "cd $$USER/${EX}/code; make clean heatgrid_1d.out; sbatch slurm/heat_1d.sh"

hpcrun/heat: push
	ssh hpc "cd $$USER/${EX}/code; make clean heatgrid_1d.out; sbatch -p hpc slurm/heat_1d.sh"

push:
	rsync -urP ../ hpc:$$USER/ --exclude out --exclude eval --exclude ex01 --exclude .git

pull:
	rsync -uPr hpc:$$USER/${EX}/code/out ./code/


########
# UPLOAD

zip:
	zip -r hpc05${EX}.zip code/src/ code/Makefile code/slurm/

########
# CLEAN

clean:
	make -C docs clean
	make -C code clean

