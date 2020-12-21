EX = $$(basename $$(pwd))

########
# EXPERIMENT

run/heat: push
	ssh hpc "cd $$USER/${EX}/code; make clean heatgrid_seq.out; sbatch slurm/heat_seq.sh"

hpcrun/heat: push
	ssh hpc "cd $$USER/${EX}/code; make clean heatgrid_seq.out; sbatch -p hpc slurm/heat_seq.sh"

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

