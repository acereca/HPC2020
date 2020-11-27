EX = $$(basename $$(pwd))

########
# EXPERIMENT

run/latency: push
	ssh hpc "cd $$USER/${EX}/code; make clean latency.out; sbatch slurm/latency1.sh; sbatch slurm/latency2.sh"

hpcrun/latency: push
	ssh hpc "cd $$USER/${EX}/code; make clean latency.out; sbatch -p hpc slurm/latency1.sh; sbatch -p hpc slurm/latency2.sh"

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

