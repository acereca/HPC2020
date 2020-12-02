EX = $$(basename $$(pwd))

########
# EXPERIMENT

run/proc: push
	ssh hpc "cd $$USER/${EX}/code; make clean mmul_para.out; sbatch slurm/proc.sh"

hpcrun/proc: push
	ssh hpc "cd $$USER/${EX}/code; make clean mmul_para.out; sbatch -p hpc slurm/proc.sh"

run/prob: push
	ssh hpc "cd $$USER/${EX}/code; make clean mmul_para.out; sbatch slurm/prob.sh"

hpcrun/prob: push
	ssh hpc "cd $$USER/${EX}/code; make clean mmul_para.out; sbatch -p hpc slurm/prob.sh"

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

