EXE = EX02

######
# C++

SRC_DIR = src
OBJ_DIR = obj

SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRC:$(SRC_DIR)/%.cpp=%.out)

CC = mpic++ --std=c++11 -fdiagnostics-color=always -O2
CPPFLAGS += -Iinclude -I../libs/argparse
CFLAGS += -Wall -Wextra -Wpedantic
LDFLAGS += -Llib
LDLIBS += -lpthread

.PHONY: all clean

all: $(OBJ) doc

%.out: $(SRC_DIR)/%.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

mpirun: $(OBJ)
	mpirun -np $$(nproc) --oversubscribe $(OBJ) -n 1000000

clean:
	$(RM) $(OBJ)
	latexrun --clean-all

########
# DOCUMENTATION

doc: document.pdf preamble.sty

%.pdf:%.tex
	latexrun $< -Wall --latex-args="-shell-escape -synctex=1"

########
# EXPERIMENT

sync:
	rsync -urP ../ hpc:$$USER/ --exclude out --exclude eval --exclude ex01 --exclude .git

run/ring:
	ssh hpc "cd $$USER/ex02; make clean ring.out; sbatch slurm_ring.sh"

pull:
	rsync -uPr hpc:$$USER/ex02/out ./
