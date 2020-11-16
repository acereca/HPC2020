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

mpirun: $(EXE)
	mpirun -np $$(nproc) --oversubscribe $(EXE).out -n 1000000

clean:
	$(RM) $(OBJ)

########
# DOCUMENTATION

doc: document.pdf

%.pdf:%.tex
	latexrun $< -Wall

########
# EXPERIMENT

sync:
	rsync -urP ../ hpc:$$USER/ --exclude out --exclude eval --exclude ex01 

run/ring:
	ssh hpc "cd $$USER/ex02; make ring.out; sbatch slurm.sh"

pull:
	rsync -uPr hpc:$$USER/ex02/out ./
