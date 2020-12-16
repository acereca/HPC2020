EX = $$(basename $$(cd ..; pwd))

######
# C++

SRC_DIR = src
OBJ_DIR = obj

SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRC:$(SRC_DIR)/%.cpp=%.out)
DBG = $(SRC:$(SRC_DIR)/%.cpp=%.dbg)

CC = mpic++ --std=c++11 -fdiagnostics-color=always -O2
CPPFLAGS += -Iinclude -I../../libs/argparse -I../../libs/export
CFLAGS += -Wall -Wextra -Wpedantic -fopenmp
LDFLAGS += -Llib
LDLIBS += -lpthread

.PHONY: all clean

all: $(OBJ) $(DBG)

%.dbg: $(SRC_DIR)/%.cpp ../../libs/mpi_debug/debug.hpp
	$(CC) $(CPPFLAGS) $(CFLAGS) -ggdb -Og -DDBG $< $(LDFLAGS) $(LDLIBS) -o $@

%.out: $(SRC_DIR)/%.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

clean:
	$(RM) $(OBJ)

test:
	echo ${EX}
