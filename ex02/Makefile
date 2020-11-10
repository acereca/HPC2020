EXE = EX02

SRC_DIR = src
OBJ_DIR = obj

SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

CC = mpic++ --std=c++11 -fdiagnostics-color=always -O2
CPPFLAGS += -Iinclude -I../libs
CFLAGS += -Wall -Wextra -Wpedantic
LDFLAGS += -Llib
LDLIBS += -lpthread

.PHONY: all clean

all: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(LDFLAGS) $^ $(LDLIBS) -o $@.out

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJ)

doc: document.pdf

%.pdf:%.tex
	latexrun $< -Wall

sync:
	rsync -uP ./src hpc:$$USER/ex02/src


