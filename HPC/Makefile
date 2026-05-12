CC = mpicc

CFLAGS_COMMON = -Ofast -flto -march=native -Wall -Wextra
CFLAGS_SERIAL = $(CFLAGS_COMMON)
CFLAGS_OMP    = $(CFLAGS_COMMON) -fopenmp

HEADERS = include
SRCS    = src
BUILD   = build
BIN     = bin

TARGET_SERIAL   = $(BIN)/stencil_serial
TARGET_PARALLEL = $(BIN)/stencil_parallel

OBJ_SERIAL   = $(BUILD)/stencil_template_serial.o
OBJ_PARALLEL = $(BUILD)/stencil_template_parallel.o

# ===== MPI PARAMETERS =====
NP ?= 4
ARGS ?=


.PHONY: all serial parallel run-serial run-parallel clean

all: serial parallel

serial: $(TARGET_SERIAL)
parallel: $(TARGET_PARALLEL)

run-serial: serial
	./$(TARGET_SERIAL) $(ARGS)

run-parallel: parallel
	mpirun -np $(NP) ./$(TARGET_PARALLEL) $(ARGS)

upload:
	eval "$(ssh-agent -s)"
	ssh-keygen -f "/home/rita_ortese/.ssh/known_hosts" -R "login.leonardo.cineca.it"
	rsync -av --exclude-from=".ignore" ./ rortese0@login.leonardo.cineca.it:./project
	chmod +x test/*

download:
	eval "$(ssh-agent -s)"
	mkdir -p output
	ssh-keygen -f "/home/rita_ortese/.ssh/known_hosts" -R "login.leonardo.cineca.it"
	rsync -avzP rortese0@login.leonardo.cineca.it:./project/output_parallel ./output/

# ===== LINK =====

$(TARGET_SERIAL): $(OBJ_SERIAL) | $(BIN)
	$(CC) $(CFLAGS_SERIAL) -I$(HEADERS) $^ -o $@

$(TARGET_PARALLEL): $(OBJ_PARALLEL) | $(BIN)
	$(CC) $(CFLAGS_OMP) -I$(HEADERS) $^ -o $@

# ===== COMPILE =====

$(OBJ_SERIAL): $(SRCS)/stencil_template_serial.c | $(BUILD)
	$(CC) $(CFLAGS_SERIAL) -I$(HEADERS) -c $< -o $@

$(OBJ_PARALLEL): $(SRCS)/stencil_template_parallel.c | $(BUILD)
	$(CC) $(CFLAGS_OMP) -I$(HEADERS) -c $< -o $@

# ===== DIRECTORIES =====

$(BUILD):
	mkdir -p $(BUILD)

$(BIN):
	mkdir -p $(BIN)

$(OUTPUT):
	mkdir -p $(OUTPUT)

# ===== CLEAN =====

clean:
	rm -rf $(BUILD) $(BIN) $(OUTPUT)