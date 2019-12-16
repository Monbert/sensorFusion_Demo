IDIR=include
SOURCE_DIR=src
TEST_DIR=test/src
CC=gcc
CFLAGS= -Wall


#CREATING BIN, BUILD AND OUTPUT FOLDERS. THESE WILL CONTAIN EXECUTABLE, .o FILES AND OUTPUT OF PROGRAM RESPECTIVELY
bin_folder := $(shell mkdir -p bin)
build_folder := $(shell mkdir -p build)
results_folder := $(shell mkdir -p data/output)

#TARGET TO COMPILE ALL .c FILES
SensorFusionAlgorithm: $(SOURCE_DIR)/algorithm.c $(SOURCE_DIR)/data.c
	$(CC) -g -c $(CFLAGS) $(SOURCE_DIR)/algorithm.c -o build/algorithm.o 
	$(CC) -g -c $(CFLAGS) $(SOURCE_DIR)/data.c -o build/data.o

main: $(SOURCE_DIR)/main.c
	$(CC) -g -c $(CFLAGS) $(SOURCE_DIR)/main.c -o build/main.o

main_test: $(TEST_DIR)/main.c
	$(CC) -g -c $(CFLAGS) $(TEST_DIR)/main.c -o build/main_test.o

executable: build/main.o build/algorithm.o build/data.o
	gcc -g -o bin/executable build/main.o build/algorithm.o build/data.o

tests: build/main_test.o build/algorithm.o build/data.o
	gcc -g -o bin/tests build/main_test.o build/algorithm.o build/data.o

#TO BE RUN TO BUILD THE OBJECT FILES AND EXECUTABLE
all: SensorFusionAlgorithm main executable

#TO BE RUN TO BUILD THE OBJECT FILES AND TEST
test: SensorFusionAlgorithm main_test tests

#RUNS THE PROGRAM
run: 
	./bin/executable

#RUN TO TEST PROGRAM
run_test: 
	./bin/tests

#USED TO CLEAR THE BUILD AND BIN DIRECTORY
.PHONY : clean
clean: 
	rm -f build/* bin/*
