CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -Iinclude -O2 -g
LDFLAGS = -lm

SRC = src/blockdiagram.c src/transfer_function.c src/signal_flow_graph.c       src/block_algebra.c src/mason.c src/block_reduction.c src/bd_graph_utils.c
OBJ = $(SRC:.c=.o)
LIB = libblockdiagram.a

TEST_SRC = tests/test_all.c
TEST_BIN = tests/test_all

EX1_SRC = examples/ex1_feedback.c
EX1_BIN = examples/ex1_feedback
EX2_SRC = examples/ex2_mason.c
EX2_BIN = examples/ex2_mason
EX3_SRC = examples/ex3_dc_motor.c
EX3_BIN = examples/ex3_dc_motor
EX4_SRC = examples/ex4_quadrotor.c
EX4_BIN = examples/ex4_quadrotor

.PHONY: all clean test examples

all: $(LIB) test examples

$(LIB): $(OBJ)
	ar rcs $@ $^

src/%.o: src/%.c include/*.h
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_BIN): $(TEST_SRC) $(LIB)
	$(CC) $(CFLAGS) $< -L. -lblockdiagram $(LDFLAGS) -o $@

$(EX1_BIN): $(EX1_SRC) $(LIB)
	$(CC) $(CFLAGS) $< -L. -lblockdiagram $(LDFLAGS) -o $@

$(EX2_BIN): $(EX2_SRC) $(LIB)
	$(CC) $(CFLAGS) $< -L. -lblockdiagram $(LDFLAGS) -o $@

$(EX3_BIN): $(EX3_SRC) $(LIB)
	$(CC) $(CFLAGS) $< -L. -lblockdiagram $(LDFLAGS) -o $@

$(EX4_BIN): $(EX4_SRC) $(LIB)
	$(CC) $(CFLAGS) $< -L. -lblockdiagram $(LDFLAGS) -o $@

test: $(TEST_BIN)
	./tests/test_all

examples: $(EX1_BIN) $(EX2_BIN) $(EX3_BIN) $(EX4_BIN)
	@echo "Examples built successfully."
	@echo "Run: ./examples/ex1_feedback"
	@echo "Run: ./examples/ex2_mason"
	@echo "Run: ./examples/ex3_dc_motor"
	@echo "Run: ./examples/ex4_quadrotor"

run-examples: examples
	./examples/ex1_feedback
	./examples/ex2_mason
	./examples/ex3_dc_motor
	./examples/ex4_quadrotor

clean:
	rm -f src/*.o $(LIB) $(TEST_BIN) $(EX1_BIN) $(EX2_BIN) $(EX3_BIN) $(EX4_BIN)
