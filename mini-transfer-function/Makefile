CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -Iinclude -O2 -g
LDFLAGS = -lm

SRC = src/transfer_function.c src/tf_polynomial.c src/tf_analysis.c \
      src/tf_interconnections.c src/tf_conversion.c src/tf_advanced.c
OBJ = $(SRC:.c=.o)
LIB = libtransferfunction.a

TEST_SRC = tests/test_tf.c
TEST_BIN = tests/test_tf

EX1_SRC = examples/ex1_dc_motor.c
EX1_BIN = examples/ex1_dc_motor
EX2_SRC = examples/ex2_mass_spring_damper.c
EX2_BIN = examples/ex2_mass_spring_damper
EX3_SRC = examples/ex3_rlc_circuit.c
EX3_BIN = examples/ex3_rlc_circuit
EX4_SRC = examples/ex4_suspension.c
EX4_BIN = examples/ex4_suspension
EX5_SRC = examples/ex5_quadrotor.c
EX5_BIN = examples/ex5_quadrotor

.PHONY: all clean test examples

all: $(LIB) test examples

$(LIB): $(OBJ)
	ar rcs $@ $^

src/%.o: src/%.c include/*.h
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_BIN): $(TEST_SRC) $(LIB)
	$(CC) $(CFLAGS) $< -L. -ltransferfunction $(LDFLAGS) -o $@

$(EX1_BIN): $(EX1_SRC) $(LIB)
	$(CC) $(CFLAGS) $< -L. -ltransferfunction $(LDFLAGS) -o $@

$(EX2_BIN): $(EX2_SRC) $(LIB)
	$(CC) $(CFLAGS) $< -L. -ltransferfunction $(LDFLAGS) -o $@

$(EX3_BIN): $(EX3_SRC) $(LIB)
	$(CC) $(CFLAGS) $< -L. -ltransferfunction $(LDFLAGS) -o $@

$(EX4_BIN): $(EX4_SRC) $(LIB)
	$(CC) $(CFLAGS) $< -L. -ltransferfunction $(LDFLAGS) -o $@

$(EX5_BIN): $(EX5_SRC) $(LIB)
	$(CC) $(CFLAGS) $< -L. -ltransferfunction $(LDFLAGS) -o $@

test: $(TEST_BIN)
	./tests/test_tf

examples: $(EX1_BIN) $(EX2_BIN) $(EX3_BIN) $(EX4_BIN) $(EX5_BIN)
	@echo "Examples built successfully."
	@echo "Run: ./examples/ex1_dc_motor"
	@echo "Run: ./examples/ex2_mass_spring_damper"
	@echo "Run: ./examples/ex3_rlc_circuit"
	@echo "Run: ./examples/ex4_suspension"
	@echo "Run: ./examples/ex5_quadrotor"

run-examples: examples
	./examples/ex1_dc_motor
	./examples/ex2_mass_spring_damper
	./examples/ex3_rlc_circuit
	./examples/ex4_suspension
	./examples/ex5_quadrotor

clean:
	rm -f src/*.o $(LIB) $(TEST_BIN) $(EX1_BIN) $(EX2_BIN) $(EX3_BIN) $(EX4_BIN) $(EX5_BIN)
