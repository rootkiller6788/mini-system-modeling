CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -Iinclude
LDFLAGS = -lm

# All source files
SRC = src/linearization_core.c \
      src/taylor_expansion.c \
      src/feedback_linearization.c \
      src/koopman_operator.c \
      src/describing_function.c \
      src/jacobian_methods.c \
      src/gain_scheduling.c \
      src/carleman_linearization.c

TESTS = tests/test_linearization.c
EXAMPLES = examples/example_pendulum.c \
           examples/example_maglev.c \
           examples/example_dc_motor_saturation.c \
           examples/example_quadrotor.c

all: test

init:
	mkdir -p build

test: init $(SRC) $(TESTS)
	$(CC) $(CFLAGS) -o build/test_runner $(TESTS) $(SRC) $(LDFLAGS)
	./build/test_runner

examples: init $(SRC)
	$(CC) $(CFLAGS) -o build/example_pendulum examples/example_pendulum.c $(SRC) $(LDFLAGS)
	$(CC) $(CFLAGS) -o build/example_maglev examples/example_maglev.c $(SRC) $(LDFLAGS)
	$(CC) $(CFLAGS) -o build/example_dc_motor examples/example_dc_motor_saturation.c $(SRC) $(LDFLAGS)
	$(CC) $(CFLAGS) -o build/example_quadrotor examples/example_quadrotor.c $(SRC) $(LDFLAGS)

run-examples: examples
	./build/example_pendulum
	./build/example_maglev
	./build/example_dc_motor
	./build/example_quadrotor

clean:
	rm -f build/test_runner build/example_pendulum build/example_maglev
	rm -f build/example_dc_motor build/example_quadrotor

.PHONY: all test examples run-examples clean init
