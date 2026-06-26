CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -Iinclude -D_GNU_SOURCE
LDFLAGS = -lm

SRC = src/electrical_elements.c \
      src/circuit_topology.c \
      src/transfer_function_electrical.c \
      src/electrical_state_space.c \
      src/two_port_network.c \
      src/electromechanical_systems.c

TESTS = tests/test_electrical_modeling.c

all: test

init:
	mkdir -p build

test: init $(SRC) $(TESTS)
	$(CC) $(CFLAGS) -o build/test_runner $(TESTS) $(SRC) $(LDFLAGS)
	./build/test_runner

examples: init $(SRC)
	$(CC) $(CFLAGS) -o build/example_rc_filter examples/example_rc_filter.c $(SRC) $(LDFLAGS)
	$(CC) $(CFLAGS) -o build/example_dc_motor examples/example_dc_motor.c $(SRC) $(LDFLAGS)
	$(CC) $(CFLAGS) -o build/example_two_port examples/example_two_port.c $(SRC) $(LDFLAGS)

run-examples: examples
	./build/example_rc_filter
	./build/example_dc_motor
	./build/example_two_port

clean:
	rm -f build/test_runner build/example_rc_filter build/example_dc_motor build/example_two_port

.PHONY: all test examples run-examples clean init
