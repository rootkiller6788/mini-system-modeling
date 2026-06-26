CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -Iinclude -D_GNU_SOURCE
LDFLAGS = -lm

SRC = src/mechanical_elements.c \
      src/mechanical_transfer_function.c \
      src/mechanical_state_space.c \
      src/rotational_systems.c \
      src/multi_dof_systems.c \
      src/mechanical_applications.c

TESTS = tests/test_mechanical_modeling.c

all: test

init:
	mkdir -p build

test: init $(SRC) $(TESTS)
	$(CC) $(CFLAGS) -o build/test_runner $(TESTS) $(SRC) $(LDFLAGS)
	./build/test_runner

examples: init $(SRC)
	$(CC) $(CFLAGS) -o build/example_sdof examples/example_sdof_vibration.c $(SRC) $(LDFLAGS)
	$(CC) $(CFLAGS) -o build/example_2dof examples/example_2dof_modal.c $(SRC) $(LDFLAGS)
	$(CC) $(CFLAGS) -o build/example_qcar examples/example_quarter_car.c $(SRC) $(LDFLAGS)

run-examples: examples
	./build/example_sdof
	./build/example_2dof
	./build/example_qcar

clean:
	rm -f build/test_runner build/example_sdof build/example_2dof build/example_qcar

.PHONY: all test examples run-examples clean init
