CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2 -Iinclude -D_GNU_SOURCE
LDFLAGS = -lm

SRC = src/mechatronic_definitions.c \
      src/electromechanical_coupling.c \
      src/mechatronic_state_space.c \
      src/actuator_modeling.c \
      src/sensor_modeling.c \
      src/mechatronic_system_design.c

TESTS = tests/test_mechatronic_modeling.c

all: test

init:
	mkdir -p build

test: init $(SRC) $(TESTS)
	$(CC) $(CFLAGS) -o build/test_runner $(TESTS) $(SRC) $(LDFLAGS)
	./build/test_runner

examples: init $(SRC)
	$(CC) $(CFLAGS) -o build/example_dc_servo examples/example_dc_servo.c $(SRC) $(LDFLAGS)
	$(CC) $(CFLAGS) -o build/example_ball_screw examples/example_ball_screw.c $(SRC) $(LDFLAGS)
	$(CC) $(CFLAGS) -o build/example_robot_joint examples/example_robot_joint.c $(SRC) $(LDFLAGS)

run-examples: examples
	./build/example_dc_servo
	./build/example_ball_screw
	./build/example_robot_joint

lean:
	@echo "Lean 4 formalization: src/mechatronic_modeling.lean"
	@echo "Run: lake build (requires Lean 4 toolchain)"

clean:
	rm -f build/test_runner build/example_dc_servo build/example_ball_screw build/example_robot_joint

.PHONY: all test examples run-examples lean clean init
