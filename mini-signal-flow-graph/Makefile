CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -Iinclude
LDFLAGS = -lm

SRCS = src/sfg_core.c src/sfg_path.c src/sfg_mason.c src/sfg_reduction.c \
       src/sfg_state_space.c src/sfg_analysis.c src/sfg_complex.c \
       src/sfg_matrix.c
OBJS = $(SRCS:.c=.o)

TEST_SRCS = tests/test_sfg.c
TEST_BINS = $(TEST_SRCS:.c=)

EXAMPLE_SRCS = examples/example_mass_spring.c examples/example_electrical.c \
               examples/example_control_tf.c examples/example_aircraft.c \
               examples/example_mimo.c
EXAMPLE_BINS = $(EXAMPLE_SRCS:.c=)

.PHONY: all clean test examples

all: $(OBJS)

test: $(TEST_BINS)
	@for t in $(TEST_BINS); do echo "Running $$t..."; ./$$t && echo "PASS: $$t" || echo "FAIL: $$t"; done

examples: $(EXAMPLE_BINS)
	@for e in $(EXAMPLE_BINS); do echo "Running $$e..."; ./$$e || echo "FAIL: $$e"; done

tests/test_sfg: tests/test_sfg.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ tests/test_sfg.c $(OBJS) $(LDFLAGS)

examples/example_mass_spring: examples/example_mass_spring.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ examples/example_mass_spring.c $(OBJS) $(LDFLAGS)

examples/example_electrical: examples/example_electrical.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ examples/example_electrical.c $(OBJS) $(LDFLAGS)

examples/example_control_tf: examples/example_control_tf.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ examples/example_control_tf.c $(OBJS) $(LDFLAGS)

examples/example_aircraft: examples/example_aircraft.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ examples/example_aircraft.c $(OBJS) $(LDFLAGS)

examples/example_mimo: examples/example_mimo.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ examples/example_mimo.c $(OBJS) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f src/*.o tests/*.o tests/*.exe $(TEST_BINS)
	rm -f examples/*.o examples/*.exe $(EXAMPLE_BINS)
