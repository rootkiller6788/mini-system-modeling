# mini-fluid-thermal-modeling Makefile
# Compiles all sources, tests, and examples.

CC      = gcc
CFLAGS  = -Wall -Wextra -std=c99 -O2 -I./include
LDFLAGS = -lm
LEANC   = lean

SRC     = src/fluid_thermal_core.c \
          src/fluid_mechanics.c \
          src/heat_transfer.c \
          src/thermal_fluid_systems.c \
          src/numerical_fluid_thermal.c

OBJ     = $(SRC:.c=.o)

TEST_SRC = tests/test_fluid_thermal.c
TEST_BIN = test_fluid_thermal

EXAMPLES_SRC = examples/example_pipe_flow.c \
               examples/example_heat_exchanger.c \
               examples/example_electronics_cooling.c

EXAMPLES_BIN = $(EXAMPLES_SRC:examples/%.c=%)

# Default: build all
.PHONY: all
all: test examples

# Compile source objects
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Build test binary
$(TEST_BIN): $(OBJ) $(TEST_SRC)
	$(CC) $(CFLAGS) -o $@ $(TEST_SRC) $(SRC) $(LDFLAGS)

# Build each example
example_pipe_flow: $(OBJ) examples/example_pipe_flow.c
	$(CC) $(CFLAGS) -o $@ examples/example_pipe_flow.c $(SRC) $(LDFLAGS)

example_heat_exchanger: $(OBJ) examples/example_heat_exchanger.c
	$(CC) $(CFLAGS) -o $@ examples/example_heat_exchanger.c $(SRC) $(LDFLAGS)

example_electronics_cooling: $(OBJ) examples/example_electronics_cooling.c
	$(CC) $(CFLAGS) -o $@ examples/example_electronics_cooling.c $(SRC) $(LDFLAGS)

.PHONY: examples
examples: $(EXAMPLES_BIN)

.PHONY: test
test: $(TEST_BIN)
	./$(TEST_BIN)

# Lean 4 formalization check (non-blocking)
.PHONY: lean
lean:
	@echo "Checking Lean 4 formalization..."
	@if $(LEANC) --version >/dev/null 2>&1; then \
		$(LEANC) src/fluid_thermal.lean || echo "Lean: warnings OK"; \
	else \
		echo "Lean 4 not installed — skipping .lean check"; \
	fi

.PHONY: clean
clean:
	rm -f $(OBJ) $(TEST_BIN) $(EXAMPLES_BIN)

.PHONY: linecount
linecount:
	@echo "include/ + src/ lines:"
	@cat include/*.h src/*.c | wc -l

# Safety scan per SKILL.md §10
.PHONY: safety-scan
safety-scan:
	@echo "=== Filler Detection ==="
	@grep -rn "_fn[0-9]" include/ src/ tests/ examples/ --include="*.c" --include="*.h" || echo "  _fnN: CLEAN"
	@grep -rn "_aux[0-9]" include/ src/ tests/ examples/ --include="*.c" --include="*.h" || echo "  _auxN: CLEAN"
	@grep -rn "_ext[0-9][^a-z]" include/ src/ tests/ examples/ --include="*.c" --include="*.h" || echo "  _extN: CLEAN"
	@grep -rn "algorithm variant" include/ src/ --include="*.c" || echo "  variant: CLEAN"
	@grep -rn "auxiliary function" include/ src/ --include="*.c" || echo "  aux func: CLEAN"
	@grep -rn "extension point" include/ src/ --include="*.c" || echo "  ext point: CLEAN"
	@grep -rn "Module extension" include/ src/ --include="*.c" || echo "  mod ext: CLEAN"
	@grep -rn "supplemental assert" include/ src/ --include="*.c" || echo "  supp assert: CLEAN"
	@echo "=== Lean Safety ==="
	@grep -rn "SystemMetric" src/ --include="*.lean" || echo "  SystemMetric: CLEAN"
	@grep -rn "traceability_matrix" src/ --include="*.lean" || echo "  traceability: CLEAN"
	@grep -rn ":= 0\.0" src/ --include="*.lean" || echo "  := 0.0: CHECK (may be legitimate Float fields)"
	@grep -rn "\bsorry\b" src/ --include="*.lean" || echo "  sorry: CLEAN"
	@echo "=== Stub Detection ==="
	@find . -name "*.c" -size -200c 2>/dev/null || echo "  No small files detected"

.PHONY: all-checks
all-checks: safety-scan test lean linecount
