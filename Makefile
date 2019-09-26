SHELL := bash
CXX := g++
CPPFLAGS := -std=c++14 -Iinclude
CXXFLAGS := -Wall -O3 -flto -fmax-errors=3 $(CPPFLAGS)
# CXXFLAGS := -Wall -g -fmax-errors=3 $(CPPFLAGS) -DDEBUG_AT
LDFLAGS :=
LDLIBS :=

BLD := .build
EXT := .cc

.PHONY: all clean

ifeq (0, $(words $(findstring $(MAKECMDGOALS), clean)))

SRCS := $(shell find -L src -type f -name '*$(EXT)')
DEPS := $(patsubst src/%$(EXT),$(BLD)/%.d,$(SRCS))

GREP_EXES := grep -rl '^ *int \+main *(' src --include='*$(EXT)'
EXES := $(patsubst src%$(EXT),bin%, $(shell $(GREP_EXES)))

TESTS_DIR := etc/test
TESTS := $(shell find -L $(TESTS_DIR) -type f | sed -n 's,/\([^.][^/]*\)$$,/.\1,p')

all: $(TESTS) $(EXES)

$(TESTS): $(TESTS_DIR)/.%: $(TESTS_DIR)/%
	@$<; touch $@

-include $(DEPS)

# -------------------------------------------------------------------
bin/bin: \
  $(BLD)/ivanp/io/mem_file.o

bin/signif: \
  $(BLD)/ivanp/linalg.o \
  $(BLD)/ivanp/wls.o
# -------------------------------------------------------------------

$(DEPS): $(BLD)/%.d: src/%$(EXT)
	@mkdir -pv $(dir $@)
	$(CXX) $(CPPFLAGS) $(C_$*) -MM -MT '$(@:.d=.o)' $< -MF $@

$(BLD)/%.o:
	@mkdir -pv $(dir $@)
	$(CXX) $(CXXFLAGS) $(C_$*) -c $(filter %$(EXT),$^) -o $@

bin/%: $(BLD)/%.o
	@mkdir -pv $(dir $@)
	$(CXX) $(LDFLAGS) $(filter %.o,$^) -o $@ $(LDLIBS) $(L_$*)

endif

clean:
	@rm -rfv $(BLD) bin

