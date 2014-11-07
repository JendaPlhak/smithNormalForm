# clang++ is a must for now
CXX     := clang++

OPT_LVL := -O0
DBG_LVL := -g3

CFLAGS  := -std=c++11 -fpic -march=native -Wall -Wextra
CFLAGS  += $(OPT_LVL) $(DBG_LVL)
CFLAGS  += -Werror-implicit-function-declaration -D_GNU_SOURCE

N_JOBS  := $(shell echo $$((2 * `grep -c "^processor" /proc/cpuinfo`)))

# it is all of -std-compile-opts -std-link-opts without -internalize
OPT_OPTS := -adce -argpromotion -basicaa -basiccg -constmerge \
	    -correlated-propagation -deadargelim -domtree -dse -early-cse \
	    -functionattrs -globaldce -globalopt -globalsmodref-aa -gvn \
	    -indvars -inline -inline-cost -instcombine -ipsccp \
	    -jump-threading -lazy-value-info -lcssa -licm -loop-deletion \
	    -loop-idiom -loop-rotate -loops -loop-simplify -loop-unroll \
	    -loop-unswitch -memcpyopt -memdep -no-aa -notti -preverify \
	    -prune-eh -reassociate -scalar-evolution -sccp -simplifycfg -sroa \
	    -strip-dead-prototypes -tailcallelim -targetlibinfo -tbaa -verify

LTO :=
LTO := -emit-llvm

BUILD := build
SRC   := src

LIBS      := -lm -larmadillo -llapack -lblas
SOURCES   := $(sort $(wildcard $(SRC)/*.cpp))
BYTECODES := $(addprefix $(BUILD)/,$(notdir $(SOURCES:.cpp=.bc)))
DEPENDS   := $(BYTECODES:.bc=.d)
BINS      := tester
LIB       := snf.so

BIN_BASES     := $(BINS)
BIN_SOURCES   := $(addprefix $(BUILD)/,$(notdir $(BIN_BASES:=.cpp)))
BIN_BYTECODES := $(BIN_SOURCES:.cpp=.bc)
LIB_SOURCES   := $(filter-out $(BIN_SOURCES),$(SOURCES))
LIB_BYTECODES := $(filter-out $(BIN_BYTECODES),$(BYTECODES))
LIB_BYTECODE  := $(BUILD)/$(LIB:.so=.bc)

# re-expand variables in subsequent rules
.SECONDEXPANSION:


# parallel make all
pall:
	make -j$(N_JOBS) all

all: $(BINS) $(LIB)

# add debug flag to flags, then compile
debug: CFLAGS += -DRRRR_INFO
debug: $(BINS)

# parallel make debug
pdebug:
	make -j$(N_JOBS) debug

# include compiler-generated dependencies, so obj files get recompiled when
# their header changes
-include $(DEPENDS)

# -MMD generates $(BUILD)/%.d dependency file (only local header dependencies,
# -MD for system-wide header dependencies)
$(BUILD)/%.bc: $(SRC)/%.cpp Makefile | $(BUILD)
	$(CXX) -c $(CFLAGS) $(LTO) -MMD -o $@ $<

$(BUILD):
	mkdir $(BUILD)

# each binary depends on its own .bc file and the library bytecode
$(BINS): $(BUILD)/$$@.bc $(LIB_BYTECODE)
	$(CXX) $(CFLAGS) $^ $(LIBS) -o $@

$(LIB): $(LIB_BYTECODE)
	$(CXX) -shared -fPIC $^ -o $@

$(LIB_BYTECODE): $(LIB_BYTECODES)
	llvm-link $^ -o $@
ifeq ($(OPT_LVL),-O3)
	opt $(OPT_OPTS) $(OPT_LVL) $@ -o $@
endif

clean:
	rm -rf $(BUILD) *.so *~ core $(BINS)
	rm -f tests/*.o tests/*~ run_tests

show:
	# $(SOURCES)
	# $(BYTECODES)
	# $(BIN_BASES)
	# $(BIN_BYTECODES)
	# $(OTHER_BYTECODES)

.PHONY: clean show all

