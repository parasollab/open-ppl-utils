# Compiler configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

CXX   := g++
FLAGS := -std=c++11 -fPIC -g3 -O2
WARN  := -Wall -Werror -Wno-parentheses
INCL  := -I.

COMPILE := $(CXX) $(FLAGS) $(WARN) $(INCL)

# Directory layout ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

OBJ_DIR  := build
SRC_DIR  := nonstd
TEST_DIR := test
DOC_DIR  := docs

$(OBJ_DIR):
	@mkdir $(OBJ_DIR)

vpath %.h        $(SRC_DIR)
vpath %.cpp      $(SRC_DIR)
vpath test_%.cpp $(TEST_DIR)

# Source file configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

SRCS      := $(shell find $(SRC_DIR) -name '*.cpp')
TEST_SRCS := $(shell find $(TEST_DIR) -name 'test_*.cpp')

# Object file configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o, $(SRCS))

$(OBJ_DIR)/%.o: %.cpp %.h | $(OBJ_DIR)
	@echo Compiling $(shell echo $< | sed 's|.*/||' )...
	@$(COMPILE) -c $< -o $@

# Libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

SHARED := libnonstd.so
STATIC := libnonstd.a

.DEFAULT_GOAL := $(STATIC)

.PHONY: shared
.PHONY: static
.PHONY: all
shared: $(SHARED)
static: $(STATIC)
all: $(SHARED) $(STATIC)

$(SHARED): $(OBJS)
	@echo Linking shared library...
	@$(CXX) -shared $^ -o $@

$(STATIC): $(OBJS)
	@echo Linking static library...
	@ar rcs $@ $^

# Test configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

TEST_OBJS := $(subst .cpp,.o, $(TEST_SRCS))

$(TEST_DIR)/%.o: %.cpp $(.DEFAULT_GOAL)
	@echo Building $(shell echo $< | sed s/.cpp// | sed 's|.*/||' )...
	@$(COMPILE) $< -L. -lnonstd -o $@

# Tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

.PHONY: tests
tests: $(TEST_OBJS)
	@echo Testing...
	@for testobj in $(TEST_OBJS); do $$testobj; done

# Docs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

HEADERS := $(shell find . -name '*.h')

.PHONY: docs
docs: $(SRCS) $(HEADERS)
	@echo Building docs...
	@doxygen $(DOC_DIR)/Doxyfile

# Clean ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

.PHONY: clean
.PHONY: cleantests
.PHONY: cleandocs
clean: cleandocs cleantests
	@rm -rf $(SHARED) $(STATIC) $(OBJ_DIR)

cleantests:
	@rm -f $(TEST_OBJS)

cleandocs:
	@rm -rf $(DOC_DIR)/html

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
