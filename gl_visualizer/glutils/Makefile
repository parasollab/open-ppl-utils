# Externals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

NONSTD_DIR  := ../nonstd
NONSTD_INCL := -I$(NONSTD_DIR)
NONSTD_LIBS := -L$(NONSTD_DIR) -lnonstd

GL_INCL     :=
GL_LIBS     := -lGLU -lGL

# Compiler configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

CXX   := g++
FLAGS := -std=c++11 -fPIC -g3 -O3
WARN  := -Wall -Werror -Wextra -Wno-parentheses

INCL  := -I. $(NONSTD_INCL) $(GL_INCL)

COMPILE := $(CXX) $(FLAGS) $(WARN) $(INCL)

# Linker configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

LIBS := $(NONSTD_LIBS) $(GL_LIBS)

# Directory layout ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

OBJ_DIR  := build
SRC_DIR  := glutils
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

NAME   := glutils
SHARED := lib$(NAME).so
STATIC := lib$(NAME).a

.DEFAULT_GOAL := $(STATIC)

.PHONY: shared
shared: $(SHARED)
$(SHARED): $(OBJS)
	@echo Linking shared library...
	@$(CXX) -shared $^ -o $@

.PHONY: static
static: $(STATIC)
$(STATIC): $(OBJS)
	@echo Linking static library...
	@ar rcs $@ $^

# Test configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

TEST_OBJS := $(subst .cpp,.o, $(TEST_SRCS))

$(TEST_DIR)/%.o: %.cpp $(.DEFAULT_GOAL)
	@echo Building $(shell echo $< | sed s/.cpp// | sed 's|.*/||' )...
	@$(COMPILE) $< -L. -l$(NAME) $(LIBS) -o $@

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
clean: cleandocs cleantests
	@rm -rf $(SHARED) $(STATIC) $(OBJ_DIR)

.PHONY: cleantests
cleantests:
	@rm -f $(TEST_OBJS)

.PHONY: cleandocs
cleandocs:
	@rm -rf $(DOC_DIR)/html

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
