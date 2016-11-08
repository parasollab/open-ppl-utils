# Externals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

NONSTD_DIR  := ../nonstd
NONSTD_INCL := -I$(NONSTD_DIR)
NONSTD_LIBS := -L$(NONSTD_DIR) -lnonstd

GL_INCL     :=
GL_LIBS     := -lGLU -lGL

# Compiler configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

CXX   := g++
FLAGS := -std=c++11 -fPIC -g3 -O2
WARN  := -Wall -Werror -Wno-parentheses

INCL  := -I. $(NONSTD_INCL) $(GL_INCL)

COMPILE := $(CXX) $(FLAGS) $(WARN) $(INCL)

# Linker configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

LIBS := $(NONSTD_LIBS) $(GL_LIBS)

# Directory layout ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

OBJ_DIR := build
SRC_DIR := glutils
DOC_DIR := docs

$(OBJ_DIR):
	@mkdir $(OBJ_DIR)

vpath %.h   $(SRC_DIR)
vpath %.cpp $(SRC_DIR)

# Source file configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

SRCS := $(shell find $(SRC_DIR) -name '*.cpp')

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

.PHONY: all
all: $(SHARED) $(STATIC)

# Docs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

HEADERS := $(shell find . -name '*.h')

.PHONY: docs
docs: $(SRCS) $(HEADERS)
	@echo Building docs...
	@doxygen $(DOC_DIR)/Doxyfile

# Clean ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

.PHONY: clean
clean: cleandocs
	@rm -rf $(SHARED) $(STATIC) $(OBJ_DIR)

.PHONY: cleandocs
cleandocs:
	@rm -rf $(DOC_DIR)/html

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
