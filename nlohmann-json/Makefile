################################################################################
# This Makefile fetches the nlohmann/json library from github and copies the
# single-include file into the install directory. The library is header-only so
# there is nothing to build. We use an install directory anyway to decouple this
# Makefile from PMPL's.
#
# This is a modern C++ json library that is vastly more usable than any other I
# could find. For reference see https://github.com/nlohmann/json.
################################################################################


# Directory Layout ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Directories for our build.
SOURCE_DIR  := $(CURDIR)/source
INSTALL_DIR := $(CURDIR)/install/nlohmann


# Config ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

TAG      := v3.7.3
GIT_REPO := https://github.com/nlohmann/json.git
GIT_OPS  := --depth 1 --branch $(TAG)


# Recipes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Build the CGAL library and generate a local 'install' here.
$(INSTALL_DIR): $(SOURCE_DIR)
	@echo Building nlohmann-json library...
	@rm -rf $(INSTALL_DIR)
	@mkdir -p $(INSTALL_DIR)
	@cp $(SOURCE_DIR)/single_include/nlohmann/json.hpp $(INSTALL_DIR)

# Fetch the CGAL git repo.
$(SOURCE_DIR):
	@echo Fetching nlohmann-json git repo...
	@git clone $(GIT_REPO) $(GIT_OPS) $(SOURCE_DIR)


# Cleanup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

.PHONY: clean
clean:
	@echo Cleaning nlohmann-json library...
	@rm -rf $(INSTALL_DIR)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
