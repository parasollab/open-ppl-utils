# This is a weird makefile that runs a bash script to set up bullet.
# Some day we should unify these.

# Script Launcher ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

BUILD_DIR := build
LIB_DIR   := $(BUILD_DIR)/install/lib
$(LIB_DIR):
	@echo Building bullet library...
	@./parasol_setup.sh

# Cleanup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

.PHONY: clean
clean:
	@echo Cleaning bullet library...
	@rm -rf $(BUILD_DIR)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
