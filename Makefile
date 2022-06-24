# SBOL Makefile
# Carol Gao <carol_gyz@gmail.com>



_UNAME_S := $(shell uname -s)

USE_DOCKER ?= yes

SHELL := bash
CWD := $(shell pwd)

TARGETS ?= python java


.PHONY: install
install:	## install compiler
	@echo "Installation done."
