# SBOL Docker Makefile
# Carol Gao <carol_gyz@gmail.com>



_UNAME_S := $(shell uname -s)

USE_DOCKER ?= yes

SHELL := bash
CWD := $(shell pwd)

TARGETS ?= python java

DOCKER_IMAGE_LATEST_TAG="sbol:latest"

.PHONY: build
build:
	@echo "Building docker image for SBOL."

	docker build --pull --rm --tag $(DOCKER_IMAGE_LATEST_TAG) .

.PHONY: dev
dev:
	@echo "Run docker image interactively with home mounted"
	docker run -it --mount type=bind,source=${HOME},target=/homeland \
	--entrypoint "/bin/bash" \
	$(DOCKER_IMAGE_LATEST_TAG)
