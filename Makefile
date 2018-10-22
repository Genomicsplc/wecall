# MAKEFILE
#
# @link        https://github.com/genomicsplc/wecall
# ------------------------------------------------------------------------------

# List special make targets that are not associated with files
.PHONY: help format clean vendor wecall test-unit env-wecall test-acceptance install package clean

PREFIX=/usr/local
BUILD=target/build

export WECALL_TEST_RESULTS=results
export WECALL_BIN=$(BUILD)

all: vendor wecall

help:
	@echo ""
	@echo "weCall Makefile."
	@echo "The following commands are available:"
	@echo ""
	@echo "    make                 : Build weCall"
	@echo "    make test-unit       : Execute unit tests"
	@echo "    make test-acceptance : Execute acceptance tests"
	@echo "    make install         : Install the executable"
	@echo "    make package         : Build DEB, RPM and TGZ packages"
	@echo "    make clean           : Remove any build artifact"
	@echo ""

vendor:
	$(MAKE) --directory=vendor 

wecall: vendor
	mkdir -p $(BUILD) \
	&& cd $(BUILD) \
	&& cmake -D CMAKE_INSTALL_PREFIX=$(PREFIX) ../../cpp
	$(MAKE) --directory=$(BUILD)
	cp vendor/samtools/samtools $(BUILD)
	cp vendor/tabix/tabix $(BUILD)
	cp vendor/tabix/bgzip $(BUILD)

test-unit: vendor wecall
	$(BUILD)/unittest
	$(BUILD)/iotest

env-wecall:
	python3 -m venv env-wecall
	bash -c "source env-wecall/bin/activate && pip install wheel"
	bash -c "source env-wecall/bin/activate && cd python && pip install ."
	bash -c "source env-wecall/bin/activate && cd test-drivers && pip install ."

test-acceptance: wecall env-wecall
	bash -c " source env-wecall/bin/activate && scripts/run-tests.sh test"

install: vendor wecall
	$(MAKE) --directory=$(BUILD) install

package: vendor wecall
	$(MAKE) --directory=$(BUILD) package

clean:
	$(MAKE) --directory=vendor clean
	-rm -rf $(BUILD)
	-rm -f cpp/src/version/version.cpp
	-rm -f doc/weCall-userguide.aux
	-rm -f doc/weCall-userguide.out
	-rm -f doc/weCall-userguide.toc
