PREFIX=/usr/local
.PHONY: vendor wecall
BUILD=build

export ECHIDNA_TEST_RESULTS=results
export ECHIDNA_BIN=$(BUILD)

all: vendor wecall

vendor:
	$(MAKE) -C vendor 

wecall: vendor
	mkdir -p $(BUILD) && cd $(BUILD) && \
	cmake -D CMAKE_INSTALL_PREFIX=$(PREFIX) ../cpp
	$(MAKE) -C $(BUILD)
	cp vendor/samtools/samtools $(BUILD)
	cp vendor/tabix/tabix $(BUILD)
	cp vendor/tabix/bgzip $(BUILD)

test-unit: vendor wecall
	build/unittest
	build/iotest

env-wecall:
	python3 -m venv env-wecall
	bash -c "source env-wecall/bin/activate  && cd python && pip install ."
	bash -c "source env-wecall/bin/activate && cd test-drivers && pip install ."

test-acceptance: wecall env-wecall
	bash -c " source env-wecall/bin/activate && scripts/run-tests.sh test"


install: vendor wecall
	$(MAKE) -C build install

package: vendor wecall
	$(MAKE) -C build package

clean:
	$(MAKE) -C vendor clean
	-rm -rf build
	-rm -f cpp/src/version/version.cpp
	-rm -f doc/weCall-userguide.aux
	-rm -f doc/weCall-userguide.out
	-rm -f doc/weCall-userguide.toc
