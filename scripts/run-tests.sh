#!/usr/bin/env bash
# All content Copyright (C) 2018 Genomics plc
set -e -u -x -o pipefail

PROJECT_SCRIPT=$( cd "$(dirname "$0")" && pwd -P)
PROJECT_HOME="${PROJECT_SCRIPT}/../"
export WECALL_TEST_RESULTS="${PROJECT_HOME}/target/test-results"
export WECALL_BIN="${PROJECT_HOME}/target/build"

mkdir -p "${WECALL_TEST_RESULTS}"

"${WECALL_BIN}/unittest" | tee >> "${WECALL_TEST_RESULTS}/unittest.log"
set +e +u
. ${PROJECT_HOME}/env-wecall/bin/activate
set -e -u

pytest ${@} \
--flakes \
--junit-xml="${WECALL_TEST_RESULTS}/acceptance-test.xml" \
--cov wecall \
--cov wecall_test_drivers \
--cov-report term:skip-covered \
--cov-report xml:"${WECALL_TEST_RESULTS}/acceptance-test-coverage.xml" \
--cov-report html:"${WECALL_TEST_RESULTS}/acceptance_test_coverage_html" \
--no-cov-on-fail
