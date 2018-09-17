#!/usr/bin/env bash
# All content Copyright (C) 2018 Genomics plc
set -e -u -x -o pipefail
PROJECT_SCRIPT=$( cd "$(dirname "$0")" && pwd -P)
PROJECT_HOME="${PROJECT_SCRIPT}/../"
export ECHIDNA_TEST_RESULTS="${PROJECT_HOME}/test-results"
export ECHIDNA_BIN="${PROJECT_HOME}/build"

mkdir -p "${ECHIDNA_TEST_RESULTS}"

"${ECHIDNA_BIN}/unittest" | tee >> "${ECHIDNA_TEST_RESULTS}/unittest.log"
set +e +u
. ${PROJECT_HOME}/env-wecall/bin/activate
set -e -u

pytest ${@} \
    --flakes \
    --junit-xml="${ECHIDNA_TEST_RESULTS}/acceptance-test.xml" \
    --cov wecall \
    --cov wecall_test_drivers \
    --cov-report term:skip-covered \
    --cov-report xml:"${ECHIDNA_TEST_RESULTS}/acceptance-test-coverage.xml" \
    --cov-report html:"${ECHIDNA_TEST_RESULTS}/acceptance_test_coverage_html" \
    --no-cov-on-fail
