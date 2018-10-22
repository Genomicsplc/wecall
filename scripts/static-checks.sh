#!/usr/bin/env bash
# All content Copyright (C) 2018 Genomics plc
set -e -u -x -o pipefail

# requires: cppcheck
if ! command -v cppcheck > /dev/null; then
    echo "Cannot find 'cppcheck', try 'sudo apt-get install cppcheck'."
    exit 1
fi

# reformatting
if python "$WECALL_SCRIPTS/clang-format-check.py"; then
    echo -e "\x1b[1;32mReformatting check passed.\x1b[0m"
else
    echo -e "\x1b[1;31mReformatting check failed.\x1b[0m"
    exit 1
fi

# check the entire codebase & record the return code
# Note: this doesn't check included files from the standard library or boost
CPPCOPTS=(-j "$(nproc)" --force --quiet --inline-suppr -UDEBUG --language=c++)
CPPCENABLES='--enable=warning,style,performance,portability,information,missingInclude'
CPPCSUPPRESSIONS=(--suppress=*:*/dependencies/samtools/include/*)
CPPCINCLUDES=(-I${WECALL_CPP_SOURCE} -I${WECALL_BUILD}/dependencies/samtools/include)
cppcheck ${CPPCOPTS[*]} "$CPPCENABLES" ${CPPCSUPPRESSIONS[*]} ${CPPCINCLUDES[*]} ${WECALL_CPP} 2>&1 | python ${WECALL_SCRIPTS}/count-errors.py
ERR=$?

if [ 0 -eq $ERR ]; then
    echo -e "\x1b[1;32mSCA passed.\x1b[0m"
else
    echo -e "\x1b[1;31mSCA failed.\x1b[0m"
fi

# propagate any errors to the test driver
exit $ERR

