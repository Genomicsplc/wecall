#!/usr/bin/env bash
# All content Copyright (C) 2018 Genomics plc
set -e -u -x -o pipefail


set -o errexit -o nounset

clang-format-3.6 -i -style=file $(find "cpp" -name '*.[hc]pp')
