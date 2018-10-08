#!/usr/bin/env bash
# All content Copyright (C) 2018 Genomics plc
set -e -u -x -o pipefail
echo $0
PROJECT_SCRIPT=$( cd "$(dirname "$0")" && pwd -P)
export PROJECT_HOME="${PROJECT_SCRIPT}/../"
WECALL_BUILD="${PROJECT_HOME}/build"

command -v pdflatex >/dev/null 2>&1 || { echo >&2 "Skipping Document generation - No pdflatex install found." ; exit 0 ; }

# generate some data
"$PROJECT_HOME/build/weCall" --help | python "$PROJECT_HOME/scripts/help_to_latex.py" >  "$PROJECT_HOME/doc/wecall-params.tex"

# make
mkdir -p "$WECALL_BUILD/doc"
cd "$PROJECT_HOME/doc"
pdflatex -interaction=nonstopmode -halt-on-error -output-directory "$WECALL_BUILD/" "$PROJECT_HOME/doc/weCall-userguide.tex"
pdflatex -interaction=nonstopmode -halt-on-error -output-directory "$WECALL_BUILD/" "$PROJECT_HOME/doc/weCall-userguide.tex"

