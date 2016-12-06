#!/usr/bin/env bash

/usr/local/Cellar/llvm/3.6.2/share/clang/tools/scan-build/scan-build -enable-checker alpha gcc -std=c99 -Wall -Werror -Wextra -pedantic fem1d.c -o fem1d
