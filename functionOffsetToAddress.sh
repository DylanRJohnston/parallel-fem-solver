#!/usr/bin/env bash

set -o pipefail
set -o nounset
set -o errexit

if [[ -z "${1:-}" || -z "${2:-}" ]]; then
    echo -e "Usage: $0 function_name offset\nGiven a function name and an offset, computes hexadecimal address for gobjdump"
    echo -e "Example: ./functionOffsetToAddress.sh main 512"
    exit 1
fi

echo "ibase=obase=16; $(echo "obase=16;ibase=10;$2" | bc) + $(nm fem1d | grep $1 | grep -oE '^\S+' | tr 'a-z' 'A-Z')" | bc | tr 'A-Z' 'a-z'
