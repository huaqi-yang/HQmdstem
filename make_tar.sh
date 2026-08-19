#!/bin/bash
# Package HQmdstemkit for migration:
#   tar -zcvf HQmdstem.tar.gz -> copy to another machine -> bash install.sh
set -e
cd "$(dirname "$0")"
tar -zcvf HQmdstem.tar.gz \
    --exclude='__pycache__' \
    --exclude='*.pyc' \
    --exclude='HQmdstem.tar.gz' \
    .
echo ""
echo "Done: $(pwd)/HQmdstem.tar.gz"
echo "On the new machine:"
echo "  tar -xzf HQmdstem.tar.gz"
echo "  cd HQmdstem && bash install.sh"
echo "  HQmdstemkit.sh"