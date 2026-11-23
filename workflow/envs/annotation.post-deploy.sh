#!/usr/bin/env bash
set -euo pipefail

# Preserve support for amino-acid contig databases from the custom Anvi'o fork.
python -m pip install --no-input --no-deps --no-build-isolation \
    git+https://github.com/kekananen/anvio@21aad5f33a2c3c2c1e48f428dda08a0ef86673c0
python -c "import anvio"
