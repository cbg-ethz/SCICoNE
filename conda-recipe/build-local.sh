#!/usr/bin/env bash
# Build the conda package from the working tree.
#
# conda-recipe/meta.yaml points at a released source tarball, which is what
# bioconda needs but is useless while iterating. This copies the recipe to a
# temporary directory, swaps the source for the local checkout, and builds that.
#
# Usage: conda-recipe/build-local.sh [extra conda-build arguments]
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo=$(cd "${here}/.." && pwd)
work=$(mktemp -d)
trap 'rm -rf "${work}"' EXIT

cp "${here}/meta.yaml" "${here}/build.sh" "${work}/"

# bioconda's global configuration defines these; a bare conda-build does not,
# and `{{ stdlib('c') }}` does not resolve without them.
cat > "${work}/conda_build_config.yaml" <<'YAML'
c_stdlib:
  - sysroot
c_stdlib_version:
  - "2.17"
YAML

python - "${work}/meta.yaml" "${repo}" <<'PY'
import re
import sys

meta_path, repo = sys.argv[1], sys.argv[2]
meta = open(meta_path).read()
meta, n = re.subn(r"^source:\n(?:[ \t]+\S.*\n)+", "source:\n  path: %s\n" % repo,
                  meta, count=1, flags=re.MULTILINE)
if n != 1:
    sys.exit("could not rewrite the source section of %s" % meta_path)
open(meta_path, "w").write(meta)
PY

conda build "${work}" -c conda-forge -c bioconda "$@"
