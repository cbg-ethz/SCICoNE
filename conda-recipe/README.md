# conda recipe

The recipe that builds SCICoNE — the C++ executables and the `scicone` Python
package — into a single conda package for [bioconda](https://bioconda.github.io/).

| File | Purpose |
| --- | --- |
| `meta.yaml` | The recipe. Source is a released tarball, which is what bioconda requires. |
| `build.sh` | Builds and installs the executables, runs the C++ unit tests, installs the Python package. |
| `build-local.sh` | Builds the recipe against the working tree instead of a release, for iterating. |

## Building it locally

```bash
conda install -c conda-forge conda-build
conda-recipe/build-local.sh
```

This does not need a release to exist: it copies the recipe to a temporary
directory, replaces the `source:` section with the local checkout, and builds
that. Everything else — the compilers, the dependency resolution, the test
phase — is what bioconda will do. To install and try the result:

```bash
conda create -n scicone-test -c local -c conda-forge -c bioconda scicone
```

## Releasing a new version

1. Bump `__version__` in `pyscicone/scicone/__init__.py`; `setup.py` reads it
   from there.
2. Tag and push: `git tag v1.0.0 && git push origin v1.0.0`, then publish the
   release on GitHub.
3. Update `meta.yaml`: set `version`, reset `build: number:` to 0, and set the
   checksum of the release tarball:

   ```bash
   curl -sL https://github.com/cbg-ethz/SCICoNE/archive/refs/tags/v1.0.0.tar.gz | sha256sum
   ```

4. Open a pull request against
   [bioconda-recipes](https://github.com/bioconda/bioconda-recipes) that copies
   `meta.yaml` and `build.sh` to `recipes/scicone/`. Bioconda's CI builds them
   for linux-64 and osx-64 and merges on green.

For a change to the recipe alone, leave `version` where it is and increment
`build: number:` instead.

## What the package installs

* `$PREFIX/bin/scicone-simulation`, `scicone-breakpoint_detection`,
  `scicone-inference`, `scicone-score`. They are prefixed because they land in a
  shared `bin/` directory, where names like `inference` and `score` would collide.
* The `scicone` Python package, which finds those executables on `PATH`.

`scicone-tests` is built and run during the build but deliberately not installed:
it resolves its input data relative to the source tree (see `scicone/Config.h.in`),
which no longer exists once the package is built.
