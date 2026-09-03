# conda recipe

The recipe that builds SCICoNE — the C++ executables and the `scicone` Python
package — into a single conda package for [bioconda](https://bioconda.github.io/).

> **The copy in bioconda-recipes is the authoritative one.**
> [`recipes/scicone/`](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/scicone)
> is what builds the package people install, and bioconda's maintainers and bots
> edit it directly — bumping the build number to rebuild against new dependency
> versions, adjusting pinnings, extending the Python version matrix. Those edits
> do not come back here.
>
> This copy exists so the recipe can be built from the working tree, locally and
> in CI, which catches a dependency or build-step problem before it reaches
> bioconda's review queue rather than days later. Treat it as the starting point
> for the next submission, not as a record of what bioconda currently has: check
> `recipes/scicone/` there before assuming the two agree.

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

4. Diff this copy against `recipes/scicone/` in bioconda-recipes and carry over
   anything they changed since the last submission, so their edits are not
   reverted.
5. Open a pull request against
   [bioconda-recipes](https://github.com/bioconda/bioconda-recipes) that copies
   `meta.yaml` and `build.sh` to `recipes/scicone/`. Bioconda's CI builds them
   for linux-64 and osx-64 and merges on green.

For a change to the recipe alone, leave `version` where it is and increment
`build: number:` instead.

Any change to what the package needs — a new import in `pyscicone`, a new build
dependency, a new executable — has to land in **both** copies. The one here
keeps CI honest; the one in bioconda-recipes is what users actually install.

## What the package installs

* `$PREFIX/bin/scicone-simulation`, `scicone-breakpoint_detection`,
  `scicone-inference`, `scicone-score`. They are prefixed because they land in a
  shared `bin/` directory, where names like `inference` and `score` would collide.
* The `scicone` Python package, which finds those executables on `PATH`.

`scicone-tests` is built and run during the build but deliberately not installed:
it resolves its input data relative to the source tree (see `scicone/Config.h.in`),
which no longer exists once the package is built.
