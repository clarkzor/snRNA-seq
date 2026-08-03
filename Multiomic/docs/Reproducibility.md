# Reproducibility

## Package installation

Package installation is separated from analysis code. See:

```r
source("scripts/install_dependencies.R")
```

Custom Xenopus EnsDb and BSgenome packages must be installed separately.

## renv

After confirming the working R environment, initialize and snapshot it from the repository root:

```r
install.packages("renv")
renv::init()
renv::snapshot()
```

Commit the generated `renv.lock`. Another machine can then restore package versions with:

```r
renv::restore()
```

At finalization the pipeline also writes `validation/package_versions.csv` and
`validation/sessionInfo.txt`.

## Current plotting-stack note

The analysis has encountered patchwork/ggplot2 operator compatibility errors in Signac
coverage plotting. Coverage is therefore optional and isolated from the core workflow.
Record the versions that successfully reproduce coverage plots before enabling them in a
frozen `renv.lock`.
