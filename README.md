# Rphenograph

R implementation of the [PhenoGraph](https://doi.org/10.1016/j.cell.2015.05.047) algorithm — a clustering method designed for high-dimensional single-cell data analysis. It works by creating a graph representing phenotypic similarities between cells by calculating the Jaccard coefficient between nearest-neighbor sets, and then identifying communities using the [Louvain method](https://sites.google.com/site/findcommunities/).

This is a maintained fork of [JinmiaoChenLab/Rphenograph](https://github.com/JinmiaoChenLab/Rphenograph) with fixes for modern R (>= 4.1) and igraph (>= 2.0).

## Installation

```r
# Install from GitHub (recommended)
pak::pkg_install("i-cyto/Rphenograph")

# Or with remotes
remotes::install_github("i-cyto/Rphenograph")
```

### Apple Silicon (M1/M2/M3/M4) notes

The package compiles cleanly on ARM64 Macs. If you hit issues:

1. **Install/update Xcode Command Line Tools:**
   ```sh
   xcode-select --install
   ```

2. **Check `~/.R/Makevars` for stale x86 flags.** If the file contains `-mtune=core2` or `-march=core2`, remove or comment out those lines. These are Intel-only flags that break ARM64 compilation.

3. **Ensure R itself is ARM64-native.** In R, run `R.version$arch` — it should say `aarch64`, not `x86_64`. If it says `x86_64`, download the native ARM64 build from <https://cran.r-project.org/bin/macosx/>.

### Windows (no compiler needed)

```r
install.packages(
  "https://github.com/i-cyto/Rphenograph/releases/download/Rphenograph_0.99.1.9004/Rphenograph_0.99.1.9004.zip",
  repos = NULL, type = "win.binary"
)
```

## Usage

```r
library(Rphenograph)
library(igraph)

iris_unique <- unique(iris)
data <- as.matrix(iris_unique[, 1:4])
out <- Rphenograph(data, k = 45)

# Cluster memberships
membership(out[[2]])

# Modularity
modularity(out[[2]])
```

## Reference

Levine JH, Simonds EF, Bendall SC, Davis KL, Amir ED, Tadmor MD, et al. Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like Cells that Correlate with Prognosis. *Cell*, 2015. <https://doi.org/10.1016/j.cell.2015.05.047>
