![GitHub Workflow Status](https://img.shields.io/github/workflow/status/bodegalab/tecount/Upload%20Python%20Package?logo=github)
[![PyPI](https://img.shields.io/pypi/v/tecount?logo=python)](https://pypi.org/project/tecount/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat&logo=anaconda)](https://bioconda.github.io/recipes/tecount/README.html)

# TEcount

A package to count reads mapping on transposable elements (TEs) subfamilies, families and classes in bulk RNA-seq.

## Features

TEcount counts high-throughput sequencing reads aligned on TEs on a reference genome.
It reports TE counts by subfamily, family and class on separate outputs.
In case of reads aligning on multiple TE loci, it counts only one alignment occurrence for each feature (i.e. subfamily, family or class).
It can use single or paired end BAM files, count in strand-specific manner and filter by a minimum read-TE overlap, among other features that can be discovered by running `TEcount --help`.

## Install

### Using conda  (recommended)
We recommend using conda, as it will install all the required packages along TEcount.
```bash
conda create -n tecount -c conda-forge -c bioconda tecount
```

### Using pip
If for any reason it's not possible or desiderable to use conda, it can be installed with pip and the following requirements must be installed manually: `python>=3.7`, `samtools>=1.14` and `bedtools>=2.30.0`.
```bash
pip install tecount
```

## Usage
Requires a BAM file sorted by coordinates and a bed6+3 file with subfamily, family and class in columns 7, 8 and 9.

```bash
TEcount -b <sorted_by_coord.bam> -r <rmsk.bed>
```

Example using test dataset:
```bash
TEcount -b ./test/Aligned.sortedByCoord.out.bam  -r ./test/rmsk.GRCm38.chr19.bed.gz
```

See all available parameters with:
```bash
TEcount --help
```

## Acknowledgements
Based on the transposable elements RNA-seq quantification strategy implemented in [1](#1).

## References
<a name="1"></a>1. Marasca, F., Sinha, S., Vadalà, R. et al. LINE1 are spliced in non-canonical transcript variants to regulate T cell quiescence and exhaustion. Nat Genet 54, 180–193 (2022). https://doi.org/10.1038/s41588-021-00989-7
