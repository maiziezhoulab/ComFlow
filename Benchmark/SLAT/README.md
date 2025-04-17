[![stars-badge](https://img.shields.io/github/stars/gao-lab/SLAT?logo=GitHub&color=yellow)](https://github.com/gao-lab/SLAT/stargazers)
[![pypi-badge](https://img.shields.io/pypi/v/scslat)](https://pypi.org/project/scslat)
[![Downloads](https://static.pepy.tech/badge/scSLAT)](https://pepy.tech/project/scSLAT)
[![build-badge](https://github.com/gao-lab/SLAT/actions/workflows/build.yml/badge.svg)](https://github.com/gao-lab/SLAT/actions/workflows/build.yml)
[![docs-badge](https://readthedocs.org/projects/slat/badge/?version=latest)](https://slat.readthedocs.io/en/latest/?badge=latest)
[![license-badge](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


# scSLAT: single cell spatial alignment tools

**scSLAT** package implements the **SLAT** (**S**patial **L**inked **A**lignment **T**ool) model to align single cell spatial omics data.

Please check our paper ***Spatial-linked alignment tool (SLAT) for aligning heterogenous slices*** on [*Nature Communications*](https://www.nature.com/articles/s41467-023-43105-5).

![Model architecture](docs/_static/Model.png)

## Directory structure

```
.
├── scSLAT/                  # Main Python package
├── env/                     # Extra environment
├── data/                    # Data files
├── evaluation/              # SLAT evaluation pipeline
├── benchmark/               # Benchmark pipeline
├── case/                    # Case studies in paper
├── docs/                    # Documentation files
├── resource/                # Other useful resource
├── pyproject.toml           # Python package metadata
└── README.md
```

## [Tutorial](https://slat.readthedocs.io/en/latest/tutorials.html)
Tutorial of `scSLAT` is [here](https://slat.readthedocs.io/en/latest/tutorials.html), if you have any question please open an issue on github

<img src='docs/_static/imgalignment.gif' width='400'>


## Installation

### PyPI
> **Note**
> Installing `scSLAT` within a new [conda environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) is recommended.

> **Warning**
> Old NVIDIA driver may raise error.

First, we create a clean environment and install `scSLAT` from PyPI. Then we also need install dependencies for `pyg` manually via `install_pyg_dependencies`.

```bash
conda create -n scSLAT python=3.11 -y && conda activate scSLAT
pip install scSLAT
install_pyg_dependencies
```

### Docker
You can pull the docker image directly from [Docker Hub](https://hub.docker.com/repository/docker/huhansan666666/slat) or refer to the [`Dockerfile`](Dockerfile) to build it.

``` bash
docker pull huhansan666666/slat:latest
```

### Development version
For development purpose, clone this repo and install:

```bash
git clone git@github.com:gao-lab/SLAT.git && cd SLAT
pip install -e ".[dev,docs]"
install_pyg_dependencies
```


## Appendix
1. The unit of runtime is seconds (s).
2. Executed on GPU.
