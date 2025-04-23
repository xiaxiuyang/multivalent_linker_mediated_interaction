# multivalent_linker_mediated_interaction
[![arXiv](https://img.shields.io/badge/arXiv-2410.23784-b31b1b.svg)](https://arxiv.org/abs/2410.23784)

Designed self-assembly of programmable colloidal atom-electron equivalents  
*Xiuyang Xia, Yuhan Peng, Ka Ki Li, Ran Ni*  
**Last updated:** 24 Apr 2025

---

## Overview
This repository contains:

1. **`pair_interaction`** – a pure-Python package that computes pairwise PAE-PAE interactions mediated by multivalent electron-equivalent (EE) linkers for *arbitrary valency*.
2. **`cpp/`** – a C++ Monte Carlo (MC) simulator capable of many-body EE-mediated interactions, suitable for large-scale off-lattice simulations.
3. Jupyter/command-line **examples** to reproduce all figures and key results in the paper.

The code is fully open-source under the MIT license.

---

## Requirements
| Component | Minimum version | Notes |
|-----------|-----------------|-------|
| Python    | 3.8             | Anaconda 3 recommended |
| NumPy     | 1.23            | |
| SciPy     | 1.10            | |
| C++       | 17              | Tested with GCC 11/Clang 15 |
| CMake     | 3.18            | For building the C++ simulator |

---

## Citation

Please cite the following paper if you find this useful in your research:
```
@misc{xia2024designedselfassemblyprogrammablecolloidal,
      title={Designed self-assembly of programmable colloidal atom-electron equivalents}, 
      author={Xiuyang Xia and Yuhan Peng and Ka Ki Li and Ran Ni},
      year={2024},
      eprint={2410.23784},
      archivePrefix={arXiv},
      primaryClass={cond-mat.soft},
      url={https://arxiv.org/abs/2410.23784}, 
}
```
