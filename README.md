# A Spectral Bundle Method for Sparse Semidefinite Programs

This repository contains an implementation of the algorithm described in the paper titled [A Spectral Bundle Method for Sparse Semidefinite Programs](https://hsmmoj.github.io/files/SpecBM-SDPs.pdf).

## Table of Contents

- [Prerequisites](#prerequisites)
- [Getting Started](#getting-started)
- [Run](#run)
- [Usage](#usage)
- [References](#references)
- [License](#license)
- [Contact](#contact)

## Prerequisites

Before you begin, ensure that you have the following software installed on your system:

- [MATLAB](https://www.mathworks.com/products/matlab.html)

## Getting Started

To get started with this project, follow these steps:

1. **Install Solvers:**
   - [Sedumi](https://sedumi.ie.lehigh.edu/)
   - [Mosek](https://www.mosek.com/)

2. **Install CDCS:**
   - Download [CDCS](https://github.com/oxfordcontrol/CDCS).
   - Install CDCS and ensure that it's correctly set up.
   - Add "+cdcs_utils" in the CDCS package to your MATLAB path.

## Run

You can run the project by executing the "Experiment1.m" script.

## Usage

This project is designed to solve standard Semidefinite Programming (SDP) problems. If your problem data is structured with a sparsity pattern, this algorithm decomposes the large semidefinite cone constraint into several smaller ones, improving efficiency.

## References

- [Link to Paper (PDF)](https://hsmmoj.github.io/files/SpecBM-SDPs.pdf)

## License

This project is licensed under the [Your License Name or Link](LICENSE). See the [LICENSE](LICENSE) file for details.

## Contact

If you have any questions, feedback, or need assistance, please don't hesitate to contact the project maintainer:

- **Your Name**
- Email: [your@email.com](mailto:your@email.com)
