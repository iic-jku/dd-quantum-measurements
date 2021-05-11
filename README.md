[![arXiv](https://img.shields.io/static/v1?label=arXiv&message=21XX.XXXXX&color=inactive&style=plastic)](https://arxiv.org/abs/21XX.XXXXX)


# Decision Diagrams for Quantum Measurements with Classical Shadows

We consider the problem of estimating quantum observables on a collection of qubits, given as a linear combination of Pauli operators, with shallow quantum circuits consisting of single-qubit rotations. 
We introduce estimators based on randomized measurements, which use decision diagrams to sample from probability distributions on measurement bases. 
This approach generalizes previously known uniform and locally-biased randomized estimators. 
The decision diagrams are constructed given target quantum operators and can be optimized considering different strategies. 
We show numerically that the estimators introduced here can produce more precise estimates on some quantum chemistry Hamiltonians, compared to previously known randomized protocols and Pauli grouping methods.   

The full paper *Decision Diagrams for Quantum Measurements with Classical Shadows* by *Stefan Hillmich, Charles Hadfield, Rudy Raymond, Antonio Mezzacapo, and Robert Wille* is available at [*ARXIV*](https://arxiv.org/abs/21XX.XXXXX).

## Installation & Usage

Clone this repository and open a shell in the cloned directory.

The following examples includes the creation of a [virtual environment](https://docs.python.org/3/library/venv.html), skip this step if you already set up an enviroment.
Prior to installing the necessaries dependencies, we recommend updating `pip`, `setuptools`, `wheel`.

```shell
dd-quantum-measurements/ $ python3 -m venv venv
dd-quantum-measurements/ $ . venv/bin/activate
(venv) dd-quantum-measurements/ $ pip install -U pip setuptools wheel
(venv) dd-quantum-measurements/ $ pip install -r requirements.txt
(venv) dd-quantum-measurements/ $ jupyter notebook measurements.ipynb
```

After the last command, a browser will open with the notebook (which contains the code alongside comments on how to use it).

## Reference

```bibtex
@misc{hillmich2021decisiondiagrams
      title={Decision Diagrams for Quantum Measurements with Classical Shadows}, 
      author={Stefan Hillmich and Charles Hadfield and Rudy Raymond and Antonio Mezzacapo and Robert Wille},
      year={2021}
}
```

