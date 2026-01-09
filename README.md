# [DEPRECATED]
This toolbox is now being developed and maintained in https://github.com/ZETA-Toolbox.

# ZETA

**Z**onotope-based **E**s**T**imation and f**A**ult diagnosis of discrete-time systems

- Check our [preprint](https://arxiv.org/abs/2504.06467) submitted to the 64th IEEE Conference on Decision and Control
- by Brenner S. Rego, Joseph K. Scott, Davide M. Raimondo, Marco H. Terra and Guilherme V. Raffo

ZETA is a MATLAB library featuring implementations of set representations based on zonotopes, namely: 
- [Zonotopes](https://doi.org/10.1007/BF02684450),
- [Constrained zonotopes](https://web.mit.edu/braatzgroup/Scott_Automatica_2016.pdf), and
- [Line zonotopes](https://arxiv.org/abs/2401.10239),
  
in addition to a basic implementation of [interval arithmetic](https://en.wikipedia.org/wiki/Interval_arithmetic). 

This library has capabilities starting from the basic set operations with
these sets, including propagations through nonlinear functions using
various approximation methods. The features of ZETA allow for
reachability analysis and state estimation of discrete-time linear,
nonlinear, and descriptor systems, in addition to active fault diagnosis
of linear systems. Efficient order reduction methods are
also implemented for the respective set representations. 

- This library requires [YALMIP](https://yalmip.github.io/) to work properly.

- ZETA supports [GUROBI](https://www.gurobi.com/) for efficient solutions of the underlying optimization problems.

- [MPT](https://www.mpt3.org/) is required in a few demonstration examples, while allowing for extended plotting capabilities.

The library is initialized by running `zeta_init.m`.

### Citation

If you have used this library as part of your work, cite it as follows:
```
@misc{Rego2025ZETA,
title={{ZETA: a library for Zonotope-based EsTimation and fAult diagnosis of discrete-time systems}}, 
author={Brenner S. Rego and Joseph K. Scott and Davide M. Raimondo and Marco H. Terra and Guilherme V. Raffo},
year={2025},
note={arXiv preprint. arXiv:2504.06467}
}
```
