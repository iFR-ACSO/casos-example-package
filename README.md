<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://github.com/iFR-ACSO/.github/blob/main/assets/logo-example-inverted.png">
  <source media="(prefers-color-scheme: light)" srcset="https://github.com/iFR-ACSO/.github/blob/main/assets/logo-example-trans.png">
  <img alt="CaΣoS example package">
</picture>

----

The **CaΣoS example package** contains tutorials, application examples, and implementations from the literature. The purpose is to provide an easy on-boarding with [CaΣoS](https://github.com/ifr-acso/casos), a nonlinear sum-of-squares optimization suite based on the symbolic framework of CasADi. New examples are continuously added.

### Requirements
The example package itself does not not need to be installed. Only a [stable version](https://github.com/ifr-acso/casos/releases/latest) of CaΣoS and a [supported conic solver](https://github.com/ifr-acso/casos/wiki#conic-solvers) is needed. Follow the instructions on the [Getting started](https://github.com/ifr-acso/casos/wiki#getting-started) page.

### Folder Structure

```text
Tutorials/						# Folder contains different tutorials
├── Basics        				# Toolbox basics
├── Pontryagin difference		# Examples for estimating inner approximations of the Pontryagin difference
├── SOS			  				# How to setup different types of SOS problems
└── Conic			  			# How to setup conic problems using sdpsol interface

Systems and Control/ 			# Examples from the systems and controls literature
├── NMPC						# Synthesis of terminal conditions for NMPC, infinitesimal-NMPC law synthesis
├── Stability					# Stability analysis, e.g., region-of-attraction estimation
└── Reachability				# Reachability analysis, e.g., inner-approximation of backwards reachable set
```



### CaΣoS Quick links

- [Getting started](https://github.com/ifr-acso/casos/wiki#getting-started)
- Available [conic solvers](https://github.com/ifr-acso/casos/wiki#conic-solvers)
- Convex and nonconvex [sum-of-squares optimization](https://github.com/ifr-acso/casos/wiki/sum%E2%80%90of%E2%80%90squares-optimization)
- Supported [vector, matrix, and polynomial cones](https://github.com/ifr-acso/casos/wiki/cones)
- Some [practical tipps](https://github.com/ifr-acso/casos/wiki/practical-sos-guide) to sum-of-squares
- [Transitioning](https://github.com/ifr-acso/casos/wiki/transitioning-from-other-toolboxes) from other toolboxes
- Example [code snippets](https://github.com/ifr-acso/casos/wiki/numerical-examples)

### Cite us

If you use CaΣoS, please cite us:

> T. Cunis and J. Olucak, ‘CaΣoS: A nonlinear sum-of-squares optimization suite’, in _2025 American Control Conference_, (Boulder, CO), pp. 1659–1666, 2025 doi: [10.23919/ACC63710.2025.11107794](https://doi.org/10.23919/ACC63710.2025.11107794).

<details>

<summary>Bibtex entry</summary>

```bibtex
@inproceedings{Cunis2025acc,
	author = {Cunis, Torbjørn and Olucak, Jan},
	title = {{CaΣoS}: {A} nonlinear sum-of-squares optimization suite},
	booktitle = {2025 American Control Conference},
	address = {Boulder, CO},
	year = {2025},
	pages = {1659--1666},
	doi = {10.23919/ACC63710.2025.11107794},
}
```

</details>

----
