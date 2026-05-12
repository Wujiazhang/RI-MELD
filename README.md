# RI-MELD

RI-MELD is a modified version of MELD for protein structure prediction. It improves the original MELD framework by introducing a rule-importance-based strategy to adaptively evaluate and update external restraints during molecular dynamics simulations.

RI-MELD is designed to improve the use of noisy, sparse, or uncertain restraint information in restraint-guided protein folding simulations.

## Method Overview

RI-MELD is based on the MELD framework, which combines molecular dynamics simulations with external restraints to enhance protein structure sampling.

In the original MELD-based adaptive scheme, the number of active restraints is sampled during the simulation, while a reward term is used to encourage the activation of additional restraints.

In RI-MELD, this reward term is adaptively updated according to the rule-importance criterion. The rule importance evaluates whether recent changes in the active restraint set improve or worsen the restraint energy. Based on this evaluation, RI-MELD dynamically adjusts the restraint activation process during sampling.

This design reduces the dependence on manually selected reward parameters and provides a more adaptive way to use external restraint information.

## Input Files

### `sequence.dat`

This file contains the amino acid sequence of the target protein.

### `ss.dat`

This file contains the predicted secondary-structure sequence of the target protein.

The secondary-structure sequence is used to construct secondary-structure-related restraints, such as helix restraints, strand restraints, and other coarse physical insights.

### `evfold.dat`

This file contains residue-residue contact information predicted from co-evolutionary signals.

The `evfold.dat` file is used in the EvFold restraint setting to provide additional contact restraints for protein folding simulations.

## Benchmarks

### CPIs

The `CPIs/` folder contains input files and running scripts for simulations using coarse physical insights-based restraints.

Files in this folder include:

```text
CPI_adapt.py
job.sh
sequence.dat
ss.dat
```

### EvFold

The `Evfold/` folder contains input files and running scripts for simulations using EvFold contact restraints together with secondary-structure information.

Files in this folder include:

```text
EVmeld_adapt.py
evfold.dat
job.sh
sequence.dat
ss.dat
```

## Running the Simulations

Example running commands are provided in the `job.sh` file under each benchmark folder.

For the CPIs benchmark, please refer to:

```bash
CPIs/job.sh
```

For the EvFold benchmark, please refer to:

```bash
Evfold/job.sh
```

Users can modify the corresponding input files and running scripts according to their own protein systems and computing environment.

## Core Modification

The core RI-MELD modification is provided in:

```text
modify.py
```

The implementation is based on the original MELD source code. The main modification is made to the `_run_param_mc()` function in `remd_runner.py`.

Compared with the original MELD adaptive sampling scheme, RI-MELD updates the reward term dynamically using the rule-importance criterion. This allows the number of active restraints to be adjusted according to the energetic contribution of recent restraint updates during the simulation.

## Original MELD Source Code

The RI-MELD implementation is modified based on the original MELD source code.

The original MELD repository is available at:

```text
https://github.com/maccallumlab/meld
```

Users should first install the original MELD package and then apply the modification provided in `modify.py` to reproduce the RI-MELD implementation.

## Notes

This repository is intended to provide the key input files, example running scripts, and core algorithmic modification used in RI-MELD.

The full simulation environment may depend on the user's local installation of MELD, OpenMM, CUDA, MPI, and other required scientific computing packages. Please make sure that the original MELD environment is correctly installed before applying the RI-MELD modification.
