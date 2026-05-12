#!/usr/bin/env python
# encoding: utf-8

import numpy as np
import meld
from meld.remd import ladder, adaptor, leader
import meld.system.montecarlo as mc
from meld import system
from meld.system import patchers
from meld import comm, vault
from meld import parse
from meld import remd
from meld.system import param_sampling
from openmm import unit as u

N_REPLICAS = 30
N_STEPS = 10000
BLOCK_SIZE = 100

hydrophobes = 'AILMFPWV'
hydrophobes_res = ['ALA', 'ILE', 'LEU', 'MET', 'PHE', 'PRO', 'TRP', 'VAL']


def gen_state(s, index):
    state = s.get_state_template()
    state.alpha = index / (N_REPLICAS - 1.0)
    return state


def make_ss_groups(subset=None):
    active = 0
    extended = 0
    sse = []
    ss = open('ss.dat', 'r').readlines()[0]
    for i, l in enumerate(ss.rstrip()):
        # print i,l
        if l not in "HE.":
            continue
        if l not in 'E' and extended:
            end = i
            sse.append((start + 1, end))
            extended = 0
        if l in 'E':
            if i + 1 in subset:
                active = active + 1
            if extended:
                continue
            else:
                start = i
                extended = 1
    print(active, ':number of E residues')
    print(sse, ':E residue ranges')
    return sse, active

def get_dist_restraints_evfold(filename, s, scaler, ramp, seq):
    dists = []
    rest_group = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    for line in lines:
        if not line:
            dists.append(s.restraints.create_restraint_group(rest_group, 1))
            rest_group = []
        else:
            cols = line.split()
            i = int(cols[0]) - 1
            name_i = cols[1]
            j = int(cols[2]) - 1
            name_j = cols[3]
            #

            rest = s.restraints.create_restraint('distance', scaler, ramp,
                                                 r1=0.0 * u.nanometer, r2=0.0 * u.nanometer, r3=0.8 * u.nanometer,
                                                 r4=1 * u.nanometer,
                                                 k=250 * u.kilojoule_per_mole / u.nanometer ** 2,
                                                 atom1=s.index.atom(i, name_i, expected_resname=seq[i][-3:]),
                                                 atom2=s.index.atom(j, name_j, expected_resname=seq[j][-3:]))
            rest_group.append(rest)
    return dists

def setup_system():
    # load the sequence
    sequence = parse.get_sequence_from_AA1(filename='sequence.dat')
    n_res = len(sequence.split())

    # build the system
    p = meld.AmberSubSystemFromSequence(sequence)
    # p = meld.AmberSubSystemFromPdbFile('TEMPLATES/minimized_protein.pdb')
    build_options = meld.AmberOptions(
        forcefield="ff14sbside",
        implicit_solvent_model='gbNeck2',
        use_big_timestep=True,
        cutoff=1.8 * u.nanometers,
        remove_com=False,
        # use_amap = False,
        enable_amap=False,
        amap_beta_bias=1.0,
    )

    builder = meld.AmberSystemBuilder(build_options)
    s = builder.build_system([p]).finalize()
    # s.temperature_scaler = meld.ConstantTemperatureScaler(300.0 * u.kelvin)
    s.temperature_scaler = meld.GeometricTemperatureScaler(0, 0.3, 300. * u.kelvin, 550. * u.kelvin)

    ramp = s.restraints.create_scaler('nonlinear_ramp', start_time=1, end_time=200,
                                      start_weight=1e-3, end_weight=1, factor=4.0)
    seq = sequence.split()
    for i in range(len(seq)):
        if seq[i][-3:] == 'HIE': seq[i] = 'HIS'
    print(seq)
    #
    # Secondary Structure
    #
    ss_scaler = s.restraints.create_scaler('constant')
    ss_rests = parse.get_secondary_structure_restraints(filename='ss.dat', system=s, scaler=ss_scaler,
                                                        ramp=ramp,
                                                        torsion_force_constant=0.01 * u.kilojoule_per_mole / u.degree ** 2,
                                                        distance_force_constant=2.5 * u.kilojoule_per_mole / u.nanometer ** 2,
                                                        quadratic_cut=2.0 * u.nanometer)
    n_ss_keep = int(len(ss_rests) * 0.85)
    s.restraints.add_selectively_active_collection(ss_rests, n_ss_keep)
    print(len(ss_rests))

    conf_scaler = s.restraints.create_scaler('constant')
    confinement_rests = []
    for index in range(n_res):
        rest = s.restraints.create_restraint('confine', conf_scaler, ramp=ramp,
                                             atom_index=s.index.atom(index, 'CA', expected_resname=seq[index][-3:]),
                                             radius=4.5 * u.nanometer,
                                             force_const=250.0 * u.kilojoule_per_mole / u.nanometer ** 2)
        confinement_rests.append(rest)
    s.restraints.add_as_always_active_list(confinement_rests)

    #
    # Setup Scaler
    #
    scaler = s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=1.0, factor=4.0)
    #

    # creates parameter sampling for strand pairing
    dists = get_dist_restraints_evfold('evfold.dat', s, scaler, ramp, seq)
    prior_ev = param_sampling.ScaledExponentialDiscretePrior(u0=1.0, temperature_scaler=s.temperature_scaler, scaler=scaler)
    sampler_ev = param_sampling.DiscreteSampler(int(1), int(1.00 * len(dists)), 1)
    param_ev = s.param_sampler.add_discrete_parameter("param_EV", int(0.8 * len(dists)), prior_ev, sampler_ev)
    s.restraints.add_selectively_active_collection(dists, param_ev)
    #s.restraints.add_selectively_active_collection(dists, int(0.8 * len(dists)))
    print(len(dists))
    # setup mcmc at startup
    movers = []
    n_atoms = s.n_atoms
    for i in range(0, n_res):
        n = s.index.atom(i, 'N', expected_resname=seq[i][-3:])
        ca = s.index.atom(i, 'CA', expected_resname=seq[i][-3:])
        c = s.index.atom(i, 'C', expected_resname=seq[i][-3:])

        atom_indxs = list(system.indexing.AtomIndex(j) for j in range(ca, n_atoms))
        mover = mc.DoubleTorsionMover(index1a=n, index1b=ca,
                                      atom_indices1=list(system.indexing.AtomIndex(i) for i in range(ca, n_atoms)),
                                      index2a=ca, index2b=c,
                                      atom_indices2=list(system.indexing.AtomIndex(j) for j in range(c, n_atoms)))

        movers.append((mover, 1))

    sched = mc.MonteCarloScheduler(movers, n_res * 60)

    # create the options
    options = meld.RunOptions(
        timesteps=14286,
        minimize_steps=20000,
        min_mc=sched,
        param_mcmc_steps=200
    )

    # create a store
    store = vault.DataStore(gen_state(s, 0), N_REPLICAS, s.get_pdb_writer(),
                            block_size=BLOCK_SIZE)  # why i need gen_state(s,0)? doubtful
    store.initialize(mode='w')
    store.save_system(s)
    store.save_run_options(options)

    # create and store the remd_runner
    l = ladder.NearestNeighborLadder(n_trials=48 * 48)
    policy_1 = adaptor.AdaptationPolicy(2.0, 50, 50)
    a = adaptor.EqualAcceptanceAdaptor(n_replicas=N_REPLICAS, adaptation_policy=policy_1, min_acc_prob=0.02)

    remd_runner = remd.leader.LeaderReplicaExchangeRunner(N_REPLICAS, max_steps=N_STEPS, ladder=l, adaptor=a)
    store.save_remd_runner(remd_runner)

    # create and store the communicator
    c = comm.MPICommunicator(s.n_atoms, N_REPLICAS, timeout=60000)
    store.save_communicator(c)

    # create and save the initial states
    states = [gen_state(s, i) for i in range(N_REPLICAS)]
    store.save_states(states, 0)

    # save data_store
    store.save_data_store()

    return s.n_atoms


setup_system()

# conda init bash
# cd 5p21evadapt
# source ~/.bashrc
# conda activate /work/env
# python EVmeld_adapt.py
# nohup mpirun -np 2 launch_remd > meld.log 2>&1 &
# chmod +x job2.sh
# chmod +x job.sh
# nohup ./job2.sh > extract.log 2>&1 &
# nohup python RMSD.py > RMSD.log 2>&1 &
# nohup python cluster2.py > cluster2.log 2>&1 &
# nohup ./job.sh > meld.log 2>&1 &
# conda create --name new_meld --clone param_meld
# CUDA_VISIBLE_DEVICES="2,3" nohup mpirun -np 2 launch_remd > meld.log 2>&1 &
