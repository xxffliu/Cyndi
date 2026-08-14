#!/usr/bin/env python
"""Round 2: combine the winning directions from round 1."""
import os

BASE = {
    'MOGA_Max_Conformers': 100,
    'MOGA_Num_Objectives': 3,
    'MOGA_Max_Generation': 300,
    'MOGA_Population_Size': 100,
    'MOGA_Max_Run': 1,
    'MOGA_Crossover_Probability': 0.85,
    'MOGA_Mutation_Probability': 0.20,
    'MOGA_SBX': '15 20',
    'MOGA_Epsilon_Quaternion': '3 0.3 0.1 2',
    'MOGA_Keep_Input_Conformer[Y/N]': 'N',
    'MOGA_Optimize_Conformer[Y/N]': 'Y',
    'MOGA_Energy_Cutoff': 60,
    'MOGA_RMSD_Scale_Factor': 0.2,
    'MOGA_Max_Opt_Iteration': 50,
    'MOGA_Split_Output[Y/N]': 'N',
    'MOGA_Random_Seed': 0.42,
}

VARIANTS = {
    't2_comb300_mut02': {},   # 300 gen + mut 0.20 (round-1 winners)
    't2_comb300_mut02_coarse': {'MOGA_Epsilon_Quaternion': '6 0.5 0.2 4'},
    't2_gen200_mut02':   {'MOGA_Max_Generation': 200},
    't2_gen500_mut02':   {'MOGA_Max_Generation': 500},
    't2_gen300_mut015':  {'MOGA_Mutation_Probability': 0.15},
    't2_gen300_mut03':   {'MOGA_Mutation_Probability': 0.30},
    't2_gen300_mut02_pop150': {'MOGA_Population_Size': 150},
    't2_gen300_mut02_pop60':  {'MOGA_Population_Size': 60},
    't2_gen300_mut02_sbx':    {'MOGA_SBX': '20 25'},
    't2_gen300_mut02_eps2':   {'MOGA_Epsilon_Quaternion': '2 0.2 0.1 2'},
}

def fmt(p):
    lines = []
    for k in p:
        v = p[k]
        lines.append('%s\t\t%s' % (k, v))
    return '\n'.join(lines) + '\n'

for name, over in VARIANTS.items():
    p = dict(BASE)
    p.update(over)
    with open(os.path.join('tune_params', name + '.in'), 'w') as fh:
        fh.write(fmt(p))
print('generated', len(VARIANTS), 'round-2 variants')
