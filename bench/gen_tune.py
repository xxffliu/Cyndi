#!/usr/bin/env python
"""Generate tuning variants from the v2 baseline param file."""
import os

BASE = {
    'MOGA_Max_Conformers': 100,
    'MOGA_Num_Objectives': 3,
    'MOGA_Max_Generation': 100,
    'MOGA_Population_Size': 100,
    'MOGA_Max_Run': 1,
    'MOGA_Crossover_Probability': 0.85,
    'MOGA_Mutation_Probability': 0.10,
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
    'tune_base': {},
    'tune_big':      {'MOGA_Population_Size': 200, 'MOGA_Max_Generation': 200},
    'tune_moregen':  {'MOGA_Max_Generation': 300},
    'tune_pop200':   {'MOGA_Population_Size': 200},
    'tune_eps_fine': {'MOGA_Epsilon_Quaternion': '1.5 0.2 0.05 1'},
    'tune_eps_coarse': {'MOGA_Epsilon_Quaternion': '6 0.5 0.2 4'},
    'tune_mut02':    {'MOGA_Mutation_Probability': 0.20},
    'tune_mut005':   {'MOGA_Mutation_Probability': 0.05},
    'tune_cutoff100':{'MOGA_Energy_Cutoff': 100},
    'tune_cutoff30': {'MOGA_Energy_Cutoff': 30},
    'tune_sbx':      {'MOGA_SBX': '20 25'},
    'tune_fast':     {'MOGA_Population_Size': 60, 'MOGA_Max_Generation': 60,
                      'MOGA_Epsilon_Quaternion': '5 0.4 0.15 3'},
}

def fmt(p):
    lines = []
    for k in sorted(p, key=lambda x: list(BASE).index(x) if x in BASE else 99):
        v = p[k]
        if k == 'MOGA_SBX':
            lines.append('%s\t\t\t%s' % (k, v))
        elif k == 'MOGA_Epsilon_Quaternion':
            lines.append('%s\t\t%s' % (k, v))
        else:
            lines.append('%s\t\t%s' % (k, v))
    return '\n'.join(lines) + '\n'

for name, over in VARIANTS.items():
    p = dict(BASE)
    p.update(over)
    with open(os.path.join('tune_params', name + '.in'), 'w') as fh:
        fh.write(fmt(p))
print('generated', len(VARIANTS), 'variants in tune_params/')
