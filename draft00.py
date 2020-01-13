import os
import numpy as np
import itertools

import niggli
import pyniggli

hf_data = lambda *x: os.path.join('data', *x)
hfe = lambda x,y,eps=1e-5:np.max(np.abs(x-y)/(np.abs(x)+np.abs(y)+eps))

def parse_data(filename):
    '''
    each row as a vector: [vec1; vec2; vec3]
    '''
    with open(filename) as fid:
        tmp0 = [x.strip() for x in fid]
    tmp0 = [x for x in tmp0[6:] if x] #remove empty
    assert len(tmp0)%4==0 and all((x+1)==int(y[2:]) for x,y in enumerate(tmp0[::4]))
    tmp0 = [([float(z) for z in y.split()]) for x,y in enumerate(tmp0) if x%4]
    ret = np.array(tmp0).reshape(-1, 3, 3)
    return ret

def get_angle(crystal_basis):
    a,b,c = np.sqrt(np.sum(crystal_basis**2, axis=1))
    ret_bc = np.arccos(np.dot(crystal_basis[1], crystal_basis[2]) / (b*c)) * (180/np.pi)
    ret_ac = np.arccos(np.dot(crystal_basis[0], crystal_basis[2]) / (a*c)) * (180/np.pi)
    ret_ab = np.arccos(np.dot(crystal_basis[0], crystal_basis[1]) / (a*b)) * (180/np.pi)
    return ret_bc, ret_ac, ret_ab

lattice_list = parse_data(hf_data('lattices.dat'))
reduced_lattice_list = parse_data(hf_data('reduced_lattices.dat'))
hf0 = lambda x: np.all(x>(90 - 1e-3)) or np.all(x<(90 + 1e-3))
assert all(hf0(np.array(get_angle(x))) for x in reduced_lattice_list)

eps = 1e-5
ret_niggli = [niggli.niggli_reduce(x.T, eps=eps).T for x in lattice_list]
ret_pyniggli = [pyniggli.niggli_reduce(x.T, eps=eps).T for x in lattice_list]
ret_pyniggli_full = [pyniggli.niggli_reduce_full(x, eps=eps)[1] @ x for x in lattice_list]

print('hfe(.dat, niggli.niggli_reduce)', max(hfe(x,y) for x,y in zip(reduced_lattice_list, ret_niggli)))
print('hfe(.dat, pyniggli.niggli_reduce)', max(hfe(x,y) for x,y in zip(reduced_lattice_list, ret_pyniggli)))
print('hfe(.dat, pyniggli.niggli_reduce_full)', max(hfe(x,y) for x,y in zip(reduced_lattice_list, ret_pyniggli_full)))
