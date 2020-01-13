import os
import numpy as np
import itertools

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


def test_reduced_data():
    reduced_lattice_list = parse_data(hf_data('reduced_lattices.dat'))
    hf0 = lambda x: np.all(x>(90 - 1e-3)) or np.all(x<(90 + 1e-3))
    assert all(hf0(np.array(get_angle(x))) for x in reduced_lattice_list)


def test_pyniggli_reduce():
    lattice_list = parse_data(hf_data('lattices.dat'))
    reduced_lattice_list = parse_data(hf_data('reduced_lattices.dat'))
    eps = 1e-5
    ret_pyniggli = [pyniggli.niggli_reduce(x.T, eps=eps).T for x in lattice_list]
    assert max(hfe(x,y) for x,y in zip(reduced_lattice_list, ret_pyniggli)) < 1e-7

def test_pyniggli_reduce_full():
    lattice_list = parse_data(hf_data('lattices.dat'))
    reduced_lattice_list = parse_data(hf_data('reduced_lattices.dat'))
    eps = 1e-5
    tmp0 = [pyniggli.niggli_reduce_full(x, eps=eps) for x in lattice_list]
    tag_success = [x[0] for x in tmp0]
    assert all(tag_success)
    ret_pyniggli_full = [x[1]@y for x,y in zip(tmp0,lattice_list)] #the second element is tmat
    assert max(hfe(x,y) for x,y in zip(reduced_lattice_list, ret_pyniggli_full)) < 1e-7


def test_compare_with_niggli():
    try:
        import niggli
    except ImportError:
        print('niggli is required to run this unittest, see https://github.com/atztogo/niggli')
    else:
        lattice_list = parse_data(hf_data('lattices.dat'))
        eps = 1e-5
        ret_niggli = [niggli.niggli_reduce(x.T, eps=eps).T for x in lattice_list]
        ret_pyniggli = [pyniggli.niggli_reduce(x.T, eps=eps).T for x in lattice_list]
        assert max(hfe(x,y) for x,y in zip(ret_niggli, ret_pyniggli)) < 1e-7
