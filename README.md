# pyniggli: niggli algorithm to reduce crystal lattice

rewrite in pure python and numpy

a drop-in replacement for the original implementation

the c-code and the formula in the paper (see link) are quite confusing for me. If anyone get a deeper insight, please help to improve the code.

## quickstart

1. download the whole repository, `git clone xxx` or download as zip file
2. at the root directory, you should see `setup.py`, `test_pyniggli.py` etc
3. install: at the root directory, run `pip install .`
4. uninstall: at the root directory, run `pip uninstall .`
5. unittest: at the root directory, run `pytest`
6. usage: `from pyniggli import niggli_reduce, niggli_reduce_full`
   * `niggli_reduce()` is mainly for backward compatilibity (used in `phonopy`)
   * `niggli_reduce_full()` recommended to use this one, see `test_pyniggli.py` for detail usage

## about license

The niggli code is released under the license (see link) (BSD-3?), but I  prefer MIT license.

The origin License file is also stored in this repo (see link) for reference, add I also add MIT license file. I am not sure whether it's legal or not, please tell me if not, thanks.

also claim: the data files in `data/` are from niggli repo.
