import re

from .re_patterns import *
from .test_files import *


def test_in_step_and_uptime():
    with open(clt_mlr_solver_fp_in, 'r') as f:
        a = f.read()
    b = re.findall(in_step_and_uptime_2, a)
    assert len(b) == 1
    assert len(b[0]) == 2
    assert b[0][0] == '1'
    assert b[0][1] == '1800'


def test_in_prt_step_and_uptime():
    with open(clt_mlr_solver_fp_in, 'r') as f:
        a = f.read()
    b = re.findall(in_prt_step_and_uptime, a)
    assert len(b) == 1
    assert len(b[0]) == 2
    assert b[0][0] == '10'
    assert b[0][1] == '1800'


def test_out_convergence_time():
    with open(max_load_solver_out, 'r') as f:
        a = re.findall(out_convergence_time, f.read())
    assert len(a) > 0
    assert (float(a[-1]) - 4920) < 1e-1
