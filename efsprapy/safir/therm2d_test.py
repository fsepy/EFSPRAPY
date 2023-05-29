import logging
import shutil
import subprocess
import tempfile
from os import path

from efsprapy.safir.test_files import output_therm1d_xml, simple_test_in
from efsprapy.safir.therm2d import PPXML, Run

logger = logging.getLogger('gui')


def test_get_line_temp_isotherm():
    """Thermal2DPPXML"""
    model = PPXML()
    with open(output_therm1d_xml, 'r') as f:
        model.xml = f.read()
    isotherm = model.get_line_temp_isotherm(100, 0, 0, 0, 0.2, 100)
    # print(isotherm[-1])
    # import matplotlib.pyplot as plt
    # plt.plot(model.t, isotherm)
    # plt.show()
    assert (isotherm[-1] - 0.02902189) <= 1e-5


def test_run():
    model = Run()
    with tempfile.TemporaryDirectory() as temp_dir:
        shutil.copy2(simple_test_in, path.join(temp_dir, path.basename(simple_test_in)))
        model.run(path.join(temp_dir, path.basename(simple_test_in)), timeout=1800, print_time=print)


def test_run_timeout():
    model = Run()
    with tempfile.TemporaryDirectory() as temp_dir:
        shutil.copy2(simple_test_in, path.join(temp_dir, path.basename(simple_test_in)))
        try:
            model.run(path.join(temp_dir, path.basename(simple_test_in)), timeout=1, print_time=print)
        except subprocess.TimeoutExpired:
            return
        raise AssertionError


if __name__ == '__main__':
    # test_get_line_temp_isotherm()
    # test_run_timeout()
    test_run()
