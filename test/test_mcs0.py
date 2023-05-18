def test_mcs0_deterministic():
    from efsprapy.mcs0.calcs import ftp_main
    from efsprapy.mcs0.inputs import EXAMPLE_INPUT_DETERMINISTIC
    *_, ftp = ftp_main(**EXAMPLE_INPUT_DETERMINISTIC)
    assert abs(ftp - 63004611.94969697) <= 1e-5


def test_mcs0():
    from efsprapy.mcs0 import MCS0, EXAMPLE_INPUT
    mcs = MCS0()
    mcs.set_inputs_dict(EXAMPLE_INPUT.copy())
    mcs.run(n_proc=1, save=True, save_archive=False)


if __name__ == "__main__":
    test_mcs0_deterministic()
    test_mcs0()
