def test_mcs0_deterministic():
    from efsprapy.mcs0.calcs import main
    from efsprapy.mcs0.inputs import EXAMPLE_INPUT_DETERMINISTIC
    *_, ftp = main(**EXAMPLE_INPUT_DETERMINISTIC)
    assert abs(ftp - 63004611.94969697) <= 1e-5


def test_mcs0():
    from efsprapy.mcs0 import MCS0, EXAMPLE_INPUT
    mcs = MCS0()
    mcs.set_inputs_dict(EXAMPLE_INPUT.copy())
    mcs.run(n_proc=1, save=True, save_archive=False)


def test_mcs0_research():
    from efsprapy.mcs0 import MCS0, EXAMPLE_INPUT

    custom_input = EXAMPLE_INPUT.copy()
    custom_input['CASE_1']['chf'] = 13_400
    custom_input['CASE_1']['ftp_index'] = 2
    custom_input['CASE_1']['ftp_target'] = 34_592
    custom_input['CASE_1']['emissivity'] = dict(dist='uniform_', lbound=0.5, ubound=1.0)
    mcs = MCS0()
    mcs.set_inputs_dict(custom_input)
    mcs.run(n_proc=1, save=True, save_archive=False)


if __name__ == "__main__":
    # test_mcs0_deterministic()
    test_mcs0()
