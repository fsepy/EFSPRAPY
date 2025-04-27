import pathlib
import shutil

from efsprapy.goal_seek import sep_parallel_any_br187
from efsprapy.mcs.parser import InputParser
from efsprapy.mcs.xlsx import dict_to_xlsx, xlsx_to_dict


def prepare_inputs_with_var_sep_dist(
        n_simulations: int, kwargs: dict, receiver_separations: list, fp_xlsx: pathlib.Path
):
    case_name = fp_xlsx.stem
    assert '-' not in case_name, f'name cannot contain -'
    kwargs_ = dict()

    for receiver_separation in receiver_separations:
        kwargs_[f'{case_name}-{receiver_separation:06.3f}'] = kwargs | dict(receiver_separation=receiver_separation)

    kwargs_input = {k_: InputParser.flatten_dict(v_) for k_, v_ in kwargs_.items()}
    dict_to_xlsx(kwargs_input, fp_xlsx.as_posix())

    input_data = xlsx_to_dict(fp_xlsx.as_posix())
    dir_cases = fp_xlsx.parents[0]
    for case_name in input_data.keys():
        if 'n_simulations' in input_data[case_name].keys():
            input_data[case_name].pop('n_simulations')
        input_parser = InputParser(input_data[case_name], n_simulations)
        try:
            shutil.rmtree(dir_cases / case_name)
        except FileNotFoundError:
            pass
        input_parser.to_cases_csv_and_json((dir_cases / case_name).as_posix())


def prepare_inputs_with_var_w_and_h(
        n_simulations: int, kwargs: dict, ws: list, hs: list, fp_xlsx: pathlib.Path
):
    kwargs_ = dict()

    for w in ws:
        for h in hs:
            kwargs_[f'{fp_xlsx.stem}-{w:06.3f}-{h:06.3f}'] = kwargs | dict(
                receiver_separation=sep_parallel_any_br187(w, h, 84, 12.6),
                room_width=w,
                room_height=h,
                opening_width=w,
                opening_height=h,
            )

    kwargs_input = {k_: InputParser.flatten_dict(v_) for k_, v_ in kwargs_.items()}
    dict_to_xlsx(kwargs_input, fp_xlsx.as_posix())

    input_data = xlsx_to_dict(fp_xlsx.as_posix())
    dir_cases = fp_xlsx.parents[0]
    for case_name in input_data.keys():
        if 'n_simulations' in input_data[case_name].keys():
            input_data[case_name].pop('n_simulations')
        input_parser = InputParser(input_data[case_name], n_simulations)
        try:
            shutil.rmtree(dir_cases / case_name)
        except FileNotFoundError:
            pass
        input_parser.to_cases_csv_and_json((dir_cases / case_name).as_posix())
