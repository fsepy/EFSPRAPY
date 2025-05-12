import multiprocessing as mp
import pathlib
import shutil
from functools import partial

from tqdm import tqdm

from efsprapy.goal_seek import sep_parallel_any_br187
from efsprapy.mcs.parser import InputParser
from efsprapy.mcs.xlsx import dict_to_xlsx, xlsx_to_dict


def __convert_input_file_xlsx_to_cases_worker(case_name, input_data, dir_cases):
    """Process a single case"""
    case_input = input_data[case_name].copy()
    n_simulations = case_input.pop('n_simulations')
    input_parser = InputParser(case_input, n_simulations)

    try:
        shutil.rmtree(dir_cases / case_name)
    except FileNotFoundError:
        pass

    InputParser(case_input, n_simulations).to_cases_csv_and_json((dir_cases / case_name).as_posix())


def convert_input_project_data_to_cases(dir_save: pathlib.Path, input_data: dict, n_proc=None):
    if n_proc is None:
        n_proc = mp.cpu_count()

    # Create a partial function with fixed parameters
    process_func = partial(__convert_input_file_xlsx_to_cases_worker, input_data=input_data, dir_cases=dir_save)

    # Create a process pool and map the function to all case names
    with mp.Pool(processes=n_proc) as pool:
        # Use imap_unordered for slightly better performance
        # and wrap with tqdm for progress tracking
        list(tqdm(pool.imap_unordered(process_func, input_data.keys()), total=len(input_data), ))


def convert_input_file_xlsx_to_cases(fp_xlsx: pathlib.Path, n_processes=None):
    """
    Convert input Excel file to case files using multiprocessing

    Args:
        fp_xlsx: Path to the Excel file
        n_processes: Number of processes to use (defaults to CPU count)
    """
    # Default to number of CPUs if not specified
    if n_processes is None:
        n_processes = mp.cpu_count()

    input_data = xlsx_to_dict(fp_xlsx.as_posix())
    dir_cases = fp_xlsx.parents[0]

    convert_input_project_data_to_cases(dir_cases, input_data, n_processes)


# def convert_input_file_xlsx_to_cases(fp_xlsx: pathlib.Path):
#     input_data = xlsx_to_dict(fp_xlsx.as_posix())
#     dir_cases = fp_xlsx.parents[0]
#     for case_name in tqdm(input_data.keys()):
#         n_simulations = input_data[case_name].pop('n_simulations')
#         input_parser = InputParser(input_data[case_name], n_simulations)
#         try:
#             shutil.rmtree(dir_cases / case_name)
#         except FileNotFoundError:
#             pass
#         input_parser.to_cases_csv_and_json((dir_cases / case_name).as_posix())


def prepare_inputs_with_var_sep_dist(
        n_simulations: int, kwargs: dict, receiver_separations: list, fp_xlsx: pathlib.Path
):
    case_name = fp_xlsx.stem
    assert '-' not in case_name, f'name cannot contain -'
    kwargs_ = dict()

    for receiver_separation in receiver_separations:
        kwargs_[f'{case_name}-{receiver_separation:06.3f}'] = kwargs | dict(
            n_simulations=n_simulations,
            receiver_separation=receiver_separation
        )

    kwargs_input = {k_: InputParser.flatten_dict(v_) for k_, v_ in kwargs_.items()}
    dict_to_xlsx(kwargs_input, fp_xlsx.as_posix())
    convert_input_file_xlsx_to_cases(fp_xlsx)


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
