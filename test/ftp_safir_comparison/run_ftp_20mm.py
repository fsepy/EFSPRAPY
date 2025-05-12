import pathlib
from os import chdir, getcwd

from efsprapy.mcs.parser import InputParser
from efsprapy.mcs.xlsx import dict_to_xlsx, xlsx_to_dict
from efsprapy.mcs1 import EXAMPLE_INPUT
from efsprapy.prepare_inputs import convert_input_file_xlsx_to_cases
from efsprapy.run_analysis import process_multiple_cases_2

if __name__ == '__main__':
    chdir(r'C:\Users\IanFu\Desktop\~1_CURRENT\efsprapy\01-analysis\ftp_safir_comparison')
    name = 'ftp_20mm'

    path_work = pathlib.Path(getcwd()) / name
    path_work.mkdir(parents=True, exist_ok=True)
    fp_xlsx = path_work / f'{name}.xlsx'

    ftp_20mm_data = dict(
        macrocarpa=dict(q_min=17, q_cr=15.4, ftp=5791, n=1.6, t_ig=394),
        beech=dict(q_min=16, q_cr=13.3, ftp=9094, n=1.7, t_ig=367),
        mdf=dict(q_min=15, q_cr=7.7, ftp=20423, n=1.7, t_ig=273),
        radiata_pine=dict(q_min=15, q_cr=11.5, ftp=8317, n=1.6, t_ig=340),
        rimu=dict(q_min=14, q_cr=8.2, ftp=45595, n=2.0, t_ig=283),
        plywood=dict(q_min=13, q_cr=8.1, ftp=27314, n=1.9, t_ig=281)
    )

    input_data_dict = dict()
    for k, v in ftp_20mm_data.items():
        input_data_dict[k] = EXAMPLE_INPUT | dict(
            ftp_chf=v['q_cr'] * 1e3, ftp_index=v['n'], ftp_target=v['ftp'],

            receiver_separation=10,
            fire_combustion_efficiency=1, receiver_ignition_temperature=-1, safir_input_file_s=None,
            fire_fuel_density_MJm2=dict(dist="gumbel_r_", lbound=10, ubound=3000, mean=780, sd=234),
            fire_hrr_density_kWm2=dict(dist="uniform_", lbound=0.15e3, ubound=0.65e3),
            fire_mode=dict(dist='constant_', lbound=0, ubound=0, values=None, weights=None),
        )

    input_data_dict_flatten = {k_: InputParser.flatten_dict(v_) for k_, v_ in input_data_dict.items()}
    dict_to_xlsx(input_data_dict_flatten, fp_xlsx.as_posix())
    convert_input_file_xlsx_to_cases(fp_xlsx)

    input_data_dict_from_xlsx = xlsx_to_dict(fp_xlsx.as_posix())
    folders = [item for item in path_work.iterdir() if item.is_dir()]
    folders.sort()
    process_multiple_cases_2(folders)
