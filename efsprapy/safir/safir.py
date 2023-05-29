import re
from os import path

import numpy as np


def out_to_pout(fp: str) -> dict:
    """Extract strain data from Safir *.out or processed output file containing strain data only and store in a dict.
    The resulting dict data structure:
    {
        list_time: [...],
        list_shell: [...],
        list_surf: [...],
        list_rebar: [...],
        list_strain: [...],
        list_strain2: [...],
    }
    All elements in the dict have the same length.
    """
    rp_time_str = re.compile(r'TIME[ ]*=[ ]+[0-9.0-9]+')
    rp_time_val = re.compile(r'[0-9.0-9]+')
    rp_shell_str = re.compile(r'SHELL:[ ]*[0-9]+')
    rp_shell_val = re.compile(r'[0-9]+')
    rp_surf_str = re.compile(r'SURF:[ ]*[0-9]+')
    rp_surf_val = re.compile(r'[0-9]+')
    rp_rebar_str = re.compile(r'REBAR:[ ]*[0-9]+')
    rp_rebar_val = re.compile(r'[0-9]+')
    rp_strain_str = re.compile(r'Total strain[ ]*:[ ]*[-0-9.]+|Strain[ ]*:[ ]*[-0-9.]+')
    rp_strain_val = re.compile(r'[-0-9.]+')
    rp_strain_str2 = re.compile(r'Stress related strain[ ]*:[ -]*[0-9.]+')
    rp_strain_val2 = re.compile(r'[-0-9.]+')

    def get_value(s, rp1, rp2):
        s1 = rp1.findall(s)
        if s1:
            return float(rp2.findall(s1[0])[0])
        else:
            return None

    time_current = 0
    list_time, list_shell, list_surf, list_rebar, list_strain, list_strain2 = [], [], [], [], [], []
    with open(fp, 'rb') as f:
        while True:
            try:
                l = f.readline().decode('utf-8')
            except UnicodeDecodeError:
                continue

            if not l:
                break
            else:
                time = get_value(l, rp_time_str, rp_time_val)
                if time:
                    time_current = time
                elif time_current and 'SHELL' in l:
                    list_time.append(time_current)
                    list_shell.append(get_value(l, rp_shell_str, rp_shell_val))
                    list_surf.append(get_value(l, rp_surf_str, rp_surf_val))
                    list_rebar.append(get_value(l, rp_rebar_str, rp_rebar_val))
                    list_strain.append(get_value(l, rp_strain_str, rp_strain_val))
                    list_strain2.append(get_value(l, rp_strain_str2, rp_strain_val2))

    list_time = np.array(list_time, float)
    list_shell = np.array(list_shell, float)
    list_surf = np.array(list_surf, float)
    list_rebar = np.array(list_rebar, float)
    list_strain = np.array(list_strain, float)
    list_strain2 = np.array(list_strain2, float)

    return dict(
        list_time=list_time, list_shell=list_shell, list_surf=list_surf,
        list_rebar=list_rebar, list_strain=list_strain, list_strain2=list_strain2,
    )


def pout_to_csv(
        fp: str,
        list_time,
        list_shell,
        list_surf,
        list_rebar,
        list_strain,
        list_strain2
):
    data = zip(list_time, list_shell, list_surf, list_rebar, list_strain, list_strain2)
    data_list = [[j for j in i] for i in data]
    data_arr = np.array(data_list, dtype=float)
    np.savetxt(fp, data_arr, delimiter=",",
               header='TIME,SHELL,SURF,REBAR,STRAIN,STRESS_STRAIN',
               fmt=['%10g', '%10g', '%10g', '%10g', '%10.7g', '%10.7g'])


def pout_plot_strain(
        unique_shell: int,
        list_time,
        list_shell,
        list_surf,
        list_rebar,
        list_strain,
        list_strain2
):
    """"""
    list_unique_surf = list(set(list_surf[list_shell == unique_shell]))
    list_unique_surf.sort()

    list_lines = []
    for unique_surf in list_unique_surf:
        list_unique_rebar = list(set(list_rebar[list_surf == unique_surf]))
        list_unique_rebar.sort()
        for unique_rebar in list_unique_rebar:
            time_ = list_time[
                (list_shell == unique_shell) & (list_surf == unique_surf) & (list_rebar == unique_rebar)
                ]
            strain_ = list_strain[
                (list_shell == unique_shell) & (list_surf == unique_surf) & (list_rebar == unique_rebar)
                ]
            label_ = f'surf {unique_surf:g} rebar {unique_rebar:g}'

            if len(strain_) > 0:
                list_lines.append(
                    dict(
                        x=time_,
                        y=strain_,
                        label=label_
                    )
                )
    return list_lines


def t0r_to_tem(fp_tem: str, fp_t0r: str = None):
    """Insert SAFIR torsion analysis output *.t0r into SAFIR thermal analysis output *.tem. The specified *.tem file
    will be overriden upon successfully completing the process.

    No official/unofficial documentation was found regarding exactly where to insert the torsion data into the tem file
    and this function is created based on "reverse-engineering" of sample files.

    :param fp_tem: file path of the *.tem file
    :param fp_t0r: file path of the *.T0R file
    :return:
    """
    # check if *.tem exist
    fp_tem = path.realpath(fp_tem)
    if not path.exists(fp_tem):
        raise FileNotFoundError(f'Unable to find tem file {fp_tem}')

    # figure not file name
    try:
        name = path.splitext(path.basename(fp_tem))[0]
    except IndexError:
        # in case of `fp_tem` without a file extension, use the file name
        name = path.basename(fp_tem)

    # check if *.t0r exist
    if fp_t0r is None:
        fp_t0r = path.join(path.dirname(fp_tem), f'{name.upper()}-1.T0R')

    if not path.exists(fp_t0r):
        raise FileNotFoundError(f'Unable to find t0r file {fp_t0r}')

    # add tor to tem
    # read t0r data
    with open(path.join(path.dirname(fp_tem), f'{name}-1.t0r'), 'r') as f:
        t0r_s = f.read()

    # obtain t0r data, tight
    gj = re.findall(r' *w[\s\-+0-9.eE]+According to the principle of virtual works,[\s]+GJ=\s*[+0-9.Ee]+', t0r_s)
    if not gj:
        raise ValueError('Unable to find torsion data in ' + path.join(path.dirname(fp_tem), f'{name}-1.t0r'))

    # read tem data
    with open(fp_tem, 'r') as f:
        tem_s = f.read()

    # check if GJ already existed in the tem data
    if re.findall('GJ', tem_s):
        raise ValueError('GJ already existed in tem file ' + fp_tem)

    # locate HOT in tem
    hot_ = re.findall(r' *HOT', tem_s)
    if not hot_:
        raise ValueError(f'Unable to find HOT in ' + fp_tem)

    # in tem data, insert t0r data
    tem_src_ = re.sub(r' *HOT', gj[0] + f'\n{hot_[0]}', tem_s)

    # save tem
    with open(fp_tem, 'r+') as f:
        f.seek(0)
        f.write(tem_src_)
        f.truncate()
