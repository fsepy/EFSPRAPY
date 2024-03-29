from os import path

from tqdm import tqdm

from efsprapy.mcs1 import MCS1

if __name__ == '__main__':
    fp_input = path.join(path.dirname(__file__), "in.xlsx")

    mcs = MCS1()
    mcs.set_inputs_file_path(fp_input)
    pbar = tqdm()
    mcs.run(
        10,
        lambda _: pbar.update(1),
        lambda _: setattr(pbar, "total", _),
        save=True,
        save_archive=False,
        concurrency_strategy=1,
    )
    pbar.close()
