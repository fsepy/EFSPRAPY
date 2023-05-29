if __name__ == '__main__':
    from efsprapy.mcs0 import MCS0

    mcs = MCS0()
    mcs.set_inputs_file_path(r"C:\Users\IanFu\Desktop\New folder (2)\efsprapy_input_file.xlsx")
    mcs.run(n_proc=6, save=True, save_archive=False)
