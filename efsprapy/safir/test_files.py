from os import path

dir_this = path.dirname(__file__)

clt_mlr_solver_fp_in = path.join(dir_this, 'tga_solver', 'input.in')
clt_mlr_solver_fp_flux = path.join(dir_this, 'tga_solver', 'flux.txt')
clt_mlr_solver_fp_tga = path.join(dir_this, 'tga_solver', 'tga.csv')

output_therm1d_xml = path.join(dir_this, 'xml', 'output.XML')

max_load_solver_in = path.join(dir_this, 'max_load_solver', 'struct3d.in')
max_load_solver_tem = path.join(dir_this, 'max_load_solver', 'therm2d.tem')
max_load_solver_load = path.join(dir_this, 'max_load_solver', 'load.fct')
max_load_solver_out = path.join(dir_this, 'max_load_solver', 'struct3d.OUT')

simple_test_in = path.join(dir_this, 'therm2d_simple', 'input.in')
simple_test_flux_txt = path.join(dir_this, 'therm2d_simple', 'flux.txt')

therm1d_hf_50 = path.join(dir_this, 'therm1d', 'hf_50.in')
