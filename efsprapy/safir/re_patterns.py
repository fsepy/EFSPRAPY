import re

__all__ = 'in_step_and_uptime_2', 'in_step_and_uptime_3', 'in_prt_step_and_uptime', 'out_convergence_time'

in_step_and_uptime_2 = re.compile(r'\n\s+TIME\s+([0-9.]+)\s+([0-9.]+)\s+END_?TIME')
in_step_and_uptime_3 = re.compile(r'\n\s+TIME\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+END_?TIME')
in_prt_step_and_uptime = re.compile(r'\n\s+TIME_?PRINT\s+([0-9.]+)\s+([0-9.]+)\s+END_?TIMEPR')
out_convergence_time = re.compile(r'\n\s+CONVERGENCE HAS BEEN OBTAINED[\s.=]*TIME[\s.=]*([0-9.]+)\n')
stdout_last_time = re.compile(r'time\s+=\s+([0-9.]+)\s+sec.')
