simple = """&HEAD VERSION = 7700, TITLE = 'CFAST Simulation' /
&MHDR NUMBER_OF_CASES = 3 /
 
!! Scenario Configuration 
&TIME SIMULATION = {t_end} PRINT = 0 SMOKEVIEW = 0 SPREADSHEET = {t_step} / 
&INIT PRESSURE = 101325 RELATIVE_HUMIDITY = 50 INTERIOR_TEMPERATURE = 20 EXTERIOR_TEMPERATURE = 20 /
 
!! Material Properties 
&MATL ID = 'CONCRETE' MATERIAL = 'Concrete Normal Weight (6 in)', 
      CONDUCTIVITY = 1.75 DENSITY = 2200 SPECIFIC_HEAT = 1, THICKNESS = 0.15 EMISSIVITY = 0.94 /
 
!! Compartments 
&COMP ID = 'ROOM'
      DEPTH = {room_depth} HEIGHT = {room_height} WIDTH = {room_width}
      CEILING_MATL_ID = 'CONCRETE' CEILING_THICKNESS = 0.15 WALL_MATL_ID = 'CONCRETE' WALL_THICKNESS = 0.15
      ORIGIN = 0, 0, 0 GRID = 50, 50, 50 /
 
!! Wall Vents
&VENT TYPE = 'WALL' ID = 'OPENING' COMP_IDS = 'ROOM' 'OUTSIDE' , BOTTOM = 0 HEIGHT = {opening_height}, WIDTH = {opening_width}
      FACE = 'FRONT'  OFFSET = 1 /
 
!! Fires 
&FIRE ID = 'Fire'  COMP_ID = 'ROOM', FIRE_ID = 'Constant Fire'  LOCATION = 1.5, 1.5 / 
&CHEM ID = 'Constant Fire' CARBON = 1 CHLORINE = 0 HYDROGEN = 4 NITROGEN = 0 OXYGEN = 0 HEAT_OF_COMBUSTION = 50000 RADIATIVE_FRACTION = 0.35 / 
{fire_hrr_curve_tabl}
 
!! Devices
&DEVC ID = 'Sprinkler_1' COMP_ID = 'ROOM' LOCATION = {sprinkler_loc_x}, {sprinkler_loc_y}, {sprinkler_loc_z} TYPE = 'SPRINKLER' SETPOINT = {sprinkler_activation_temperature}, RTI = {sprinkler_rti} SPRAY_DENSITY = 7E-05 /
 
&TAIL /
"""
