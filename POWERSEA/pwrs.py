###################################################
#    POWERSEA automation script developed by      #
#         Hardik Patil, Ph.D. Candidate           #
#               hardikyp@umich.edu                #
#   Civil Engineering, University of Michigan     #       
###################################################

import re
import time
import numpy as np
from pywinauto.application import Application

###################################################
######      READ GEOMETRY DATA FROM FILE      #####
###################################################
print('Reading geometry data from file...')
geom_prop_path = "D:\\RESEARCH\\Navy_Boat_Project\\Curved-crease-origami-hull\\CAD_Files\\Geom_Prop.mpr"
with open(geom_prop_path) as geom_file:
        line = geom_file.readline()
        while line:
            if re.search('Centroid', line) is not None:
                num = re.search('\d*[.]\d*', line)
                LCG = '%.5f' % (float(num.group())/1000)
                line = geom_file.readline()
                line = geom_file.readline()
                num = re.search('\d*[.]\d*', line)
                VCG = '%.5f' % (float(num.group())/1000)
            elif re.search('Radii of gyration', line) is not None:
                line = geom_file.readline()
                num = re.search('\d*[.]\d*', line)
                radius_of_gyration = '%.5f' % (float(num.group())/1000)
            elif re.search('aftPerp', line) is not None:
                num = re.search('\d*[.]\d*', line)
                aft_perp = '%.5f' % (float(num.group())/1000)
            elif re.search('hullLength', line) is not None:
                num = re.search('\d*[.]\d*', line)
                hull_length = '%.5f' % (float(num.group())/1000)
            line = geom_file.readline()

print('Setting up POWERSEA simulation for calm water...')
input_speed = '5.05' # m/s

###################################################
#####           APPLICATION STARTUP           #####
###################################################
app = Application(backend="uia").start("C:\\Program Files (x86)\\Ship Motion Associates\\PWRS\\pwrs.exe")
time.sleep(3)
pwrs = app.window(title="POWERSEA Simulator - POWERSEA1")

###################################################
#####       IMPORT GEOMETRY IN POWERSEA       #####
###################################################
file_menu = pwrs.child_window(title="Application", auto_id="MenuBar", control_type="MenuBar").child_window(title="File", control_type="MenuItem")
import_menu = file_menu.child_window(title="File", control_type="Window").child_window(title="File", control_type="Menu").child_window(title="Import", control_type="MenuItem")
geometry = import_menu.child_window(title="Import", control_type="Window").child_window(title="Import", control_type="Menu").child_window(title="Geometry...", auto_id="32889", control_type="MenuItem")

# File -> Import -> Geometry...
file_menu.click_input()
import_menu.click_input()
geometry.click_input()
# Import design data
pwrs.child_window(title="Import Design Data", control_type="Window").child_window(title="Append to Existing Design", auto_id="1062", control_type="RadioButton").click_input()
pwrs.child_window(title="Import Design Data", control_type="Window").child_window(title="OK", auto_id="1", control_type="Button").click_input()
# Import Geometry Data 
import_geometry_data = pwrs.child_window(title="Import Geometry Data", control_type="Window")
file_name = import_geometry_data.child_window(title="File name:", auto_id="1152", control_type="Edit")
open_geometry = import_geometry_data.child_window(title="Open", auto_id="1", control_type="Button")
# File name (path) -> File type (select igs from dropdown) -> open geometry
file_name.click_input()
geom_path = '"D:\\RESEARCH\\Navy_Boat_Project\\Curved-crease-origami-hull\\CAD_Files\\Planing_Hull_CCO.igs"'
import_geometry_data.type_keys(geom_path, with_spaces=True)
open_geometry.click_input()

###################################################
#####          SET BOAT PARAMETERS            #####
###################################################
boat_menu = pwrs.child_window(title="Application", auto_id="MenuBar", control_type="MenuBar").child_window(title="Boat", control_type="MenuItem")
units = boat_menu.child_window(title="Boat", control_type="Window").child_window(title="Boat", control_type="Menu").child_window(title="Units...", auto_id="32799", control_type="MenuItem")
water_prop = boat_menu.child_window(title="Boat", control_type="Window").child_window(title="Boat", control_type="Menu").child_window(title="Water Properties", auto_id="32797", control_type="MenuItem")
curves = boat_menu.child_window(title="Boat", control_type="Window").child_window(title="Boat", control_type="Menu").child_window(title="Choose Curves...", auto_id="32888", control_type="MenuItem")
vess_params = boat_menu.child_window(title="Boat", control_type="Window").child_window(title="Boat", control_type="Menu").child_window(title="Vessel Params...", auto_id="32774", control_type="MenuItem")
vess_coefficients = boat_menu.child_window(title="Boat", control_type="Window").child_window(title="Boat", control_type="Menu").child_window(title="Vessel Coefficients...", auto_id="32977", control_type="MenuItem")
propulsion = boat_menu.child_window(title="Boat", control_type="Window").child_window(title="Boat", control_type="Menu").child_window(title="Propulsion...", auto_id="32919", control_type="MenuItem")
simulation_properties = boat_menu.child_window(title="Boat", control_type="Window").child_window(title="Boat", control_type="Menu").child_window(title="Simulation Properties...", auto_id="32956", control_type="MenuItem")

# Boat -> Units...
boat_menu.click_input()
units.click_input()

# SET UNITS
pwrs.child_window(title="Units", control_type="Window").child_window(title="Metric", auto_id="1009", control_type="RadioButton").click_input()
pwrs.child_window(title="Units", control_type="Window").child_window(title="Convert Design Data", auto_id="1177", control_type="CheckBox").click_input()
pwrs.child_window(title="Units", control_type="Window").child_window(title="OK", auto_id="1", control_type="Button").click_input()

# Boat -> Water properties -> Select default
boat_menu.click_input()
water_prop.click_input()
pwrs.child_window(title="Water Properties", control_type="Window").child_window(title="OK", auto_id="1", control_type="Button").click_input()

# CHOOSE CURVES
boat_menu.click_input()
curves.click_input()
choose_curves_window = pwrs.child_window(title="Choose Curves", control_type="Window")

# SET KEEL AND CHINE CURVES
choose_curves_window.type_keys('{TAB 2}{DOWN}{TAB}{DOWN 2}{ENTER}')

# VESSEL PARAMETERS
boat_menu.click_input()
vess_params.click_input()
vessel_characteristics_window = pwrs.child_window(title="Vessel Characteristics", control_type="Window")
vessel_characteristics_window.type_keys('{TAB 3}')
vessel_characteristics_window.type_keys(LCG)
vessel_characteristics_window.type_keys('{TAB}')
vessel_characteristics_window.child_window(title="pwrs", control_type="Window").child_window(title="OK", auto_id="2", control_type="Button").click_input()
vessel_characteristics_window.type_keys(VCG)
vessel_characteristics_window.type_keys('{TAB}')
vessel_characteristics_window.type_keys(radius_of_gyration)
vessel_characteristics_window.type_keys('{TAB}')
vessel_characteristics_window.type_keys('89.3')
vessel_characteristics_window.type_keys('{TAB 3}')
vessel_characteristics_window.type_keys('0.01')
vessel_characteristics_window.type_keys('{TAB}')
vessel_characteristics_window.type_keys(aft_perp)
vessel_characteristics_window.child_window(title="OK", auto_id="1", control_type="Button").click_input()

# VESSEL COEFFICIENTS
boat_menu.click_input()
vess_coefficients.click_input()
vessel_coeff_window = pwrs.child_window(title="Vessel Coefficients", control_type="Window")
vessel_coeff_window.type_keys('{TAB 3}')
vessel_coeff_window.type_keys(input_speed)
vessel_coeff_window.type_keys('{TAB}')
vessel_coeff_window.child_window(title="pwrs", control_type="Window").child_window(title="OK", auto_id="2", control_type="Button").click_input()
vessel_coeff_window.child_window(title="Do Not Use Residual Forces", auto_id="1329", control_type="RadioButton").click_input()
vessel_coeff_window.child_window(title="Automatic", auto_id="1452", control_type="RadioButton").click_input()
vessel_coeff_window.child_window(title="OK", auto_id="1", control_type="Button").click_input()

# PROPULSION PROPERTIES
boat_menu.click_input()
propulsion.click_input()

propulsion_window = pwrs.child_window(title="Propulsion", control_type="Window")
propulsion_window.type_keys('{TAB 4}')
propulsion_window.type_keys(hull_length)
propulsion_window.type_keys('{TAB}0')
propulsion_window.child_window(title="OK", auto_id="1", control_type="Button").click_input()

# SIMULATION PROPERTIES
boat_menu.click_input()
simulation_properties.click_input()
simulation_properties_window = pwrs.child_window(title="Simulation Properties", control_type="Window")
simulation_properties_window.child_window(title="OK", auto_id="1", control_type="Button").click_input()

# CHECK HYDRO SECTIONS
pwrs.type_keys('{F2}')

###################################################
#### SET INITIAL CONDITIONS FROM HYDROSTATICS #####
###################################################
analysis_menu = pwrs.child_window(title="Application", auto_id="MenuBar", control_type="MenuBar").child_window(title="Analysis", control_type="MenuItem")
run_hydrostatics = analysis_menu.child_window(title="Analysis", control_type="Window").child_window(title="Analysis", control_type="Menu").child_window(title="Simple Hydrostatics", auto_id="32924", control_type="MenuItem")

analysis_menu.click_input()
run_hydrostatics.click_input()
pwrs.child_window(title="pwrs", control_type="Window").child_window(title="Yes", auto_id="6", control_type="Button").click_input()

###################################################
#####        SET ANALYSIS CONDITIONS          #####
###################################################
conditions_menu = pwrs.child_window(title="Application", auto_id="MenuBar", control_type="MenuBar").child_window(title="Conditions", control_type="MenuItem")
initial_conditions = conditions_menu.child_window(title="Conditions", control_type="Window").child_window(title="Conditions", control_type="Menu").child_window(title="Initial Conditions...", auto_id="32776", control_type="MenuItem")
incident_waves = conditions_menu.child_window(title="Conditions", control_type="Window").child_window(title="Conditions", control_type="Menu").child_window(title="Incident Waves...", auto_id="32785", control_type="MenuItem")
propulsion_conditions = conditions_menu.child_window(title="Conditions", control_type="Window").child_window(title="Conditions", control_type="Menu").child_window(title="Propulsion Conditions...", auto_id="32920", control_type="MenuItem")

# INITIAL CONDITIONS
conditions_menu.click_input()
initial_conditions.click_input()
initial_conditions_window = pwrs.child_window(title="Initial Conditions", control_type="Window")
initial_conditions_window.type_keys('{TAB 5}')
initial_conditions_window.type_keys(input_speed)
initial_conditions_window.type_keys('{TAB}')
initial_conditions_window.type_keys('0')
initial_conditions_window.child_window(title="OK", auto_id="1", control_type="Button").click_input()

# INCIDENT WAVES
conditions_menu.click_input()
incident_waves.click_input()
incident_waves_window = pwrs.child_window(title="Incident Waves", control_type="Window")
incident_waves_window.child_window(title="OK", auto_id="1", control_type="Button").click_input()

# PROPULSION CONDITIONS
conditions_menu.click_input()
propulsion_conditions.click_input()
propulsion_conditions_window = pwrs.child_window(title="Thrust/Propulsion Conditions", control_type="Window")
propulsion_conditions_window.child_window(title="Constant:", auto_id="1228", control_type="RadioButton").click_input()
propulsion_conditions_window.type_keys('{TAB}')
propulsion_conditions_window.type_keys(input_speed)
propulsion_conditions_window.child_window(title="OK", auto_id="1", control_type="Button").click_input()

###################################################
#####        SET SIMULATION CONTROLS          #####
###################################################
# OPEN RUN CONTROL
pwrs.type_keys('{F4}')
run_control_window = pwrs.child_window(title="Run Control", control_type="Window")
run_control_window.type_keys('{TAB 3}0.001') # set time step
run_control_window.type_keys('{TAB}30{TAB 3}') # set stop time
run_control_window.type_keys(LCG)
run_control_window.type_keys('{TAB}')
run_control_window.type_keys(VCG)
run_control_window.child_window(title="OK", auto_id="1", control_type="Button").click_input()

###################################################
#####             RUN SIMULATION              #####
###################################################
start_time = time.time()
pwrs.type_keys('{F5}')
status_message = pwrs.child_window(auto_id="59393", control_type="StatusBar").Static2.window_text()
status = status_message == 'POWERSEA: Ready'
if not status:
    print('Simulating calm water motion for planing speed =', input_speed, 'm/s')
while (not status):
    status_message = pwrs.child_window(auto_id="59393", control_type="StatusBar").Static2.window_text()
    status = status_message == 'POWERSEA: Ready'
end_time = time.time()
print('Elapsed time:', (end_time - start_time), 'sec')
print('Calm water simulation complete!')

###################################################
#####             EXPORT SIM DATA             #####
###################################################
results_menu = pwrs.child_window(title="Application", auto_id="MenuBar", control_type="MenuBar").child_window(title="Results", control_type="MenuItem")
save_data = results_menu.child_window(title="Results", control_type="Window").child_window(title="Results", control_type="Menu").child_window(title="Save Data...", auto_id="32796", control_type="MenuItem")
results_menu.click_input()
save_data.click_input()
save_results_window = pwrs.child_window(title="Save Simulation Results", control_type="Window")

# CLEAR DATA AND SELECT NEW PARAMTERS TO SAVE
save_results_window.child_window(title="Remove All", auto_id="1119", control_type="Button").click_input()
save_results_window.type_keys('{TAB 7}{SPACE}{DOWN 3}{SPACE}{DOWN}{SPACE}{DOWN}{SPACE}{DOWN 10}{SPACE}{DOWN}{SPACE}{DOWN 7}{SPACE}')
save_results_window.child_window(title="Add", auto_id="1114", control_type="Button").click_input()

# SET FORMAT OF DATA
save_results_window.child_window(title="Format", auto_id="1201", control_type="Button").click_input()
save_results_window.child_window(title="Report Format", control_type="Window").child_window(title="Comma", auto_id="1112", control_type="RadioButton").click_input()
save_results_window.child_window(title="Report Format", control_type="Window").child_window(title="OK", auto_id="1", control_type="Button").click_input()
save_results_window.child_window(title="OK", auto_id="1", control_type="Button").click_input()

# SET PATH AND FILE NAME
save_output_window = pwrs.child_window(title="Save 'Output' Data", control_type="Window")
output_path = '"D:\\RESEARCH\\Navy_Boat_Project\\Curved-crease-origami-hull\\POWERSEA\\Data\\Calm_speed_' + input_speed + '.txt"'
save_output_window.type_keys(output_path, with_spaces=True)
save_output_window.child_window(title="Save", auto_id="1", control_type="Button").click_input()

##################################################
####       RUN LOOP FOR VARYING WAVES        #####
##################################################

print('Setting up wavy simulations...')
input_wavelength = str(2 * float(hull_length))

# CHANGE SPEED IN VESSEL COEFFICIENTS
boat_menu.click_input()
vess_coefficients.click_input()
vessel_coeff_window = pwrs.child_window(title="Vessel Coefficients", control_type="Window")
vessel_coeff_window.type_keys('{TAB 3}')
vessel_coeff_window.type_keys(input_speed)
vessel_coeff_window.type_keys('{TAB}')
vessel_coeff_window.child_window(title="Calculate Residual Forces", auto_id="1334", control_type="RadioButton").click_input()
vessel_coeff_window.child_window(title="Automatic", auto_id="1452", control_type="RadioButton").click_input()
vessel_coeff_window.child_window(title="OK", auto_id="1", control_type="Button").click_input()

# CHANGE SPEED IN INITIAL CONDITIONS
conditions_menu.click_input()
initial_conditions.click_input()
initial_conditions_window = pwrs.child_window(title="Initial Conditions", control_type="Window")
initial_conditions_window.type_keys('{TAB 5}')
initial_conditions_window.type_keys(input_speed)
initial_conditions_window.type_keys('{TAB}')
initial_conditions_window.type_keys('0')
initial_conditions_window.child_window(title="OK", auto_id="1", control_type="Button").click_input()

# CHANGE SPEED IN PROPULSION CONDITIONS
conditions_menu.click_input()
propulsion_conditions.click_input()
propulsion_conditions_window = pwrs.child_window(title="Thrust/Propulsion Conditions", control_type="Window")
propulsion_conditions_window.child_window(title="Constant:", auto_id="1228", control_type="RadioButton").click_input()
propulsion_conditions_window.type_keys('{TAB}')
propulsion_conditions_window.type_keys(input_speed)
propulsion_conditions_window.child_window(title="OK", auto_id="1", control_type="Button").click_input()

# INCIDENT WAVES
conditions_menu.click_input()
incident_waves.click_input()
incident_waves_window = pwrs.child_window(title="Incident Waves", control_type="Window")
incident_waves_window.type_keys('{TAB 4}{DOWN 3}{TAB}0.02{TAB}')
incident_waves_window.type_keys(input_wavelength)
incident_waves_window.child_window(title="OK", auto_id="1", control_type="Button").click_input()

# RUN SIMULATION
start_time = time.time()    
pwrs.type_keys('{F5}')
pwrs.child_window(title="pwrs", control_type="Window").type_keys('{ENTER}')
pwrs.child_window(title="pwrs", control_type="Window").type_keys('{ENTER}')
status_message = pwrs.child_window(auto_id="59393", control_type="StatusBar").Static2.window_text()
status = status_message == 'POWERSEA: Ready'
if not status:
    print('Running simulation for wavelength = ', input_wavelength, ', planing speed = ', input_speed, 'm/s')
while (not status):
    status_message = pwrs.child_window(auto_id="59393", control_type="StatusBar").Static2.window_text()
    status = status_message == 'POWERSEA: Ready'
end_time = time.time()
print('Elapsed time:', (end_time - start_time), 'sec')
print('Simulation complete!')

# SAVE DATA
results_menu.click_input()
save_data.click_input()
save_results_window = pwrs.child_window(title="Save Simulation Results", control_type="Window")

# CLEAR DATA AND SELECT NEW PARAMTERS TO SAVE
save_results_window.child_window(title="Remove All", auto_id="1119", control_type="Button").click_input()
save_results_window.type_keys('{TAB 5}')
save_results_window.type_keys(LCG)
save_results_window.type_keys('{TAB}')
save_results_window.type_keys(VCG)
save_results_window.type_keys('{TAB}{SPACE}{DOWN 3}{SPACE}{DOWN}{SPACE}{DOWN}{SPACE}{DOWN}{SPACE}{DOWN}{SPACE}{DOWN 8}{SPACE}{DOWN}{SPACE}{DOWN}{SPACE}{DOWN 2}{SPACE}{DOWN 4}{SPACE}')
save_results_window.child_window(title="Add", auto_id="1114", control_type="Button").click_input()

# SET FORMAT OF DATA
save_results_window.child_window(title="Format", auto_id="1201", control_type="Button").click_input()
save_results_window.child_window(title="Report Format", control_type="Window").child_window(title="Comma", auto_id="1112", control_type="RadioButton").click_input()
save_results_window.child_window(title="Report Format", control_type="Window").child_window(title="OK", auto_id="1", control_type="Button").click_input()
save_results_window.child_window(title="OK", auto_id="1", control_type="Button").click_input()

# SET PATH AND FILE NAME
save_output_window = pwrs.child_window(title="Save 'Output' Data", control_type="Window")
output_path = '"D:\\RESEARCH\\Navy_Boat_Project\\Curved-crease-origami-hull\\POWERSEA\\Data\\Wavy_Speed_' + input_speed + '_wavelength_' + input_wavelength + '.txt"'
save_output_window.type_keys(output_path, with_spaces=True)
save_output_window.child_window(title="Save", auto_id="1", control_type="Button").click_input()

###################################################
#####            CLOSING SEQUENCE             #####
###################################################
pwrs.type_keys('%{F4}')
pwrs.child_window(title="pwrs", control_type="Window").child_window(title="Yes", auto_id="6", control_type="Button").click_input()
save_path = '"D:\\RESEARCH\\Navy_Boat_Project\\Curved-crease-origami-hull\\POWERSEA\\POWERSEA1.pws"'
pwrs.child_window(title="Save As", control_type="Window").type_keys(save_path, with_spaces=True)
pwrs.child_window(title="Save As", control_type="Window").type_keys('{ENTER}')
time.sleep(3)

print('ALL SIMULATIONS COMPLETE! :)')