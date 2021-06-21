import pandas as pd
import numpy as np
import time
import cantera as ct
from matplotlib import pyplot as plt
import csv
import math
import os
import sys
import re
import itertools
import logging
from collections import defaultdict
import git
import json

from rmgpy.molecule import Molecule
from rmgpy.data.base import Database

def save_pictures(git_path="", species_path="", overwrite=False):
    """
    Save a folder full of molecule pictures, needed for the pretty dot files.

    Saves them in the results directory, in a subfolder "species_pictures".
    Unless you set overwrite=True, it'll leave alone files that are
    already there.
    """
    dictionary_filename = git_path + "/chemkin/species_dictionary.txt"
    specs = Database().get_species(dictionary_filename, resonance=False)

    images_dir = os.path.join(species_path)
    os.makedirs(images_dir, exist_ok=True)
    for name, species in specs.items():
        filepath = os.path.join(images_dir, name + ".png")
        if not overwrite and os.path.exists(filepath):
            continue
        species.molecule[0].draw(filepath)

def prettydot(species_path, dotfilepath, strip_line_labels=False):
    """
    Make a prettier version of the dot file (flux diagram)

    Assumes the species pictures are stored in a directory
    called 'species_pictures' alongside the dot file.
    """
    import os, sys, re
    import subprocess
    pictures_directory = f"{species_path}/"

    if strip_line_labels:
        print("stripping edge (line) labels")

    reSize = re.compile('size="5,6"\;page="5,6"')
    #reNode = re.compile('(?P<node>s\d+)\ \[\ fontname="Helvetica",\ label="(?P<label>[^"]*)"\]\;')
    reNode = re.compile('(?P<node>s\d+)\ \[\ fontname=\"Helvetica\",\ label=\"(?P<label>[^\"]*)\"\]\;')
    #rePicture = re.compile("(?P<smiles>.+?)\((?P<id>\d+)\)\.png")
    rePicture = re.compile('(?P<smiles>.+?)\((?P<id>\d+)\)\.png')
    #reLabel = re.compile("(?P<name>.+?)\((?P<id>\d+)\)$")
    reLabel = re.compile('(?P<name>.+?)\((?P<id>\d+)\)$')

    species_pictures = dict()
    for picturefile in os.listdir(pictures_directory):
        match = rePicture.match(picturefile)
        if match:
            species_pictures[match.group("id")] = picturefile
        else:
            pass
            # print(picturefile, "didn't look like a picture")

    filepath = dotfilepath

    if not open(filepath).readline().startswith("digraph"):
        raise ValueError("{0} - not a digraph".format(filepath))

    infile = open(filepath)
    prettypath = filepath.replace(".dot", "", 1) + "-pretty.dot"
    outfile = open(prettypath, "w")

    for line in infile:
        (line, changed_size) = reSize.subn('size="12,12";page="12,12"', line)
        match = reNode.search(line)
        if match:
            label = match.group("label")
            idmatch = reLabel.match(label)
            if idmatch:
                idnumber = idmatch.group("id")
                if idnumber in species_pictures:
                    line = (f'%s [ image="{pictures_directory}%s" label="" width="0.5" height="0.5" imagescale=false fixedsize=false color="none" ];\n'% (match.group("node"), species_pictures[idnumber]))

        # rankdir="LR" to make graph go left>right instead of top>bottom
        if strip_line_labels:
            line = re.sub('label\s*=\s*"\s*[\d.]+"', 'label=""', line)

        # change colours
        line = re.sub('color="0.7,\ (.*?),\ 0.9"', r'color="1.0, \1, 0.7*\1"', line)

        outfile.write(line)

    outfile.close()
    infile.close()
    print(f"Graph saved to: {prettypath}")
    os.system(f'dot {prettypath} -Tpng -o{prettypath.replace(".dot", "", 1) + ".png"} -Gdpi=200')
    return prettypath

def show_flux_diagrams(self, suffix="", embed=False):
    """
    Shows the flux diagrams in the notebook.
    Loads them from disk.
    Does not embed them, to keep the .ipynb file small,
    unless embed=True. Use embed=True if you might over-write the files,
    eg. you want to show flux at different points.
    """
    import IPython

    for element in "CHONX":
        for phase_object in (self.gas, self.surf):
            phase = phase_object.name
            img_file = (
                f"reaction_path_{element}_{phase}{'_' if suffix else ''}{suffix}.png"
            )
            display(IPython.display.HTML(f"<hr><h2>{element} {phase}</h2>"))
            if embed:
                display(IPython.display.Image(filename=img_file, width=400, embed=True))
            else:
                display(IPython.display.Image(url=img_file, width=400, embed=False))

        # Now do the combined
        img_file = f"reaction_path_mass{'_' if suffix else ''}{suffix}.png"
        display(IPython.display.HTML(f"<hr><h2>Combined mass</h2>"))
        if embed:
            display(IPython.display.Image(filename=img_file, width=400, embed=True))
        else:
            display(IPython.display.Image(url=img_file, width=400, embed=False))

def save_flux_diagrams(*phases, suffix="", timepoint="", species_path=""):
    """
    Saves the flux diagrams. The filenames have a suffix if provided,
    so you can keep them separate and not over-write.
    """
    for element in "CHONX":
        for phase_object in phases:
            phase = phase_object.name

            diagram = ct.ReactionPathDiagram(phase_object, element)
            diagram.title = f"Reaction path diagram following {element} in {phase}"
            diagram.label_threshold = 0.001

            dot_file = f"{suffix}/reaction_path_{element}_{phase}_{timepoint}.dot"
            img_file = f"{suffix}/reaction_path_{element}_{phase}_{timepoint}.png"
            dot_bin_path = (
                "/Users/lee.ting/Code/anaconda3/pkgs/graphviz-2.40.1-hefbbd9a_2/bin/dot" 
            )
            img_path = os.path.join(os.getcwd(), img_file)
            diagram.write_dot(dot_file)

            #also make a prettydot file
            prettydot(species_path, dot_file, strip_line_labels=False)

            # print(diagram.get_data())

            print(
                f"Wrote graphviz input file to '{os.path.join(os.getcwd(), dot_file)}'."
            )
            os.system(f"dot {dot_file} -Tpng -o{img_file} -Gdpi=200")
            print(f"Wrote graphviz output file to '{img_path}'.")

def run_reactor(
    cti_file,
    t_array=[700],
    surf_t_array=[700],
    p_array=[1],
    v_array=[1e-6],
    o2_array=[0.02],
    nh3_array=[0.01],
    h2o_array=[0.05],
    he_array=[0.91],
    rtol=1.0e-11,
    atol=1.0e-22,
    reactor_type=1,
    energy="off",
    sensitivity=False,
    sensatol=1e-6,
    sensrtol=1e-6,
    reactime=1e5,
):


    try:
        array_i = int(os.getenv("SLURM_ARRAY_TASK_ID"))
    except TypeError:
        array_i = 0

    # get git commit hash and message
    rmg_model_path = "../../../Rh211"
    repo = git.Repo(rmg_model_path)
    date = time.localtime(repo.head.commit.committed_date)
    git_date = f"{date[0]}_{date[1]}_{date[2]}_{date[3]}{date[4]}"
    git_sha = str(repo.head.commit)[0:6]
    git_msg = str(repo.head.commit.message)[0:50].replace(" ", "_").replace("'", "_").replace("\n", "")
    git_file_string = f"{git_date}_{git_sha}_{git_msg}"

    # set sensitivity string for file path name
    if sensitivity:
        sensitivity_str = "on"
    else: 
        sensitivity_str = "off"

    # this should probably be outside of function
    settings = list(itertools.product(t_array, surf_t_array, p_array, v_array, o2_array, nh3_array, h2o_array))

    # constants
    pi = math.pi

    # set initial temps, pressures, concentrations
    temp = settings[array_i][0]  # kelvin
    temp_str = str(temp)[0:100]
    
    surf_temp = settings[array_i][1]
    surf_temp_str = str(surf_temp)[0:10]
    
    pressure = settings[array_i][2] * ct.one_atm  # Pascals

    X_o2 = settings[array_i][4]
    x_O2_str = str(X_o2)[0:10].replace(".", "_")
    
    X_nh3 = (settings[array_i][5])
    x_NH3_str = str(X_nh3)[0:10].replace(".", "_")
    
    X_h2o = (settings[array_i][6])
    x_H2O_str = str(X_h2o)[0:10].replace(".", "_")

    X_he = 1 - X_o2 - X_nh3 - X_h2o
    
    mw_nh3 = 17.0306e-3  # [kg/mol]
    mw_o2 = 31.999e-3  # [kg/mol]
    mw_h2o = 18.015e-3  # [kg/mol]
    mw_he = 4.003e-3  # [kg/mol]
    o2_ratio = X_nh3 / X_o2
    # O2/NH3/H2O/He: typical is
    concentrations_rmg = {"O2(2)": X_o2, "NH3(6)": X_nh3, "H2O(3)": X_h2o, "He": X_he}
 
    # initialize cantera gas and surface
    gas = ct.Solution(cti_file, "gas")
    surf = ct.Interface(cti_file, "surface1", [gas])

    # initialize T and P
    gas.TPX = temp, pressure, concentrations_rmg
    surf.TP = temp, pressure # change this to surf_temp when we want a different starting temperature for the surface

    # if a mistake is made with the input, 
    # cantera will normalize the mole fractions. 
    # make sure that we are reporting/using 
    # the normalized values
    X_o2 = float(gas["O2(2)"].X)
    X_nh3 = float(gas["NH3(6)"].X)
    X_h2o = float(gas["H2O(3)"].X)
    X_he = float(gas["He"].X)


    # create gas inlet
    inlet = ct.Reservoir(gas)

    # create gas outlet
    exhaust = ct.Reservoir(gas)

    # Reactor volume
    rvol = 1e-6 # 0.00014711#1.4711 #m3

    # Catalyst Surface Area
    site_density = (
        surf.site_density * 1000
    )  # [mol/m^2]cantera uses kmol/m^2, convert to mol/m^2
    cat_area = 1e-2#100*0.00065461 #6.5461  #1.387*1.5023*3.1416=6.5461[m^3] 
    
    # reactor initialization
    if reactor_type == 0:
        r = ct.Reactor(gas, energy=energy)
        reactor_type_str = "Reactor"
    elif reactor_type == 1:
        r = ct.IdealGasReactor(gas, energy=energy)
        reactor_type_str = "IdealGasReactor"
    elif reactor_type == 2:
        r = ct.ConstPressureReactor(gas, energy=energy)
        reactor_type_str = "ConstPressureReactor"
    elif reactor_type == 3:
        r = ct.IdealGasConstPressureReactor(gas, energy=energy)
        reactor_type_str = "IdealGasConstPressureReactor"

    # calculate the available catalyst area in a differential reactor
    rsurf = ct.ReactorSurface(surf, r, A=cat_area)
    r.volume = rvol
    surf.coverages = "X(1):1.0"

    # flow controllers
    one_atm = ct.one_atm
    FC_temp = 293.15
    volume_flow = settings[array_i][3]  # [m^3/s]
    molar_flow = volume_flow * one_atm / (8.3145 * FC_temp)  # [mol/s]
    mass_flow = molar_flow * (X_nh3 * mw_nh3 + X_o2 * mw_o2 + X_h2o * mw_h2o+ X_he * mw_he)  # [kg/s]
    mfc = ct.MassFlowController(inlet, r, mdot=mass_flow)

    # A PressureController has a baseline mass flow rate matching the 'master'
    # MassFlowController, with an additional pressure-dependent term. By explicitly
    # including the upstream mass flow rate, the pressure is kept constant without
    # needing to use a large value for 'K', which can introduce undesired stiffness.
    outlet_mfc = ct.PressureController(r, exhaust, master=mfc, K=0.01)

    # initialize reactor network
    sim = ct.ReactorNet([r])

    # set relative and absolute tolerances on the simulation
    sim.rtol = 1.0e-11
    sim.atol = 1.0e-22
    #################################################
    # Run single reactor
    #################################################

    # round numbers so they're easier to read
    # temp_str = '%s' % '%.3g' % tempn
    cat_area_str = "%s" % "%.3g" % cat_area
    species_path = (
        os.path.dirname(os.path.abspath(__file__))
        + f"/results/{git_file_string}/species_pictures"
    )
    print(species_path)
    results_path = (
        os.path.dirname(os.path.abspath(__file__))
       + f"/results/{git_file_string}/{reactor_type_str}/energy_{energy}/sensitivity_{sensitivity_str}/results"
    )
    logging.warning(f"Saving results in {results_path}, the file's name is _temp_{temp}_O2_{x_O2_str}_NH3_{x_NH3_str}.csv")
    flux_path = (
        os.path.dirname(os.path.abspath(__file__))
        + f"/results/{git_file_string}/{reactor_type_str}/energy_{energy}/sensitivity_{sensitivity_str}/flux_diagrams"
    )
    # create species folder for species pictures if it does not already exist
    try:
        os.makedirs(species_path, exist_ok=True)
        save_pictures(git_path=rmg_model_path, species_path=species_path)
    except OSError as error:
        print(error)

    try:
        os.makedirs(results_path, exist_ok=True)
    except OSError as error:
        print(error)

    try:
        os.makedirs(flux_path, exist_ok=True)
    except OSError as error:
        print(error)

    gas_ROP_str = [i + " ROP [kmol/m^3 s]" for i in gas.species_names]

    # surface ROP reports gas and surface ROP. these values are not redundant
    gas_surf_ROP_str = [i + " surface ROP [kmol/m^2 s]" for i in gas.species_names]
    surf_ROP_str = [i + " ROP [kmol/m^2 s]" for i in surf.species_names]

    # gasrxn_ROP_str = [i + " ROP [kmol/m^3 s]" for i in gas.reaction_equations()]
    surfrxn_ROP_str = [i + " ROP [kmol/m^2 s]" for i in surf.reaction_equations()]

    output_filename = (
        results_path
        + f"/Spinning_basket_area_{cat_area_str}_energy_{energy}"
        + f"_temp_{temp}_O2_{x_O2_str}_NH3_{x_NH3_str}.csv"
    )

    outfile = open(output_filename, "w")
    writer = csv.writer(outfile)

    # Sensitivity atol, rtol, and strings for gas and surface reactions if selected
    # slows down script by a lot
    if sensitivity:
        sim.rtol_sensitivity = sensrtol
        sim.atol_sensitivity = sensatol
        sens_species = ["NH3(6)", "O2(2)", "N2(4)", "NO(5)", "N2O(7)"]  #change THIS to your species, can add "," and other species

        # turn on sensitive reactions/species
        # for i in range(gas.n_reactions):
        #     r.add_sensitivity_reaction(i)

        for i in range(surf.n_reactions):
            rsurf.add_sensitivity_reaction(i)

        # for i in range(gas.n_species):
        #     r.add_sensitivity_species_enthalpy(i)

        # for i in range(surf.n_species):
        #     rsurf.add_sensitivity_species_enthalpy(i)

        for j in sens_species:
            # gasrxn_sens_str = [
            #     j + " sensitivity to " + i for i in gas.reaction_equations()
            # ]
            surfrxn_sens_str = [
                j + " sensitivity to " + i for i in surf.reaction_equations()
            ]
            # gastherm_sens_str = [j + " thermo sensitivity to " + i for i in gas.species_names]
            # surftherm_sens_str = [j + " thermo sensitivity to " + i for i in surf.species_names]
            sens_list = gasrxn_sens_str + surfrxn_sens_str  # + gastherm_sens_str

        writer.writerow(
            [
                "T (K)",
                "P (Pa)",
                "V (M^3/s)",
                "X_nh3 initial",
                "X_o2 initial",
                "X_h2o initial",
                "X_he initial",
                "(NH3/O2)",
                "T (K) final",
                "Rtol",
                "Atol",
                "reactor type",
                "energy on?"
            ]
            + gas.species_names
            + surf.species_names
            + gas_ROP_str
            + gas_surf_ROP_str
            + surf_ROP_str
            # + gasrxn_ROP_str
            + surfrxn_ROP_str
            + sens_list
        )

    else:
        writer.writerow(
            [
                "T (K)",
                "P (Pa)",
                "V (M^3/s)",
                "X_nh3 initial",
                "X_o2 initial",
                "X_h2o initial",
                "X_he initial",
                "(NH3/O2)",
                "T (K) final",
                "Rtol",
                "Atol",
                "reactor type",
                "energy on?"
            ]
            + gas.species_names
            + surf.species_names
            + gas_ROP_str
            + gas_surf_ROP_str
            + surf_ROP_str
            # + gasrxn_ROP_str
            + surfrxn_ROP_str
        )

    t = 0.0
    dt = 1e4
    iter_ct = 0
    # run the simulation
    first_run = True

    while t < reactime:
        # save flux diagrams at beginning of run
        if first_run == True:
            save_flux_diagrams(gas, suffix=flux_path, timepoint="beginning", species_path=species_path)
            save_flux_diagrams(surf, suffix=flux_path, timepoint="beginning", species_path=species_path)

            first_run = False
        t += dt
        sim.advance(t)
#         if t % 10 < 0.01:

        if sensitivity:
            # get sensitivity for sensitive species i (e.g. methanol) in reaction j
            for i in sens_species:
                g_nrxn = gas.n_reactions
                s_nrxn = surf.n_reactions
                # g_nspec = gas.n_species
                # s_nspec = surf.n_species

                gas_sensitivities = [sim.sensitivity(i, j) for j in range(g_nrxn)]
                surf_sensitivities = [
                    sim.sensitivity(i, j) for j in range(g_nrxn, g_nrxn + s_nrxn)
                ]
                # gas_therm_sensitivities = [sim.sensitivity(i,j)
                # for j in range(g_nrxn+s_nrxn,g_nrxn+s_nrxn+g_nspec)]
                # surf_therm_sensitivities = [sim.sensitivity(i,j)
                # for j in range(g_nrxn+s_nrxn+g_nspec,g_nrxn+s_nrxn+g_nspec+s_nspec)]

                sensitivities_all = (
                    gas_sensitivities
                    + surf_sensitivities
                    # + gas_therm_sensitivities
                )

            writer.writerow(
                [
                    temp,
                    pressure,
                    volume_flow,
                    X_nh3,
                    X_o2,
                    X_h2o,
                    X_he,
                    o2_ratio,
                    gas.T,
                    sim.rtol,
                    sim.atol,
                    reactor_type_str,
                    energy,
                ]
                + list(gas.X)
                + list(surf.X)
                + list(gas.net_production_rates)
                + list(surf.net_production_rates)
                # + list(gas.net_rates_of_progress)
                + list(surf.net_rates_of_progress)
                + sensitivities_all
            )

        else:
            writer.writerow(
                [
                    temp,
                    pressure,
                    volume_flow,
                    X_nh3,
                    X_o2,
                    X_h2o,
                    X_he,
                    o2_ratio,
                    gas.T,
                    sim.rtol,
                    sim.atol,
                    reactor_type_str,
                    energy,
                ]
                + list(gas.X)
                + list(surf.X)
                + list(gas.net_production_rates)
                + list(surf.net_production_rates)
                # + list(gas.net_rates_of_progress)
                + list(surf.net_rates_of_progress)
            )


        iter_ct += 1

    outfile.close()

    # save flux diagrams at the end of the run
    save_flux_diagrams(gas, suffix=flux_path, timepoint="end", species_path=species_path)
    save_flux_diagrams(surf, suffix=flux_path, timepoint="end", species_path=species_path)
    return

#######################################################################
# Input Parameters for combustor
#######################################################################

# filepath for writing files
git_repo = "../../../Rh211/"
cti_file = git_repo + "cantera/chem_annotated.cti"

# Reactor settings arrays for run
Temps = np.linspace(400,1300,100)
#Temps = [300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400]
Pressures = [1] # 1 bar
volume_flows = [6e-5] #10 ml/min = 1.6666666666667E-7 m^3/s

O2_fraction = [0.02] #O2 partial pressure(atm)
NH3_fraction = [0.001]
H2O_fraction = [0.05]

# reaction time
reactime = 1e8

# sensitivity settings
sensitivity = False
sensatol = 1e-6
sensrtol = 1e-6

run_reactor(
    cti_file=cti_file,
    t_array=Temps,
    reactor_type=1,
    v_array=volume_flows,
    o2_array=O2_fraction,
    nh3_array=NH3_fraction,
    h2o_array=H2O_fraction,
    energy="off",
    sensitivity=sensitivity,
    sensatol=sensatol,
    sensrtol=sensrtol,
    reactime=reactime,
)
