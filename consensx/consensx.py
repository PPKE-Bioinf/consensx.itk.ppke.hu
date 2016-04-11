#!/usr/bin/env python3

"""
This is the Django webapp implemetation of
              _____      _   _  _____ ______          __   __
             / ____|    | \ | |/ ____|  ____|         \ \ / /
            | |     ___ |  \| | (___ | |__   _ __  ___ \ V /
            | |    / _ \| . ` |\___ \|  __| | '_ \/ __| > <
            | |___| (_) | |\  |____) | |____| | | \__ \/ . \
             \_____\___/|_| \_|_____/|______|_| |_|___/_/ \_\

     Compliance of NMR-derived Structural Ensembles with experimental data

Authors: Zoltán Gáspári, Dániel Dudola
Fork of: https://github.com/derPuntigamer/CoNSEnsX
"""

# standard modules
import os
import time

# own modules
import consensx.csx_libs.calc    as csx_calc
import consensx.csx_libs.methods as csx_func
import consensx.csx_libs.objects as csx_obj

# Django server
from django.shortcuts import render
from .models import CSX_upload

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def run_calculation(request, calc_id):
    ts = time.time()

    abspath = os.path.abspath(__file__)
    dirname = os.path.dirname(abspath)
    csx_func.check_3rd_party(dirname)

    my_id = calc_id
    print("ID IS:", my_id)

    from django.conf import settings
    my_path = os.path.join(settings.MEDIA_ROOT, my_id) + '/'
    DB_entry = CSX_upload.objects.get(id_code=calc_id)


    #-----------------  Setting up working directory and files  ---------------#
    csx_obj.CSV_buffer.working_dir = my_path

    my_PDB = my_path + DB_entry.PDB_file

    csx_func.pdb_cleaner(my_PDB)                   # bringing PDB to format
    model_count = csx_func.pdb_splitter(my_path, my_PDB)

    csx_func.get_model_list(my_PDB, model_count)

    pdb_models = []                                # list of models (PDB)
    for file in os.listdir(my_path):
        if file.startswith("model_") and file.endswith(".pdb"):
            pdb_models.append(file)

    pdb_models = csx_func.natural_sort(pdb_models)
    data_found = False

    #----------------------  Read  and parse NOE file   -----------------------#
    if DB_entry.NOE_file:
        NOE_name = DB_entry.NOE_file
        my_NOE = my_path + DB_entry.NOE_file
        saveShifts = csx_func.getNOE(my_NOE)
        NOE_violations = csx_calc.calcNOEviolations(my_NOE, saveShifts,
                                                    my_path, DB_entry.r3average)
        NOE_violations = str(NOE_violations) + " distance restraints found"
        PRIDE_data = csx_calc.calcNMR_Pride(pdb_models, my_path)

        NOE_PRIDE_data = {
            "NOE_violations": NOE_violations,
            "NOE_hist": os.path.join(settings.MEDIA_ROOT, my_id,
                                     "NOE_hist.svg"),
            "best_score": PRIDE_data[0],
            "worst_score": PRIDE_data[1],
            "average_score": '{0:.3f}'.format(PRIDE_data[2]),
            "deviation": '{0:.3f}'.format(PRIDE_data[3]),
            "PRIDE_hist": str(os.path.join(settings.MEDIA_ROOT, my_id,
                                                 "PRIDE-NMR_score.svg"))
        }


        print(str(os.path.join(settings.MEDIA_ROOT, my_id,
                                                 "PRIDE-NMR_score.svg")))
        data_found = True
    else:
        NOE_name = "[NOT PRESENT]"
        NOE_violations = ""
        NOE_PRIDE_data = None


    #----------------------  Read  and parse STR file   -----------------------#
    if DB_entry.STR_file:
        STR_name = DB_entry.STR_file
        my_STR = my_path + DB_entry.STR_file
        try:
            parsed = csx_func.parseSTR(my_STR)
        except Exception as e:
            return e

        #-------------------------  RDC calculation  --------------------------#
        RDC_lists = csx_func.get_RDC_lists(parsed.value)

        if RDC_lists:
            SVD_enabled = DB_entry.svd_enable
            lc_model = DB_entry.rdc_lc
            csx_calc.calcRDC(RDC_lists, pdb_models, my_path,
                             SVD_enabled, lc_model)
            data_found = True

        #-----------------------------  S2 calc  ------------------------------#
        S2_dict = csx_func.parseS2_STR(parsed.value)

        if S2_dict:
            csx_calc.calcS2(
                S2_dict, my_path,
                fit=DB_entry.superimpose,
                fit_range=DB_entry.fit_range
            )
            data_found = True

        S2_sidechain = csx_func.parse_sidechain_S2_STR(parsed.value)

        if S2_sidechain:
            csx_calc.calcS2_sidechain(S2_sidechain, my_path,
                                      fit=DB_entry.superimpose)
            data_found = True

        #-------------------------  J-coupling calc  --------------------------#
        Jcoup_dict = csx_func.parseJcoup_STR(parsed.value)

        if Jcoup_dict:
            csx_calc.calcJCouplings(DB_entry.karplus, Jcoup_dict,
                                    my_PDB, my_path)
            data_found = True

        #------------------------  Chemical shift calc  -----------------------#
        ChemShift_lists = csx_func.parseChemShift_STR(parsed.value)

        if ChemShift_lists:
            csx_calc.calcChemShifts(ChemShift_lists, pdb_models, my_path)
            data_found = True
    else:
        STR_name = "[NOT PRESENT]"


    csx_obj.CSV_buffer.writeCSV()

    te = time.time()

    if data_found:
        return render(request, "consensx/calculation.html", {
            "my_id": my_id,
            "my_PDB": DB_entry.PDB_file,
            "n_model": model_count,
            "my_NOE": NOE_name,
            "n_NOE" : NOE_violations,
            "my_STR": STR_name,
            "NOE_PRIDE_data": NOE_PRIDE_data
        })
    else:
        return "NO DATA FOUND IN STAR-NMR FILE"
