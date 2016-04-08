#!/usr/bin/env python3

"""
              _____      _   _  _____ ______          __   __
             / ____|    | \ | |/ ____|  ____|         \ \ / /
            | |     ___ |  \| | (___ | |__   _ __  ___ \ V /
            | |    / _ \| . ` |\___ \|  __| | '_ \/ __| > <
            | |___| (_) | |\  |____) | |____| | | \__ \/ . \
             \_____\___/|_| \_|_____/|______|_| |_|___/_/ \_\

     Compliance of NMR-derived Structural Ensembles with experimental data

Authors:    Zoltán Gáspári, Dániel Dudola
Fork me at: https://github.com/derPuntigamer/CoNSEnsX
"""

# standard modules
import os
import sys
import time

# own modules
import consensx.csx_libs.calc    as csx_calc
import consensx.csx_libs.parser  as csx_parser
import consensx.csx_libs.methods as csx_func
import consensx.csx_libs.output  as csx_out
import consensx.csx_libs.objects as csx_obj

from .models import CSX_upload

version = 1.9

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def run_calculation(calc_id):
    ts = time.time()

    abspath = os.path.abspath(__file__)
    dirname = os.path.dirname(abspath)
    csx_func.check_3rd_party(dirname)

    my_id = calc_id

    print("ID IS:", my_id)

    my_path = os.path.join(BASE_DIR, 'static', 'calculations', my_id) + '/'
    DB_entry = CSX_upload.objects.get(id_code=calc_id)


    #--------------  Setting up working directory and results HTML  -----------#
    csx_obj.CSV_buffer.working_dir = my_path
    csx_out.writeHeaderHTML(my_path, version)

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


    #----------------------  Read  and parse STR file   -----------------------#
    if DB_entry.STR_file:
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

    return "ALL OK!"


if __name__ == '__main__':


    #------------------------  Read  and parse STR file   -------------------------#
    if args.NOE_file:
        saveShifts = csx_func.getNOE(args.NOE_file)
        NOE_violations = csx_calc.calcNOEviolations(args, saveShifts,
                                                    my_path, args.r3_averaging)
        PRIDE_data = csx_calc.calcNMR_Pride(pdb_models, my_path)
        csx_out.writeFileTable(my_path, args, my_PDB, my_id,
                               model_count, args.NOE_file,
                               csx_obj.Restraint_Record.getRestraintCount())
        csx_out.write_bottom_table(my_path, NOE_violations, PRIDE_data)

        data_found = True
    else:
        csx_out.writeFileTable(my_path, args, my_PDB, my_id, model_count)


    #---------------------------------  S2 calc  ----------------------------------#
        S2_dict = csx_func.parseS2_STR(parsed.value)
        S2_dump = [S2_dict, args.fit, args.fit_range]
        S2_dict_path = my_path + "/S2_dict.pickle"
        pickle.dump(S2_dump, open(S2_dict_path, "wb"))

        if S2_dict:
            csx_calc.calcS2(S2_dict, my_path, fit=args.fit, fit_range=args.fit_range)
            data_found = True

        S2_sidechain = csx_func.parse_sidechain_S2_STR(parsed.value)

        if S2_sidechain:
            csx_calc.calcS2_sidechain(S2_sidechain, my_path, fit=args.fit)
            data_found = True


    #-----------------------------  J-coupling calc  ------------------------------#
        Jcoup_dict = csx_func.parseJcoup_STR(parsed.value)
        Jcoup_dict_path = my_path + "/Jcoup_dict.pickle"
        pickle.dump(Jcoup_dict, open(Jcoup_dict_path, "wb"))

        if Jcoup_dict:
            csx_calc.calcJCouplings(args.d, Jcoup_dict, my_PDB, my_path)
            data_found = True


    #----------------------------  Chemical shift calc  ---------------------------#
        ChemShift_lists = csx_func.parseChemShift_STR(parsed.value)
        ChemShift_lists_path = my_path + "/ChemShift_lists.pickle"
        pickle.dump(ChemShift_lists, open(ChemShift_lists_path, "wb"))

        if ChemShift_lists:
            csx_calc.calcChemShifts(ChemShift_lists, pdb_models, my_path)
            data_found = True

    csx_obj.CSV_buffer.writeCSV()
    csx_out.close_HTML(my_path)
    csx_out.add_PHP_variables(my_path, csx_obj.PHP_variables.PHP_dict, model_count)

    te = time.time()

    if data_found:
        print("total runtime", te-ts)
        print("Your ID was: \033[0;35m" + my_id + "\033[0m")
    else:
        print("ERROR, no data for back-calculation found!")
