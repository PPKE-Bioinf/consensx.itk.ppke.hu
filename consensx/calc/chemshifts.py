import pickle
import os
import multiprocessing
from multiprocessing.pool import ThreadPool
import subprocess
from time import perf_counter

import consensx.graph as graph

from consensx import thirdparty
from .measure import correlation, q_value, rmsd
from consensx.parse import ChemShiftRecord
from consensx.misc.natural_sort import natural_sort


chemshift_types = ["HA", "CA", "CB", "N", "H", "C"]


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def shiftx_runner(file_list):
    for config in file_list:
        subprocess.call(
            [
                thirdparty.ThirdParty.shiftx,
                "1",
                config["pdb_file"],
                config["out_name"]
            ]
        )


def call_shiftx_on(my_path, pdb_files, bme_weights=None):
    """Call ShiftX on PDB models. Each output is appended to 'out_name'"""
    t1_start = perf_counter()
    shiftx_runs = []

    for i, pdb_file in enumerate(pdb_files):
        shiftx_runs.append({
            "pdb_file": my_path + pdb_file,
            "out_name": my_path + "/modell_" + str(i + 1) + ".out",
        })

    num_threads = int(multiprocessing.cpu_count() / 2)
    print(f"Running shiftx on {num_threads} threads")
    thread_pool = ThreadPool(num_threads)
    combination_chunks_list = list(chunks(shiftx_runs, num_threads))
    thread_pool.map(shiftx_runner, combination_chunks_list)

    t1_stop = perf_counter()
    print("[chemshift] Call Shiftx in seconds:", t1_stop - t1_start)

    shiftx_output_files = []
    average_ha, average_h, average_n = {}, {}, {}
    average_ca, average_cb, average_c = {}, {}, {}
    mod_ha, mod_h, mod_n, mod_ca, mod_cb, mod_c = {}, {}, {}, {}, {}, {}
    model_data_list = []
    model_weight = 1
    model_weight_sum = len(pdb_files)

    if bme_weights:
        model_weight_sum = sum(bme_weights)

    for some_file in natural_sort(os.listdir(my_path)):
        if some_file.startswith("modell") and some_file.endswith(".out"):
            shiftx_output_files.append(some_file)

    t1_start = perf_counter()

    for i, shiftx_output in enumerate(shiftx_output_files):
        if bme_weights:
            model_weight = bme_weights[i]

        out_file = open(my_path + shiftx_output)
        part = 0

        for line in out_file:
            if line.strip().startswith("NUM"):
                part += 1

                if part == 2:
                    break

                continue

            if line.strip().startswith("---"):
                continue

            if line.strip():
                line_values = line.split()

                try:
                    resnum = int(line_values[0])
                except ValueError:
                    resnum = int(line_values[0][1:])

                HA = float(line_values[2])
                H = float(line_values[3])
                N = float(line_values[4])
                CA = float(line_values[5])
                CB = float(line_values[6])
                C = float(line_values[7])

                mod_ha[resnum] = HA * model_weight
                mod_h[resnum] = H * model_weight
                mod_n[resnum] = N * model_weight
                mod_ca[resnum] = CA * model_weight
                mod_cb[resnum] = CB * model_weight
                mod_c[resnum] = C * model_weight

                if resnum in list(average_ha.keys()):
                    average_ha[resnum] += HA * model_weight
                else:
                    average_ha[resnum] = HA * model_weight

                if resnum in list(average_h.keys()):
                    average_h[resnum] += H * model_weight
                else:
                    average_h[resnum] = H * model_weight

                if resnum in list(average_n.keys()):
                    average_n[resnum] += N * model_weight
                else:
                    average_n[resnum] = N * model_weight

                if resnum in list(average_ca.keys()):
                    average_ca[resnum] += CA * model_weight
                else:
                    average_ca[resnum] = CA * model_weight

                if resnum in list(average_cb.keys()):
                    average_cb[resnum] += CB * model_weight
                else:
                    average_cb[resnum] = CB * model_weight

                if resnum in list(average_c.keys()):
                    average_c[resnum] += C * model_weight
                else:
                    average_c[resnum] = C * model_weight

        out_file.close()

        model_data_list.append(
            {
                "HA": mod_ha,
                "H": mod_h,
                "N": mod_n,
                "CA": mod_ca,
                "CB": mod_cb,
                "C": mod_c,
            }
        )

        mod_ha, mod_h, mod_n, mod_ca, mod_cb, mod_c = {}, {}, {}, {}, {}, {}

    t1_stop = perf_counter()
    print("[chemshift] Getting averages from Pales outputs in seconds:", t1_stop - t1_start)

    averages = [
        average_ha,
        average_h,
        average_n,
        average_ca,
        average_cb,
        average_c,
    ]

    t1_start = perf_counter()

    for avg_dict in averages:
        for key in avg_dict:
            avg_dict[key] /= model_weight_sum

    t1_stop = perf_counter()
    print("[chemshift] Calculating averages in seconds:", t1_stop - t1_start)

    return (
        {
            "HA": average_ha,
            "H": average_h,
            "N": average_n,
            "CA": average_ca,
            "CB": average_cb,
            "C": average_c,
        },
        model_data_list,
    )

chemshift_corrections_prev = {
    "ALA": {"H": 0.088,  "HA": -0.046, "C": 0.068,  "CA": 0.161,  "CB": -0.135, "N": -0.639},
    "ARG": {"H": 0.251,  "HA": -0.004, "C": -0.246, "CA": 0.029,  "CB": -0.038, "N": 1.909},
    "ASP": {"H": 0.081,  "HA": -0.024, "C": 0.059,  "CA": 0.290,  "CB": -0.182, "N": 0.679},
    "ASN": {"H": 0.130,  "HA": -0.004, "C": -0.052, "CA": 0.293,  "CB": -0.205, "N": 0.699},
    "CYS": {"H": 0.293,  "HA": 0.025,  "C": -0.417, "CA": -0.006, "CB": 0.008,  "N": 3.100},
    "GLU": {"H": 0.174,  "HA": -0.048, "C": -0.032, "CA": 0.171,  "CB": -0.119, "N": 1.226},
    "GLN": {"H": 0.225,  "HA": -0.009, "C": -0.153, "CA": 0.125,  "CB": -0.064, "N": 1.700},
    "GLY": {"H": 0.000,  "HA": 0.000,  "C": 0.000,  "CA": 0.000,  "CB": 0.000,  "N": 0.000},
    "HIS": {"H": 0.141,  "HA": -0.019, "C": -0.243, "CA": 0.069,  "CB": -0.093, "N": 1.761},
    "ILE": {"H": 0.250,  "HA": 0.010,  "C": -0.045, "CA": -0.087, "CB": -0.111, "N": 4.205},
    "LEU": {"H": 0.105,  "HA": -0.042, "C": -0.115, "CA": 0.083,  "CB": -0.143, "N": 0.835},
    "LYS": {"H": 0.174,  "HA": -0.031, "C": -0.117, "CA": 0.195,  "CB": -0.102, "N": 1.311},
    "MET": {"H": 0.205,  "HA": -0.018, "C": -0.118, "CA": 0.028,  "CB": -0.139, "N": 1.407},
    "PHE": {"H": 0.079,  "HA": -0.040, "C": -0.672, "CA": -0.163, "CB": 0.059,  "N": 2.238},
    "PRO": {"H": 0.292,  "HA": -0.032, "C": -0.040, "CA": 0.021,  "CB": -0.151, "N": 0.569},
    "SER": {"H": 0.230,  "HA": 0.006,  "C": -0.128, "CA": 0.155,  "CB": -0.153, "N": 2.297},
    "THR": {"H": 0.211,  "HA": -0.005, "C": -0.130, "CA": 0.143,  "CB": -0.072, "N": 2.680},
    "TRP": {"H": -0.555, "HA": -0.168, "C": -0.465, "CA": 0.014,  "CB": -0.046, "N": 0.795},
    "TYR": {"H": -0.052, "HA": -0.035, "C": -0.633, "CA": -0.245, "CB": 0.017,  "N": 2.729},
    "VAL": {"H": 0.259,  "HA": -0.007, "C": -0.198, "CA": 0.075,  "CB": -0.147, "N": 4.507},
}

chemshift_corrections_next = {
    "ALA": {"H": -0.056, "HA": -0.042, "C": -0.816, "CA": -0.187, "CB": 0.018,  "N": -0.004},
    "ARG": {"H": -0.039, "HA": -0.017, "C": -0.625, "CA": -0.258, "CB": 0.046,  "N": -0.065},
    "ASP": {"H": -0.003, "HA": -0.009, "C": -0.853, "CA": -0.095, "CB": 0.079,  "N": -0.373},
    "ASN": {"H": -0.023, "HA": -0.033, "C": -0.849, "CA": -0.067, "CB": -0.013, "N": -0.219},
    "CYS": {"H": 0.375,  "HA": -0.008, "C": -0.610, "CA": 0.047,  "CB": -0.170, "N": 0.785},
    "GLU": {"H": -0.014, "HA": -0.005, "C": -0.503, "CA": -0.071, "CB": 0.065,  "N": -0.100},
    "GLN": {"H": -0.016, "HA": -0.039, "C": -0.551, "CA": -0.098, "CB": 0.004,  "N": -0.087},
    "GLY": {"H": 0.000,  "HA": 0.000,  "C": 0.000,  "CA": 0.000,  "CB": 0.000,  "N": 0.000},
    "HIS": {"H": -0.136, "HA": -0.059, "C": -0.677, "CA": -0.248, "CB": 0.041,  "N": -0.225},
    "ILE": {"H": -0.071, "HA": 0.000,  "C": -0.559, "CA": -0.129, "CB": 0.093,  "N": -0.167},
    "LEU": {"H": -0.086, "HA": -0.021, "C": -0.612, "CA": -0.184, "CB": -0.053, "N": -0.397},
    "LYS": {"H": -0.053, "HA": -0.011, "C": -0.520, "CA": -0.142, "CB": 0.029,  "N": 0.068},
    "MET": {"H": -0.066, "HA": -0.018, "C": -0.271, "CA": -0.061, "CB": -0.205, "N": -0.145},
    "PHE": {"H": -0.133, "HA": -0.052, "C": -0.891, "CA": -0.247, "CB": 0.023,  "N": -0.484},
    "PRO": {"H": -0.028, "HA": 0.276,  "C": -2.874, "CA": -2.423, "CB": -0.426, "N": 1.072},
    "SER": {"H": -0.002, "HA": 0.030,  "C": -0.483, "CA": -0.151, "CB": 0.087,  "N": -0.053},
    "THR": {"H": 0.042,  "HA": 0.065,  "C": -0.327, "CA": -0.122, "CB": 0.086,  "N": 0.076},
    "TRP": {"H": 0.032,  "HA": -0.082, "C": -0.694, "CA": -0.201, "CB": 0.069,  "N": -0.273},
    "TYR": {"H": -0.093, "HA": -0.044, "C": -0.848, "CA": -0.148, "CB": 0.019,  "N": -0.387},
    "VAL": {"H": -0.052, "HA": 0.021,  "C": -0.662, "CA": -0.214, "CB": 0.044,  "N": -0.047},
}

chemshift_corrections = {
    "ALA": {"H": 8.158, "HA": 4.224, "C": 178.418, "CA": 52.599, "CB": 19.102, "N": 123.906},
    "ARG": {"H": 8.232, "HA": 4.239, "C": 176.821, "CA": 56.088, "CB": 30.691, "N": 121.288},
    "ASP": {"H": 8.217, "HA": 4.537, "C": 176.987, "CA": 54.331, "CB": 41.089, "N": 120.207},
    "ASN": {"H": 8.366, "HA": 4.632, "C": 175.825, "CA": 53.231, "CB": 38.790, "N": 118.668},
    "CYS": {"H": 8.410, "HA": 4.447, "C": 174.927, "CA": 58.327, "CB": 28.085, "N": 119.068},
    "GLU": {"H": 8.304, "HA": 4.222, "C": 177.125, "CA": 56.650, "CB": 30.225, "N": 120.769},
    "GLN": {"H": 8.258, "HA": 4.254, "C": 176.510, "CA": 55.840, "CB": 29.509, "N": 120.224},
    "GLY": {"H": 8.307, "HA": 3.980, "C": 174.630, "CA": 45.236, "CB": 0.00,   "N": 108.783},
    "HIS": {"H": 8.310, "HA": 4.585, "C": 175.349, "CA": 55.964, "CB": 29.719, "N": 118.930},
    "ILE": {"H": 7.963, "HA": 4.076, "C": 176.897, "CA": 61.247, "CB": 38.563, "N": 120.512},
    "LEU": {"H": 8.088, "HA": 4.260, "C": 178.037, "CA": 55.260, "CB": 42.212, "N": 121.877},
    "LYS": {"H": 8.221, "HA": 4.237, "C": 177.224, "CA": 56.412, "CB": 32.921, "N": 121.353},
    "MET": {"H": 8.209, "HA": 4.425, "C": 176.953, "CA": 55.591, "CB": 32.690, "N": 120.002},
    "PHE": {"H": 8.107, "HA": 4.573, "C": 176.368, "CA": 57.934, "CB": 39.660, "N": 120.138},
    "PRO": {"H": 0.00,  "HA": 4.339, "C": 177.542, "CA": 63.180, "CB": 32.072, "N": 136.612},
    "SER": {"H": 8.215, "HA": 4.392, "C": 175.236, "CA": 58.352, "CB": 63.766, "N": 115.935},
    "THR": {"H": 8.047, "HA": 4.252, "C": 175.122, "CA": 61.926, "CB": 69.794, "N": 114.024},
    "TRP": {"H": 7.725, "HA": 4.567, "C": 174.549, "CA": 57.500, "CB": 29.380, "N": 120.733},
    "TYR": {"H": 8.026, "HA": 4.504, "C": 176.284, "CA": 57.761, "CB": 38.750, "N": 120.228},
    "VAL": {"H": 8.037, "HA": 4.009, "C": 176.772, "CA": 62.347, "CB": 32.674, "N": 120.403},
}


def chemshifts(
    csv_buffer,
    calced_data_storage,
    chem_shift_lists,
    pdb_models,
    my_path,
    bme_weights,
):
    """Back calculate chemical shifts from given chemical shift list and PDB
       models"""
    cs_data = []
    cs_calced, model_data = call_shiftx_on(my_path, pdb_models, bme_weights)

    for n, cs_list in enumerate(chem_shift_lists):
        bme_exp_filename = "chemshift_" + str(n) + "_exp.dat"
        bme_calc_filename = "chemshift_" + str(n) + "_calc.dat"

        with open(my_path + bme_exp_filename, "w") as exp_dat_file:
            exp_dat_file.write("# DATA=CS PRIOR=GAUSS\n")

            for atom_type in chemshift_types:
                try:
                    for i in cs_list[atom_type]:
                        exp_dat_file.write(
                            str(i.resnum)
                            + "_"
                            + str(i.atom_name)
                            + "\t"
                            + str(i.value)
                            + "\t0.1\n"
                        )
                except KeyError:
                    continue

        with open(my_path + bme_calc_filename, "w") as calc_dat_file:
            for index, model in enumerate(model_data):
                calc_dat_file.write(str(index) + " ")

                for atom_type in chemshift_types:
                    try:
                        for i in cs_list[atom_type]:
                            calc_dat_file.write(
                                str(model[atom_type][i.resnum]) + " "
                            )
                    except KeyError:
                        continue

                calc_dat_file.write("\n")

        for corrected in [False, True]:
            for CS_type in sorted(list(cs_list.keys())):
                calc_dict = {}

                if corrected:
                    cs_list[CS_type + "_corrected"] = []

                for i, record in enumerate(cs_list[CS_type]):
                    if corrected:
                        prev_record = None
                        if i != 0:
                            prev_record = cs_list[CS_type][i - 1]

                        try:
                            next_record = cs_list[CS_type][i + 1]
                        except IndexError:
                            next_record = None

                        calc_dict[record.resnum] = (
                            cs_calced[CS_type][record.resnum] -
                            chemshift_corrections[record.res_name][record.atom_name]
                        )

                        corrected_record = ChemShiftRecord(
                            record.resnum, record.res_name, record.atom_name, record.value
                        )

                        corrected_record.value = (
                                record.value -
                                chemshift_corrections[record.res_name][record.atom_name]
                        )

                        if prev_record:
                            calc_dict[record.resnum] -= chemshift_corrections_prev[prev_record.res_name][prev_record.atom_name]
                            corrected_record.value -= chemshift_corrections_prev[prev_record.res_name][prev_record.atom_name]

                        if next_record:
                            calc_dict[record.resnum] -= chemshift_corrections_next[next_record.res_name][next_record.atom_name]
                            corrected_record.value -= chemshift_corrections_next[next_record.res_name][next_record.atom_name]

                        cs_list[CS_type + "_corrected"].append(corrected_record)
                    else:
                        calc_dict[record.resnum] = cs_calced[CS_type][record.resnum]

                # if corrected:
                #     CS_type = CS_type + "_corrected"

                model_corrs = []

                for model_num, model in enumerate(model_data):
                    inner_exp = {}

                    for i, record in enumerate(cs_list[CS_type]):
                        if corrected:
                            prev_record = None
                            if i != 0:
                                prev_record = cs_list[CS_type][i - 1]

                            try:
                                next_record = cs_list[CS_type][i + 1]
                            except IndexError:
                                next_record = None

                            inner_exp[record.resnum] = model[CS_type][record.resnum] - chemshift_corrections[record.res_name][record.atom_name]

                            if prev_record:
                                inner_exp[record.resnum] -= chemshift_corrections_prev[prev_record.res_name][prev_record.atom_name]

                            if next_record:
                                inner_exp[record.resnum] -= chemshift_corrections_next[next_record.res_name][next_record.atom_name]

                            # TODO _corrected types have less values as normal types
                            model_data[model_num][CS_type + "_corrected"] = inner_exp
                        else:
                            inner_exp[record.resnum] = model[CS_type][record.resnum]

                    if corrected:
                        model_corrs.append(
                            correlation(inner_exp, cs_list[CS_type + "_corrected"])
                        )
                    else:
                        model_corrs.append(
                            correlation(inner_exp, cs_list[CS_type])
                        )

                avg_model_corr = sum(model_corrs) / len(model_corrs)

                # REMOVE ME
                if corrected:
                    CS_type = CS_type + "_corrected"

                my_correl = correlation(calc_dict, cs_list[CS_type])
                my_q_value = q_value(calc_dict, cs_list[CS_type])
                my_rmsd = rmsd(calc_dict, cs_list[CS_type])

                corr_key = "CS_" + CS_type + "_corr"
                qval_key = "CS_" + CS_type + "_qval"
                rmsd_key = "CS_" + CS_type + "_rmsd"

                calced_data_storage.update(
                    {
                        corr_key: "{0}".format("{0:.3f}".format(my_correl)),
                        qval_key: "{0}".format("{0:.3f}".format(my_q_value)),
                        rmsd_key: "{0}".format("{0:.3f}".format(my_rmsd)),
                    }
                )

                csv_buffer.add_data(
                    {
                        "name": "ChemShifts (" + CS_type + ")",
                        "calced": calc_dict,
                        "experimental": cs_list[CS_type],
                    }
                )

                print("CHEM SHIFT", CS_type)
                print("Correl: ", my_correl)
                print("Q-val:  ", my_q_value)
                print("RMSD:   ", my_rmsd)
                print()

                graph_name = str(n + 1) + "_CS_" + CS_type + ".svg"
                graph.values_graph(my_path, calc_dict, cs_list[CS_type], graph_name)

                corr_graph_name = str(n + 1) + "_CS_corr_" + CS_type + ".svg"
                graph.correl_graph(
                    my_path, calc_dict, cs_list[CS_type], corr_graph_name
                )

                mod_corr_graph_name = "CS_mod_corr_" + CS_type + ".svg"
                graph.mod_correl_graph(
                    my_path,
                    my_correl,
                    avg_model_corr,
                    model_corrs,
                    mod_corr_graph_name,
                )

                my_id = my_path.split("/")[-2] + "/"
                CS_type_name = CS_type

                if corrected:
                    CS_type_name = CS_type + " (corrected)"

                cs_data.append(
                    {
                        "CS_type": CS_type_name,
                        "CS_model_n": len(cs_list[CS_type]),
                        "correlation": "{0:.3f}".format(my_correl),
                        "q_value": "{0:.3f}".format(my_q_value),
                        "rmsd": "{0:.3f}".format(my_rmsd),
                        "corr_graph_name": my_id + corr_graph_name,
                        "graph_name": my_id + graph_name,
                        "mod_corr_graph_name": my_id + mod_corr_graph_name,
                        "input_id": "CS_" + CS_type,
                    }
                )

        # TODO move this to a place where the corrected values exist
        cs_model_data_path = my_path + "/ChemShift_model_data.pickle"
        pickle.dump(model_data, open(cs_model_data_path, "wb"))

    chem_shift_lists_path = my_path + "/ChemShift_lists.pickle"
    pickle.dump(chem_shift_lists, open(chem_shift_lists_path, "wb"))
    return cs_data
