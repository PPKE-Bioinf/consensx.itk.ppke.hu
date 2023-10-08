import os
import sys
import subprocess

from consensx import thirdparty
import consensx.graph as graph

def saxs(my_path, pdb_file_path, saxs_data_file_path):
    fit_file_path = f"{my_path}saxs.fit"

    subprocess.call(
        [
            thirdparty.ThirdParty.pepsiSAXS,
            pdb_file_path,
            saxs_data_file_path,
            "-o",
            fit_file_path
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT
    )

    chi2 = -1

    with open(fit_file_path) as fit_file:
        for line in fit_file:
            if "Chi2" in line:
                chi2 = float(line.split(':')[-1].strip())
                break

    graph.saxs_graph(my_path, fit_file_path)

    return chi2