import subprocess

from consensx import thirdparty
import consensx.graph as graph
from consensx.selection import get_pdb_models, averageChi2_on

def saxs(calced_data_storage, my_path, saxs_data_file_path):
    pdb_models = get_pdb_models(my_path)
    pdb_model_nums = []

    for model_num, _ in enumerate(pdb_models):
        fit_file_path = f"{my_path}/model_{model_num + 1}.fit"
        pdb_file_path = f"{my_path}/model_{model_num + 1}.pdb"

        pdb_model_nums.append(model_num)

        try:
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
        except PermissionError:
            pass

    chi2 = averageChi2_on(my_path, pdb_model_nums, write_fit_file=True)

    calced_data_storage.update(
        {
            "saxs_chi2": "{0}".format("{0:.3f}".format(chi2)),
        }
    )

    graph.saxs_graph(my_path, f"{my_path}my_fit.fit")

    return chi2
