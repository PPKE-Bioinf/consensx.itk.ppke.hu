import math
import subprocess
import os
import matplotlib.pyplot as plt

from consensx import thirdparty


plt.switch_backend("Agg")


def nmr_pride(pdb_models, my_path, noe_restraints):
    """Calculate NMR-PRIDE score on given PDB models"""
    pwd = os.getcwd()
    os.chdir(my_path)

    # write model list text file
    model_list = open("model_list.txt", "w")

    for model in pdb_models:
        model_list.write(model + "\n")

    model_list.write("END\n")
    model_list.close()

    # write distance dict to text file
    restraints = noe_restraints.get_pride_restraints()
    pride_input = open("pride_input.txt", "w")

    pride_input.write("HEADER\n")

    prime_distances = list(restraints.keys())
    prime_distances.sort()

    for distance in prime_distances:
        pride_input.write(
            str(distance) + " " + str(restraints[distance]) + "\n"
        )

    pride_input.write("END\n")
    pride_input.close()

    # create binary database for PRIDE-NMR
    devnull = open(os.devnull, "w")
    hhdb_log = open("hhdb.log", "w")
    model_list = open("model_list.txt", "r")
    subprocess.call(
        [thirdparty.ThirdParty.prideDB, "-D", "HHDB"],  # model list
        stdin=model_list,
        stdout=devnull,
        stderr=hhdb_log,
    )

    hhdb_log.close()
    devnull.close()
    model_list.close()

    # run PRIDE-NMR
    devnull = open(os.devnull, "w")
    pride_input = open("pride_input.txt", "r")
    pride_output = open("pride_output.txt", "w")
    subprocess.call(
        [
            thirdparty.ThirdParty.prideNMR,
            "-D",
            "HHDB",
            "-d",
            str(56),
            "-b",
            str(len(pdb_models)),
            "-m",
            str(3),
        ],
        stdin=pride_input,
        stdout=pride_output,
        stderr=devnull,
    )

    pride_input.close()
    pride_output.close()
    devnull.close()

    pride_scores = {}
    pride_output = open("pride_output.txt", "r")
    for line in pride_output:
        if line.startswith("PRIDENMR:"):
            model_num = int(line.split()[-1])
            model_score = float(line.split()[1])
            pride_scores[model_num] = model_score

    scores = list(pride_scores.values())
    avg = sum(scores) * 1.0 / len(scores)
    variance = [(x - avg) ** 2 for x in scores]
    standard_deviation = math.sqrt(sum(variance) * 1.0 / len(variance))
    pride_data = []

    print("PRIDE-NMR calculation")
    print("MAX: ", max(pride_scores, key=pride_scores.get))
    pride_data.append(max(pride_scores, key=pride_scores.get))
    print("MIN: ", min(pride_scores, key=pride_scores.get))
    pride_data.append(min(pride_scores, key=pride_scores.get))
    print("AVG: ", avg)
    pride_data.append(avg)
    print("DEV: ", standard_deviation, "\n")
    pride_data.append(standard_deviation)

    os.chdir(pwd)
    make_pride_graph(my_path, scores, avg)

    return pride_data


def make_pride_graph(my_path, graph_data, avg_score):
    graph_data.sort()

    plt.figure(figsize=(6, 5), dpi=80)
    plt.plot(
        graph_data,
        linewidth=2.0,
        color="blue",
        label="Model scores",
        alpha=0.7,
    )
    plt.plot(
        list(range(0, len(graph_data))),
        [avg_score] * len(graph_data),
        linewidth=2.0,
        color="green",
        label="Average score",
        alpha=0.7,
    )
    plt.axis([-1, len(graph_data), 0, 1])
    plt.xlabel("models by score (worse to best)")
    plt.ylabel("PRIDE-NMR score")
    plt.title("PRIDE-NMR scores")
    plt.tight_layout()
    plt.legend(loc="lower left")
    ax = plt.axes()
    ax.yaxis.grid()
    plt.savefig(my_path + "/PRIDE-NMR_score.svg", format="svg")
    plt.close()
