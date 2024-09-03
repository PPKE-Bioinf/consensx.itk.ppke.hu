import matplotlib.pyplot as plt
import numpy as np

plt.switch_backend('Agg')

def saxs_graph(my_path, fit_file_path):
    x_values = []
    y_exp = []
    y_calc = []

    with open(fit_file_path) as fitfile:
        for line in fitfile:

            if line[0] == '#':
                continue

            line_list = line.split()
            x_values.append(float(line_list[0]))
            y_exp.append(float(line_list[1]))
            y_calc.append(float(line_list[3]))

    plt.figure(figsize=(10, 5), dpi=80)
    plt.yscale("log") 

    plt.plot(
        x_values, y_exp,
        color='#FD6C6C', marker='x', label='experimental', alpha=.7
    )

    plt.plot(
        x_values, y_calc,
        color='#027A8B', marker='x', label='pepsi-SAXS', alpha=.7
    )

    plt.legend(loc='upper right')
    plt.xlabel('s, 1/Å')
    plt.ylabel('I, relative')
    plt.tight_layout(pad=1.08)
    plt.savefig(f"{my_path}fit_graph.svg", format="svg", transparent=False)
    plt.close()
