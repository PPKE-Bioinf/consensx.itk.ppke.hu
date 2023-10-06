import matplotlib.pyplot as plt
import numpy as np

plt.switch_backend('Agg')


x_values = []
y_exp = []
y_calc = []

with open("sah2_1_PHOS2_0035_sub_Drbn_sah_5_99p_secondmin_1778.fit") as fitfile:
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
plt.savefig("fit_graph.svg", format="svg", transparent=False)
plt.close()
