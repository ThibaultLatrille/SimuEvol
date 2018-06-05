import glob
import csv
import matplotlib.pyplot as plt

folder = "data_sfs"
folder_path = "/home/thibault/SimuEvol/{0}".format(folder)

for tsv_path in sorted(glob.glob("{0}/np*.tsv".format(folder_path))):
    with open(tsv_path, 'r') as tsvfile:
        tsvin = csv.reader(tsvfile, delimiter='\t')
        data = [row for row in tsvin if (len(row) > 3 and row[3].count(" ") == 0)]
        y_axis = list(set([row[2] for row in data]))
        print(y_axis)
        y_values = [list() for _ in y_axis]
        x_values = [list() for _ in y_axis]
        for row in data:
            if row[2] in y_axis:
                index = y_axis.index(row[2])
                x_values[index].append(int(row[1]))
                y_values[index].append(float(row[3]))

        for index, y_label in enumerate(y_axis):
            my_dpi = 196
            plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
            plt.scatter(x_values[index], y_values[index], label=y_label)
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            plt.xlabel('$t$', fontsize=24)
            plt.ylabel(y_label, fontsize=24)
            if min(y_values[index]) > 0:
                plt.yscale("log")
            plt.legend()
            plt.tight_layout()
            plt.savefig("../figures/{0}_{1}".format(tsv_path.split("/")[-1].replace("tsv", "png"), y_label), format="png")
            plt.show()
            plt.clf()
            plt.close('all')
