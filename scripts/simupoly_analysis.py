import glob
import csv
import matplotlib.pyplot as plt
from codons import get_color

folder = "data_sfs"
folder_path = "/home/thibault/SimuEvol/{0}".format(folder)

for tsv_path in sorted(glob.glob("{0}/np*.tsv".format(folder_path))):
    with open(tsv_path, 'r') as tsvfile:
        tsvin = csv.reader(tsvfile, delimiter='\t')
        y_axis = ["Ps", "Pn", "HapNbr", "CpxSites"]
        y_values = [list() for _ in y_axis]
        x_values = [list() for _ in y_axis]
        for row in tsvin:
            if len(row) >= 2 and (row[1] in y_axis):
                index = y_axis.index(row[1])
                x_values[index].append(int(row[0]))
                y_values[index].append(int(row[2]))

        my_dpi = 196
        plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
        for index, y_label in enumerate(y_axis):
            plt.plot(x_values[index], y_values[index], color=get_color(index), linewidth=2, label=y_label)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)
        plt.xlabel('$t$', fontsize=24)
        plt.ylabel('Number of sites', fontsize=24)
        plt.legend()
        plt.tight_layout()
        plt.savefig("../figures/{0}".format(tsv_path.split("/")[-1].replace("tsv", "png")), format="png")
        plt.show()
        plt.clf()
        plt.close('all')
