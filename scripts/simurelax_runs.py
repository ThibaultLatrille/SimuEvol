import os
from subprocess import run

cluster = False
nbr_cpu = 1

if cluster:
    current_dir = "/panhome/tlatrill/SimuEvol"
    cmd = "qsub"
else:
    current_dir = "/home/thibault/SimuEvol"
    cmd = "sh"

os.makedirs(current_dir, exist_ok=True)
os.makedirs(current_dir + "/qsub", exist_ok=True)
data_path = "{0}/data_relax".format(current_dir)
os.makedirs(data_path, exist_ok=True)


def runs(pop_size, k, mu, r, n, m, a, q):
    id_prefix = "{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}".format(pop_size, k, mu, r, n, m, a, q)
    qsub_path = "{0}/qsub/{1}.pbs".format(current_dir, id_prefix)

    qsub_str = "#!/bin/bash\n"
    qsub_str += "#\n"
    qsub_str += "#PBS -q q1day\n"
    qsub_str += "#PBS -l nodes=1:ppn={0},mem=1gb\n".format(nbr_cpu)
    qsub_str += "#PBS -o /pandata/tlatrill/read/read_out_{0}\n".format(id_prefix)
    qsub_str += "#PBS -e /pandata/tlatrill/read/read_err_{0}\n".format(id_prefix)
    qsub_str += "#PBS -j oe\n"
    qsub_str += "#PBS -W umask=022\n"
    qsub_str += "#PBS -r n\n"
    qsub_str += "#PBS -r n\n"

    simu_evol_cmd = current_dir
    simu_evol_cmd += "/SimuRelax --pop_size={0} --k={1} --mu={2} --r={3} --n={4} --m={5} --a={6} --q={7} --dir={8}"
    simu_evol_cmd = simu_evol_cmd.format(pop_size, k, mu, r, n, m, a, q, data_path)
    qsub_str += simu_evol_cmd

    qsub = open(qsub_path, 'w')
    qsub.write(qsub_str)
    qsub.close()

    print("Running " + qsub_path)
    run("{0} {1}".format(cmd, qsub_path), shell=True)
    print("Finished running")


pop_size = 100  # The population size
k = 23  # The number of pair of chromosomes
mu = 10  # The mean number of mutations
r = 1e-3  # The effect of a mutation
n = 3  # The phenotype complexity (dimension)
m = 3  # The pleiotropy of mutations (<=n)
a = 1  # The landscape flatness
q = 1  # The landscape flatness
for n in [3, 25, 200]:
    runs(pop_size, k, mu, r, n, m, a, q)
