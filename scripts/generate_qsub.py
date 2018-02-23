import numpy as np
import os

current_dir = "/home/thibault/SimuEvol"

nbr_points = 20
alpha = 0.5
nbr_sites = 200

protein = "np"
data_path = "{0}/data_hyphy/{1}_{2}".format(current_dir, protein, alpha)
os.makedirs(data_path, exist_ok=True)

os.makedirs("{0}/qsub".format(current_dir), exist_ok=True)
os.makedirs("{0}/qsub_running".format(current_dir), exist_ok=True)

newick_path = "{0}/data_trees/{1}.newick".format(current_dir, protein)

prefs_path = "{0}/data_prefs/{1}_{2}.txt".format(current_dir, protein, alpha)
preferences = np.random.dirichlet(alpha*np.ones(20))
prefs_file = open(prefs_path, 'w')
prefs_file.write("# POSITION WT SITE_ENTROPY PI_A PI_C PI_D PI_E PI_F PI_G PI_H PI_I PI_K PI_L PI_M PI_N PI_P PI_Q PI_R PI_S PI_T PI_V PI_W PI_Y\n")
for i in range(1, nbr_sites + 1):
    prefs_file.write("{0} A {0} ".format(i, alpha) + " ".join([str(i) for i in preferences]) + "\n")
prefs_file.close()

for qsub_id, mut_bias in enumerate(np.logspace(-1, 1, nbr_points)):
    qsub_name = "{0}_{1}_id{2}_m{3:.3f}".format(protein, alpha, qsub_id, mut_bias)
    qsub_path = "{0}/qsub/{1}.pbs".format(current_dir, qsub_name)
    qsub = open(qsub_path, 'w')
    qsub.write("#!/bin/bash\n")
    qsub.write("#\n")
    qsub.write("#PBS -q q1day\n")
    qsub.write("#PBS -l nodes=1:ppn=8,mem=8gb\n")
    qsub.write("#PBS -o /pandata/tlatrill/read/read_out_{0}\n".format(qsub_name))
    qsub.write("#PBS -e /pandata/tlatrill/read/read_err_{0}\n".format(qsub_name))
    qsub.write("#PBS -j oe\n")
    qsub.write("#PBS -W umask=022\n")
    qsub.write("#PBS -r n\n")
    qsub.write("#PBS -r n\n")

    # Run SimuEvol with the given mutional bias, and amino-acid preferences file and the tree
    ali_path = "{0}/data_alignment/{1}".format(current_dir, qsub_name)
    simu_evol_result = "{0}.txt".format(ali_path)
    mu = 20 / (1 + mut_bias)
    simu_evol_cmd = "SimuEvol --preferences={0} --newick={1} --output={2} --mu={3} --lambda={4} > {2}.txt\n"
    qsub.write(simu_evol_cmd.format(prefs_path, newick_path, ali_path, mu, mut_bias))

    # Run SimuEvol with the given mutional bias, and amino-acid preferences file and the tree
    fasta_path = ali_path + ".fasta"
    scripts_dir = "{0}/scripts".format(current_dir)
    hyphy_batch_path = "{0}/id{1}_m{2:.3f}.bf".format(data_path, qsub_id, mut_bias)
    batchfile_cmd = "python3 {0}/hyphy_jinja.py {0} {1} {2} {3}\n"
    qsub.write(batchfile_cmd.format(scripts_dir, hyphy_batch_path, fasta_path, newick_path))

    qsub.write("HYPHYMP {0}\n".format(hyphy_batch_path))

    qsub.write("rm -f {0}\n".format(qsub_path))
    qsub.close()
    print(qsub_name)
