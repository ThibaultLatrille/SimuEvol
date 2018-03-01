import numpy as np
import os
from subprocess import run

current_dir = "/panhome/tlatrill/SimuEvol"

nbr_sites = 300
protein = "lactamase"
mixture = True
alpha = 0.1
chain = 1
id_prefix = "{0}_{1}_{2}_{3}_{4}".format(nbr_sites, protein, mixture, alpha, chain)

nbr_points = 20

data_path = "{0}/data_hyphy/{1}".format(current_dir, id_prefix)
os.makedirs(data_path, exist_ok=True)

os.makedirs("{0}/qsub".format(current_dir), exist_ok=True)

newick_path = "{0}/data_trees/{1}.newick".format(current_dir, protein)

prefs_path = "{0}/data_prefs/{1}.txt".format(current_dir, id_prefix)
preferences = np.random.dirichlet(alpha*np.ones(20))
prefs_file = open(prefs_path, 'w')
prefs_file.write("# POSITION WT SITE_ENTROPY PI_A PI_C PI_D PI_E PI_F PI_G PI_H PI_I PI_K PI_L PI_M PI_N PI_P PI_Q PI_R PI_S PI_T PI_V PI_W PI_Y\n")
for i in range(1, nbr_sites + 1):
    if mixture:
        preferences = np.random.dirichlet(alpha * np.ones(20))
    prefs_file.write("{0} A {0} ".format(i, alpha) + " ".join([str(i) for i in preferences]) + "\n")
prefs_file.close()

for qsub_id, mut_bias in enumerate(np.logspace(-1, 1, nbr_points)):
    id_sufix = "id{0}_m{1:.5f}".format(qsub_id, mut_bias)
    qsub_name = "{0}_{1}".format(id_prefix, id_sufix)

    # Run SimuEvol with the given mutional bias, and amino-acid preferences file and the tree
    ali_path = "{0}/data_alignment/{1}".format(current_dir, qsub_name)
    simu_evol_result = "{0}.txt".format(ali_path)
    mu = 10 / (1 + mut_bias)
    simu_evol_cmd = current_dir + "/SimuEvol --preferences={0} --newick={1} --output={2} --mu={3} --lambda={4}"
    print("Running SimuEvol with mu={0} and lambda={1}".format(mu, mut_bias))
    run(simu_evol_cmd.format(prefs_path, newick_path, ali_path, mu, mut_bias), shell=True)
    print("SimuEvol completed !")
    
    for model, nbr_cpu in [("hyphy", 8), ("bayescode", 1)]:
        qsub_path = "{0}/qsub/{1}_{2}.pbs".format(current_dir, qsub_name, model)

        qsub_str = "#!/bin/bash\n"
        qsub_str += "#\n"
        qsub_str += "#PBS -q q1day\n"
        qsub_str += "#PBS -l nodes=1:ppn={0},mem=8gb\n".format(nbr_cpu)
        qsub_str += "#PBS -o /pandata/tlatrill/read/read_out_{0}_{1}\n".format(qsub_name, model)
        qsub_str += "#PBS -e /pandata/tlatrill/read/read_err_{0}_{1}\n".format(qsub_name, model)
        qsub_str += "#PBS -j oe\n"
        qsub_str += "#PBS -W umask=022\n"
        qsub_str += "#PBS -r n\n"
        qsub_str += "#PBS -r n\n"
        qsub_str += "TMP=/tmp/tlatrill$RANDOM\n"
        qsub_str += "export TMPDIR=$TMP\n"

        if model == "hyphy":
            # Run SimuEvol with the given mutional bias, and amino-acid preferences file and the tree
            fasta_path = ali_path + ".fasta"
            scripts_dir = "{0}/scripts".format(current_dir)
            hyphy_batch_path = "{0}/{1}.bf".format(data_path, id_sufix)
            batchfile_cmd = "python3 {0}/hyphy_jinja.py -d {0} -b {1} -f {2} -t {3} -p '0-0_1_3-1_4' \n"

            qsub_str += batchfile_cmd.format(scripts_dir, hyphy_batch_path, fasta_path, newick_path)
            qsub_str += "HYPHYMP {0} CPU={1}\n".format(hyphy_batch_path, nbr_cpu)
        else:
            bayescode_path = "{0}/data_bayescode/{1}".format(current_dir, qsub_name)
            bayescode_cmd = "globom -d {0} -t {1} -x 1 400 {2}\n".format(ali_path, newick_path, bayescode_path)

            qsub_str += bayescode_cmd

        qsub_str += "rm -f {0}\n".format(qsub_path)
        
        qsub = open(qsub_path, 'w')
        qsub.write(qsub_str)
        qsub.close()

        run("qsub {0}".format(qsub_path), shell=True)
        print(qsub_path)
