import os

nbr_cpu = 4
current_dir = "/panhome/tlatrill/SimuEvol"

qsub_id = 0
for file in os.listdir(current_dir + "/data_prefs"):
    prefix = file.strip().replace(".txt", "")
    if "500_" in prefix:
        prefs_path = "{0}/data_prefs/{1}".format(current_dir, file)
        param = prefix.split("_")
        protein = param[1]
        newick_path = "{0}/data_trees/{1}.newick".format(current_dir, protein)

        hyphy_path = "{0}/data_hyphy/{1}".format(current_dir, prefix)
        if os.path.isdir(hyphy_path):
            batchfiles = [batch[:-3] for batch in os.listdir(hyphy_path) if batch[-3:] == ".bf"]

            for batch in batchfiles:
                mut_bias = float(batch.split("_m")[-1])
                ali_path = "{0}/data_alignment/{1}_{2}.ali".format(current_dir, prefix, batch)

                id_sufix = "id{0}_m{1:.5f}".format(qsub_id, mut_bias)
                qsub_name = "{0}_{1}".format(prefix, id_sufix)
                qsub_path = "{0}/qsub/{1}.pbs".format(current_dir, qsub_name)
                qsub = open(qsub_path, 'w')
                qsub.write("#!/bin/bash\n")
                qsub.write("#\n")
                qsub.write("#PBS -q q1day\n")
                qsub.write("#PBS -l nodes=1:ppn={0},mem=8gb\n".format(nbr_cpu))
                qsub.write("#PBS -o /pandata/tlatrill/read/read_out_{0}\n".format(qsub_name))
                qsub.write("#PBS -e /pandata/tlatrill/read/read_err_{0}\n".format(qsub_name))
                qsub.write("#PBS -j oe\n")
                qsub.write("#PBS -W umask=022\n")
                qsub.write("#PBS -r n\n")
                qsub.write("#PBS -r n\n")

                # Run BayesCode with the given file
                bayescode_path = "{0}/data_bayescode/{1}".format(current_dir, qsub_name)
                bayescode_cmd = "mpirun -n {0} globom -d {1} -T {2} -x 1 400 {3}\n"
                qsub.write(bayescode_cmd.format(nbr_cpu, ali_path, newick_path, bayescode_path))
                qsub.write("rm -f {0}\n".format(qsub_path))

                qsub.close()
                print(qsub_name)
