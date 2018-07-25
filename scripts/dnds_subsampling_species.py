from ete3 import Tree
import os
from subprocess import run

cluster = False
nbr_cpu = 8

if cluster:
    current_dir = "/panhome/tlatrill/SimuEvol"
    pb_mpi_path = "/panhome/tlatrill/pbmpi2/data"

else:
    current_dir = "/home/thibault/SimuEvol"
    pb_mpi_path = "/home/thibault/Tools/pbmpi2/data"

protein = "np"

# Run SimuEvol for the 3 conditions
newick_path = "{0}/data_trees/{1}.newick".format(current_dir, protein)
prefs_path = "{0}/data_prefs/{1}.txt".format(current_dir, protein)
data_path = "{0}/data_pb".format(current_dir)

os.makedirs(current_dir + "/qsub", exist_ok=True)
os.makedirs(data_path, exist_ok=True)

for s, p in [(0.0, 0.0), (10.0, 0.0), (0.0, 1.0)]:
    output_path = "{0}/data_alignment/{1}_{2}_{3}".format(current_dir, protein, s, p)
    ali_path = output_path + ".ali"

    simu_evol_cmd = current_dir + "/SimuEvol --preferences={0} --newick={1} --output={2} --mu={3} --lambda={4}"
    simu_evol_cmd = simu_evol_cmd.format(prefs_path, newick_path, output_path, 2.5, 1.0)
    simu_evol_cmd += " --s={0} --p={1} --a=False".format(s, p)

    print("Running SimuEvol")
    print(simu_evol_cmd)
    run(simu_evol_cmd, shell=True)
    print("Finished running SimuEvol")

    t = Tree(newick_path)
    min_n = 125

    for index, node in enumerate(t.traverse("levelorder")):
        n = len(node.get_leaves())
        if n > min_n:
            print(index)
            tree_path = "{0}/data_trees/{1}_{2}.newick".format(current_dir, protein, index)
            node.write(format=1, outfile=tree_path)

            # Create and submit .pbs for the subsampled trees
            for option, opt_name in [("-siteomega", "siteomega"),
                                     ("-mutsel -dp", "mutsel")]:
                for chain in ["1", "2"]:
                    qsub_id = "_".join([protein, index, opt_name, chain])

                    qsub_path = "{0}/qsub/{1}.pbs".format(current_dir, qsub_id)

                    qsub_str = "#!/bin/bash\n"
                    qsub_str += "#\n"
                    qsub_str += "#PBS -q q1day\n"
                    qsub_str += "#PBS -l nodes=1:ppn={0},mem=4gb\n".format(nbr_cpu)
                    qsub_str += "#PBS -o /pandata/tlatrill/out_err/out_{0}\n".format(qsub_id)
                    qsub_str += "#PBS -e /pandata/tlatrill/out_err/err_{0}\n".format(qsub_id)
                    qsub_str += "#PBS -j oe\n"
                    qsub_str += "#PBS -W umask=022\n"
                    qsub_str += "#PBS -r n\n"
                    qsub_str += "TMP=/tmp/tlatrill$RANDOM\n"
                    qsub_str += "export TMPDIR=$TMP\n"

                    qsub_str += "mpirun -n {0} {1}/pb_mpi -f -s -x 1 400 {2}".format(nbr_cpu, pb_mpi_path, option)
                    qsub_str += " -d {0}".format(ali_path)
                    qsub_str += " -T {0}".format(tree_path)
                    qsub_str += " {0}/{1}\n".format(data_path, qsub_id)

                    if opt_name == "mutsel":
                        qsub_str += "{0}/readpb_mpi -x 100 -ss {1}/{2}\n".format(pb_mpi_path, data_path, qsub_id)
                        qsub_str += "{0}/readpb_mpi -x 100 -om {1}/{2}\n".format(pb_mpi_path, data_path, qsub_id)

                    qsub_str += "rm -rf $TMP\n"
                    qsub_str += "rm {0}\n".format(qsub_path)

                    qsub = open(qsub_path, 'w')
                    qsub.write(qsub_str)
                    qsub.close()

                    if cluster:
                        print("Submitting " + qsub_path + " to the cluster")
                        run("qsub {0}".format(qsub_path), shell=True)
                    else:
                        print("Running " + qsub_path)
                        run("sh {0}".format(qsub_path), shell=True)

print('Job completed')
