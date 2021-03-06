#!python3
import os
import argparse
import pandas as pd
import numpy as np


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--protein', required=False, type=str, default='np.prefs',
                        dest="protein", metavar="<protein>", help="")
    parser.add_argument('-c', '--ncat', required=False, type=int, default=0,
                        dest="ncat", metavar="<ncat>", help="")
    parser.add_argument('-n', '--nsite', required=False, type=int, default=500,
                        dest="nsite", metavar="<nsite>", help="")
    args = parser.parse_args()
    file = os.getcwd() + "/" + args.protein
    preferences = pd.read_csv(file, sep=",")
    if args.ncat > len(preferences) or args.ncat == 0:
        args.ncat = len(preferences)
    name = file.replace(".prefs", "") + "_{0}_{1}.prefs".format(args.nsite, args.ncat)
    subset = preferences.sample(n=args.ncat).values
    df = pd.DataFrame(subset[np.random.randint(subset.shape[0], size=args.nsite), :], columns=preferences.columns)
    df["site"] = range(1, args.nsite + 1)
    df.to_csv(name, sep=',', encoding='utf-8', index=False)
