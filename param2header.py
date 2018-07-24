#!/usr/bin/env python3


import os
import sys


if __name__ == "__main__":

    filter_prefix = ""
    try:
        filter_prefix = sys.argv[2]
    except:
        pass


    lines = []
    with open(sys.argv[1], "r") as pf:
        _lines = pf.readlines()
        lines = [l.strip() for l in _lines]

    params = []
    for line in lines:

        if  (len(line) <= 0 or line.startswith("#")):
            # empty line
            continue

        param = line.split()[0]
        # print("%20s <== %s" % (param,line))
        if (param.find("(") > 0):
            # this line needs to be repeated as it contains an array
            param = line.split("(")[0].strip()
            n_times = int(line.split("(")[1].split(")")[0].strip())
            for i in range(n_times):
                params.append(param)
        else:
            params.append(param)

    # add all filter prefixes
    params = ["%s%s" % (filter_prefix, p) for p in params]

    print("\n".join(params))
