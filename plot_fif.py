#!/usr/bin/env python3

import os
import sys
import re
import json
import matplotlib
import numpy as np
import scipy
from matplotlib import rc
from matplotlib import cm
from docopt import docopt
import matplotlib.pyplot as plt

doc = """Plot FiF

Usage:
    plot_fif.py <file>... [options]

Options:
    -o OUTPUT           Save the plot to an output file
"""
args = docopt(doc)

parse_log_fname = re.compile(".*mini-cinema-(\d+)n-([^0-9]+)-(\d+)\.txt")
match_config = re.compile("Rendering Config: (.*)")
match_fif = re.compile("FRAMES_IN_FLIGHT = (\d+)")
match_frame_time = re.compile("Frame (\d+) took (\d+)ms")
match_renders_done = re.compile("All renders completed in (\d+)ms")
match_all_done = re.compile("Full run.*completed in (\d+)ms")

scaling_runs = []
for filename in args["<file>"]:
    m_config = parse_log_fname.match(filename)
    if not m_config:
        print("Unrecognized filename pattern {}".format(filename))
        sys.exit(1)

    print(m_config.groups())
    with open(filename, "r") as f:
        fif = -1
        all_done = 0
        renders_done = 0
        for l in f:
            m = match_config.match(l)
            if m:
                config = json.loads(m.group(1))
                continue
            m = match_fif.search(l)
            if m:
                fif = int(m.group(1))
                continue
            m = match_renders_done.search(l)
            if m:
                renders_done = float(m.group(1)) / 1000.0
                continue
            m = match_all_done.search(l)
            if m:
                all_done = float(m.group(1)) / 1000.0
                continue
        if not config:
            print("[error]: Did not find config line for {}!".format(filename))
        elif fif == -1:
            print("[error]: Did not find FIF line for {}!".format(filename))
        else:
            print("fif: {}, renders: {}, total: {}".format(fif, renders_done, all_done))
            scaling_runs.append({
                "fif": fif,
                "render": renders_done,
                "all": all_done
            })

scaling_runs.sort(key=lambda r: r["fif"])

if args["-o"] and os.path.splitext(args["-o"])[1] == ".pdf":
    #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    ## for Palatino and other serif fonts use:
    rc('font',**{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True)

fig, ax = plt.subplots()
ax.set_xscale("log", basex=2, nonposx="clip")

x = [r["fif"] for r in scaling_runs]
y = [r["render"] for r in scaling_runs]
plt.plot(x, y, linewidth=2)

ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.yaxis.set_ticks_position("left")
ax.set_ylim(bottom=0)
ax.xaxis.set_ticks_position("bottom")
ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%d"))

ax.set_ylabel("Total Time (s)")
ax.set_xlabel("Frames in Flight")

if args["-o"]:
    if os.path.splitext(args["-o"])[1] == ".png":
        plt.savefig(args["-o"], dpi=150, bbox_inches="tight")
    else:
        plt.savefig(args["-o"], bbox_inches="tight")
    print("saved to {}".format(args["-o"]))
else:
    plt.show()

