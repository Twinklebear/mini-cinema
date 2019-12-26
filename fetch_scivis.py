#!/usr/bin/env python3

import os
import sys
import requests
import json
from urllib.parse import urlparse

if len(sys.argv) != 2:
    print("Usage: ./fetch_scivis.py <dataset name>")
    print("The dataset name should be all lowercase and without spaces to match the website")
    sys.exit(1)

r = requests.get("http://sci.utah.edu/~klacansky/cdn/open-scivis-datasets/{}/{}.json".format(sys.argv[1], sys.argv[1]))

script_path = os.path.dirname(os.path.realpath(__file__))
work_path = os.getcwd()

meta = r.json()
volume_path = urlparse(meta["url"])
meta["volume"] = os.path.normpath(work_path + "/" + os.path.basename(volume_path.path))

meta["spp"] = 1
meta["camera"] = [
        {
            "pos": [meta["size"][0] / 2, meta["size"][1] / 2, -meta["size"][2]],
            "dir": [0, 0, 1],
            "up": [0, 1, 0]
        },
        {
            "orbit": 1
        }]
meta["colormap"] = [os.path.normpath(script_path + "/configs/paraview_cool_warm.png")]
meta["isosurface_color"] = [0.2, 0.5, 1.0]
meta["background_color"] = [0, 0, 0]
meta["image_size"] = [512, 512]

print("Data set Information")
print(json.dumps(meta, indent=4))

with open(sys.argv[1] + "_cinema.json", "w") as f:
    f.write(json.dumps(meta, indent=4))

if not os.path.isfile(meta["volume"]):
    print("Fetching volume from {}".format(meta["url"]))
    r = requests.get(meta["url"])
    data = bytes(r.content)
    with open(os.path.basename(meta["volume"]), "wb") as f:
        f.write(data)
else:
    print("File {} already exists, not re-downloading".format(meta["volume"]))

