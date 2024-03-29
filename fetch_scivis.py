#!/usr/bin/env python3

import os
import sys
import requests
import json
from urllib.parse import urlparse
from argparse import ArgumentParser

sizeof = {'uint8' : 1, 'uint16' : 2, 'int16' : 2, 'float32' : 4, 'float64' : 8} # bytes

parser = ArgumentParser(description='Fetch a scientific visualization dataset')
group  = parser.add_mutually_exclusive_group()
group.add_argument('-l', '--list', dest='do_list', action='store_true', help='List available datasets')
group.add_argument('-d', '--dataset', type=str, help='Name of the dataset to fetch')
args = parser.parse_args()

r = requests.get("https://klacansky.com/open-scivis-datasets/datasets.json")
index = r.json()

if args.do_list:
  print('Available datasets:')
  for x in index:
    name = x['name'].lower().replace(' ', '_')
    dims = [int(d) for d in x['size']]
    size = dims[0] * dims[1] * dims[2] * sizeof[x['type']] / 1048576
    print(f'  {name}', end=' ')
    print(f'({size/1024:.0f} GB)' if size > 1024 else f'({size:.0f} MB)' if size > 1 else '(<1 MB)')
  sys.exit(0)

meta = [x for x in index if x["name"].lower().replace(' ', '_') == args.dataset][0]
print(json.dumps(meta, indent=4))

script_path = os.path.dirname(os.path.realpath(__file__))
work_path = os.getcwd()

volume_path = urlparse(meta["url"])
meta["volume"] = os.path.normpath(work_path + "/" + os.path.basename(volume_path.path))

meta["spp"] = 1
meta["camera"] = [
        {
            "pos": [meta["size"][0] / 2 * meta["spacing"][0],
                    meta["size"][1] / 2 * meta["spacing"][1],
                    -meta["size"][2]],
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

with open(args.dataset + "_cinema.json", "w") as f:
    f.write(json.dumps(meta, indent=4))

if not os.path.isfile(meta["volume"]):
    with open(os.path.basename(meta["volume"]), "wb") as f:
        print("Fetching volume from {}".format(meta["url"]))
        r = requests.get(meta["url"], stream=True)
        file_size = r.headers.get('content-length')

        if file_size is None:
            data = bytes(r.content)
            f.write(data)
        else:
            downloaded = 0
            file_size = int(file_size)
            for data in r.iter_content(chunk_size=4096):
                downloaded += len(data)
                f.write(bytes(data))
                progress = downloaded / file_size
                bar = int(50 * progress)
                sys.stdout.write(f'\r[{"#"*bar}{" "*(50-bar)}] {100*progress:.1f}%')
                sys.stdout.flush()
else:
    print("File {} already exists, not re-downloading".format(meta["volume"]))

