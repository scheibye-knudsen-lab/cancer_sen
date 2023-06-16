
import csv
import sys
import imageio

#
# script to identify tissue region for senescent nuclei.
# the location of each nucleus is checked using region mask images
#


# determine if each "centers" point is located in region_img mask
def scan_nuc(nuc_mask_key, centers, region_img, out_writer):
    if len(centers) == 0:
        print("No blobs, skipping...")
        return None

    for idx, xy in enumerate(centers):
        out_key = nuc_mask_key + "_{}_{}".format(xy[0], xy[1])
        x, y = xy
        ishere = region_img[y, x] > 0
        out_writer.writerow([out_key, ishere])


#
# Main
#
# nuc_sen_path: path of senescence score file, where each row includes a key
# region_path: base path where for regions tiles.  region tiles are 1024x1024 patches,
#   encoded with >0 where there are relevant features.
#
_, nuc_sen_path, region_path, img_name_pattern, region_ext, out_path = sys.argv

csvfile = open(nuc_sen_path)
csv_reader = csv.reader(csvfile)

# read nuclei sen keys, building map from image tiles to xy coordinates
nuc_map = {}
for row in csv_reader:
    key = row[0]

    # extract x,y coordinates from key
    p = key.rfind("_")
    y = int(key[p + 1:])
    key = key[0:p]
    p = key.rfind("_")
    x = int(key[p + 1:])
    key = key[0:p]

    # key has coords removed, now identifying tiles

    if key not in nuc_map:
        nuc_map[key] = []
    nuc_map[key].append([x, y])

print(list(nuc_map.keys())[0:10])


csv_filename = out_path

csvfile = open(csv_filename, "w")
out_writer = csv.writer(csvfile)

# iterate through tiles
for nuc_mask_fn, centers in nuc_map.items():

    # generate region path using nuclei keys
    region_key = nuc_mask_fn + "." + region_ext
    region_fn = region_path + "/" + region_key

    try:
        region_img = imageio.imread(region_fn)
    except FileNotFoundError as e:
        print("skipping due to region loading: ", region_fn)
        continue

    # map each nuclei center for this region
    scan_nuc(nuc_mask_fn, centers, region_img, out_writer)

csvfile.close()
