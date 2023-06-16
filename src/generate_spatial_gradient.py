
MODEL = "ir" # ir, rs, aad, uni

SEN_PATH = f".../coded-sen-{MODEL}-nuc.csv"
OUT_FILE = f".../spatial-grad-by30-{MODEL}.csv"

# reduced tile size
TILE_H, TILE_W = 128, 128

import re
import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix
import csv

# distance buckets, from 0 to 600 pixels, in steps of 30
NEAR_DISTANCES = np.arange(0, 601, 30)[1:]

# data coded by slide, code => [[x,y,tissue,sen]]
coded_data = {}

tissues = []

# for a given nucleus key, find all other nuclei in each distance bucket,
# and calculate metrics on it, including mean, rank, norm
def analyze_nucleus(key, tissue, score, score_rank, score_norm, dmr,
                    scorer, scorer_ranks, scorer_norms, log_writer):
    last_dist = 0
    for dist_idx, dist in enumerate(NEAR_DISTANCES):
        # same tissue and in dist range
        qsel = (dmr > last_dist) & (dmr <= dist)
        last_dist = dist
        mscore = np.mean(scorer[qsel])
        mscore_rank = np.mean(scorer_ranks[qsel])
        mscore_norm = np.mean(scorer_norms[qsel])
        count = len(scorer[qsel])

        log_writer.writerow([key, tissue, score, score_rank, score_norm,
                            dist_idx, count, mscore, mscore_rank, mscore_norm])


# analyze all nuclei per individual
def analyze_individual(points, log_writer):
    pts, tis, scores = [], [], []
    for point in points:
        pts.append(point[0:2])
        tis.append(point[2])
        scores.append(point[3])

    tis = np.asarray(tis)  # so can eval tis==tissue for each elem
    scores = np.asarray(scores)

    sdf = pd.Series(scores)
    if sdf.max() - sdf.min() == 0:
        return

    scores_rank = sdf.rank(pct=True)
    scores_norm = sdf.apply(lambda x: (
        x - sdf.min()) / (sdf.max() - sdf.min()))
    scores_rank = scores_rank.to_numpy()
    scores_norm = scores_norm.to_numpy()

    dm = distance_matrix(pts, pts)

    # for each nucleus
    for pt_idx, point in enumerate(points):
        tissue = point[2]
        score = point[3]
        key = point[4]
        dmr = dm[pt_idx, :]

        # same cell type
        tsel = (tis == tissue)
        scorer = scores[tsel]
        scorer_ranks = scores_rank[tsel]
        scorer_norms = scores_norm[tsel]
        dmr = dmr[tsel]
        analyze_nucleus(key, tissue, score, scores_rank[pt_idx], scores_norm[pt_idx],
                        dmr, scorer, scorer_ranks, scorer_norms, log_writer)


#
# Main
#
# load all nuclei with sen scores, determine tile row, col and nucleus x,y.
# group by "code" slide, and then analyze spatial per individual
#

# load sen score file
pdf1 = pd.read_csv(SEN_PATH)
for idx, row in pdf1.iterrows():
    key, code, tissue = row['key'], row['code'], row['tissue']
    sen = float(row['sen'])

    p = key.rfind("_")
    p = key.rfind("_", 0, p - 1)
    slide = key[0:p]

    m = re.search("_(\d+)x(\d+)_(\d+)_(\d+)$", key)
    row, col, x, y = int(m.group(1)), int(
        m.group(2)), int(m.group(3)), int(m.group(4))

    if code not in coded_data:
        coded_data[code] = []

    if tissue not in tissues:
        tissues.append(tissue)

    oy = row * TILE_H + round(y * TILE_H / 1024)
    ox = col * TILE_W + round(x * TILE_W / 1024)

    coded_data[code].append([ox, oy, tissue, sen, key])


csvfile = open(OUT_FILE, "w")
log_writer = csv.writer(csvfile)

for code in coded_data:
    points = coded_data[code]
    analyze_individual(points, log_writer)

csvfile.close()
