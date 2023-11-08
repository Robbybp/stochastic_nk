import argparse 
import pandas as pd
import numpy as np
from pipe import select, where
from k_means_constrained import KMeansConstrained
import math

parser = argparse.ArgumentParser(description='Python script to cluster the buses')
parser.add_argument('-c','--case', help='case file name', default='CATS_buses.csv')
parser.add_argument('-n', '--num_clusters', help='number of clusters', type=int, default=25)
args = parser.parse_args()

data = pd.read_csv(args.case)
lat = list(list(data.columns) | where(lambda x: 'lat' in str.casefold(x)))[0]
lon = list(list(data.columns) | where(lambda x: 'lon' in str.casefold(x)))[0]
bus = list(list(data.columns) | where(lambda x: 'bus' in str.casefold(x)))[0]
features = data[[lat, lon]]
x = np.array(features)
count = len(data)

clf = KMeansConstrained(
    n_clusters = args.num_clusters, 
    size_min = math.ceil(count/args.num_clusters) - 5, 
    size_max = math.ceil(count/args.num_clusters) + 5,
    random_state = 0)

clf.fit_predict(x)
labels = clf.labels_
data['cluster'] = labels

case_name = args.case.strip('buses.csv')
data[[bus, 'cluster']].to_csv(case_name + 'cluster.csv', index = False, header = ['bus_id', 'cluster_id'])