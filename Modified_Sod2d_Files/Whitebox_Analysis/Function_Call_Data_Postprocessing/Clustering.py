import json
import pandas as pd
from alive_progress import alive_bar
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min

def create_clusters(store_path: str, function_names: list, num_clusters: int = 30):
    centroinds = {}
    with alive_bar(len(function_names)) as bar:
        for name in function_names:
            dataset = pd.read_csv(f"{store_path}/{name}.csv")

            # Min-Max Normalization
            # min_max_norm_dataset = dataset.copy(deep=True)
            # for column in min_max_norm_dataset.columns:
            #     min_max_norm_dataset[column] = (min_max_norm_dataset[column] - min_max_norm_dataset[column].min()) / (min_max_norm_dataset[column].max() - min_max_norm_dataset[column].min())    

            # Mean-Normalization
            # mean_norm_dataset = dataset.copy(deep=True)
            # for column in min_max_norm_dataset.columns:
            #     mean_norm_dataset[column] = (mean_norm_dataset[column]-mean_norm_dataset[column].mean())/mean_norm_dataset[column].std()

            # Z-Score Normalization
            # z_score_norm_dataset = dataset.copy(deep=True)
            # for column in z_score_norm_dataset.columns:
            #   z_score_norm_dataset[column] = (z_score_norm_dataset[column] - z_score_norm_dataset[column].mean()) / z_score_norm_dataset[column].std()    

            # Max Normalization
            max_norm_dataset = dataset.copy(deep=True)
            for column in max_norm_dataset.columns:
                max_norm_dataset[column] = max_norm_dataset[column] / max_norm_dataset[column].max()  

            data = max_norm_dataset.copy(deep=True)
            data = data.dropna(axis=1)

            n_clusters = 30
            kmeans = KMeans(n_clusters=n_clusters, init='k-means++', n_init='auto').fit(data[data.columns])
            y_predict = kmeans.predict(data[data.columns])  

            data_clustered = data.copy(deep=True)
            data_clustered['Cluster'] = y_predict

            closest, _ = pairwise_distances_argmin_min(kmeans.cluster_centers_, data[data.columns])
            closest.sort()
            centroinds[name] = closest
            bar()

    with open(f"{store_path}/representative_calls.json", "w") as outfile:
        json.dump(centroinds, outfile)