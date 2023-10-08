import json
import pandas as pd
from alive_progress import alive_bar
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min

def create_clusters(store_path: str, function_names: list, num_clusters: int = 30):
    """
    Create clusters for each function's data and store representative indices in a JSON file.

    Args:
        store_path (str): The path where data files are stored.
        function_names (list): List of function names.
        num_clusters (int): Number of clusters to create (default is 30).

    Returns:
        None
    """
    # Dictionary to store representative indices for each function
    centroids = {}
    
    # Initialize progress bar
    with alive_bar(len(function_names)) as bar:
        for name in function_names:
            # Read data from Excel file
            dataset = pd.read_excel(f'{store_path}/{name}.xlsx')

            # Max Normalization
            max_norm_dataset = dataset.copy(deep=True)
            for column in max_norm_dataset.columns:
                max_norm_dataset[column] = max_norm_dataset[column] / max_norm_dataset[column].max()

            data = max_norm_dataset.copy(deep=True)
            data = data.dropna(axis=1)

            # Perform K-Means clustering
            kmeans = KMeans(n_clusters=num_clusters, init='k-means++', n_init='auto').fit(data[data.columns])
            y_predict = kmeans.predict(data[data.columns])

            data_clustered = data.copy(deep=True)
            data_clustered['Cluster'] = y_predict

            # Find closest data points to cluster centers
            closest, _ = pairwise_distances_argmin_min(kmeans.cluster_centers_, data[data.columns])
            closest.sort()
            
            # Store representative indices for the function
            centroids[name] = closest
            bar()

    # Save the representative indices to a JSON file
    with open(f'{store_path}/representative_calls.json', 'w') as outfile:
        json.dump(centroids, outfile)
