# Import packages
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA
from sklearn.manifold import LocallyLinearEmbedding as LLE
from sklearn.manifold import Isomap
from sklearn.manifold import MDS
from sklearn.decomposition import FastICA
from sklearn.decomposition import KernelPCA
from sklearn.decomposition import SparsePCA
from sklearn.decomposition import TruncatedSVD
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from yellowbrick.classifier import PrecisionRecallCurve
from sklearn.multiclass import OneVsRestClassifier
import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
from sklearn.model_selection import train_test_split

# Perform dimensionality reduction on training data
def dimensionality_reduce(data, label_data, n_dimension, method_key):

    df = pd.read_csv(data, index_col=0).T

    # TSNE
    if method_key == 1:
        tsne = TSNE(
            n_components=3,
            random_state=0,
            learning_rate=4500,
            perplexity=28,
            early_exaggeration=14,
        )
        X_np = tsne.fit_transform(df)

    # PCA
    elif method_key == 2:
        pca = PCA(n_components=n_dimension)
        X_np = pca.fit_transform(df)

    # LLE
    elif method_key == 3:
        embed_lle = LLE(
            n_neighbors=110,  # default=5, number of neighbors to consider for each point.
            n_components=n_dimension,  # default=2, number of dimensions of the new space
            reg=0.001,  # default=1e-3, regularization constant, multiplies the trace of the local covariance matrix of the distances.
            eigen_solver="auto",  # {‘auto’, ‘arpack’, ‘dense’}, default=’auto’, auto : algorithm will attempt to choose the best method for input data
            max_iter=100000,  # default=100, maximum number of iterations for the arpack solver. Not used if eigen_solver==’dense’.
            method="modified",  # {‘standard’, ‘hessian’, ‘modified’, ‘ltsa’}, default=’standard’
            modified_tol=1e-12,  # default=1e-12, Tolerance for modified LLE method. Only used if method == 'modified'
            neighbors_algorithm="auto",  # {‘auto’, ‘brute’, ‘kd_tree’, ‘ball_tree’}, default=’auto’, algorithm to use for nearest neighbors search, passed to neighbors.NearestNeighbors instance
            random_state=42,  # default=None, Determines the random number generator when eigen_solver == ‘arpack’. Pass an int for reproducible results across multiple function calls.
            n_jobs=-1,  # default=None, The number of parallel jobs to run. -1 means using all processors.
        )

        X_np = embed_lle.fit_transform(df)

    # Isomap
    elif method_key == 4:
        embed_isomap = Isomap(n_neighbors=70, n_components=n_dimension, n_jobs=-1)

        X_np = embed_isomap.fit_transform(df)

    # Multidimensional scaling (MDS)
    elif method_key == 5:
        mds = MDS(n_components=n_dimension, metric=True, random_state=2)
        X_np = mds.fit_transform(df)

    # ICA
    elif method_key == 6:
        ica_transformer = FastICA(
            n_components=n_dimension, random_state=0, whiten="unit-variance"
        )
        X_np = ica_transformer.fit_transform(df)

    # KernelPCA
    elif method_key == 7:
        kernel_pca = KernelPCA(n_components=n_dimension, kernel="linear")
        X_np = kernel_pca.fit_transform(df)

    # SparsePCA
    elif method_key == 8:
        sparse_pca = SparsePCA(n_components=n_dimension, random_state=0)
        X_np = sparse_pca.fit_transform(df)

    # LDA
    elif method_key == 9:
        lda = LDA(n_components=n_dimension)
        X_np = lda.fit_transform(df, get_label(label_data))
        
        # Generate loadings matrix (coefficients)
        loadings = pd.DataFrame(lda.coef_)
        loadings.to_csv('LDA_loadings.csv')

    # Truncated SVD
    else:
        svd = TruncatedSVD(n_components=n_dimension)
        X_np = svd.fit_transform(df)

    # Generate cell and PC index
    X = pd.DataFrame(X_np)
    cell = list(df.index.values)
    X.index = cell
    X.rename(columns=lambda x: "PC_" + str(x + 1), inplace=True)
    X.to_csv("reduced_data.csv")
    return X


# Get the class label from meta data
def get_label(path_to_meta_file):
    df_label = pd.read_csv(path_to_meta_file, index_col=0)
    y = np.ravel(df_label["cell_type"])
    return y


# Train model and output score
def run_model(X, y, test_size):
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=42
    )
    model = LogisticRegression()
    model.fit(X_train, y_train)
    return model.score(X_test, y_test)


# Run the entire pipeline
def run_pipeline(data, label, output_path, n, key):
    X = dimensionality_reduce(
        data=data, label_data=label, n_dimension=n, method_key=key
    )
    X.to_csv(output_path)
    y = get_label(label)
    score = run_model(X, y, 0.1)
    print("Accuracy:", score)
    print(X)
    print("----------------------------")


# Testing
# run_pipeline('scaled_training_data_sample_official.csv', 'scaled_training_label_sample_official.csv', 'Truncated_SVD_output.csv', 100, 10)

# Train model on one-vs-all classification and graph PRC curve
def run_model_prc(pc_file, label_file, graph_output):
    X = pd.read_csv(pc_file, index_col=0)
    df_label = pd.read_csv(label_file, index_col=0)
    y = np.ravel(df_label["cell_type"])
    from collections import Counter

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.1, random_state=42
    )
    print(Counter(y_train))
    print("-------------------------")
    print(Counter(y_test))
    model = LogisticRegression()
    model_ova = OneVsRestClassifier(model)
    model_ova.fit(X_train, y_train)

    prc_lg = PrecisionRecallCurve(
        model_ova,
        classes=model_ova.classes_,
        iso_f1_curves=True,
        per_class=True,
        micro=False,
        size=(1000, 800),
    )
    prc_lg.fit(X_train, y_train)
    prc_lg.score(X_test, y_test)
    prc_lg.show(outpath=graph_output)
    print(str(pc_file))
    print("Score :", model_ova.score(X_test, y_test))


# run_model_prc('LDA_output.csv', 'scaled_training_label_sample_official.csv', 'PRC_LDA.png')
