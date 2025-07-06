import numpy as np
import array


def numpy_elem(df):
    # Convert Python built-in array.array to numpy array
    df = df.map(lambda x: np.array(x) if isinstance(x, array.array) else x)
    return df


def cat_elem(df):
    # Convert each element in the dataframe to a numpy array
    col_list = [np.concatenate(df[col].values) for col in df.columns]
    # Concatenate all columns
    array = np.stack(col_list, axis=1)
    return array

def zscore(X, axis=0):
    X = np.copy(X)
    X = X - np.nanmean(X, axis=axis)[None,:]
    X = X / np.nanstd(X, axis=axis)[None,:]
    X[np.where(np.isnan(X))] = 0
    return X