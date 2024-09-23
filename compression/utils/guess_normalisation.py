import anndata
import numpy as np


def robust_z_score(data, number):
    """
    Calculate the robust z-score for a given number relative to a dataset.
    
    Parameters:
    data (array-like): The dataset to calculate the z-score from.
    number (float or int): The number to calculate the z-score for.

    Returns:
    float: The robust z-score, which is based on the median and median absolute deviation (MAD).
    
    This function uses the median and MAD (as opposed to mean and standard deviation) to handle 
    datasets that may contain outliers, making it more robust for data with irregular distributions.
    """
    data = np.array(data)
    median = np.median(data)
    mad = np.median(np.abs(data - median))
    if mad == 0:
        return 0
    robust_z = (number - median) / mad
    return robust_z


def guess_unit_and_log(adata):
    """
    Guess the unit and log-transformation status of single-cell data based on the matrix values.
    
    Parameters:
    adata (AnnData): The AnnData object containing the single-cell dataset.
    
    Returns:
    tuple: A tuple where the first value is the guessed unit ('raw', 'cptt', 'cpm') and the second value 
           is a boolean indicating whether the data is log-transformed or not.
    
    This function tries to determine if the data is untransformed ('raw'), normalized to counts per ten thousand ('cptt'),
    or normalized to counts per million ('cpm'). It also checks if the data is log-transformed.
    """
    log_transformed = False
    threshold = 1
    X = adata.X
    is_integer = (np.floor(X.data[:30]) == np.ceil(X.data[:30])).all()
    if is_integer:
        return "raw", log_transformed
    row_sums = [i[0, 0] for i in X.sum(axis=1)]
    if abs(robust_z_score(row_sums, 10000)) < threshold:
        return "cptt", log_transformed
    if abs(robust_z_score(row_sums, 1000000)) < threshold:
        return "cpm", log_transformed
    
    Xe = np.expm1(adata.X)
    log_transformed = True
    is_integer = (np.floor(Xe.data[:30]) == np.ceil(Xe.data[:30])).all()
    if is_integer:
        return "raw", log_transformed
    row_sums = [i[0, 0] for i in Xe.sum(axis=1)]
    if abs(robust_z_score(row_sums, 10000)) < threshold:
        return "cptt", log_transformed
    if abs(robust_z_score(row_sums, 1000000)) < threshold:
        return "cpm", log_transformed
    
    return None, False


def guess_unit_and_log_relaxed_method(adata):
    """
    A relaxed method to guess the normalization and log-transformation status of single-cell data.
    
    Parameters:
    adata (AnnData): The AnnData object containing the single-cell dataset.
    
    Returns:
    tuple: A tuple where the first value is the guessed unit ('raw', 'cptt', 'cpm') and the second value 
           is a boolean indicating whether the data is log-transformed or not.
    
    This method takes a more relaxed approach by analyzing a sample of 100 cells (or fewer if there are less than 100),
    and checks for common normalization ranges, such as counts per ten thousand ('cptt') or counts per million ('cpm').
    """
    
    if adata.n_obs == 0:
        raise ValueError(
            "Single cell data set must have at least one observation/cell." ""
        )

    nsample = min(100, adata.n_obs)
    
    sum0 = int(adata.X[:nsample].sum() / nsample)
    if 9000 < sum0 < 11000:
        return "cptt", False
    if 900000 < sum0 < 1100000:
        return "cpm", False

    is_integer = (np.floor(adata.X.data[:30]) == np.ceil(adata.X.data[:30])).all()
    if is_integer:
        return "raw", False

    Xe = np.expm1(adata.X)
    adatae = anndata.AnnData(X=Xe)
    sum0 = int(adatae.X[:nsample].sum() / nsample)
    if 9000 < sum0 < 11000:
        return "cptt", True
    if 900000 < sum0 < 1100000:
        return "cpm", True

    is_integer = (np.floor(adata.X.data[:30]) == np.ceil(adata.X.data[:30])).all()
    if is_integer:
        return "raw", True
    
    return None, False