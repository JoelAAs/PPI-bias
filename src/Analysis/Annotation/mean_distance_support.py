import numpy as np
import pandas as pd


def get_cumulative_sum(df, value_column, cumulative_columns, greater=True, min_samples=200):
    """
    Get mean of value column above/below limit column
    :param df: (pandas Dataframe) Dataframe with sum column and limit colum
    :param value_column: (string) column name of value columns
    :param cumulative_columns: (list) column names of columns mean value above/below value column
    :param greater: (bool) is greater slide descend and sum above among cumulative columns, vice versa otherwise
    :param min_samples: (integer) min observations before starting recording mean
    :return: (pandas DataFrame) containing sliding mean of each cumulative columns above/below value column
    """
    bins = df[value_column].unique()
    bins.sort()

    values = df[value_column].values
    idx_val = values.argsort()
    if greater:
        idx_val = idx_val[::-1]
        bins = bins[::-1]
    values = values[idx_val]


    if not greater:
        bins = -bins
        values = -values

    measurement_matrix = df[cumulative_columns].iloc[idx_val]
    na_matrix = measurement_matrix.isna().to_numpy()
    measurement_matrix = measurement_matrix.to_numpy()
    cumulative_sum = np.zeros(len(cumulative_columns),dtype=np.float64)
    cumulative_na = np.zeros(len(cumulative_columns),dtype=int)

    previous = 0
    i = 0
    j = 0
    rows = [{}] * len(bins)
    for threshold in bins:
        while (i < len(values) and threshold <= values[i]) or i < min_samples:
            i += 1
        if previous != i:
            cumulative_sum += np.nansum(measurement_matrix[previous:i], axis=0, dtype=np.float64)
            cumulative_na += np.nansum(na_matrix[previous:i], axis=0, dtype=int)
            previous = i
            rows[j] = {
                "limit": value_column,
                **{
                    "value": (threshold if greater else -threshold)
                 },
                **{
                    f"sum_{c}": v for c,v in zip(cumulative_columns, cumulative_sum)
                },
                **{
                    f"non_na_pairs_{c}": i-v for c, v in zip(cumulative_columns, cumulative_na)
                }
            }
            j += 1
    return pd.DataFrame(rows).dropna()