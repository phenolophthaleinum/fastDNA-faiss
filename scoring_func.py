import pandas as pd
import numpy as np


def meanmax_harmonic(df: pd.DataFrame):
    mean_p = df.apply(lambda x: np.mean(x))
    max_p = df.apply(lambda x: max(x))
    df_out = ((mean_p.pow(-1) + max_p.pow(-1)) / 2) ** -1
    return df_out


def meanmax_geometric(df: pd.DataFrame):
    mean_p = df.apply(lambda x: np.mean(x))
    max_p = df.apply(lambda x: max(x))
    df_out = (mean_p * max_p) ** 1 / 2
    return df_out


def meanmax_arithmetic(df: pd.DataFrame):
    mean_p = df.apply(lambda x: np.mean(x))
    max_p = df.apply(lambda x: max(x))
    df_out = (mean_p + max_p) / 2
    return df_out


def hit_polynomial(df: pd.DataFrame, exponent: float = 2):
    df_out = df.apply(lambda x: x ** exponent).apply(lambda x: np.mean(x))
    return df_out


def hit_max(df: pd.DataFrame):
    df_out = df.apply(lambda x: max(x))
    return df_out


def hit_amean(df: pd.DataFrame):
    df_out = df.apply(lambda x: np.mean(x))
    return df_out


def hit_hmean(df: pd.DataFrame):
    df_out = df.apply(lambda x: len(x) / np.sum(1.0 / x))
    return df_out


def hit_gmean(df: pd.DataFrame):
    df_out = df.apply(lambda x: x.prod() ** (1.0 / len(x)))
    return df_out


def hit_decile1(df: pd.DataFrame):
    df_out = df.apply(lambda x: np.quantile(x, 0.10))
    return df_out


def hit_decile2(df: pd.DataFrame):
    df_out = df.apply(lambda x: np.quantile(x, 0.20))
    return df_out


def hit_decile3(df: pd.DataFrame):
    df_out = df.apply(lambda x: np.quantile(x, 0.30))
    return df_out


def hit_decile4(df: pd.DataFrame):
    df_out = df.apply(lambda x: np.quantile(x, 0.40))
    return df_out


def hit_decile5(df: pd.DataFrame):
    df_out = df.apply(lambda x: np.quantile(x, 0.50))
    return df_out


def hit_decile6(df: pd.DataFrame):
    df_out = df.apply(lambda x: np.quantile(x, 0.60))
    return df_out


def hit_decile7(df: pd.DataFrame):
    df_out = df.apply(lambda x: np.quantile(x, 0.70))
    return df_out


def hit_decile8(df: pd.DataFrame):
    df_out = df.apply(lambda x: np.quantile(x, 0.80))
    return df_out


def hit_decile9(df: pd.DataFrame):
    df_out = df.apply(lambda x: np.quantile(x, 0.90))
    return df_out


def rank_median(df: pd.DataFrame):
    ranks = df.rank(ascending=False, axis=1)
    df_out = ranks.apply(lambda x: (len(x) - np.median(x)) / len(x))
    return df_out


def rank_mean(df: pd.DataFrame):
    ranks = df.rank(ascending=False, axis=1)
    df_out = ranks.apply(lambda x: (len(x) - np.mean(x)) / len(x))
    return df_out


def tomato_index(df: pd.DataFrame):
    max_vals = df.apply(lambda x: max(x), axis=1)
    filtered_fragments = df[df.max(axis=1) > np.median(max_vals)]
    ranks = filtered_fragments.rank(ascending=False, axis=1)
    df_out = ranks.apply(lambda x: (len(x) - np.mean(x)) / len(x))
    return df_out


def potato_index(df: pd.DataFrame):
    max_vals = df.apply(lambda x: max(x), axis=1)
    filtered_fragments = df[df.max(axis=1) < np.median(max_vals)]
    ranks = filtered_fragments.rank(ascending=False, axis=1)
    df_out = ranks.apply(lambda x: (len(x) - np.mean(x)) / len(x))
    return df_out


def potato_harmony(df: pd.DataFrame):
    max_vals = df.apply(lambda x: max(x), axis=1)
    filtered_fragments = df[df.max(axis=1) < np.median(max_vals)]
    df_out = filtered_fragments.apply(lambda x: len(x) / np.sum(1.0 / x))
    return df_out


def potato_mean(df: pd.DataFrame):
    max_vals = df.apply(lambda x: max(x), axis=1)
    filtered_fragments = df[df.max(axis=1) < np.median(max_vals)]
    df_out = filtered_fragments.apply(lambda x: np.mean(x))
    return df_out


def potato_geometry(df: pd.DataFrame):
    max_vals = df.apply(lambda x: max(x), axis=1)
    filtered_fragments = df[df.max(axis=1) < np.median(max_vals)]
    df_out = filtered_fragments.apply(lambda x: x.prod() ** (1.0 / len(x)))
    return df_out


def beet_harmony(df: pd.DataFrame):
    s = df.apply(lambda x: sum(x), axis=1)
    filtered_fragments = df[df.sum(axis=1) < np.median(s)]
    df_out = filtered_fragments.apply(lambda x: len(x) / np.sum(1.0 / x))
    return df_out


def beet_mean(df: pd.DataFrame):
    s = df.apply(lambda x: sum(x), axis=1)
    filtered_fragments = df[df.sum(axis=1) < np.median(s)]
    df_out = filtered_fragments.apply(lambda x: np.mean(x))
    return df_out


def beet_geometry(df: pd.DataFrame):
    s = df.apply(lambda x: sum(x), axis=1)
    filtered_fragments = df[df.sum(axis=1) < np.median(s)]
    df_out = filtered_fragments.apply(lambda x: x.prod() ** (1.0 / len(x)))
    return df_out


def hit_rank_adjusted_row(df: pd.DataFrame):
    rank_penalties = df.rank(ascending=False, axis=1) ** -1
    df_out = (df * rank_penalties).apply(lambda x: np.mean(x))
    return df_out


def hit_rank_adjusted_col(df: pd.DataFrame):
    rank_penalties = df.rank(ascending=False, axis=0) ** -1
    df_out = (df * rank_penalties).apply(lambda x: np.mean(x))
    return df_out
