import statistics as stat
import pandas as pd


def avg(df: pd.DataFrame):
    df_out = df.apply(lambda x: stat.mean(x))
    return df_out


def harmonic(df: pd.DataFrame):
    mean_p = df.apply(lambda x: stat.mean(x))
    max_p = df.apply(lambda x: max(x))
    df_out = ((mean_p.pow(-1) + max_p.pow(-1)) / 2) ** -1
    return df_out


def polynomial(df: pd.DataFrame, exponent: float = 2):
    df_out = df.apply(lambda x: x ** exponent).apply(lambda x: stat.mean(x))
    return df_out


def choose_max(df: pd.DataFrame):
    df_out = df.apply(lambda x: max(x))
    return df_out


def rank_adjusted_row(df: pd.DataFrame):
    rank_penalties = df.rank(ascending=False, axis=1) ** -1
    df_out = (df * rank_penalties).apply(lambda x: stat.mean(x))
    return df_out


def rank_adjusted_col(df: pd.DataFrame):
    rank_penalties = df.rank(ascending=False, axis=0) ** -1
    df_out = (df * rank_penalties).apply(lambda x: stat.mean(x))
    return df_out
