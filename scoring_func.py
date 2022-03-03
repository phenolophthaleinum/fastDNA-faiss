import statistics as stat


def avg(df):
    df_out = df.apply(lambda x: stat.mean(x))
    return df_out


def harmonic(df):
    df_out = df.apply(lambda x: stat.harmonic_mean(x))
    return df_out


def choose_max(df):
    df_out = df.apply(lambda x: max(x))
    return df_out
