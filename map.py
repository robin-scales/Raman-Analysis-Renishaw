from generalLoad import load
import scipy.interpolate
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from line_profile import LineProfile

def normal_interp(x, y, a, xi, yi):
    rbf = scipy.interpolate.Rbf(x, y, a)
    ai = rbf(xi, yi)
    return ai

def rescaled_interp(x, y, a, xi, yi):
    a_rescaled = (a - a.min()) / np.ptp(a)
    ai = normal_interp(x, y, a_rescaled, xi, yi)
    ai = np.ptp(a) * ai + a.min()
    return ai

def pMAP_plot(df):
        x = df["XList"]
        y = df["YList"]
        z = df["IList"]
        # print(x); print(y); print(z)

        xi, yi = np.mgrid[x.min():x.max():500j, y.min():y.max():500j]
        z_rescale = rescaled_interp(x, y, z, xi, yi)

        # plot
        fig, ax = plt.subplots()

        im = ax.imshow(z_rescale.T, origin='lower',
                        extent=[x.min(), x.max(), y.min(), y.max()]) # , vmin=566.3399156718697, vmax=566.5873836107614
        # ax.scatter(x, y, c=a)

        ax.scatter(x, y, s=15, c=z, cmap="grey") #  , vmin=566.3399156718697, vmax=566.5873836107614
        ax.set_ylim(ax.get_ylim()[::-1])

        # ax.set(xlim=(0, 8), xticks=np.arange(1, 8),
        #        ylim=(0, 8), yticks=np.arange(1, 8))

        print(np.min(z_rescale))
        print(np.max(z_rescale))

        cb = plt.colorbar(im)
        plt.show()


filepath = r"xxx"
print(f'Working on {filepath}')
df, scanType = load(filepath=filepath)

match scanType:
    case "uMAP":
        x = df["X"]
        y = df["Y"]
        wave = df["Wave"]
        intensity = df["Intensity"]
        unique_combinations = df[['X', 'Y']].drop_duplicates()
        print(unique_combinations)
        line_profiles = []
        IList = np.empty((len(unique_combinations.index),))
        IList[:] = np.nan
        print(IList)
        i = 0
        for index, row in unique_combinations.iterrows():
            valid_comb = (row['X'], row['Y'])
            df2 = df.set_index(['X', 'Y']).loc[valid_comb] # .reset_index()
            l = LineProfile(x=df2['Wave'].to_numpy(), y_raw=df2['Intensity'].to_numpy())
            l.fit()
            # line_profiles.append(l)
            df_l = l.fit_result['df_table']
            peak_wave = df_l.loc[df_l['Peak Index'] == 2, 'Center Grvty'].iloc[0]
            IList[i] = peak_wave
            i += 1
        XList = unique_combinations['X'].to_list()
        YList = unique_combinations['Y'].to_list()

        s1 = pd.Series(XList, name='XList')
        s2 = pd.Series(YList, name='YList')
        s3 = pd.Series(IList, name='IList')
        df_out = pd.concat([s1, s2, s3], axis=1)
        print(df_out)
        df_out.to_csv("output.txt", sep="\t", index=None)
        pMAP_plot(df=df_out)

    case "pMAP":
        pMAP_plot(df=df)