from generalLoad import load
import scipy.interpolate
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from line_profile import LineProfile
from input_file import InputFile

def normal_interp(x, y, a, xi, yi):
    rbf = scipy.interpolate.Rbf(x, y, a)
    ai = rbf(xi, yi)
    return ai

def rescaled_interp(x, y, a, xi, yi):
    a_rescaled = (a - a.min()) / np.ptp(a)
    ai = normal_interp(x, y, a_rescaled, xi, yi)
    ai = np.ptp(a) * ai + a.min()
    return ai

def pMAP_plot(df, plot_as_image:bool=True, plot_as_scatter:bool=False, interp_resolution:complex=500j, cmap:str|list[str] = 'viridis', plot_title:str = ""):
        x = df["XList"]
        y = df["YList"]
        z = df["IList"]

        xi, yi = np.mgrid[x.min():x.max():500j, y.min():y.max():500j]
        z_rescale = rescaled_interp(x, y, z, xi, yi)

        fig, ax = plt.subplots()

        if type(cmap) is list:
            cmap_image = cmap[0]
            cmap_scatter = cmap[1]
        else:
            cmap_image = cmap
            cmap_scatter = cmap

        if plot_as_image:
            im = ax.imshow(z_rescale.T, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap_image) # , vmin=566.3399156718697, vmax=566.5873836107614
        if plot_as_scatter:
            ax.scatter(x, y, s=15, c=z, cmap=cmap_scatter) #  , vmin=566.3399156718697, vmax=566.5873836107614
        ax.set_ylim(ax.get_ylim()[::-1])

        print(f'Data range is {np.min(z_rescale)} to {np.max(z_rescale)}...')

        plt.title(plot_title)

        cb = plt.colorbar(im)
        plt.show()


filepath = r"xxx"
inputfile = InputFile
peak_number:int = 2 # peak_number starts from 1
peak_name:str = inputfile.peak_names[peak_number-1] # -1 as peak_number starts from 1
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
            if i == 0:
                 print(df_l)
            peak_wave = df_l.loc[df_l['Peak Index'] == peak_number, 'Center Grvty'].iloc[0]
            IList[i] = peak_wave
            i += 1
        XList = unique_combinations['X'].to_list()
        YList = unique_combinations['Y'].to_list()

        s1 = pd.Series(XList, name='XList')
        s2 = pd.Series(YList, name='YList')
        s3 = pd.Series(IList, name='IList')
        df_out = pd.concat([s1, s2, s3], axis=1)
        print(df_out)
        
        # df_out.to_csv("output.txt", sep="\t", index=None)

        plot_title = f"Peak {peak_name}"

        pMAP_plot(df=df_out, cmap=['grey', 'viridis'], plot_title=plot_title, plot_as_scatter=True)

    case "pMAP":
        pMAP_plot(df=df, cmap='grey')

print('Finished map.py')