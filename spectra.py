from generalLoad import load
import scipy.interpolate
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from line_profile import LineProfile
from input_file import InputFile
import os

def normal_interp(x, y, a, xi, yi):
    rbf = scipy.interpolate.Rbf(x, y, a)
    ai = rbf(xi, yi)
    return ai

def rescaled_interp(x, y, a, xi, yi):
    a_rescaled = (a - a.min()) / np.ptp(a)
    ai = normal_interp(x, y, a_rescaled, xi, yi)
    ai = np.ptp(a) * ai + a.min()
    return ai

def pMAP_plot(df, plot_as_image:bool=True, plot_as_scatter:bool=False, interp_resolution:complex=500j, cmap:str|list[str] = 'viridis', plot_title:str = "", plot_variable:str = "wave"):
        inputfile = InputFile

        plot_variable = plot_variable.lower()

        x = df["XList"]
        y = df["YList"]
        z = df["IList"]

        peak_number = 1

        reference_method = "initial"
        match reference_method.lower():
            case "mean":
                referenceWave = np.mean(z)
            case "initial":
                referenceWave = inputfile.initial_centers[peak_number-1] # np.min(z)
        
        shiftWave = z - referenceWave

        if shift_per_GPa is not None:
            stress = (-shiftWave/shift_per_GPa) # Postitive shiftWave means compressive and neagtive means tensile, see https://doi.org/10.1002/sia.1134

        match plot_variable:
            case "wave":
                variable2plot = z
                colorbar_label = r"Peak Wavenumber $[cm{-1}]$"
            case "shift":
                variable2plot = shiftWave
                colorbar_label = f"x-{referenceWave:.2f}"+ r" $[cm^{-1}]$"
            case "stress":
                variable2plot = 1000*stress # Into MPa
                colorbar_label = f"Stress (rel. to {referenceWave:.2f}) [MPa]"
            case "strain":
                variable2plot = 100*(stress/inputfile.GPa_per_strain[peak_number-1]) # From GPa to strain
                colorbar_label = f"Strain (rel. to {referenceWave:.2f}) [%]"

        xi, yi = np.mgrid[x.min():x.max():500j, y.min():y.max():500j]
        z_rescale = rescaled_interp(x, y, variable2plot, xi, yi)

        fig, ax = plt.subplots()

        if type(cmap) is list:
            cmap_image = cmap[0]
            cmap_scatter = cmap[1]
        else:
            cmap_image = cmap
            cmap_scatter = cmap

        if plot_as_image:
            im = ax.imshow(z_rescale.T, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap_image, vmin=794, vmax=797) # 
        if plot_as_scatter:
            ax.scatter(x, y, s=15, c=z, cmap=cmap_scatter) #  , vmin=566.3399156718697, vmax=566.5873836107614
        ax.set_ylim(ax.get_ylim()[::-1])

        print(f'Data range is {np.min(z_rescale)} to {np.max(z_rescale)}...')

        plt.title(plot_title)

        cb = plt.colorbar(im)
        cb.set_label(colorbar_label)

        plt.show()


filepath = r"C:\Users\mans3428\OneDrive - Nexus365\Postdoc Lilly Liu\Data\Raman\30102025\spectra.csv"
normaliseY = "off" # "off", "max", "norm"
inputfile = InputFile
peak_number:int = 1 # peak_number starts from 1
peak_name:str = inputfile.peak_names[peak_number-1] # -1 as peak_number starts from 1
shift_per_GPa = inputfile.shift_per_GPa[peak_number-1] # Either float or None. If None does not calculate, otherwise, the option to plot stress is available.
print(f'Working on {filepath}')
# df, scanType = load(filepath=filepath)

files = []
files_names = []

if os.path.isfile(filepath): # Does bob.txt exist?  Is it a file, or a directory?
    mode = "file"
    head, tail = os.path.split(filepath)
    if os.path.basename(tail) == "spectra.csv":
        print(os.path.basename(tail))   
        df0 = pd.read_csv(filepath, header=0, engine='python')
        for index, row in df0.iterrows():
            # if row['plot'] == 0:
            #     continue
            if row['magnification'] != 20:
                continue
            files.append(row['filepath'])
            files_names.append(row['name'])
    else:
        files.append(filepath)
        files_names.append(filepath)
elif os.path.isdir(filepath):
    mode = "files"
    for file in os.listdir(filepath):
        if file.endswith(".txt"):
            fullfilepath = os.path.join(filepath, file)
            files.append(fullfilepath)
            files_names.append(file)
    
print("Files loaded")
print(files_names)

plt.figure()

for i, filepath in enumerate(files):
    df, scanType = load(filepath=filepath)
    # print(scanType)
    if scanType != "SP":
        continue

    x = df['Wave'].copy()
    y = df['Intensity'].copy()
    maxY = np.max(y)
    minY = np.min(y)
    match normaliseY:
        case "max":
            ylabel = "Counts/max(Counts) [1]"
            y = y/maxY
        case "norm":
            ylabel = "(Counts-min)/(max-min) [1]"
            y = (y-minY)/(maxY-minY)
        case _:
            ylabel = "Counts [1]"

    # if maxY < 450:
    #     print(f"Skipping {files_names[i]} as below max Y limit")
    #     continue

    plt.plot(x, y, label=files_names[i])

# match mode:
#     case "file":
#         x = df["X"]
#         y = df["Y"]

#     case "files":
plt.ylabel(ylabel)
plt.xlabel(r"Wavenumber [$cm^-2$]")
plt.xlim(left=0)
plt.legend()
plt.show()

print('Finished spectra.py')