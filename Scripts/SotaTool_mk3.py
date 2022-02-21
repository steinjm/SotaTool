
# Copyright © Jeroen Stein, <steinjeroen@mgmail.com>
# Leiden University Medical Center, 2022
# This work is free. You can redistribute it and/or modify it under the
# terms of the Do What The Fuck You Want To Public License, Version 2,
# as published by Sam Hocevar. See http://www.wtfpl.net/ for more details.


# Compile .exe with:
# >>> pyinstaller --hidden-import pandas.plotting._matplotlib C:/Users/stein/Documents/Python/RvH_Project/SotaTool_mk3.py --noconfirm

# First libraries here to decrease start-up time of interface
from tkinter import filedialog
from tkinter import *

def browse():
    input_dir = filedialog.askdirectory()
    folder.delete(0, END)
    folder.insert(0, input_dir)
    folder.grid(row=3, column=2)

def help1():
    help1 = Tk()
    help1.geometry('400x100')
    help1.title('Help...')
    help1label = Label(help1, text='Background subtraction might be needed when image quality is low, \n'
                                   'but is also recommended in high quality images. \n'
                                   'Default ball size is automatically calculated as \n'
                                   'round(resolution*2) (About 2 \u03bcm)')
    help1label.place(relx=.5, rely=0.4, anchor="center")
    ok1 = Button(help1, text='Ok', command=help1.destroy).place(relx=.5, rely=0.8, anchor="center")

def help2():
    help1 = Tk()
    help1.geometry('400x180')
    help1.title('Help...')
    help1label = Label(help1, text='Enter an integer here, like 4. The raw \n'
                                   'images will be tiled in a 4x4 fashion and saved in a subdirectory. \n'
                                   'This increases the ability of the program to find\n'
                                   'the sarcomere length, but computing time is greatly increased.\n'
                                   'Any pre-segmented images found in the subfolder\n'
                                   'will be deleted to make place for new segments. ')
    help1label.place(relx=.5, rely=0.4, anchor="center")
    ok1 = Button(help1, text='Ok', command=help1.destroy).place(relx=.5, rely=0.8, anchor="center")

def help3():
    help1 = Tk()
    help1.geometry('400x100')
    help1.title('Help...')
    help1label = Label(help1, text='Enter the maximum offset in \u03bcm. Default is \n'
                                   '4*resolution. Use an offset bigger then the maximum \n'
                                   'sarcomere length. Reducing the offset here has the \n'
                                   'decreases computing time')
    help1label.place(relx=.5, rely=0.4, anchor="center")
    ok1 = Button(help1, text='Ok', command=help1.destroy).place(relx=.5, rely=0.8, anchor="center")

def help4():
    help1 = Tk()
    help1.geometry('400x100')
    help1.title('Help...')
    help1label = Label(help1, text='Override the default rolling ball size of the \n'
                                   'background subtraction. Default is 2 \u03bcm. ')
    help1label.place(relx=.5, rely=0.4, anchor="center")
    ok1 = Button(help1, text='Ok', command=help1.destroy).place(relx=.5, rely=0.8, anchor="center")

def help5():
    help1 = Tk()
    help1.geometry('400x100')
    help1.title('Help...')
    help1label = Label(help1, text='Disregard the lowest bin of the GLCM. This will decrease \n'
                                   'background noise, and increase the sarcomere scores. \n'
                                   'Default is 3. This will affect the organization score \n'
                                   'but not the sarcomere length')
    help1label.place(relx=.5, rely=0.4, anchor="center")
    ok1 = Button(help1, text='Ok', command=help1.destroy).place(relx=.5, rely=0.8, anchor="center")

def help6():
    help1 = Tk()
    help1.geometry('400x100')
    help1.title('Help...')
    help1label = Label(help1, text='Enter a maximum distance in um for the peak detection. \n'
                                   'This might be needed when there is a lot of noise in the image. ')
    help1label.place(relx=.5, rely=0.4, anchor="center")
    ok1 = Button(help1, text='Ok', command=help1.destroy).place(relx=.5, rely=0.8, anchor="center")

def help7():
    help1 = Tk()
    help1.geometry('400x100')
    help1.title('Help...')
    help1label = Label(help1, text='Check this box if you want no graph output. This is not \n'
                                   'recommended, as the graph is the only way to see if the \n'
                                   'analysis has been performed correctly.')
    help1label.place(relx=.5, rely=0.4, anchor="center")
    ok1 = Button(help1, text='Ok', command=help1.destroy).place(relx=.5, rely=0.8, anchor="center")

def help8():
    help1 = Tk()
    help1.geometry('400x100')
    help1.title('Help...')
    help1label = Label(help1, text='Check this box if you also want the interpolated \n'
                                   'curve data saved as a csv. ')
    help1label.place(relx=.5, rely=0.4, anchor="center")
    ok1 = Button(help1, text='Ok', command=help1.destroy).place(relx=.5, rely=0.8, anchor="center")

def adv():
    root.geometry('350x370')

    label4 = Label(root, text='Segmentation')
    label4.grid(row=6, column=1, sticky=E)
    Checkbutton(root, variable=var2, command=tileboxes).grid(row=6, column=2, sticky=W)

    label6 = Label(root, text='Enter offset distance (\u03bcm)').grid(row=7, column=1, sticky=E)
    offset = Entry(root, textvariable=offset_val).grid(row=7, column=2)

    label7 = Label(root, text='Enter rolling ball size (\u03bcm)')
    label7.grid(row=8, column=1, sticky=E)
    ballsize = Entry(root, textvariable=ballsize_val)
    ballsize.grid(row=8, column=2)

    label9 = Label(root, text='Disregard lowest GLCM-bin')
    label9.grid(row=9, column=1, sticky=E)
    noise_red.grid(row=9, column=2)

    label8 = Label(root, text='Max sarcomere length (\u03bcm)')
    label8.grid(row=10, column=1, sticky=E)
    max_length = Entry(root, textvariable=max_length_val)
    max_length.grid(row=10, column=2)

    label8 = Label(root, text='Skip graphs to save time')
    label8.grid(row=11, column=1, sticky=E)
    Checkbutton(root, variable=var3).grid(row=11, column=2, sticky=W)

    label10 = Label(root, text='Output interpolated data')
    label10.grid(row=12, column=1, sticky=E)
    Checkbutton(root, variable=var4).grid(row=12, column=2, sticky=W)

    button2 = Button(root, text='?', command=help2).grid(row=6, column=4)
    button3 = Button(root, text='?', command=help3).grid(row=7, column=4)
    button4 = Button(root, text='?', command=help4).grid(row=8, column=4)
    button5 = Button(root, text='?', command=help5).grid(row=9, column=4)
    button6 = Button(root, text='?', command=help6).grid(row=10, column=4)
    button7 = Button(root, text='?', command=help7).grid(row=11, column=4)
    button8 = Button(root, text='?', command=help8).grid(row=12, column=4)

def tileboxes():
    label4 = Label(root, text='Segmentation')
    label4.grid(row=6, column=1, sticky=E)
    labeltx = Label(root, text='x').place(x=200, y=138)
    tileboxx.place(x=210, y=140)
    labelty = Label(root, text='y').place(x=240, y=138)
    tileboxy.place(x=250, y=140)

def analyze():
    input_dir = input1.get()
    tileval_x = tilevalx.get()
    tileval_y = tilevaly.get()
    if input_dir == '...' and voxelvalue.get() == '':
        stop1 = Tk()
        stop1.geometry('400x100')
        stop1.title('Error')
        stop1label = Label(stop1, text='You must enter an input folder and resolution [px/\u03bcm]')
        stop1label.place(relx=.5, rely=0.4, anchor="center")
        ok1 = Button(stop1, text='Ok', command=stop1.destroy)
        ok1.place(relx=.5, rely=0.8, anchor="center")
    if not input_dir == '...' and voxelvalue.get() == '':
        stop1 = Tk()
        stop1.geometry('400x100')
        stop1.title('Error')
        stop1label = Label(stop1, text='You must enter the resolution [px/\u03bcm]')
        stop1label.place(relx=.5, rely=0.4, anchor="center")
        ok1 = Button(stop1, text='Ok', command=stop1.destroy)
        ok1.place(relx=.5, rely=0.8, anchor="center")
    if input_dir == '...' and not voxelvalue.get() == '':
        stop1 = Tk()
        stop1.geometry('400x100')
        stop1.title('Error')
        stop1label = Label(stop1, text='You must enter input folder')
        stop1label.place(relx=.5, rely=0.4, anchor="center")
        ok1 = Button(stop1, text='Ok', command=stop1.destroy)
        ok1.place(relx=.5, rely=0.8, anchor="center")
    if int(tileval_x) > int(8) or int(tileval_y) > int(8):
        stop1 = Tk()
        stop1.geometry('400x100')
        stop1.title('Error')
        stop1label = Label(stop1, text='Thats a lot of tiles')
        stop1label.place(relx=.5, rely=0.4, anchor="center")
        ok1 = Button(stop1, text='Ok', command=stop1.destroy)
        ok1.place(relx=.5, rely=0.8, anchor="center")
    if not input_dir == '...' and not voxelvalue.get() == '':
        root.destroy()

# starting the interface
root = Tk()
root.geometry('350x200')
root.title('SotaTool')

# Title
label1 = Label(root, text='Sarcomere Organization and \nTexture Analysis Tool')
label1.place(relx=0.5, y=30, anchor="center")
label1.config(font=(20))

# Define default advanced settings here, to avoid confusion if adv is skipped
tilevalx = StringVar()
tileboxx = Entry(root, text=tilevalx, width=4)
tileboxx.insert(0, '4')
tilevaly = StringVar()
tileboxy = Entry(root, text=tilevaly, width=4)
tileboxy.insert(0, '4')

ballsize_val = StringVar()
max_length_val = StringVar()
offset_val = StringVar()
var2 = IntVar()
var3 = IntVar()
var4 = IntVar()
noise_red_val = StringVar()
noise_red = Entry(root, text=noise_red_val)
noise_red.insert(0, '3')

# choose dir
label2 = Label(root, text='Choose input directory')
label2.grid(row=3, column=1, sticky=E)
labelempty = Label(root, text='').grid(row=0, column=0, sticky=E)
labelempty1 = Label(root, text='').grid(row=1, column=0, sticky=E)
labelempty2 = Label(root, text='').grid(row=2, column=0, sticky=E)
labelempty3 = Label(root, text='').grid(row=13, column=0, sticky=E)

input1 = StringVar()
input_dir = '...'
folder = Entry(root, textvariable=input1)
folder.insert(0, input_dir)
folder.grid(row=3, column=2)

button = Button(root, text='Browse', command=browse)
button.grid(row=3, column=4)

# Enter voxel size
label5 = Label(root, text='    Enter resolution (pixel/\u03bcm)')
label5.grid(row=4, column=1, sticky=E)

voxelvalue = StringVar()
voxel = Entry(root, textvariable=voxelvalue)
voxel.grid(row=4, column=2)

# Other settings
label3 = Label(root, text='Background subtraction')
label3.grid(row=5, column=1, sticky=E)
var1 = IntVar(value=1)
Checkbutton(root, variable=var1).grid(row=5, column=2, sticky=W)

# Help screens
button1 = Button(root, text='?', command=help1).grid(row=5, column=4)

button_adv = Button(root, text='Advanced', command=adv).grid(row=14, column=1, sticky=W)

Button(root, text='Start analysis', command=analyze).grid(row=14, column=2, sticky=W)

root.mainloop()

# rest of the libraries
from skimage import (restoration)
import warnings
import numpy as np
from skimage.feature import greycomatrix, greycoprops
import cv2
import os
import matplotlib.pyplot as plt
import pandas as pd
import glob
import imutils
from scipy.signal import find_peaks
from scipy.interpolate import CubicSpline
import matplotlib
from tqdm import tqdm
import seaborn as sns
import multiprocessing
from joblib import Parallel, delayed

matplotlib.use("Agg")
warnings.filterwarnings("ignore")  # I know, I know..
num_cores = multiprocessing.cpu_count()

# CoV-formula, input from user translations and output folder creation
cv = lambda x: np.std(x, ddof=1) / np.mean(x) * 100
input_dir = input1.get()
output_dir = input_dir + '/output'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
folderLen = len(input_dir)
Summary = pd.DataFrame()
voxel = float(voxelvalue.get())
noise_red_val = int(noise_red_val.get())

# If boxes are empty, use default values for ballsize and maxlength
if len(offset_val.get()) == 0:  # 4 um default offset size
    prompt = round(voxel * 4)
else:
    prompt = round(float(offset_val.get()) * voxel)
if len(ballsize_val.get()) == 0:  # rolling ball size, default 2um
    ballsize_val = round(voxel * 2)
else:
    ballsize_val = round(float(ballsize_val.get()) * voxel)
if len(max_length_val.get()) == 0:
    max_length_val = prompt
else:
    max_length_val = round(float(max_length_val.get()) * voxel)

x = np.arange(1, prompt, 1)
x = x * (1 / voxel)
x_new = np.linspace((1 / voxel), max(x), prompt * 100)
distances = [*range(1, prompt, 1)]
distances2 = [*range(0, prompt - 1, 1)]
angles = [*range(1, 181, 1)]

# Tiling
if var2.get() == 1:
    tiling_dir = input_dir + '/segments'
    if not os.path.exists(tiling_dir):
        os.makedirs(tiling_dir)
    files = glob.glob(tiling_dir + '/*')
    for f in files:
        os.remove(f)
    tileval_y = tilevaly.get()
    tileval_x = tilevalx.get()
    for img in glob.glob(input_dir + "/*.tif"):
        image = cv2.imread(img)
        image2 = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
        numrows = int(tileval_x)
        numcols = int(tileval_y)
        height = int(image2.shape[0] / numrows)
        width = int(image2.shape[1] / numcols)
        for row in range(numrows):
            for col in range(numcols):
                y0 = row * height
                y1 = y0 + height
                x0 = col * width
                x1 = x0 + width
                cv2.imwrite(tiling_dir + '/' + img[folderLen:].replace('.tif', '') + '_tile_%d%d.tif' % (row, col), image2[y0:y1, x0:x1])

def algorithm(iAngle):
    rotateI = imutils.rotate_bound(gray_image, iAngle)
    glcm = greycomatrix(rotateI,
                        distances=distances,
                        angles=[0],
                        symmetric=False,
                        normed=False)[noise_red_val:, noise_red_val:, :, :]

    correlation = greycoprops(glcm, 'correlation')
    df = pd.DataFrame(correlation)
    df2[iAngle] = df
    asm = greycoprops(glcm[:, :, 0:1, :], 'ASM').mean()
    homogeneity = greycoprops(glcm[:, :, 0:1, :], 'homogeneity').mean()


# List files:
if var2.get() == 1:
    files = glob.glob(input_dir + "/*.tif") + glob.glob(tiling_dir + "/*.tif")
else:
    files = glob.glob(input_dir + "/*.tif")

# Analysis GLCM and sarcomere length, loop over folder
for img in tqdm(files):
    image = cv2.imread(img)
    image_ID = '/' + img[folderLen:].replace('.tif', '').replace('\\', '').replace('/segments', '')
    # image = cv2.imread('C:/Users/stein/Documents/Python/Raw data for article/Matthew cells/TID/Image005_TID_ch01.tif')

    gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

    # Background subtraction
    if var1.get() == 1:
        background = restoration.rolling_ball(gray_image, radius=int(ballsize_val))
        gray_image = gray_image - background

    # prepare dataframes for loop
    df2 = pd.DataFrame(index=distances2, columns=angles)
    df2 = df2.astype(np.float)

    # multithread function of GLCM and save df2 to .csv
    Parallel(n_jobs=num_cores, prefer="threads")(delayed(algorithm)(iAngle) for iAngle in angles)

    df2.to_csv(output_dir + image_ID + '_corr' + ".csv")

    # Peak detection in and amplitude calculation in all curves
    # maximum sarcomere length taken into account here
    peaks = []
    for column in df2:
        ser = pd.Series(df2[column][:max_length_val])
        indices = find_peaks(ser)[0]
        flag = not np.any(indices)
        if flag:
            indices = np.array([1, 1])
        peak = pd.Series(ser[indices[0]] - min(ser[0:indices[0]]))
        peaks.append(peak)

    peak2 = {'angle': pd.Series(range(1, 181)), 'SOS': [i[0] for i in peaks]}
    peaklist = pd.DataFrame(peak2)
    SOS_max = peaklist.iloc[peaklist.SOS.argmax()]

    # So now we have the angle and SOS of the highest curve. Now we need to extract the sarcomere length from df2
    # using the splined data
    max_curve = pd.Series(df2[SOS_max[0]])

    y = np.array(max_curve.values.tolist())

    f = CubicSpline(x, y, bc_type='natural')
    y_new = f(x_new)

    scaled_maxcurve = pd.DataFrame(pd.Series(x_new), columns=['Offset'])
    scaled_maxcurve['Correlation'] = y_new.tolist()
    if var4.get() == 1:
        scaled_maxcurve.to_csv(output_dir + image_ID + 'int' + ".csv")

    peak_s_length = find_peaks(y_new)[0]
    flag = not np.any(peak_s_length)
    if flag:  # these parameters should be 0 if there are no peaks in the splined graph
        s_length = 'NA'
        angle_val = 'NA'
        SOS = 'NA'
    else:
        s_length = round(scaled_maxcurve.iloc[peak_s_length[0]][0], 3)
        angle_val = int(SOS_max[0])
        SOS = SOS_max[1]

    # Sure, just calculate the Coefficient of Variation real quick
    CoV_offset = find_peaks(np.array(df2[SOS_max[0]]))
    if not np.any(CoV_offset[0]):
        CoV = 'NA'
        CoV_offset = np.array([[0.0], [0]], dtype=np.int64)
    else:
        CoV = cv(df2.iloc[[CoV_offset[0][0] + 1]].squeeze()[1:])

    # Create epic graph which took way too long to design, fml
    # Melt df, order hues on CoV-offset in a list, choose palette, plot all lines with seaborn
    if var3.get() == 0:
        df2.insert(loc=0, column='Offset', value=x)
        df2 = df2.drop([90], 1)
        df2 = df2.drop([180], 1)
        angles2 = [e for e in angles if e not in (90, 180)]
        df3 = pd.melt(df2, 'Offset')
        hue_order = df2.iloc[[CoV_offset[0][0] + 1]].drop(['Offset'], axis=1).transpose().assign(angle=angles2)
        hue_order = hue_order.sort_values(by=hue_order.columns[0])
        hue_order = hue_order['angle'].values.tolist()

        colors = sns.color_palette('light:b', n_colors=len(hue_order))

        plt.close('all')
        plt.scatter(x, y, marker='x', c='black')  # Add raw points of best graph
        ax = sns.lineplot(x='Offset', y='value', hue='variable', data=df3, alpha=0.3, palette=colors, legend=False,
                          hue_order=hue_order)  # Add all angles
        ax.grid(b=True, which='major', color='black', linewidth=0.075)
        if s_length == 'NA':  # Create abline, but on 0 if s_length not detected
            plt.axvline(x=0, c='black')
        else:
            plt.axvline(x=s_length, c='black')
        plt.plot(scaled_maxcurve['Offset'], scaled_maxcurve['Correlation'], 'b')  # add splined graph
        plt.xlabel('Offset [um]')
        plt.ylabel('Correlation')
        fig = plt.gcf()
        fig.savefig(output_dir + image_ID + '_output' + '.tif')

    # Create new GLCM with furthest pixel offset for ASM and Homogeneity
    glcm = greycomatrix(gray_image,
                        distances=[max(distances2)],
                        angles=[0],
                        symmetric=False,
                        normed=False)[noise_red_val:, noise_red_val:, :, :]
    asm = greycoprops(glcm[:, :, 0:1, :], 'ASM').mean()
    homogeneity = greycoprops(glcm[:, :, 0:1, :], 'homogeneity').mean()

    ### Compile useful output in summary excel file
    data = [
        [img[folderLen:].replace('.tif', '').replace('\\', '').replace('segments', ''), s_length, angle_val, SOS, CoV,
         asm, homogeneity]]
    final = pd.DataFrame(data, columns=['Name', 'Sarcomere length', 'Angle', 'Sarcomere Organization Score',
                                        'Alignment index', 'Uniformity', 'Homogeneity'])
    Summary = Summary.append(final)
    Summary.to_csv(output_dir + "/Summary" + ".csv")


