import matplotlib
print(matplotlib.get_backend())
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import glob
import sys
import re

from skimage import measure,io
from skimage.filters import threshold_otsu
from skimage.segmentation import clear_border
from skimage.morphology import closing, square
from skimage.color import label2rgb
#############Define clicking function###############
def click_on_cell(event):
    print('leave_figure', event.canvas.figure)
    event.canvas.figure.patch.set_facecolor('grey')
    event.canvas.draw()
#########################
#Read data and plot it before the clicking begins
########################
images_path = r'Z:\Pers_Amalia\Screening\clicking_test\images'
masks_path = r'Z:\Pers_Amalia\Screening\clicking_test\masks'

#Label outlines and get coordinates and shit
home = os.getenv("HOME")
#masks_img = glob.glob(masks_path + '//*//*.tiff')
masks_img = [r'/mnt/quantex/Pers_Amalia/Screening/clicking_test/masks/dir_194/F3--W00063--P00007--Z00000--T00003--Blank.tif']

all_dfs = []

for img in masks_img:
    print(img)
    mask=io.imread(img)
    thresh = threshold_otsu(mask)
    bw = closing(mask > thresh, square(3))
    cleared = clear_border(bw)
    label_image = measure.label(cleared)
    contours = measure.find_contours(label_image, 0.8)#Measure contours

    area = []
    length = []
    width = []
    solidity = []

    for region in measure.regionprops(label_image):
        area.append(region.area)
        length.append(region.major_axis_length)
        width.append(region.minor_axis_length)
        solidity.append(region.solidity)

    image_name = os.path.basename(img)

    image_df = pd.DataFrame(
        {'Image_name': image_name,
         'Area': area,
         'Major_axis_length': length,
         'Minor_axis_length': width,
         'Solidity': solidity
        })
    all_dfs.append(image_df)

    proper_image = re.sub('masks','images', img)
    input_image = io.imread(proper_image)
    #Plot contours on top of image
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.imshow(input_image,cmap=plt.cm.gray) #Need cmap

    for n, contour in enumerate(contours):
        ax.plot(contour[:, 1], contour[:, 0], linewidth=1)

    ax.set_axis_off()
    plt.tight_layout()

    fig.canvas.mpl_connect('button_press_event', click_on_cell)

    plt.show()
