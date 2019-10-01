import opts as opt
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import glob
import sys

from skimage import measure,io
from skimage.filters import threshold_otsu
from skimage.segmentation import clear_border
from skimage.morphology import closing, square
from skimage.color import label2rgb

def plot_outlines(input_name, contours):
    input_img = io.imread(opt.input_directory + '/'+ input_name)
    fig, ax = plt.subplots()
    ax.imshow(input_img, cmap=plt.cm.gray)

    for n, contours in enumerate(contours):
        ax.plot(contours[:, 1], contours[:, 0], linewidth=0.5)

    ax.axis('image')
    ax.set_xticks([])
    ax.set_yticks([])
    input_no_ext = os.path.splitext(input_name)[0]
    plt.savefig(opt.output_directory + 'outlines/' + input_no_ext + '.pdf')
    plt.close()


path = sys.argv[1]
input_images = os.listdir(path+'/masks')
all_dfs = []
for image_name in input_images:

	mask=io.imread(path +'/masks/'+ image_name)
	thresh = threshold_otsu(mask)
	bw = closing(mask > thresh, square(3))

	cleared = clear_border(bw)

	label_image = measure.label(cleared)
	area = []
	length = []
	width = []
	solidity = []
	for region in measure.regionprops(label_image):
		area.append(region.area)
		length.append(region.major_axis_length)
		width.append(region.minor_axis_length)
		solidity.append(region.solidity)
		centroid.apped(region.centroid)

	image_df = pd.DataFrame(
		{'Image_name':image_name,
		'Area':area,
		'Major_axis_length':length,
		'Minor_axis_length':width,
		'Solidity':solidity
		})
	all_dfs.append(image_df)

all_data = pd.concat(all_dfs)
all_data.to_csv(path.rstrip('/') + '_all_data.csv')
print(path.rstrip('/') + '_all_data.csv')
