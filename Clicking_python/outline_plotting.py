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
from skimage.color import label2rgb,rgb2gray

def plot_outlines(input_name, contours, file_path):
	input_img = io.imread(path + '/'+ input_name)
	fig, ax = plt.subplots()
	ax.imshow(input_img, cmap=plt.cm.gray)

	for n, contours in enumerate(contours):
		ax.plot(contours[:, 1], contours[:, 0], linewidth=0.5)

	ax.axis('image')
	ax.set_xticks([])
	ax.set_yticks([])
	input_no_ext = os.path.splitext(input_name)[0]
	if not os.path.exists(path + '/outlines/'):
		os.makedirs(path + '/outlines/')

	plt.savefig(path + '/outlines/' + input_no_ext + '.jpg')
	plt.close()

path = sys.argv[1]
input_images =  sorted([f for f in os.listdir(path) if f.endswith('.tif')])
print(input_images)
mask_images = sorted([f for f in os.listdir(path) if f.endswith('.png')])
all_dfs = []

print('Paths are correct')
for image_name, mask_name in zip(input_images, mask_images) :

	mask= io.imread(path +'/' + mask_name)
	mask_bw = rgb2gray(mask)
	thresh = threshold_otsu(mask_bw)
	bw = closing(mask_bw > thresh, square(3))

	cleared = clear_border(bw)

	label_image = measure.label(cleared)
	area = []
	length = []
	width = []
	solidity = []
	centroid = []
	for region in measure.regionprops(label_image):
		area.append(region.area)
		length.append(region.major_axis_length)
		width.append(region.minor_axis_length)
		solidity.append(region.solidity)
		centroid.append(region.centroid)

	plot_outlines(input_name = image_name,
				contours = measure.find_contours(mask_bw, 0.9),
				file_path = path
	)
	print('Print Outlines')
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
