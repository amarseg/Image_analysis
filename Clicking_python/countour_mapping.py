import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import glob
import sys

from scipy import ndimage as ndi
from skimage import measure, io
from skimage.filters import threshold_otsu
from skimage.segmentation import clear_border, mark_boundaries
from skimage.morphology import closing, square,watershed
from skimage.color import label2rgb,rgb2gray
from skimage.util import invert,img_as_float64
from skimage.feature import canny, peak_local_max
from skimage.morphology import erosion, dilation, opening, closing, white_tophat


path = sys.argv[1]
mask_images = sorted([f for f in os.listdir(path) if f.endswith('.npy')])[0]
t = os.path.splitext(mask_images)[0].split('_P')[0]


all_dfs = []
#Get just one to test
#for image_name, mask_name in zip(input_images, mask_images):
mask = np.load(path + '/' + mask_images)
edges_canny = canny(mask[:,:,1])
edges_otsu = threshold_otsu(mask[:,:,1])

# fig, (ax1,ax2) = plt.subplots(nrows=1, ncols=2,
#                                     sharex=True, sharey=True)
#
# ax1.imshow(edges_canny, cmap=plt.cm.gray)
# ax1.axis('off')
# ax1.set_title('Canny filter, $\sigma=1$', fontsize=20)
#
# ax2.imshow(mask[:,:,1] < edges_otsu, cmap='gray', interpolation='nearest')
# ax2.axis('off')
# ax2.set_title('Otsu threshold', fontsize=20)
#
# plt.show()

bw = closing(mask[:,:,1] > edges_otsu, square(3))

cleared = clear_border(bw)

label_image = measure.label(cleared)

image = io.imread(path +'/' + t + '.tif')

image_label_overlay = label2rgb(label_image, image=image)

# fig, ax = plt.subplots(figsize=(10, 6))
# ax.imshow(image_label_overlay)
#
#
#
# ax.set_axis_off()
# plt.tight_layout()
# plt.show()

##############Watershed###############################
distance = ndi.distance_transform_edt(mask[:,:,1])
local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((3, 3)),
                            labels=mask[:,:,1])
markers = ndi.label(local_maxi)[0]
labels = watershed(-distance, markers, mask=mask[:,:,1])

fig, axes = plt.subplots(ncols=3, figsize=(9, 3), sharex=True, sharey=True)
ax = axes.ravel()

ax[0].imshow(mask[:,:,1], cmap=plt.cm.gray)
ax[0].set_title('Overlapping objects')
ax[1].imshow(-distance, cmap=plt.cm.gray)
ax[1].set_title('Distances')
ax[2].imshow(labels, cmap=plt.cm.nipy_spectral)
ax[2].set_title('Separated objects')

for a in ax:
    a.set_axis_off()

fig.tight_layout()
plt.show()
