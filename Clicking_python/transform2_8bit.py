import os
import numpy as np
from skimage import measure,io
import sys

path = sys.argv[1]
input_images =  sorted([f for f in os.listdir(path) if f.endswith('.tif')])
print(input_images)

for  img_name in input_images:
    img16 = io.imread(path + '/' +img_name)
    img8 = (img16/256).astype('uint8')
    io.imsave(path + '/' + img_name, img8)
    print('8 bit saved')
