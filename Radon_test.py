import numpy as np
import matplotlib.pyplot as plt

from skimage.io import imread
from skimage import data_dir
from skimage.transform import radon, rescale
from skimage.color import rgb2gray

image = imread("img0070.png")
image=rgb2gray(image)
image = rescale(image, scale=0.4, mode='reflect')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4.5))
ax1.set_title("Original")
ax1.imshow(image, cmap=plt.cm.Greys_r)

theta = np.linspace(0., 180., max(image.shape), endpoint=False)
sinogram = radon(image, theta=theta)

ax2.set_title("Radon transform\n(Sinogram)")
ax2.set_xlabel("Projection angle (deg)")
ax2.set_ylabel("Projection position (pixels)")
#ax2.imshow(sinogram, extent=(0, 180, 0, sinogram.shape[0]), aspect='auto')



image_const = np.ones(image.shape)
sinogram_const = radon(image_const, theta=theta)
data = np.mean(sinogram**2, axis=0)
data_const = np.mean(sinogram_const**2, axis=0)
lons=np.linspace(0,180,sinogram.shape[1])
ax2.plot(lons,data/data_const)
#ax2.plot(data)
fig.tight_layout()

plt.show()
