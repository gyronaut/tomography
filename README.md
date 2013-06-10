tomography
==========

This contains c++ code for functions associated with tomographic reconstruction,
primarily of reconstruction of parallel beam tomographic images, and with image
filtering algorithms, such as ring removal algorithms, and code to find the center
of rotation based on a set of images. All code is set up to be compiled on 64-bit
linux clusters, though they should be fairly OS independent.


 Ring Removal Algrothim (Recon Domain)
------------------------------------------
This code filters out ring artifacts, and works from the reconstructed image.  It uses
classes that perform certain standard filters and image transforms that should be set up
so that they can be used for a variety of applications. It's based on the paper "Comparison of ring artifact correction methods for flat-detector CT." by Prell, Kyriakou, Kalender.

The basic algorithm is as follows:
+ Transform the image into polar coordinates, using bilinear interpolation (with oversampling to limit data loss).
+ Perform 1D median filter in the Radial direction with a kernel size of the maximum width of the ring artifacts
+ Subtract this filtered image from the original image, and apply a thresholding to ignore possible structures
+ Perform 1D mean filter in the Azimuthal direction with a kernel size corresponding to the maximum angular arc you expect the rings to have.
+ Transform this filtered image back into Cartesian coordinates, leaving an image that only contains the ring artifacts.
+ subtract this ring artifact image from the original, leaving just the corrected image.

Currently, the slowest part of the algorithm is the median filter.  I've made efforts to speed it up, but it's still slower than many algorithms out there (it currently uses a rolling median filter version, where for every row/col, it does a median filter centered on the first element, then subtracts the last element added, adds the next element, and "bubbles" that element into place - basically adding a single element to an already ordered array, which takes much less time than just brute forcing every pixel of the image). There

Note: This code currently uses the TIFF library for file i/o, and relies on the input image being a 32 bit floating point tiff image.  However, the actual filter/transform functions could easily be changed to work on a different image datatype if necessary (the image data is stored in a double pointer of the form image[row][col], so any i/o that reads and writes the data to/from that format should work fine).
