Ring Removal Algorithm in Recon Domain
======================================

This is code for a ring removal algorithm that
works in the reconstructed domain.  It's based
on the paper "Comparison of ring artifact correction
methods for flat-detector CT" by Prell, et al. It works by:

* Performing a polar transform on the reconstructed image.
* Doing a median filter in the radial direction (with a 
kernel width that depends on the radius).
* Taking the difference between the filtered and original 
polar image.
* Performing a mean filter in the azimuthal direction 
(again with a kernel width that depends on the radius).
* Doing an inverse polar transformation (giving an image
of just the rings).
* Subtracting the ring image from the original image.
