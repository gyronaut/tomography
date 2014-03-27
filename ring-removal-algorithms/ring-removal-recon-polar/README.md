Recon Ring Remover v. 1.1.0
======================================

This is code for a ring artifact removal algorithm that
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

To Build on Ubuntu:
---------------------

You'll need the latest version of libtiff and associated files (which can be obtained using "sudo apt-get install libtiff-dev"), as well as g++ or a similar compiler.  Once you have those, to compile simply run:

    g++ -w -o recon_ring_remover.1.1.0 main.cpp image_filters.cpp image_transforms.cpp tiff_io.cpp -ltiff

Usage:
---------------

(calling the program without any arguments will also print the usage to stdout)

ring_remover_recon [input path] [output path] [input root] [output root] [first file num] [last file num] [center x y] [max ring width] [thresh min max] [ring threshold] [angular min] [verbose]

                       [input path]    Path to the folder containing the input images (32 bit floating point tiffs).

                      [output path]    Path to the folder to hold the filtered images.

                  [input file root]    Root name of the image to filter, with no file

                                       number or extension.

                 [output file root]    Root of the output file name (number and file extension

                                       will be appended to this).

                   [first file num]    Number of first file to filter.

                    [last file num]    Number of last file to filter.

                     [center: x, y]    X and Y value of the center of the rings (floats), with origin at the
		     
		     	               top left corner of the image. If 0 is entered for both X and Y, the
				       
				       exact middle of the image is assumed (standard use case)

                   [max ring width]    maximum width of the rings to be filtered in pixels (int).

                   [thresh min max]    min and max values for portion of image to filter (floats).

                   [ring threshold]    Rings are treated as offsets to the corrected image. This

                                       threshold is the maximum value of an offset due to a

                                       ring artifact (float).

                      [angular min]    minimum angle in degrees (int) to be considered ring artifact

                          [verbose]    0 = No Output messages (default), 1 = Output messages

An example call on reconstructed images "recon_slice_10.tif" through "recon_slice_25.tif", inside folder /path/to/recon, with the center of the rings in the center of the image, a max ring with of 30, min/max values of -300 300, threshold of 300, and an angular minimum of 30 degrees:

    ./recon_ring_remover.x /path/to/recon/ /path/to/recon/corrected/ recon_slice_ recon_slice_corrected_ 10 25 0 0 30 -300 300 300 30 1

*Notes on Standard Usage:*

While the algorithm is set-up to allow for thresholding (so that only the areas of the image within the min/max range are filtered), for our use cases (samples with a wide range of materials/densities) the artifact removal works best if the entire image is considered, so the min and max are set past the actual min and max pixel values (-300 and 300 work fine).

Likewise, we find that our ring intensity can vary significantly, so a ring threshold of 300 is often used to consider any pixel as possibly part of an artifact.

By ignoring the thresholding component, it becomes easier to incidentally filter out parts of the sample as well, so a ring width of 30 pixels, and an angular minimum of 30 degrees is often best for removing rings and leaving sample structures intact.



