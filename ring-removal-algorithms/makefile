objects = main.cpp image_filters.cpp image_transforms.cpp tiff_io.cpp
SHELL = /bin/csh
ring_remover_make_test.x : $(objects)
	module load tiff; \
	g++ -w -o ring_remover_make_test.x $(objects) -L$$TIFF_DIR/lib -I$$TIFF_DIR/include -ltiff
