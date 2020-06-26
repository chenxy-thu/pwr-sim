This folder contains the dataset used in manuscript "Enhanced reconstruction of structured illumination microscopy on polarized specimen"

The raw file of sample lambdaDNA and in vitro actin is saved in the corresponding folder.

The file "psf_dir_inVitroActin.mat" and "psf_dir_lambdaDNA.mat" contains the presaved direction PSF, which is calculated based on the SIM modulation parameter.

The file "actin_AF488.mat", "actin_atto633.mat", "inVitroActin.mat", "lambdaDNA.mat" contains the presaved direction PSF, polarization components and conventional SIM image. Due to the file limitation of 25MB in Github, the data is saved in "single" format, which should be retrieved to "double" format. Since the dynamic range of the sample images are not so high, this compression should not matter a lot.

These dataset are supposed to be used as the input for "PWRSIM_deconvolution_experiment.m" to run the spatial-angular deconvolution.