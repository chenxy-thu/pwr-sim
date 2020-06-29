Code and example dataset for manuscript "Enhanced reconstruction of structured illumination microscopy on polarized specimen"

The SIM reconstruction frame work is based on open-source software fairSIM (www.fairSIM.org), and re-written in Matlab. All the functions related to SIM reconstruction is located in the Folder "functions"

The function related to spatial-angular deconvolution is located in folder "deconvolution"

The file "PWRSIM_deconvolution_experiment.m" is to apply spatial-angular on the input dataset. The dataset is required to have polarization components, directional PSF, and optional SIM image. The example dataset is pre-saved in folder "dataset".

The file "PWRSIM_simulation_demo.m" is the simulation of PWR-SIM.

The file "polarization_components_generation.m" is to generated the polarization components from raw SIM images and prepared the dataset for spatial-angular deconvolution.

The directional PSF and polarization components can also be generated from the fairSIM plugin in ImageJ. The code and detailed instruction is in file "polarization_components_load_from_fairSIM.m".

For any question, please contact Xingye Chen, chenxy16@mails.tsinghua.edu.cn