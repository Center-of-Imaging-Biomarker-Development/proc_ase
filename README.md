# proc_ase
Processing tool for quantitative modeling of asymmetric spin echo (ASE) signal data. ASE is an MRI pulse sequence in which the 180° refocusing pulse in incrementally shifted by time 𝜏 (ms). Modeling of the resulting signal data allows for voxel-wise quantification of R2', venous cerebral blood volume (vCBV), and oxygen extraction fraction (OEF). For more information on the sequence itself, please visit our [lab website](https://www.vumc.org/donahue-lab/methods). 

# Installation
## Prerequisites
To run this processing pipeline, you'll need:
- Python 3.9+
- MRI files in NIFTI format

# Usage
## Inputs
- 4D ASE data in NIFTI format (first echo)
>[!Note]
>The ASE MR pulse sequence may include multiple echoes; and thus the acquired data is 5-dimensional (nx, ny, nz, ndyn, nechoes) which is not supported in NIFTI file format. When converted to NIFTI, the signal data will be split into multiple 4D files (nx, ny, nz ,ndyn) dependent on the number of echoes.
- ASE sequence and modeling parameters in YAML format
>[!Note]
>The ASE sequence and modeling parameters must be provided in a YAML file that will be read by the pipeline. An example of the formatting can be found in config/ASE_HARMON.yml
## Outputs
### Voxel-wise
- R2prime.nii.gz
- vCBV.nii.gz
- rOEF.nii.gz
- Rsquared.nii.gz
- voxelwise_report.pdf

## Command options
- hct: Hemocrit level (range: 0-1; default=0.41)
- qa: Save intermediary files for quality assurance
- e2_fpath: Filepath to 2nd echo of ASE data in NIFTI fromat if using dual-echo sequence
- b0_fpath: Filepath to B0 map in NIFTI format if collected (must be same dimensions as ASE)
- kern_sz: Kernel size for Gaussian smoothing (in voxels)
- mot_cor: Motion correction using ANTs registration
- rm_tau: Tau shifts to remove during modeling of signal data (i.e. 17 17.5 18)

## Example usage
```
python proc_ase.py --in_fpath [/PATH/TO/IMAGE/] --out_fpath [/PATH/TO/SAVE/FILES] --img_seq ASE_HARMON --hct 0.41 --qa
```

# Version history
- 0.0.1
  - Initial release

# Authors
- Alex Song @alexksong
