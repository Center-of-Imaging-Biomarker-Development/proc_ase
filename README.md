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

## Outputs
### Voxel-wise
- R2prime.nii.gz
- vCBV.nii.gz
- rOEF.nii.gz
- Rsquared.nii.gz
- voxelwise_report.pdf

## Command options


## Example usage
```
python proc_ase.py --in_fpath [/PATH/TO/IMAGE/] --out_fpath [/PATH/TO/SAVE/FILES] --img_seq ASE_HARMON --hct 0.41 --qa
```

# Version history
- 0.0.1
  - Initial release

# Authors
- Alex Song @alexksong
