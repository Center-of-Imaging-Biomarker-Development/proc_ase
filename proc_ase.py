# Main script for processing ASE data
# Generic Python import
import argparse
import os

# Project-specific import
import utils.asedata as asedata
import utils.preprocess as preprocess
import utils.voxelfit as voxelfit
import utils.asereport as asereport

def proc_ase(img_fpath, img_seq, hct, out_fpath, qa=False, e2_fpath=None, b0_fpath=None, yml_fpath=None, rm_tau=None, kernsz=3, motcor=False):
    '''
    ======= Processing Framework =======
    1) Load image data and imaging parameters
    2) Pre-process image data - includes re-order, gauss smooth, and motion correction
    3) Voxel-wise processing
    4) Preliminary gray and white matter results
    5) Report
    ====================================
    '''
    print('Starting...')
    # ----- 1) Load image data and imaging parameters ----- #
    if img_fpath.lower().endswith('.nii'):
        subj_ase = asedata.AseDataNifti1(img_fpath=img_fpath, img_seq=img_seq, hct=hct, out_fpath=out_fpath, e2_fpath=e2_fpath,
                                         rm_tau=rm_tau, b0_fpath=b0_fpath, yml_fpath=yml_fpath)
    elif img_fpath.lower().endswith('.dcm'):
        pass # Comment: Will need to add this later on - aks 2/9/27

    # ----- 2) Pre-process image data - includes re-order, gauss smooth, and/or motion correction ----- #
    if not motcor:
        data_prep = preprocess.Preprocess()
    else:
        data_prep = preprocess.PreprocessANTs()
    data_prep.run(data=subj_ase, kernsz=kernsz, qa=qa, out_fpath=out_fpath+'/'+subj_ase.params['mrID'])

    subj_ase.set_fitdata()
    
    # ----- 3) Voxel-wise processing ----- #
    vout_fpath = out_fpath+'/'+subj_ase.params['mrID']+'/voxel'
    os.makedirs(vout_fpath, exist_ok=True)
    voxel_fitting = voxelfit.VoxelFit(data=subj_ase)
    voxel_fitting.run(data=subj_ase, out_fpath=vout_fpath)

    # ----- 5) Create report for voxel-wise results ----- #
    in_str = asereport._gen_input(img_fpath, img_seq, hct, out_fpath, qa, e2_fpath, b0_fpath, yml_fpath, rm_tau, kernsz, motcor)
    asereport.generate_outputs(vout_fpath, subj_ase.params['mrID'], subj_ase, in_str, vout_fpath+'/voxelmaps_qa.png', vout_fpath+'/voxelfit_qa.png')

if __name__ == "__main__":
    # Parsing arguments
    parser = argparse.ArgumentParser(description='Quantify oxygen extraction fraction (OEF) from Asymmetric Spin Echo (ASE) data using Yablonskiy-Haacke 1-compartment model')
    parser.add_argument('--img_fpath', type=str, required=True, help='Filepath to ASE image (echo1 if dual-echo)')
    parser.add_argument('--img_seq', type=str, required=True, help='ASE sequence variant. Current options:\n1)ASE_HARMON\n2)Other')
    parser.add_argument('--yml_fpath', type=str, required=False, help='If ASE sequence variant is Other, include filepath to config.yml')
    parser.add_argument('--hct', type=float, required=False, default=0.41, help='Hematocrit level (0-1; default=0.41)')
    parser.add_argument('--out_fpath', type=str, required=True, help='Filepath to save outputs')
    parser.add_argument('--qa', action='store_true', help='Save intermediary files for quality assurance')

    parser.add_argument('--e2_fpath', type=str, required=False, default=None, help='Filepath to 2nd echo of ASE image (NIFTI1 formats)')
    parser.add_argument('--b0_fpath', type=str, required=False, default=None, help='Filepath to B0 image if acquired (must be same dimensions as ASE image); default=3.0T')
    parser.add_argument('--rm_tau', type=float, required=False, default=None, nargs='+', help='Tau values to be excluded from fitting (ex: --rm_tau 17 17.5 18)')
    parser.add_argument('--kern_sz', type=int, required=False, default=3, help='Gaussian smoothing kernel size (default=3x3x3 voxels)')
    parser.add_argument('--mot_cor', action='store_true', help='Apply ANTs motion correction (reference frame = total frames / 2)')

    args = parser.parse_args()

    # Argument Parser checks
    if (args.img_seq == 'Other') and (args.yml_fpath is None):
        parser.error("If img_seq is Other, --yml_fpath required")

    # Required arguments
    img_fpath = args.img_fpath
    img_seq = args.img_seq
    out_fpath = args.out_fpath

    # Optional arguments
    yml_fpath = args.yml_fpath
    hct = args.hct
    qa=args.qa
    e2_fpath=args.e2_fpath
    b0_fpath = args.b0_fpath
    rm_tau = args.rm_tau
    kernsz = args.kern_sz
    motcor = args.mot_cor

    #
    proc_ase(img_fpath, img_seq, hct, out_fpath, qa=qa, e2_fpath=e2_fpath, b0_fpath=b0_fpath,
               yml_fpath=yml_fpath, rm_tau=rm_tau, kernsz=kernsz, motcor=motcor)