#! /usr/bin/env python
# -*- coding: utf-8 -*-




import os
import argparse
import nibabel as nib
import numpy as np

DESCRIPTION =   'Joint correction of eddy currents and gradient non linearity correction using FSL tools. Cornelius Eichner 2018'

def buildArgsParser():
    p = argparse.ArgumentParser(description=DESCRIPTION)
    p.add_argument('--in', dest='data', action='store', type=str,
                            help='Path of the input volume (nifti format)')

    p.add_argument('--bvec', dest='bvec', action='store', type=str,
                            help='Path of the bvec file (default \'bvec\')')

    p.add_argument('--bval', dest='bval', action='store', type=str,
                           help='Path of the bval file (default \'bval\')')

    p.add_argument('--mask', dest='mask', action='store', type=str,
                            help='Path of the brain mask')

    p.add_argument('--acqp', dest='acqp', action='store', type=str,
                        help='Path of eddy acquisition parameter file (default \'acq_param.txt\')')
    
    p.add_argument('--index', dest='index', action='store', type=str,
                            help='Path of the eddy index file (default \'index.txt\')')

    p.add_argument('--topup', dest='topup', action='store', type=str,
                           help='Base name of topup file structure (default \'topup/topup\')')

    p.add_argument('--out', dest='out', action='store', type=str,
                            help='Path of the output volume')

    p.add_argument('--mp', dest='openmp', action='store', type=int, default='1', 
                            help='Optional: OpenMP parallelization factor (default 8)')

    return p


def main():
    parser = buildArgsParser()
    args = parser.parse_args()

    if args.bvec is None:
        BVEC = 'bvec'
    else:
        BVEC = os.path.realpath(args.bvec)

    if args.bval is None:
        BVAL = 'bval'
    else:
        BVAL = os.path.realpath(args.bval)

    if args.acqp is None:
        ACQP = 'acq_param.txt'
    else:
        ACQP = os.path.realpath(args.acqp)

    if args.index is None:
        INDEX = 'index.txt'
    else:
        INDEX = os.path.realpath(args.index)

    if args.topup is None:
        TOPUP = 'topup/topup'
    else:
        TOPUP = os.path.realpath(args.topup)

    if args.acqp is None:
        OPENMP_THREADS = 8
    else:
        OPENMP_THREADS = args.openmp


    PATH = os.path.dirname(os.path.realpath(args.data)) + '/'
    os.chdir(PATH)

    DATA = os.path.realpath(args.data)
    MASK = os.path.realpath(args.mask)
    
    OUT = args.out


    ######################
    # Calculate Gradient Non Linearities

    # Extract first data volume for correction
    cmd = 'fslroi ' + DATA + ' ' + PATH + 'single.nii.gz 0 1'
    os.system(cmd)

    cmd = '~/.local/bin/gradient_unwarp.py '+ PATH + 'single.nii.gz ' + PATH + 'nlgc.nii.gz siemens -g /home/raid/ceichner/Software/gradunwarp_data/coeffConnectom.grad -n'
    os.system(cmd)
    
    cmd = 'mv fullWarp_abs.nii.gz ' + PATH + 'nlgc_warp.nii.gz'
    os.system(cmd)

    """
    Resulting output files:
    nlgc.nii.gz           first volume, corrected for gradient non linearities, will not be further employed
    fullWarp_abs.nii.gz   warp field for gradient correction
    """


    ######################
    # Calculate eddy displacement fields

    # Set number of parallel processing threads
    cmd = 'set OMP_NUM_THREADS = ' + str(OPENMP_THREADS)
    os.system(cmd)

    # Run eddy from home directory installation
    cmd =\
    '/home/raid/ceichner/Software/eddy/eddy_openmp \
        --imain='+ DATA + ' \
        --mask=' + MASK + ' \
        --index=' + INDEX + ' \
        --acqp=' + ACQP + ' \
        --bvecs=' + BVEC + ' \
        --bvals=' + BVAL + ' \
        --out=' + PATH + 'eddy \
        --topup=' + TOPUP + ' \
        --residuals \
        --repol \
        --dfields \
        --data_is_shelled \
        -v'
    os.system(cmd)

    # Move all displacement field files in subdirectory
    os.system('mkdir -p dfield')
    os.system('mv *displacement_fields* dfield/')
    DFIELD_FILES = sorted(os.listdir('dfield/'))


    ######################
    # Combine both warp fields 

    os.system('mkdir -p comb_warp')

    # This code combines both warp fields by simple addition.
    '''
    for i in DFIELD_FILES:
        cmd = 'fslmaths nlgc_warp.nii.gz -add dfield/' + i + ' comb_warp/nlcg.' + i
        os.system(cmd)
    '''

    # Combination of warpfields by concatenation

    for i in DFIELD_FILES:
        cmd = 'convertwarp \
                    -o comb_warp/nlcg.' + i + '\
                    -r nlgc_warp.nii.gz \
                    --warp1=dfield/' + i + ' \
                    --warp2=nlgc_warp.nii.gz \
                    -v'
        os.system(cmd)


    COMB_WARP_FILES = sorted(os.listdir('comb_warp/'))


    ######################
    # Apply warp for using combined field

    os.system('mkdir -p split_data')
    os.system('mkdir -p corr_data')

    cmd = 'fslsplit ' + PATH + '/eddy.eddy_outlier_free_data.nii.gz split_data/data -t'
    os.system(cmd)

    SPLIT_DATA_FILES = sorted(os.listdir('split_data/'))


    for i in xrange(len(os.listdir('split_data/'))):
        cmd = 'applywarp \
                -i split_data/' + SPLIT_DATA_FILES[i] + ' \
                -r split_data/' + SPLIT_DATA_FILES[0] + ' \
                -o corr_data/' + SPLIT_DATA_FILES[i] + ' \
                -w comb_warp/' + COMB_WARP_FILES[i] + ' \
                --interp=spline \
                --datatype=float'
        os.system(cmd)

    # Wait for al warps to be done
    while len(os.listdir('corr_data/')) < len(os.listdir('split_data/')):
        pass

    cmd = 'fslmerge -t ' + OUT + ' corr_data/*'
    os.system(cmd)



if __name__ == '__main__':
    main()

