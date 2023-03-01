# Script to run COMMIT and COMMIT2 on CATO output.
# This script is based on the COMMIT wiki: https://github.com/daducci/COMMIT/wiki/COMMIT2

import numpy as np
import os
import commit
from commit import trk2dictionary
import argparse

argParser = argparse.ArgumentParser()

argParser.add_argument("--dwiProcessedFile", help="preprocessed DWI file", required=True)
argParser.add_argument("--dwiSchemeFile", help="DWI scheme file", required=True)
argParser.add_argument("--fiberFile", help="Input fiber file", required=True)
argParser.add_argument("--outputCommitDir", help="Output directory for COMMIT.", required=True)
argParser.add_argument("--intermediateConnectomeFile", help="Connectome file with number of streamline per connection (for internal use)", required=True)
argParser.add_argument("--subjectDir", help="Subject directory", required=True)
args = argParser.parse_args()

# Import tractogram
# TODO: check fiber_shift
os.chdir(args.subjectDir)

trk2dictionary.run(
    filename_tractogram = args.fiberFile,
    path_out = 'DWI_processed_test/COMMIT',
    fiber_shift         = 0
)

# load the data
mit = commit.Evaluation( '.', '.' )
mit.load_data(args.dwiProcessedFile, args.dwiSchemeFile)

# use a forward-model with 1 Stick for the streamlines and 2 Balls for all the rest
mit.set_model( 'StickZeppelinBall' )
d_par       = 1.7E-3             # Parallel diffusivity [mm^2/s]
d_perps_zep = []                 # Perpendicular diffusivity(s) [mm^2/s]
d_isos      = [ 1.7E-3, 3.0E-3 ] # Isotropic diffusivity(s) [mm^2/s]
mit.model.set( d_par, d_perps_zep, d_isos )

mit.generate_kernels( regenerate=True )
mit.load_kernels()

# create the sparse data structures to handle the matrix A
mit.load_dictionary( 'DWI_processed_test/COMMIT' )
mit.set_threads()
mit.build_operator()

# perform the fit
mit.fit( tol_fun=1e-3, max_iter=1000, verbose=False )
mit.save_results( path_suffix="_COMMIT1" )

# Retrieve the streamline contributions estimated by COMMIT for later use in COMMIT2:
x_nnls, _, _ = mit.get_coeffs( get_normalized=False )

# Preparing the anatomical prior on bundles
C = np.loadtxt( args.intermediateConnectomeFile, delimiter=',' )
C = np.triu( C ) # be sure to get only the upper-triangular part of the matrix
group_size = C[C>0].astype(np.int32)

tmp = np.insert( np.cumsum(group_size), 0 , 0)
group_idx = np.array( [np.arange(tmp[i],tmp[i+1]) for i in range(len(tmp)-1)], dtype=np.object_ )

group_w = np.empty_like( group_size, dtype=np.float64 )
for k in range(group_size.size) :
    group_w[k] = np.sqrt(group_size[k]) / ( np.linalg.norm(x_nnls[group_idx[k]]) + 1e-12 )

# TODO: examine best lambda setting.
reg_lambda = 5e-4 # change to suit your needs

prior_on_bundles = commit.solvers.init_regularisation(
    mit,
    regnorms    = [commit.solvers.group_sparsity, commit.solvers.non_negative, commit.solvers.non_negative],
    structureIC = group_idx,
    weightsIC   = group_w,
    lambdas     = [reg_lambda, 0.0, 0.0]
)

mit.fit( tol_fun=1e-3, max_iter=1000, regularisation=prior_on_bundles, verbose=False )
mit.save_results( path_suffix="_COMMIT2" )


