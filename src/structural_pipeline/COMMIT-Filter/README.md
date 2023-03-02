# Documentation COMMIT-Filter

The COMMIT-Filter add-on adds a post-processing pipeline step to the structural pipeline, which filters streamlines in the fiber cloud using [COMMIT2](https://github.com/daducci/COMMIT) and reconstructs connectivity matrices from the resulting filtered fiber cloud. The current integration (in `COMMIT_script.py`) is based on the [COMMIT2 example](https://github.com/daducci/COMMIT/wiki/COMMIT2) provided on the COMMIT wiki.

## Installation
Before using this add-on, ensure that COMMIT and its dependencies are installed. You can install COMMIT by following the instructions on the  [COMMIT wiki](https://github.com/daducci/COMMIT/wiki/Installation).


## Usage

To use the COMMIT-Filter add-on, follow these steps:


1. Add the following items to your configuration file `CATO.conf` to filter the fibercloud constructed using constrained spherical deconvolution (`csd`) using information about fiber bundles form the Desikan-Killiany (`aparc`) parcellation.

```
"commit_filter":{
  "reconstructionMethods":[
     "csd"
  ], 
  "templates":[
     "aparc"
  ], 
  "outputCommitDir":"OUTPUTDIR/commit", 
  "filteredFiberFile":"OUTPUTCOMMITDIR/SUBJECT_fibers_commit_METHOD.trk", 
  "filteredFiberPropertiesFile":"OUTPUTCOMMITDIR/SUBJECT_fiber_properties_commit_METHOD_TEMPLATE.mat", 
  "schemeFile":"OUTPUTCOMMITDIR/dwi.scheme", 
  "intermediateConnectomeFile":"OUTPUTCOMMITDIR/connectome.csv", 
  "pythonInterpreter":"/usr/bin/python", 
  "commitScriptFile":"TOOLBOXDIR/structural_pipeline/COMMIT-filter/COMMIT_script.py", 
  "fiberWeightsFile":"OUTPUTCOMMITDIR/Results_StickZeppelinBall_COMMIT2/streamline_weights.txt", 
  "filteredConnectivityMatrixFile":"OUTPUTCOMMITDIR/SUBJECT_connectivity_commit_METHOD_TEMPLATE.mat", 
  "lambda":5E-4
}
```

**Note:** If you installed COMMIT and its dependencies in a specific Andaconda environment, update the `pythonInterpreter` parameter. For example, if you are using Anaconda and the environment is called COMMIT use:


```
"pythonInterpreter":"/Applications/anaconda3/envs/COMMIT/bin/python", 
```

**Note:** The suitable value of the regularization term `lambda` depends on the dataset. Update this parameter in accordance with your data.

2. Add the `commit_filter` step to the pipeline steps by editing the `reconstructionSteps` parameter in the configuration file:
```
"reconstructionSteps":{ 
  "structural_preprocessing": true,
  "parcellation": true,
  "collect_region_properties": true,
  "reconstruction_diffusion": true,
  "reconstruction_fibers": true,
  "reconstruction_fiber_properties": true,
  "reconstruction_network": true,
  "commit_filter": true
},
```

3. Run the structural pipeline with the additional `commit_filter` step:


```
subjectDir = '/Volumes/Example/0001';
configurationFile = '/Volumes/Example/CATO.conf';
structural_pipeline(subjectDir, ...
    'configurationFile', configurationFile);
```

This will execute COMMIT and provide the COMMIT output, filtered fiber clouds and filtered connectivity matrices in the directory `OUTPUTDIR/commit`.












