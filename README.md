# Computational Modelling for Diffusion MRI

Diffusion MRI is a special type of MRI where the image intensity is sensitive to dispersion of water molecules. 
This method is highly beneficial in research as the cellular architecture has a specific water dispersion patterns. By modelling
the diffusion patterns from images, we can gain insights into the tissue structure. In this project, I investigated different 
parametric computational models for understanding diffusion patterns of water in brain MRIs, which could be used to describe 
the nerve fibre architecture in the brain.

## Dataset

This project used the High Angular Resolution Imaging (HARDI) dataset which was
acquired as part of the Human Connectome Project (http://www.humanconnectome.org/). It has 90 diffusion weighted MRIs
with b-value = 1000 s/mm<sup>2</sup> and 18 MRIs with b-value=0.

## Structure

This project is executed in MATLAB R2019b. The repo has the following structure:

```bash
├── Figures_for_report.pdf
├── images
├── parametric_models.mlx
├── README.md
├── Report.pdf
└── scripts
    ├── BallStickEval.m
    ├── BallStickSSD.m
    ├── BallStickSSD_transformation.m
    ├── createParameterMaps.m
    ├── findBestFit.m
    ├── generateDataSet.m
    ├── MCMC.m
    ├── PBootstrap.m
    ├── ZeppelinStickSSD.m
    └── ZeppelinStickTortuositySSD.m
```
This main code is in ``parametric_models.mlx``. The indivitual functions for the models 
can be accessed in the `scripts/` directory. A final report of what I did and my results
can be found in `Report.pdf` and `Figures_for_report.pdf`.

## Models Investigated

Experiments were performed with 4 types of models.

1. Diffusion tensor model - This is a simple model where a 7-dimensional parameter vector is computed, given the 
   voxel intensities in the image volume, the b-values describing a diffusion weighting factor 
   and the q-hat vectors which describe the gradient directions in the diffusion image. The problem is essentially a least
   squares problem of the form ``A = Gx`` - A is the voxel intensities, x is the parameter vector and G is a design matrix
   of the values of b and qhat.
   
2. Ball and Stick model - This is two-compartment model where I model water diffusion outside the nerve fibres and inside the
   nerves separately (hence two compartments).Water molecules outside the nerve bundles are slightly more free hence we call this
   a "ball" component while inside the axons, the water movement is more restricted, which can described by a
   "stick" model.
   
2. Zeppelin-stick model - This is another two-compartment model like a Ball and Stick, except now instead of a 
   Ball component, a Zeppelin component is used. The zeppelin component has the shape of an oblate spheroid, which 
   has a primary eigenvector describing the principal direction and two other minor eigenvectors with equal eigenvalues.
   
3. Zeppelin-stick and Tortuosity - A variation of the Zeppelin-stick where the principal eigenvalue is a function of the second eigenvalue.

## Results

The models above were compared using two metrics - Akaike's information content (AIC) and Bayesian information criterion (BIC).
These metrics use information about the model such as the number of parameters N (encodes the complexity of the model) and the log likelihood of the data under these 
model parameters (encodes the errors of the model). The results are as follows:

Model | AIC | BIC
---|---|---|
Diffusion tensor | 3.8727e+03 | 1.4164e+04
Ball and Stick  | -1.9771e+04 | -9.4913e+03
Zeppelin Stick | -2.0975e+04 | -1.0690e+04
Zeppelin Stick and Tortuosity | -2.0723e+04| -1.0444e+04

The Zeppelin stick model has the smallest AIC and BIC, suggesting that this model best describes water diffusion within the brain.

## Extensions

Possible extensions to this project would be to model diffusion patterns using machine learning models such as Deep neural networks or
Convolutional neural networks. 
