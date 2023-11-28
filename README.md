# OrganDeformationModels
MATLAB-implementation of different statistical models for organ motion and deformation commonly used in radiotherapy research. These models are based on organ contours given as surface point clouds, and require that deformable contour registration has already been performed. 

The four implemented methods are:

Patient specific model:
The patient specific organ deformation based on principal component analysis was first presented in 
Söhn M, Birkner M, Yan D, Alber M. Modelling individual geometric variation based on dominant eigenmodes of organ deformation: implementation and evaluation. Phys Med Biol. 2005;50(24):5893-5908. doi:10.1088/0031-9155/50/24/009

Population model:
The population model was first presented in
Budiarto E, Keijzer M, Storchi PRM, Heemink AW, Breedveld S, Heijmen BJM. Computation of mean and variance of the radiotherapy dose for PCA-modeled random shape and position variations of the target. Phys Med Biol. 2013;59(2):289. doi:10.1088/0031-9155/59/2/289

And finally, two Bayesian models:
The Normal-Inverse-Wishart model and the Variational Bayes model were both introduced in 
Rørtveit ØL, Hysing LB, Stordal AS, Pilskog S. An organ deformation model using Bayesian inference to combine population and patient-specific data. Phys Med Biol. 2023;68(5):055009. doi:10.1088/1361-6560/acb8fc

The data set that was used to generate the results in the latter paper is [available for download from dataverse.NO here](https://dataverse.no/dataset.xhtml?persistentId=doi:10.18710/DKVPIJ).


Note 1: The parameter "nu" for the two Bayesian models is interpreted slightly different from its use in the paper: Here, nu needs to be greater than P+1, where P is the number of (intra-patient) principal components used. In the paper, this was implicit, so setting nu = P+n for some n in this code is equivalent to setting nu = n in the equations in the paper.

Note 2: The "psi" matrix, or scale matrix for the inverse Wishart distribution, is scaled such that it is in the same scale as the covariance matrix; i.e. psi here is the expected value of the inverse Wishart distribution. This is consistent throughout the code, and goes for both the NIW and variational Bayes model. This differs from the standarad representation of psi by a factor of sqrt(1/(nu-P-1)).

Note 3: There are different ways to generate samples from a patient; either through point estimates of mu and R (as in the paper) or through the posterior predictive distribution. The difference is discussed in the document [On random sampling of patient data from Bayesian models.pdf](https://github.com/oyvindlr/OrganDeformationModels/blob/main/On%20random%20sampling%20of%20patient%20data%20from%20Bayesian%20models.pdf).

Note 4: The code includes functions to estimate the parameters nu and kappa for the NIW-model using maximum likelihood, both when using posterior predictive and point estimates. For the variational Bayes model, no such code is available. 

See example.m for an example of how to use the code.

