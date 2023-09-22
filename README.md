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

See example.m for an example of how to use the code.

