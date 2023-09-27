# Summer-project-on-time-dependent-infectivity-and-detectability

The result is summarized in 3 different jupyter notebook files. two of which concerning the change in infectibility with time, and one concerning change in detectability. 

The jupyter notebook "approximating logit" outlines the method inwhich logistic disease progression can be approximated with multiple box model, and show that a different set of parameter can be fitted from the epidemic thus generated to create an simpler SEI model
where the prediction is not visually different from the actual, more complex model.

The note book titled "forward fitting approach" shows that when a disease progression is known, the resulting epidemic is not significantly different between the multiple box model with accurate approximation of the logistic curve and the cruder approximation using
an monomolecular curve. This contradicts previous finding that I have presented before. This is because previously I have used the same step function heuristic for the fitting of the simpler two box model. However, this is not required, and direct optimization is possible. Therefore, in this notebook, a better approximation with monomolecular curve is achieved, and the prediction are not visibly different between the simpler two box (SEI) model and the more complex multiple box model.

Therefore, monomolecular model is probably sufficient for modelling changes with infectivity through time. However, this will require the rate parameter for the movement between boxes to be fitted with the whole course of disease progression, instead of the observation of non infective period. For example, for Huang-Long-Bing (Bassanezi & Bassanezi 2008), where the severity and hence presumably infectivity increases after the onset of symptoms, instead of letting the rate parameter to be the inverse of latent period, it should be fitted to the entire curve, which includes the non-infectious initial period and the later gradual increase in infectivity. This thus represent a conceptual change from the E boxes being a "latent period" box, it should be conceptualized as merely the first part of a two box infectivity model. In the two notebooks, I have only looked at the case when the "real" change is logistic, where more complicated curves may be obtained from experiments, in which case it should be assessed on a case by case basis.

The notebook titled "different culling radius" concerns the change of detectability with time and their implication on controlling method. The original idea of using an truncated logistic distribution for the latent period is discarded since the expected latent period becomes variable. The notebook thus includes two models, first with SCI model with latent period gamma distributed with a fixed mean and variable scale parameter. The second model have the probability of detection of an individual gradually increasing, following an logistic curve with fixed inflection point but variable rate of increase. In both cases, the optimal strategy remains in the same region, but the expected final epidemic size (the number of removed tree when all infected individual is removed) is expected to be higher when the variance of time of detection is lower, represented in the first case by smaller rate parameter (variance of gamma distribution = mean * rate) and in the second case a sharper increase in detectability around inflection point. 


Among the repositories there are also julia codes that are affiliated with the notebook "different culling radius" and some graphs plotted from it.


citation
 Bassanezi, R.B. & Bassanezi, R.C. 2008 An approach to model the impact of Huanglongbing on citrus yield. In Proceedings of the International Research Conference on Huanglongbing (eds. T.R. Gottwald & J.H. Graham). Orlando, Florida, Plant Management Network. https://swfrec.ifas.ufl.edu/hlb/database/pdf/22_IRCHLB_08.pdf 
