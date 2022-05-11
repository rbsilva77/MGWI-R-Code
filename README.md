# MGWI-R-Code

This repository contains R codes for generating and fitting stationary and non-stationary MGWI processes.

$\bullet$ For the non-stationary case, we consider trend and seasonal covariates in the regression structure as follows:

$$\mu_t = \exp(\beta_0 + \beta_1 t/n + \beta_2 \cos(2 \pi t / 12)) \quad \text{and} \quad \alpha_t = \exp(\gamma_0 + \gamma_1 t/n),$$ for $t = 1, \ldots, n$.

$\bullet$ The above setup aims to mimic realistic situations when dealing with epidemic diseases, for example. It can be edited as long as the rest of the code is accordingly adapted.  
