# MGWI-R-Code

This repository contains R codes for generating and fitting stationary and non-stationary MGWI processes.

* For the non-stationary case, we consider trend and seasonal covariates in the regression structure as follows:

![formula](https://render.githubusercontent.com/render/math?math=\mu_t%20=%20\exp(\beta_0%20%2B%20\beta_1%20t/n%20%2B%20\beta_2%20\cos(2%20\pi%20t%20/%2012))%20\quad%20\text{and}%20\quad%20\alpha_t%20=%20\exp(\gamma_0%20%2B%20\gamma_1%20t/n))

* The above setup aims to mimic realistic situations when dealing with epidemic diseases, for example. It can be edited as long as the rest of the code is accordingly adapted.  
