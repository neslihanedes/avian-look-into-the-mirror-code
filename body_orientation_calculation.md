The body orientation was analyzed with a time-series modeling technique called state-space model (Commandeur, & Koopman, 2007). 
We focused on dynamic fluctuation of body orientation. The state space modeling technique is a general scheme of time-series analysis (Commandeur, & Koopman, 2007; Matsui et al., 2018). 
The aim of time-series modeling is to deal with the dependency of each data point. 
The ordinary state-space model accommodates dependency by predicting the data at _t + 1_ depending on what they are at t:

> μ<sub>t+1</sub> ∼ Normal(μ<sub>t</sub>, σ<sub>μ</sub>)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(1)

where μ is the estimated average of the data. Intuitively, equation (1) estimates the overall temporal fluctuation across a session. The actual data includes the noise, which is postulated using an observation noise term,

> Y<sub>t</sub> ∼ Normal(μ<sub>t</sub>, σ<sub>Y</sub>)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(2)

where Y is the data (i.e., body orientation). 
Equations (1) and (2) are called “state model” and “observation model”, respectively. 
We assumed that the pigeons from all conditions follow the common state model. 
It was also assumed that the effect of the presentation of the stranger or mirror would start from a certain time point, if any. 
This assumption was necessary because all pigeons were involved in feeding, which produced highly aligned poses toward the feeder that did not largely differ. 
We modeled this assumption by modifying equation (2):

> Y<sub>t</sub> ∼ Normal(μ<sub>t</sub>, σ<sub>Y</sub>)(t ≤ t<sub>cp</sub>) Y<sub>t</sub>  ∼ Normal(μ<sub>t</sub> + δ<sub>condition</sub>, σ<sub>Y</sub>)(t > t<sub>cp</sub>)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(3)

where δ<sub>condition</sub> is the effect of the presentation of the stranger or mirror, and hence δ<sub>stranger</sub> or δ<sub>mirror</sub> in practice. 
The t<sub>cp</sub> is a change point, which is a starting time of difference of δ<sub>condition</sub>. 
We estimated these parameters and the difference was evaluated by checking whether Bayesian credible intervals overlapped with 0. 
The parameters were estimated via MCMC sampling software STAN (Carpenter et al., 2017) with a configuration of iteration 400000, burn-in periods 350000, 4 chains, and 100 thinning between samples.


## References
*  Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B., Betancourt, M., ... & Riddell, A. (2017). Stan: a probabilistic programming language. Grantee Submission, 76(1), 1-32.
*  Commandeur, J. J., & Koopman, S. J. (2007). An introduction to state space time series analysis. Oxford University Press.
*  Matsui, H., Yamada, K., Sakagami, T., & Tanno, T. (2018). Modeling bout–pause response patterns in variable-ratio and variable-interval schedules using hierarchical Bayesian methodology. Behavioural processes, 157, 346-353.

<img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1">