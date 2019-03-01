# 3HDM-charged-Higgs
light charge Higgs decay products calculation codes
Also PDF file of paper

invariant.py : All possible invariant values and parameters which can be used for other codes.
>> Need to put with others togehter. 

eehh.py :  Code to produce charged Higgs from e+ e- collider approach based on LEP2 results. 
> With b-tagging which can produce signal / sqrt(background) of specific jets cases(2 or 4 jets)
>> simply by choose e_c (miss tagging of c quark from b quark) single or list.

exercise.py : Code to produce BR( H+> cb, cs, taunv ) based on Yukawa couplings X and Y and Z. 
> Additional, to select 4 fundamnetal parameters of Mixing matrix by different types of 3HDM and plot the branching ratio of charged Higgs decay.
> The signal / sqrt(background) of 2jets (cb/cs + taunv) and 4 jets (cb/cs + cs/cb ) cases  plots are included.
>> This is connected with eehh.py, so after select e_c as single, then a choice for specific parameter for specific type 3HDM is required.

b_sgamma.py : Code to produce Branching ratio of B_bar > Xs + gamma(photon) based on NLO Wilson coefficients from matching scale (mu_w) to lower scale (mu_b).
> Functions are based on paper PHYSICAL REVIEW D 58 074004. 
>> It generates basic functions which requires to put i = Y^2 and j = XY^* for each function.

gammaeffective.py : anomalous dimension matrix. It includes LO and NLO together.
>> Requires to put this with b>sgamma.py together.



