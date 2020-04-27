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

bsgamma.py : Code to produce Branching ratio of B > Xs + gamma(photon) based on NLO Wilson coefficients from matching scale (mu_w) to lower scale (mu_b).
> Functions are based on paper PHYSICAL REVIEW D 58 074004 and formulas (3,4,5,6,7,8) from paper: arXiv:hep-ph/9803368.
>> It generates basic functions which requires to put i_{1,2} = Y and j_{1,2} = XY^* for each function.
>> require Y and XY^* arrays from exercise.py if possible, and load to A_cp(s2,s1,mass1,mass2,i1,j1,i2,i2) function to create final CP-asymmetry limitation.
> The final function produces BR(B > Xs + gamma(photon)) and CP-asymmetry. Both functions require 8 input parameters:
> s2: low scale. (ie. mb); s1: high scale.(ie. mw); mass1: first MH+ range; mass2: second MH+ range.
> i1,j1: i1 = Y array , j1 = XY^* array for the first MH+ contribution.
> i2,j2: i2 = Y array , j2 = XY^* array for the second MH+ contribution.
>>> Set a parameter called: mass_differ = first MH+ and second MH+ mass difference.

bsgonly.py : Plot functions or codes to show dependence between input parameters of H1+,H2+ masses, scale values(mu_i, i = mw/mb/mhch), mixing matrix parameters (theta, tanbeta,tangamma,delta) against BR(B_bar > X_s + gamma(photon)), CP-asymmetry (B_bar > X_s gamma, X_d gamma), CP-asymmetry difference and Untagg-asymmetry formula (B_bar > X_{s+d} gamma) results.
> several plots based on papers: arXiv:hep-ph/9803368 and arXiv:1605.05881.
> read b_sgamma.py functions and import exercise.py to plot dependence between (theta, tanbeta,tangamma,delta) and BR(B_bar > X_s + gamma(photon)).
>> Do same procedure in exercise.py first, and input specific type-3HDM needed, then b_sgammaplot.py will read through exercise.py to get X,Y,Z_{2,3}, complexfunction (X_2Y_2^* ), complex3function (X_3Y_3^* ) arrays. Using array4() function to obtain input (theta, tanbeta,tangamma,delta) to get result of BR(B_bar > X_s + gamma(photon)) in function called: 
plot_under_Heatherbasis().


gammaeffective.py : anomalous dimension matrix. It includes LO and NLO together.
>> Requires to put this with b>sgamma.py together.

CEPCeehh.py: Codes for CEPC (more energy available range. ie. Centre of mass energy = 240 GeV) to produce charged Higgs from eehh.py program.
> Choose specific integrated luminosity value as input.
>> combine with CEPC_plot.py (produces graphs at CEPC centre of mass energy).
>> Function to plot significances of H+ against e_c tagging values and invriant mass cut values.

CEPC_plot.py: Codes to produce graphs of branching ratios of charged Higgs decay (Leptonic and Hadronic).
>> CEPCeehh.py produces charged Higgs events based on chosen luminoscities. Using the events that code produced can plot for parameter space against branching ratios and other stuffs.

neutronedm.py: Codes to produce 3HDM charged Higgs contribution for Neutron EDM limit based on the (X_iY_i^* ) where i = 2,3 reading values from exercise.py. 
>> Calculation taken from mu_tH(m_t) scale to mu_hadron(hadron = 1GeV) scale.
>> Two contributions are mainly useful, (C)3HDM and weinberg operator. Whole formulas are based on the paper :  arXiv:1308.6283.

edmonly.py: Codes to plot Neutron-EDM and Electron-EDM based on the neutronedm.py file. The effective parameters are from charged Higgs masses H+1, H+2, and four mixing parameters (theta, tanbeta,tangamma,delta ).

Hplus3decay.py : Codes for H+> AW* / hW* 3 body decay. Based on the SUSY Higgs particle paper: https://arxiv.org/pdf/hep-ph/9511342.pdf. 
>> Results produce the partial width of charged Higgs H+ decay to Pseudoscalar and off-shell W boson which will decays further to 2 decays (Thus, 3-body decay released from charged Higgs).
