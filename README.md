# wound-calcium
The files in this repository are Mathematica notebooks (.nb) and package files (.m). These notebooks were written to model wound-induced calcuim responses in epithelia that follow from an upstream protease cleavage of a pro-ligand to release active ligand that then binds to a G-protein coupled receptor. The models are described in full in the publication </BR></BR>
  "A protease-initiated model of wound detection"</BR>
  James T. O’Connor, Aaron C. Stevens, Erica K. Shannon, Fabiha Bushra Akbar, Kimberly S. LaFever, Neil Narayanan, M. Shane Hutson, Andrea Page-McCaw</BR>
  https://www.biorxiv.org/content/10.1101/2020.12.08.415554v1
<BR/>
<BR/>
Main Files<BR/>
<UL>
<LI>SecondExpansionData.m: Data that was used for fitting</LI>
<LI>LRModelRedo_RandomInitialConditions.m: Notebook (it is an m file in order to run on ACCRE) that was used to randomly determine sets of initial guesses for the fitting</LI>
<LI>initialGuesses_*.m: Initial guesses used for fitting data set *.</LI>
<LI>LRModelRedo_ParallelFitting_15.m: Notebook (it is an m file in order to run on ACCRE) that did the fitting. This was specifically for dataset 15. Fits to other sets are practically identical, so I did not include them</LI>
<LI>Analysis_AllFits_06.nb: Notebook that did the analysis of the fits</LI>
<LI>SplitExpression.nb: Notebook to analyze the split-expression system model (Mthl10)</LI>
</UL>

Supplementary Files<BR/>

<UL>
<LI>LRModelRedo_ParallelFitting_SameParameters.m: Notebook (it is an m file in order to run on ACCRE) that attempted simultaneous fitting to all data sets with the same parameter set</LI>
<LI>inits_SameParameters.m: Initial guesses for the simultaneous-same parameter fitting</LI>
<LI>LRModelRedo_ParallelFitting_woundScaled.m: Notebook (it is an m file in order to run on ACCRE) that attempted simultaneous fitting to all data sets with a similar set of parameters with wound-dependent parameters that scaled with wound size</LI>
<LI>inits_woundScaled.m: Initial guesses for the simultaneous wound scaled fitting</LI>
</UL>
