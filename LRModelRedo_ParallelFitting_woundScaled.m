(* ::Package:: *)

(* ::Input::Initialization:: *)
(*Run to Launch all Kernels for parallel computation in fitting*)
LaunchKernels[16];

rd=10.;(* \[Mu]m. spatial scale.*)
t0=47; (*s. median time for 2nd expansion signal to start from Erica's paper. *)
DL=260; (*\[Mu]m^2/s. Estimated free diffusion constant of GBP*)
LRth =.5; (*Dimensionless. Ratio of [L.R] Subscript[to [R], T] needed for signal to be "on"*)
dr=.001; (*"small" r since we cannot use the origin in polar coordinates, \[Mu]m*)
\[Rho]Far=1000/rd (*r that is "infinity" for numeric solving. scaled*);
(*Differential Equations*)
s=\[Gamma]d/\[Gamma]w^2*Exp[-(r-dr)^2/(2*\[Gamma]w^2)-\[Gamma]d*t];
dx=D[x[r,t],t]==(1+(\[Gamma]pL*yT[r,t])/(1+\[Gamma]P*x[r,t])^2)^-1*(\[Gamma]DP*Laplacian[x[r,t],{r,\[Theta]},"Polar"]+s+\[Gamma]c*\[Gamma]pL/\[Gamma]P*((\[Gamma]P*x[r,t])/(1+\[Gamma]P*x[r,t]))^2*yT[r,t]);
dyT=D[yT[r,t],t]==(-\[Gamma]c*\[Gamma]P*x[r,t]*yT[r,t])/(1+\[Gamma]P*x[r,t]);
dz=D[z[r,t],t]==(1+\[Gamma]R/(1+\[Gamma]L*z[r,t])^2)^-1*(\[Gamma]DL*Laplacian[z[r,t],{r,\[Theta]},"Polar"]+(\[Gamma]c*\[Gamma]P*x[r,t]*yT[r,t])/(1+\[Gamma]P*x[r,t]));


(*Initial Conditions and Boundary Conditions*)
initx=x[r,0]==0.; (*Initially no active protease*)
bc1x=(D[x[r,t],r]==0)/.r->dr; (*Radially symmetric --> no flux through origin*)
bc2x=x[\[Rho]Far,t]==0.; (*Concentration to 0 far away*)

inityT=yT[r,0]==1.; (*proligand initially exists at all points in space*)

initz=z[r,0]==0; (*Initially no ligand*)
bc1z=(D[z[r,t],r]==0)/.r->dr; (*Radially symmetric --> no flux through origin*)
bc2z=z[\[Rho]Far,t]==0.; (*Concentration to 0 far away*)

(*Put it all together*)
deqnsXY={dx,dyT,initx,bc1x,bc2x,inityT};
deqnsZ={dz,initz,bc1z,bc2z};
deqns=Join[deqnsXY,deqnsZ];

(*secondExpansionData is a list of elements for each data set {calcium radius at all times, calcium radius during second expansion, wound sizes*)
data2ndExp=Import["SecondExpansionData.m"][[{1,15,17,19},2]];
woundSize=Import["SecondExpansionData.m"][[{1,15,17,19},3]];
scaledData=data2ndExp/.{t_?NumericQ,r_?NumericQ}->{t/t0,r/rd};

(*Add to data extra points for fitting*)
maxRad=MaximalBy[#,Last][[1,2]]&/@scaledData;
\[CapitalDelta]t=(#[[2,1]]-#[[1,1]])&/@scaledData;
scaledDataFitting=Table[Join[scaledData[[set]],Table[{scaledData[[set,-1,1]]+i*\[CapitalDelta]t[[set]],maxRad[[set]]},{i,Length@scaledData[[set]]}]],{set,Length@scaledData}];

(*Set max time to just be a little past this final data point. No need to go farther in time for fitting. Can look at farther time if inspecting particular parameter sets*)
\[Tau]max=(1.1*#[[;;,1]]//Last)&/@scaledDataFitting;
times=scaledDataFitting[[;;,;;,1]];(*Only need to pick out time points corresponding to the data*)

wScale=woundSize/woundSize[[1]];

getOA[gen_,numRuns_,numLevels_]:=
Block[{genRows=Log[numLevels,numRuns],ga,x},
ga=Table[RotateRight[gen,i],{i,0,genRows-1}];
x=Tuples[Range[0,numLevels-1],genRows];
Return[Mod[#.ga&/@x,numLevels]]
]

getOA[gen_,numRuns_,numLevels_,const_]:=
Block[{genRows=Log[numLevels,numRuns],ga,x},
ga=Table[Append[RotateRight[gen,i],const],{i,0,genRows-1}];
x=Tuples[Range[0,numLevels-1],genRows];
Return[Mod[#.ga&/@x,numLevels]]
]

testOA[oa_,numLevels_,strength_]:=
Block[{columnCheck,digitCheck,maxN=numLevels^strength-1},
columnCheck=Select[Subsets[Range[Length[oa[[1]]]],{strength}],OrderedQ];
digitCheck=Flatten[Table[ConstantArray[Range[0,maxN][[i]],Length[oa]/numLevels^strength],{i,maxN+1}]];

Return[Table[Sort[FromDigits[#,numLevels]&/@oa[[;;,columnCheck[[i]]]]]==digitCheck,{i,Length[columnCheck]}]]
]
oaInit=getOA[{2,2,1,2,2,2,1,1,1,2,1},3^6,3,1];
And@@testOA[oaInit,3,5];
oa=(oaInit/1.)[[;;,;;9]]/.{(2.->-1.)}; (*Divide by 1. so that data type changes. Take first 9 columns for 9 parameters. Change 2. to -1.*)

getParams[{\[CapitalGamma]DP_,\[CapitalGamma]DL_,\[CapitalGamma]d_,\[CapitalGamma]c_,\[CapitalGamma]P_,\[CapitalGamma]pL_,\[CapitalGamma]L_,\[CapitalGamma]R_,\[CapitalGamma]w_}]:=
{
\[Gamma]DP->12.22*10^#[[1]] ,(*Estimated from scales and DP 12.22 \[TildeTilde] .1DL*)
\[Gamma]DL->t0/rd^2*DL*#[[2]], (*known value*)
\[Gamma]d->10^#[[3]] ,
\[Gamma]c->10^#[[4]],
\[Gamma]P->10^#[[5]],
\[Gamma]pL->10^#[[6]],
\[Gamma]L->10^#[[7]],
\[Gamma]R->10^#[[8]],
\[Gamma]w->#[[9]]
}&[{\[CapitalGamma]DP,\[CapitalGamma]DL,\[CapitalGamma]d,\[CapitalGamma]c,\[CapitalGamma]P,\[CapitalGamma]pL,\[CapitalGamma]L,\[CapitalGamma]R,\[CapitalGamma]w}];
undoParams[rules_]:=
{
Log[10,\[Gamma]DP/12.22],
rd^2/(t0*DL)*\[Gamma]DL,
Log[10,\[Gamma]d],
Log[10,\[Gamma]c],
Log[10,\[Gamma]P],
Log[10,\[Gamma]pL],
Log[10,\[Gamma]L],
Log[10,\[Gamma]R],
\[Gamma]w
}/.rules;

(*Pick random guesses from predetermined initial guesses*)
initialGuess=Import["bestInits.m"][[ToExpression[Environment["SLURM_ARRAY_TASK_ID"]]]];

getRadFit//ClearAll;
getRadFit[{\[CapitalGamma]DP_,\[CapitalGamma]DL_,\[CapitalGamma]d_,\[CapitalGamma]c_,\[CapitalGamma]P_,\[CapitalGamma]pL_,\[CapitalGamma]L_,\[CapitalGamma]R_,\[CapitalGamma]w_},set_]:=getRadFit[{\[CapitalGamma]DP,\[CapitalGamma]DL,\[CapitalGamma]d,\[CapitalGamma]c,\[CapitalGamma]P,\[CapitalGamma]pL,\[CapitalGamma]L,\[CapitalGamma]R,\[CapitalGamma]w},set]=
Block[{params=getParams[{\[CapitalGamma]DP,\[CapitalGamma]DL,\[CapitalGamma]d,\[CapitalGamma]c,\[CapitalGamma]P,\[CapitalGamma]pL,\[CapitalGamma]L,\[CapitalGamma]R,\[CapitalGamma]w}],soln,rStep,points,interps,rads,maxes},

soln=NDSolveValue[deqns/.{\[Gamma]w->\[Gamma]w*wScale[[set]],\[Gamma]P->\[Gamma]P*wScale[[set]]^2}/.params,z,{r,dr,\[Rho]Far},{t,0,\[Tau]max[[set]]},Method->{"IndexReduction"->Automatic,"EquationSimplification"->"Residual","PDEDiscretization"->{"MethodOfLines","SpatialDiscretization"->{"TensorProductGrid","MinPoints"->200,"MaxPoints"->500,"DifferenceOrder"->2}}}];

(*Check if maximum ever gets above threshold*)
maxes=((\[Gamma]L*soln[dr,#])/(1+\[Gamma]L*soln[dr,#]))/.params&/@times[[set]];

(*Find Ligand radius v. time points. Replace negative/ no solution cases with 0, and use a better (but longer) solving method if points are bad (just in case)*)
rads=
Table[If[maxes[[i]]<0.5,0.,r/.FindRoot[(((\[Gamma]L*soln[r,times[[set,i]]])/(1+\[Gamma]L*soln[r,times[[set,i]]]))/.params)==.5,{r,dr,\[Rho]Far}]],
{i,Length@times[[set]]}];

points=Table[{times[[set,i]],If[rads[[i]]<maxRad[[set]]&&i>(Length@times[[set]])/2,maxRad[[set]],rads[[i]]]},{i,Length@times[[set]]}];

Return[points];
];

(*Function to get sum of square distance between corresponding data and simulation points*)
getSumSquare[sim_]:=Total[Flatten[Table[(sim[[set,;;,2]]-scaledDataFitting[[set,;;,2]])^2,{set,Length@scaledDataFitting}]]];

(*Function that modifies runs before testing them by modifiying columns of the OA so that we do not go out of bounds of the region we want to test in parameter space*)
runBound[trial_,{\[CapitalGamma]DPB_,\[CapitalGamma]DLB_,\[CapitalGamma]dB_,\[CapitalGamma]cB_,\[CapitalGamma]PB_,\[CapitalGamma]pLB_,\[CapitalGamma]LB_,\[CapitalGamma]RB_,\[CapitalGamma]wB_}]:=
Block[{row=trial, bounds={\[CapitalGamma]DPB,\[CapitalGamma]DLB,\[CapitalGamma]dB,\[CapitalGamma]cB,\[CapitalGamma]PB,\[CapitalGamma]pLB,\[CapitalGamma]LB,\[CapitalGamma]RB,\[CapitalGamma]wB}},
Do[
If[row[[param]]>bounds[[param,2]],row[[param]]=bounds[[param,2]]];
If[row[[param]]<bounds[[param,1]],row[[param]]=bounds[[param,1]]];
,{param,Length[row]}];
Return[row];
]

findFitOA[{\[CapitalGamma]DP_?NumberQ,\[CapitalGamma]DL_?NumberQ,\[CapitalGamma]d_?NumberQ,\[CapitalGamma]c_?NumberQ,\[CapitalGamma]P_?NumberQ,\[CapitalGamma]pL_?NumberQ,\[CapitalGamma]L_?NumberQ,\[CapitalGamma]R_?NumberQ,\[CapitalGamma]w_?NumberQ},{tPStart_,tPEnd_}]:=
Reap[
Block[
{bestFit={\[CapitalGamma]DP,\[CapitalGamma]DL,\[CapitalGamma]d,\[CapitalGamma]c,\[CapitalGamma]P,\[CapitalGamma]pL,\[CapitalGamma]L,\[CapitalGamma]R,\[CapitalGamma]w},
bounds=
{{-4.,0.},
{.01,1.1},
{-4.,4.},
{-4.,4.},
{-4.,4.},
{-4,4.},
{-4.,3.},
{-4.,4.},
{.1,5.}},
foundBest=False,
runs,
sds},


Do[
While[!foundBest,
(*Distribute everything (need function definitions) to make parallel computing faster. This might already be taken care of in the ParallelTable function, but just in case*)
DistributeDefinitions["Global`"];

runs=ParallelTable[getRadFit[runBound[bestFit+(oa/2.^\[CapitalDelta]d)[[run]],bounds],dataSet],{run,Length@oa},{dataSet,Length@scaledDataFitting}];

(*Obtain sum of square distances*)
sds=getSumSquare[#]&/@runs;
bestFit=runBound[bestFit+(oa/2.^\[CapitalDelta]d)[[Ordering[sds][[1]]]],bounds];
Sow[bestFit,fitLabel];
Sow[Min[sds],sdsLabel];

(*Check if best fit is the original parameter list (OA all 0's row)*)
If[Ordering[sds][[1]]==1,foundBest=True];
];

foundBest=False;
,{\[CapitalDelta]d,tPStart,tPEnd}];
bestFit
]
]

ParallelEvaluate[Off[NDSolveValue::eerr]];
ParallelEvaluate[Off[General::munfl]];
ParallelEvaluate[Off[FindRoot::brmp]];

timing=AbsoluteTiming[fit=findFitOA[initialGuess,{1.,4.}];];

finalPlot=Show[ListPlot[scaledDataFitting,PlotRange->All],ListLinePlot[Table[getRadFit[fit[[1]],set],{set,Length@scaledDataFitting}],PlotRange->All]]

Export["finalPlot_"<>Environment["SLURM_ARRAY_TASK_ID"]<>".png",finalPlot];
Export["fit_"<>Environment["SLURM_ARRAY_TASK_ID"]<>".m",fit];
Export["timing_"<>Environment["SLURM_ARRAY_TASK_ID"]<>".m",timing];

CloseKernels[];
