(* ::Package:: *)

(* ::Input::Initialization:: *)
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
{dataAllAll,data2ndExpAll,woundSizesAll}=Transpose@(Import["SecondExpansionData.m"]);
scaledDataAll=data2ndExpAll/.{t_?NumericQ,r_?NumericQ}->{t/t0,r/rd};
scaledDataAllAll=dataAllAll/.{t_?NumericQ,r_?NumericQ}->{t/t0,r/rd};
data2ndExp=data2ndExpAll[[ToExpression[Environment["SLURM_ARRAY_TASK_ID"]]]];
scaledData=scaledDataAll[[ToExpression[Environment["SLURM_ARRAY_TASK_ID"]]]];

(*Add to data extra points for fitting*)
maxRad=MaximalBy[scaledData,Last][[1,2]];
\[CapitalDelta]t=(scaledData[[2,1]]-scaledData[[1,1]]);
scaledDataFitting=Join[scaledData,Table[{scaledData[[-1,1]]+i*\[CapitalDelta]t,maxRad},{i,Length@scaledData}]];

(*Set max time to just be a little past this final data point. No need to go farther in time for fitting. Can look at farther time if inspecting particular parameter sets*)
\[Tau]max=(1.1*scaledDataFitting[[;;,1]]//Last);
times=scaledDataFitting[[;;,1]];(*Only need to pick out time points corresponding to the data*)

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
failedExpression=scaledDataFitting/.{x_,y_}->{x,0};
getRadFit//ClearAll
getRadFit[{\[CapitalGamma]DP_,\[CapitalGamma]DL_,\[CapitalGamma]d_,\[CapitalGamma]c_,\[CapitalGamma]P_,\[CapitalGamma]pL_,\[CapitalGamma]L_,\[CapitalGamma]R_,\[CapitalGamma]w_}]:=
Catch[
TimeConstrained[
Block[{params,soln,rStep,points,interps,rads,rads1,maxes},
params=getParams[{\[CapitalGamma]DP,\[CapitalGamma]DL,\[CapitalGamma]d,\[CapitalGamma]c,\[CapitalGamma]P,\[CapitalGamma]pL,\[CapitalGamma]L,\[CapitalGamma]R,\[CapitalGamma]w}];

soln=Check[

NDSolveValue[deqns/.params,z,{r,dr,\[Rho]Far},{t,0,\[Tau]max},Method->{"IndexReduction"->Automatic,"EquationSimplification"->"Residual","PDEDiscretization"->{"MethodOfLines","SpatialDiscretization"->{"TensorProductGrid","MinPoints"->200,"MaxPoints"->500,"DifferenceOrder"->2}}}],

Throw[failedExpression],

{NDSolveValue::eerr,NDSolveValue::ndsz}
];

(*Check if maximum ever gets above threshold*)
maxes=Check[((\[Gamma]L*soln[dr,#])/(1+\[Gamma]L*soln[dr,#]))/.params&/@times,Throw[failedExpression],{Power::infy}];

(*Find Ligand radius v. time points. Replace negative/ no solution cases with 0, and use a better (but longer) solving method if points are bad (just in case)*)
rads1=
Table[
If[maxes[[i]]<0.5,0.,r/.Check[FindRoot[(((\[Gamma]L*soln[r,times[[i]]])/(1+\[Gamma]L*soln[r,times[[i]]]))/.params)==.5,{r,dr,\[Rho]Far}],Throw[failedExpression],{FindRoot::cvmit,FindRoot::nlnum,FindRoot::brmp,Power::infy}]],
{i,Length@times}];

(*Replace all points in the latter half with the constraint data value if they are below it*)
rads=Table[
{times[[i]],If[i>Length[times]/2&&rads1[[i]]<=maxRad,maxRad,rads1[[i]]]}
,{i,Length@times}];

Return[rads];
],
5,
Throw[failedExpression]
]
]
getSSD[sim_]:=Total[(sim-scaledDataFitting[[;;,2]])^2]
Quiet@
Block[{count=0.,randomGuess,rads},

initialSets=Reap[
While[count<32,

randomGuess=RandomReal/@{
{-3.,1.},
{.01,1.0},
{-3.,3.},
{-3.,3.},
{-3.,3.},
{-3.,3.},
{-4.,3.},
{-4.,1.},
{.1,5.}
};

rads=getRadFit[randomGuess][[;;,2]];

If[getSSD[rads]<=getSSD[failedExpression[[;;,2]]]/10.&&AllTrue[rads[[(Length@times)/2+1;;Length@times]],#<=maxRad&],
Sow[randomGuess];
count=count+1;
];

]
][[2,1]]
];

Export["initialGuesses_"<>Environment["SLURM_ARRAY_TASK_ID"]<>".m",initialSets];

