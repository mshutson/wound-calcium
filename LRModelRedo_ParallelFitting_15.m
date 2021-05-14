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
{dataAll,data2ndExp,woundSizes}=Import["SecondExpansionData.m"][[15]];
scaledData=data2ndExp/.{t_?NumericQ,r_?NumericQ}->{t/t0,r/rd};
(*Add to data extra points for fitting*)
maxRads=MaximalBy[scaledData,Last][[1,2]];
\[CapitalDelta]t=scaledData[[2,1]]-scaledData[[1,1]];
scaledDataFitting=Join[scaledData,Table[{scaledData[[-1,1]]+i*\[CapitalDelta]t,maxRads},{i,Length@scaledData}]];
(*Set max time to just be a little past this final data point. No need to go farther in time for fitting. Can look at farther time if inspecting particular parameter sets*)
\[Tau]maxes=(1.1*scaledDataFitting[[;;,1]]//Last);
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
getRad//ClearAll
getRad[{\[CapitalGamma]DP_,\[CapitalGamma]DL_,\[CapitalGamma]d_,\[CapitalGamma]c_,\[CapitalGamma]P_,\[CapitalGamma]pL_,\[CapitalGamma]L_,\[CapitalGamma]R_,\[CapitalGamma]w_}]:=getRad[{\[CapitalGamma]DP,\[CapitalGamma]DL,\[CapitalGamma]d,\[CapitalGamma]c,\[CapitalGamma]P,\[CapitalGamma]pL,\[CapitalGamma]L,\[CapitalGamma]R,\[CapitalGamma]w}]=
Block[{params,soln,rStep,points,interps,rads,maxes},
params=getParams[{\[CapitalGamma]DP,\[CapitalGamma]DL,\[CapitalGamma]d,\[CapitalGamma]c,\[CapitalGamma]P,\[CapitalGamma]pL,\[CapitalGamma]L,\[CapitalGamma]R,\[CapitalGamma]w}];

soln=NDSolveValue[deqns/.params,z,{r,dr,\[Rho]Far},{t,0,\[Tau]maxes},Method->{"IndexReduction"->Automatic,"EquationSimplification"->"Residual","PDEDiscretization"->{"MethodOfLines","SpatialDiscretization"->{"TensorProductGrid","MinPoints"->200,"MaxPoints"->500,"DifferenceOrder"->2}}}];

(*Check if maximum ever gets above threshold*)
maxes=((\[Gamma]L*soln[dr,#])/(1+\[Gamma]L*soln[dr,#]))/.params&/@times;

(*Find Ligand radius v. time points. Replace negative/ no solution cases with 0, and use a better (but longer) solving method if points are bad (just in case)*)
rads=
Table[
{times[[i]],If[maxes[[i]]<0.5,0.,r/.FindRoot[(((\[Gamma]L*soln[r,times[[i]]])/(1+\[Gamma]L*soln[r,times[[i]]]))/.params)==.5,{r,dr,\[Rho]Far}]]},
{i,Length@times}];

Return[rads];
]
(*Look at initial guess. Chose a best fit parameter set to data set 1*)
initialGuess=Import["initialGuesses_15.m"][[ToExpression[Environment["SLURM_ARRAY_TASK_ID"]]]];
initialGuessFitting=Table[{#[[i]],initialGuess[[i]]},{i,Length@#}]&@{\[CapitalGamma]DP,\[CapitalGamma]DL,\[CapitalGamma]d,\[CapitalGamma]c,\[CapitalGamma]P,\[CapitalGamma]pL,\[CapitalGamma]L,\[CapitalGamma]R,\[CapitalGamma]w};
(*Failed expression for when errors are thrown*)
failedExpression=scaledDataFitting/.{x_,y_}->{x,0};
(*getRadFit function is similar to getRad function, for the later time points it sets the output equal to the constraint data if it is below it so that there is no peanalty for going below the constraint*)
getRadFit//ClearAll
getRadFit[{\[CapitalGamma]DP_,\[CapitalGamma]DL_,\[CapitalGamma]d_,\[CapitalGamma]c_,\[CapitalGamma]P_,\[CapitalGamma]pL_,\[CapitalGamma]L_,\[CapitalGamma]R_,\[CapitalGamma]w_}]:=getRadFit[{\[CapitalGamma]DP,\[CapitalGamma]DL,\[CapitalGamma]d,\[CapitalGamma]c,\[CapitalGamma]P,\[CapitalGamma]pL,\[CapitalGamma]L,\[CapitalGamma]R,\[CapitalGamma]w}]=
Catch[
TimeConstrained[
Block[{params,soln,rStep,points,interps,rads,rads1,maxes},
params=getParams[{\[CapitalGamma]DP,\[CapitalGamma]DL,\[CapitalGamma]d,\[CapitalGamma]c,\[CapitalGamma]P,\[CapitalGamma]pL,\[CapitalGamma]L,\[CapitalGamma]R,\[CapitalGamma]w}];

soln=Check[

NDSolveValue[deqns/.params,z,{r,dr,\[Rho]Far},{t,0,\[Tau]maxes},Method->{"IndexReduction"->Automatic,"EquationSimplification"->"Residual","PDEDiscretization"->{"MethodOfLines","SpatialDiscretization"->{"TensorProductGrid","MinPoints"->200,"MaxPoints"->500,"DifferenceOrder"->2}}}],

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
{times[[i]],If[i>Length[times]/2&&rads1[[i]]<=maxRads,maxRads,rads1[[i]]]}
,{i,Length@times}];

Return[rads];
],
4,
Throw[failedExpression]
]
]

fittingForm//ClearAll
fittingForm[{\[CapitalGamma]DP_?NumericQ,\[CapitalGamma]DL_?NumericQ,\[CapitalGamma]d_?NumericQ,\[CapitalGamma]c_?NumericQ,\[CapitalGamma]P_?NumericQ,\[CapitalGamma]pL_?NumericQ,\[CapitalGamma]L_?NumericQ,\[CapitalGamma]R_?NumericQ,\[CapitalGamma]w_?NumericQ},time_?NumericQ]:=getRadFit[{\[CapitalGamma]DP,\[CapitalGamma]DL,\[CapitalGamma]d,\[CapitalGamma]c,\[CapitalGamma]P,\[CapitalGamma]pL,\[CapitalGamma]L,\[CapitalGamma]R,\[CapitalGamma]w}][[IntegerPart[(time-First@times)/\[CapitalDelta]t+1],2]]
timing=AbsoluteTiming[
fit=
NonlinearModelFit[
scaledDataFitting,

{
fittingForm[{\[CapitalGamma]DP,\[CapitalGamma]DL,\[CapitalGamma]d,\[CapitalGamma]c,\[CapitalGamma]P,\[CapitalGamma]pL,\[CapitalGamma]L,\[CapitalGamma]R,\[CapitalGamma]w},time],
-3.<=\[CapitalGamma]DP<=3.,
.01<=\[CapitalGamma]DL<=1.1,
-4.<=\[CapitalGamma]d<=4.,
-4.<=\[CapitalGamma]c<=4.,
-4.<=\[CapitalGamma]P<=4.,
-4.<=\[CapitalGamma]pL<=4.,
-4.<=\[CapitalGamma]L<=3.,
-4.<=\[CapitalGamma]R<=4.,
1.<=\[CapitalGamma]w<=5.
},

initialGuessFitting,

time,

Method->{NMinimize,Method->"NelderMead"}
];
];
bestFit={\[CapitalGamma]DP,\[CapitalGamma]DL,\[CapitalGamma]d,\[CapitalGamma]c,\[CapitalGamma]P,\[CapitalGamma]pL,\[CapitalGamma]L,\[CapitalGamma]R,\[CapitalGamma]w}/.fit["BestFitParameters"];

finalPlot=Show[ListPlot[scaledDataFitting,PlotRange->All],ListLinePlot[getRadFit[bestFit],PlotRange->All,PlotStyle->Red]];

Export["finalPlot_"<>Environment["SLURM_ARRAY_TASK_ID"]<>".png",finalPlot];
Export["fit_"<>Environment["SLURM_ARRAY_TASK_ID"]<>".m",fit];
Export["bestParams_"<>Environment["SLURM_ARRAY_TASK_ID"]<>".m",bestFit];
Export["timing_"<>Environment["SLURM_ARRAY_TASK_ID"]<>".m",timing];
