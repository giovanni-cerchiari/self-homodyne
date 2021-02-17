(* ::Package:: *)

(*
Copyright: Giovanni Cerchiari, Lorenzo Dania, Dmitry Bykov
e-mail: giovanni.cerchiari@uibk.ac.at
date : 02/2021
*)
(*This file is a free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This file is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the other repository files.
    If not, it can be found at <https://www.gnu.org/licenses/>.
*)
(*
-----------------------------------------------------------------------
This file contains the script-program that has been used to evaluat
the results presented in the article:

G. Cerchiari, L. Dania, D. S. Bykov, R. Blatt, T. Northup
"Position measurement of a dipolar scatterer via self-homodyne detection"
arXiv

If you wish to cite this work, we prepared a citation file "selfhomodyne.bib"
in bibtex format in the repository.
*)
(*
-----------------------------------------------------------------------
The script is written for Mathematica 11.3 .
It executes in approximately 20 min on a machine with the following specifications
- operating system : Windows 10
- processor : Intel(R) Core(TM) i5-7300U CPU @ 2.60GHz 2.71 GHz
- RAM : 16 GB
-----------------------------------------------------------------------
*)
(*---------------------------------------------------------------------------------------------*)
(*---------------------------------------------------------------------------------------------*)
(*---------------------------------------------------------------------------------------------*)
ClearAll["Global`*"]
SetCoordinates[Spherical]
(*Generic direction in space*)
nr[\[Theta]_,\[Phi]_]:={Cos[\[Phi]]*Sin[\[Theta]],Sin[\[Phi]]*Sin[\[Theta]],Cos[\[Theta]]};
Print["Generic direction in space (\[Theta],\[Phi]) = ", MatrixForm[nr[\[Theta],\[Phi]]]]
(*---------------------------------------------------------------------------------------------*)
(*---------------------------------------------------------------------------------------------*)
(*RADIATED POWER PER UNIT SOLID ANGLE*)
(*differential scattered power into the solid andgle d\[CapitalOmega] = Sin[\[Theta]]d\[Theta]d\[Phi] *)
(*eq. 6 of the main text*)
(*Power radiated by a dipole:
- Pdip integral power of the dipole
- \[Theta]0,\[Phi]0 direction of linear polarization
- \[Theta],\[Phi] direction scattered light*)
(*Please remember the Sin[\[Theta]] missing factor at integration*)
(*linear polarization general case*)
dpdipd\[CapitalOmega][\[Theta]_,\[Phi]_,\[Theta]0_,\[Phi]0_]:=3/(8*\[Pi])*Pdip*(1-(nr[\[Theta],\[Phi]].nr[\[Theta]0,\[Phi]0])^2);
Print["(dpdipd\[CapitalOmega]) differential power radiate by a dipole
 with linear polarization: ",dpdipd\[CapitalOmega][\[Theta],\[Phi],\[Theta]p,\[Phi]p]]


(*----------------------------------------------------------------------------------------------*)
(*----------------------------------------------------------------------------------------------*)
(*MEASUREMENT IMPRECISION: IDEAL AND SELF-HOMODYNE. LINEAR POLARIZATION*)
(*The label "ideal" in the names of the variables refers to the ideal configuration described in
--
Felix Tebbenjohanns, Martin Frimmer, and Lukas Novotny
"Optimal position detection of a dipolar scatterer in a focused field"
Phys. Rev. A 100, 043821 \[Dash] Published 14 October 2019
url: https://link.aps.org/doi/10.1103/PhysRevA.100.043821
--
In this configuration the reference field for homodyne detection is a second dipolar
field located in the origin*)
(*----------------------------------------------------------------------------------------------*)
(*The label "mirror" in the names of the variables refers to the self-homodyne method*)
(*----------------------------------------------------------------------------------------------*)
(*The imprecision s and S are described in the article Eqs 10 and 11: 
s = differential power spectral density of the imprecision noise associated with a differential 
    detector under direction (\[Theta], \[Phi])
S = integral of s in the angular integral \[Phi]=[0, 2\[Pi]] and \[Theta]=[0,\[Theta]]. Note the free parameter \[Theta].
With SS and ss natural units are introduced*)
(*Assumptions: shot noise dominated by the constant term*)
(*-------------------*)
(*differential imprecision for the ideal configuration*)
sideal[\[Theta]_,\[Phi]_,\[Theta]p_,\[Phi]p_,\[Theta]0_,\[Phi]0_]:= Evaluate[Simplify[
	(hbar*c/\[Lambda])*\[Gamma]^2/( (2*\[Gamma]*(2*\[Pi]/\[Lambda])*(nr[\[Theta]0,\[Phi]0].nr[\[Theta],\[Phi]]))^2 *dpdipd\[CapitalOmega][\[Theta],\[Phi],\[Theta]p,\[Phi]p]),
	Assumptions->{0<=\[Theta]<=\[Pi], 0<=\[Phi]<2*\[Pi]}]];
Print["sideal = ", sideal[\[Theta],\[Phi],\[Theta]p,\[Phi]p,\[Theta]0,\[Phi]0]]
(*differential imprecision for the self-homodyne configuration*)
smirror[\[Theta]_,\[Phi]_,\[Theta]p_,\[Phi]p_,\[Theta]0_,\[Phi]0_]:= Evaluate[Simplify[
	(hbar*c/\[Lambda])*(1+R^2)/( (4*R*(2*\[Pi]/\[Lambda])*(nr[\[Theta]0,\[Phi]0].nr[\[Theta],\[Phi]]))^2 *dpdipd\[CapitalOmega][\[Theta],\[Phi],\[Theta]p,\[Phi]p]),
	Assumptions->{0<=\[Theta]<=\[Pi], 0<=\[Phi]<2*\[Pi]}]];
Print["smirror = ", smirror[\[Theta],\[Phi],\[Theta]p,\[Phi]p,\[Theta]0,\[Phi]0]]
(*Imprecision for the ideal configuration*)
Sideal[\[Theta]_,\[Theta]p_,\[Phi]p_,\[Theta]0_,\[Phi]0_]:= Evaluate[Simplify[
	(Integrate[Sin[\[Theta]\[Theta]]/sideal[\[Theta]\[Theta],\[Phi],\[Theta]p,\[Phi]p,\[Theta]0,\[Phi]0],{\[Theta]\[Theta],0,\[Theta]},{\[Phi],0,2*\[Pi]}])^-1,
	Assumptions->{0<=\[Theta]<=\[Pi]/2}]];
Print["Sideal = ", Sideal[\[Theta],\[Theta]p,\[Phi]p,\[Theta]0,\[Phi]0]]
(*Imprecision for the self-homodyne configuration*)
Smirror[\[Theta]_,\[Theta]p_,\[Phi]p_,\[Theta]0_,\[Phi]0_]:=Evaluate[Simplify[
	(Integrate[Sin[\[Theta]\[Theta]]/smirror[\[Theta]\[Theta],\[Phi],\[Theta]p,\[Phi]p,\[Theta]0,\[Phi]0],{\[Theta]\[Theta],0,\[Theta]},{\[Phi],0,2*\[Pi]}])^-1,
	 Assumptions->{0<=\[Theta]<=\[Pi]/2}]];
Print["Smirror = ", Smirror[\[Theta],\[Theta]p,\[Phi]p,\[Theta]0,\[Phi]0]]
(*Enforce reference system with light polarization in the xz-plane and introducing natural units.*)
SSmirror[\[Theta]_,\[Theta]p_,\[Theta]0_,\[Phi]0_]:= Evaluate[Simplify[Smirror[\[Theta],\[Theta]p,0,\[Theta]0,\[Phi]0],
	Assumptions->{0<=\[Theta]<=\[Pi]/2, 0<=\[Phi]<2*\[Pi], 0<=\[Theta]0<=\[Pi], 0<=\[Phi]0<2*\[Pi], hbar==1, \[Lambda]==1, c==1, R==1, Pdip==1}]];
ssmirror[\[Theta]_,\[Phi]_,\[Theta]p_,\[Theta]0_,\[Phi]0_]:= Evaluate[FullSimplify[smirror[\[Theta],\[Phi],\[Theta]p,0,\[Theta]0,\[Phi]0],
	Assumptions->{0<=\[Theta]<=\[Pi]/2, 0<=\[Phi]<2*\[Pi], 0<=\[Theta]0<=\[Pi], 0<=\[Phi]0<2*\[Pi], hbar==1, \[Lambda]==1, c==1, R==1, Pdip==1}]];


(*-------------------------------------------------------------------------------------------------*)
(*-------------------------------------------------------------------------------------------------*)
(*EQUIVALENCE SELF-HOMODYNE AND IDEAL CONFIGURATIONS*)
sidealovers = Simplify[sideal[\[Theta],\[Phi],\[Theta]p,\[Phi]p,\[Theta]0,\[Phi]0]/smirror[\[Theta],\[Phi],\[Theta]p,\[Phi]p,\[Theta]0,\[Phi]0]];
Print["s(ideal)/s(self-homodyne) = ", sidealovers];
SidealoverS= Simplify[Sideal[\[Pi],\[Theta]p,\[Phi]p,\[Theta]0,\[Phi]0]/Smirror[\[Pi]/2,\[Theta]p,\[Phi]p,\[Theta]0,\[Phi]0]];
Print["S(ideal)/S(self-homodyne) = ", SidealoverS];
(*-------------------------------------------------------------------------------------------------*)
(*-------------------------------------------------------------------------------------------------*)
(*SYMMETRY OF IMPRECISION*)
Print["Is the total imprecision symmetric in direction: s[\[Theta]\[Epsilon],\[Theta]0,\[Phi]0]-s[\[Theta]\[Epsilon],\[Pi]-\[Theta]0,\[Pi]+\[Phi]0]=? (yes == 0): ", Simplify[ssmirror[\[Theta],\[Phi],\[Theta]p,\[Theta]0,\[Phi]0]-ssmirror[\[Theta],\[Phi],\[Theta]p,\[Pi]-\[Theta]0,\[Pi]+\[Phi]0]]]
Print["Is the total imprecision symmetric in direction: s[\[Theta]\[Epsilon],\[Theta]0,\[Phi]0]-s[\[Pi]-\[Theta]\[Epsilon],\[Theta]0,\[Pi]+\[Phi]0]=? (yes == 0): ", Simplify[ssmirror[\[Theta],\[Phi],\[Theta]p,\[Theta]0,\[Phi]0]-ssmirror[\[Theta],\[Phi],\[Pi]-\[Theta]p,\[Theta]0,\[Pi]+\[Phi]0], Assumptions->{0<=\[Theta]<=\[Pi], 0<=\[Phi]<=2*\[Pi]}]]
Print["Is the total imprecision symmetric in direction: S[\[Theta]\[Epsilon],\[Theta]0,\[Phi]0]-S[\[Theta]\[Epsilon],\[Pi]-\[Theta]0,\[Pi]+\[Phi]0]=? (yes == 0): ", Simplify[SSmirror[\[Theta],\[Theta]p,\[Theta]0,\[Phi]0]-SSmirror[\[Theta],\[Theta]p,\[Pi]-\[Theta]0,\[Pi]+\[Phi]0]]]
Print["Is the total imprecision symmetric in direction: S[\[Theta]\[Epsilon],\[Theta]0,\[Phi]0]-S[\[Pi]-\[Theta]\[Epsilon],\[Theta]0,\[Pi]+\[Phi]0]=? (yes == 0): ", Simplify[SSmirror[\[Theta],\[Theta]p,\[Theta]0,\[Phi]0]-SSmirror[\[Theta],\[Pi]-\[Theta]p,\[Theta]0,\[Pi]+\[Phi]0]]]
Print["Is the total imprecision symmetric in direction: S[\[Pi]/2,\[Theta]0,\[Phi]0]-S[\[Pi]/2,\[Pi]-\[Theta]0,\[Phi]0]=? (yes == 0): ", Simplify[SSmirror[\[Theta],\[Pi]/2,\[Theta]0,\[Phi]0]-SSmirror[\[Theta],\[Pi]/2,\[Pi]-\[Theta]0,\[Phi]0]]]
Print["Is the total imprecision symmetric in direction: S[\[Pi]/2,\[Theta]0,\[Phi]0]-S[\[Pi]/2,\[Theta]0,\[Pi]+\[Phi]0]=? (yes == 0): ", Simplify[SSmirror[\[Theta],\[Pi]/2,\[Theta]0,\[Phi]0]-SSmirror[\[Theta],\[Pi]/2,\[Pi]-\[Theta]0,\[Phi]0]]]
Print["Is the total imprecision symmetric in direction: S[\[Pi]/2,\[Theta]0,\[Phi]0]-S[\[Pi]/2,\[Theta]0,-\[Phi]0]=? (yes == 0): ", Simplify[SSmirror[\[Theta],\[Pi]/2,\[Theta]0,\[Phi]0]-SSmirror[\[Theta],\[Pi]/2,\[Theta]0,-\[Phi]0]]]
Print["Is the total imprecision symmetric in direction: S[\[Pi]/2,\[Theta]0,\[Phi]0]-S[\[Pi]/2,\[Theta]0,\[Pi]-\[Phi]0]=? (yes == 0): ", Simplify[SSmirror[\[Theta],\[Pi]/2,\[Theta]0,\[Phi]0]-SSmirror[\[Theta],\[Pi]/2,\[Theta]0,\[Pi]-\[Phi]0]]]
Print["Is the total imprecision symmetric in direction: S[\[Pi]/2,\[Theta]0,\[Phi]0]-S[\[Pi]/2,\[Theta]0,\[Pi]+\[Phi]0]=? (yes == 0): ", Simplify[SSmirror[\[Theta],\[Pi]/2,\[Theta]0,\[Phi]0]-SSmirror[\[Theta],\[Pi]/2,\[Theta]0,\[Pi]+\[Phi]0]]]
(*-------------------------------------------------------------------------------------------------*)
(*-------------------------------------------------------------------------------------------------*)
(*CONFIGURATION OF MINIMAL DIFFERENTIAL IMPRECISION*)
(*Maximizing the inverse imprecision is equivalent to finding the
 minimum of imprecision, but easier to solve*)
sol\[Theta]=Maximize[1/ssmirror[0,0,\[Theta]p,\[Theta]0,\[Phi]0],{\[Theta]p,\[Theta]0,\[Phi]0}]
\[Theta]pmax = Simplify[Mod[\[Theta]p/.sol\[Theta][[2]],\[Pi]]];
\[Theta]0max = Simplify[Mod[\[Theta]0/.sol\[Theta][[2]],\[Pi]]];
\[Phi]0max = Simplify[Mod[\[Phi]0/.sol\[Theta][[2]],2*\[Pi]]];
Print["Minimal imprecision conditions:"]
Print["S[\[Theta]pmax, \[Theta]0max, \[Phi]0max] = ", sol\[Theta][[1]]]
Print["\[Theta]pmax = ", \[Theta]pmax, " = ", N[\[Theta]pmax]*(180/\[Pi])," deg"]
Print["{\[Theta]0max, \[Phi]0max} = ", {\[Theta]0max, \[Phi]0max}," rad = ", N[{\[Theta]0max, \[Phi]0max}*(180/\[Pi])]]
Print["{\[Theta]0max, \[Phi]0max} = ", Simplify[{Mod[\[Pi]-\[Theta]0max,\[Pi]], Mod[\[Pi]+\[Phi]0max,2*\[Pi]]}],
	" rad = ", N[{Mod[\[Pi]-\[Theta]0max,\[Pi]], Mod[\[Pi]+\[Phi]0max,2*\[Pi]]}*(180/\[Pi])]]


(*----------------------------------------------------------------------------------------------*)
(*----------------------------------------------------------------------------------------------*)
(*MEASUREMENT IMPRECISION: SELF-HOMODYNE WITH CIRCULARLY POLARIZED LIGHT*)
(*Circular polarization corresponds to the atomic \[Sigma] tranistions.*)
(*The expression of the polarization vector is x+-iz*)
(*Atomic \[Pi] transitions are already described by the linearly polarize calculation*)
(*The naming convention follows from the previous part of the script with the 
addition of the label "circ"*)
(*---------------------------------------------------------*)
(*circular polarization*)
polcirc[\[Theta]_,\[Phi]_]:= (nr[\[Pi]/2,\[Pi]/2]+I*nr[0,0])/Sqrt[2];
(*dipole scalar product for circularly polarized light*)
dipolecirc[\[Theta]_,\[Phi]_]:=Simplify[nr[\[Theta],\[Phi]].polcirc[\[Theta],\[Phi]]*nr[\[Theta],\[Phi]].Conjugate[polcirc[\[Theta],\[Phi]]],
 Assumptions->{0<=\[Theta]<=\[Pi], 0<=\[Phi]<2*\[Pi]}];
(*differential radiated power from circularly polarized light*)
dpdipd\[CapitalOmega]circ[\[Theta]_,\[Phi]_]:=Evaluate[Simplify[3/(8*\[Pi])*Pdip*(1-dipolecirc[\[Theta],\[Phi]]),
 Assumptions->{0<=\[Theta]<=\[Pi], 0<=\[Phi]<2*\[Pi]}]];
Print["(dpdipd\[CapitalOmega]) differential power radiate by a dipole with polarization \[Epsilon]=x+iz: ",dpdipd\[CapitalOmega]circ[\[Theta],\[Phi]]]
(*differential imprecision for self-homodyne with circular polarization*)
smirrorcirc[\[Theta]_,\[Phi]_,\[Theta]0_,\[Phi]0_]:=Evaluate[Simplify[
	(hbar*c/\[Lambda])*(1+R^2)/( (4*R*(2*\[Pi]/\[Lambda])*(nr[\[Theta]0,\[Phi]0].nr[\[Theta],\[Phi]]))^2 *dpdipd\[CapitalOmega]circ[\[Theta],\[Phi]]),
	Assumptions->{0<=\[Theta]<=\[Pi]/2, 0<=\[Phi]<2*\[Pi]}]];
Print["smirrorcir = ", smirrorcirc[\[Theta],\[Phi],\[Theta]0,\[Phi]0]]
(*Total imprecision for self-homodyne with circular polarization*)
Smirrorcirc[\[Theta]_,\[Theta]0_,\[Phi]0_]:=Evaluate[Simplify[
	(Integrate[Sin[\[Theta]\[Theta]]/smirrorcirc[\[Theta]\[Theta],\[Phi],\[Theta]0,\[Phi]0],{\[Theta]\[Theta],0,\[Theta]},{\[Phi],0,2*\[Pi]}])^-1,
	Assumptions->{0<=\[Theta]<=\[Pi]/2, 0<=\[Phi]<2*\[Pi]}]];
Print["Smirrorcirc = ", Smirrorcirc[\[Theta],\[Theta]0,\[Phi]0]]
(*Enforce reference system and introducing natural units*)
SSmirrorcirc[\[Theta]_,\[Theta]0_,\[Phi]0_]:= Evaluate[FullSimplify[Smirrorcirc[\[Theta],\[Theta]0,\[Phi]0],
	Assumptions->{0<=\[Theta]<=\[Pi]/2, 0<=\[Phi]<2*\[Pi], 0<=\[Theta]0<=\[Pi], 0<=\[Phi]0<2*\[Pi], hbar==1, \[Lambda]==1, c==1, R==1, Pdip==1}]];
ssmirrorcirc[\[Theta]_,\[Phi]_,\[Theta]0_,\[Phi]0_]:= Evaluate[Simplify[smirrorcirc[\[Theta],\[Phi],\[Theta]0,\[Phi]0],
	Assumptions->{0<=\[Theta]<=\[Pi], 0<=\[Phi]<2*\[Pi], 0<=\[Theta]0<=\[Pi], 0<=\[Phi]0<2*\[Pi], hbar==1, \[Lambda]==1, c==1, R==1, Pdip==1}]];


(*----------------------------------------------------------------------------------------------*)
(*----------------------------------------------------------------------------------------------*)
(*MEASUREMENT IMPRECISION: UNPOLARIZED*)
(*Atom integral emission considering \[Pi] and \[Sigma] transitions via the Wigner-Eckart theorem*)
(*These variables follows the naming convention already introduced and can be identified
 by the label "atom"*)
(*----------------------------------------------*)
(*differential radiated power*)
dpdipd\[CapitalOmega]atom[\[Theta]_,\[Phi]_]:=Evaluate[Simplify[3/(8*\[Pi])*Pdip*2/3, Assumptions->{0<=\[Theta]<=\[Pi], 0<=\[Phi]<2*\[Pi]}]];
Print["(dpdipd\[CapitalOmega]) differential power radiate by a dipole with polarization \[Epsilon]=x+iz: ",dpdipd\[CapitalOmega]atom[\[Theta],\[Phi]]]
(*differential imprecision for \[Pi] and \[Sigma] transitions*)
smirroratom[\[Theta]_,\[Phi]_,\[Theta]0_,\[Phi]0_]:=Evaluate[Simplify[
	(hbar*c/\[Lambda])*(1+R^2)/( (4*R*(2*\[Pi]/\[Lambda])*(nr[\[Theta]0,\[Phi]0].nr[\[Theta],\[Phi]]))^2 *dpdipd\[CapitalOmega]atom[\[Theta],\[Phi]]),
	Assumptions->{0<=\[Theta]<=\[Pi]/2, 0<=\[Phi]<2*\[Pi]}]];
Print["smirroratom = ", smirroratom[\[Theta],\[Phi],\[Theta]0,\[Phi]0]]
(*Total imprecision for \[Pi] and \[Sigma] transitions*)
Smirroratom[\[Theta]_,\[Theta]0_,\[Phi]0_]:=Evaluate[Simplify[
	(Integrate[Sin[\[Theta]\[Theta]]/smirroratom[\[Theta]\[Theta],\[Phi],\[Theta]0,\[Phi]0],{\[Theta]\[Theta],0,\[Theta]},{\[Phi],0,2*\[Pi]}])^-1,
	Assumptions->{0<=\[Theta]<=\[Pi]/2, 0<=\[Phi]<2*\[Pi]}]];
Print["Smirroratom = ", Smirroratom[\[Theta],\[Theta]0,\[Phi]0]]
(*Enforce reference system and introducing natural units*)
SSmirroratom[\[Theta]_,\[Theta]0_,\[Phi]0_]:= Evaluate[FullSimplify[Smirroratom[\[Theta],\[Theta]0,\[Phi]0],
	Assumptions->{0<=\[Theta]<=\[Pi]/2, 0<=\[Phi]<2*\[Pi], 0<=\[Theta]0<=\[Pi], 0<=\[Phi]0<2*\[Pi], hbar==1, \[Lambda]==1, c==1, R==1,Pdip==1}]];
ssmirroratom[\[Theta]_,\[Phi]_,\[Theta]0_,\[Phi]0_]:= Evaluate[Simplify[smirroratom[\[Theta],\[Phi],\[Theta]0,\[Phi]0],
	Assumptions->{0<=\[Theta]<=\[Pi], 0<=\[Phi]<2*\[Pi], 0<=\[Theta]0<=\[Pi], 0<=\[Phi]0<2*\[Pi], hbar==1, \[Lambda]==1, c==1, R==1,Pdip==1}]];


(*----------------------------------------------------------------------------------------------*)
(*----------------------------------------------------------------------------------------------*)
(*SELF-HOMDYNE DETECTION WITH A QUADRANT PHOTO DETECTOR*)
(*Calculating the inverse imprecision of the self-homodyne method obtained with a quadrant photo detector (QPD).
The derivation of the imprecision SmirrorQPD is described in the Appendix D in the main text:
\[Beta] = differential sensitivity to displacements along the direction (\[Theta]0, \[Phi]0) associated 
	with a differential detector under direction (\[Theta], \[Phi]) (eq. 8 in the main text)
B = integral of the differential sensitivity \[Beta] in the angular interval corresponding to a quadrant sector
	of the QPD. The integral angles are \[Phi]=[0, \[Pi]/2] and \[Theta]=[0,\[Theta]]. Note the free parameter \[Theta].
\[Sigma] = differential noise power spectral density of dpdipd\[CapitalOmega] at the differential detector 
	under direction (\[Theta], \[Phi]) (eq. 9 in the main text)
\[CapitalStigma] = integral of the noise power spectral density of dpdipd\[CapitalOmega] for the integral angles \[Phi]=[0, 2\[Pi]] and \[Theta]=[0,\[Theta]].
	Note the free parameter \[Theta].
SmirrorQPD = power spectral density of the imprecision noise detected with a QPD for particle displacements
			along the direction (\[Theta]0, \[Phi]0) and for variable \[Theta] . (eq. A.17 in the main text)*)
(*Assumptions: light polarization is fixed to (\[Theta]p, \[Phi]p) \[Rule] (\[Theta]pmax, 0) = (\[Pi]/2,0)*)
(*-------------------*)
(*Differential displacement sensitivity*)
\[Beta][\[Theta]_,\[Phi]_,\[Theta]0_,\[Phi]0_]:= Evaluate[Simplify[2*Sqrt[R]*(4*\[Pi]/\[Lambda])*(nr[\[Theta]0,\[Phi]0].nr[\[Theta],\[Phi]])*dpdipd\[CapitalOmega][\[Theta],\[Phi],\[Theta]pmax,0],Assumptions->{0<=\[Theta]<=\[Pi], 0<=\[Phi]<=Pi/2}]];
Print["Differential sensitivity to x displacemnts \[Beta][\[Theta],\[Phi],Pi/2,0] = ", \[Beta][\[Theta],\[Phi],Pi/2,0] ]
Print["Differential sensitivity to y displacemnts \[Beta][\[Theta],\[Phi],Pi/2,Pi/2] = ", \[Beta][\[Theta],\[Phi],Pi/2,Pi/2] ]
Print["Differential sensitivity to z displacemnts \[Beta][\[Theta],\[Phi],0,0] = ", \[Beta][\[Theta],\[Phi],0,0] ]
(*displacement sensitivity on a single QPD quadrant*)
B[\[Theta]_,\[Theta]0_,\[Phi]0_]:= Evaluate[Simplify[Integrate[\[Beta][\[Theta]\[Theta],\[Phi]\[Phi],\[Theta]0,\[Phi]0]*Sin[\[Theta]\[Theta]],{\[Theta]\[Theta],0,\[Theta]},{\[Phi]\[Phi],0,\[Pi]/2},Assumptions->{0<=\[Theta]<=\[Pi]/2}]]];
Print["QPD sensitivity to x displacemnts 4*B[\[Theta],Pi/2,0] = ", 4*B[\[Theta],Pi/2,0]]
Print["QPD Sensitivity to y displacemnts 4*B[\[Theta],Pi/2,Pi/2] = ", 4*B[\[Theta],Pi/2,Pi/2]]
Print["QPD Sensitivity to z displacemnts 4*B[\[Theta],0,0] = ", 4*B[\[Theta],0,0]]
(*Differential shot noise*)
\[Sigma][\[Theta]_,\[Phi]_]:=Evaluate[Simplify[hbar*c/\[Lambda]*(R+1)*dpdipd\[CapitalOmega][\[Theta],\[Phi],\[Theta]pmax,0],Assumptions->{0<=\[Theta]<=\[Pi], 0<=\[Phi]<2*Pi}]]; 
Print["Differential shot noise PSD \[Sigma][\[Theta],\[Phi]] = ", \[Sigma][\[Theta],\[Phi]] ]
(*Shot noise on the full QPD area*)
\[CapitalStigma][\[Theta]_]:=Evaluate[Simplify[Integrate[\[Sigma][\[Theta]\[Theta],\[Phi]\[Phi]]*Sin[\[Theta]\[Theta]],{\[Theta]\[Theta],0,\[Theta]},{\[Phi]\[Phi],0,2*Pi},Assumptions->{0<=\[Theta]<=\[Pi]/2}]]];
Print["Integral shot noise PSD on the qpd  \[CapitalStigma][\[Theta]]= ", \[CapitalStigma][\[Theta]]]
(*Imprecision for the self-homodyne configuration with a QPD *)
SmirrorQPD[\[Theta]_,\[Theta]0_,\[Phi]0_]:=Evaluate[FullSimplify[1/(16*B[\[Theta],\[Theta]0,\[Phi]0]^2)*\[CapitalStigma][\[Theta]]]];
Print["SmirrorQPD =", SmirrorQPD[\[Theta],\[Theta]0,\[Phi]0]]
(*introducing natural units*)
SSmirrorQPD[\[Theta]_,\[Theta]0_,\[Phi]0_]:=Evaluate[FullSimplify[SmirrorQPD[\[Theta],\[Theta]0,\[Phi]0],Assumptions->{hbar==1,\[Lambda]==1, c==1, R==1, Pdip==1}]];


(*----------------------------------------------------------------------------------------------*)
(*----------------------------------------------------------------------------------------------*)
(*GENERAL FORMULAS FOR SELF-HOMDYNE DETECTION*)
Print[""]
(*Maximum and minimum imprecision*)
minS = FullSimplify[SSmirror[ArcSin[1],\[Theta]pmax,\[Theta]0max,\[Phi]0max]];
maxS = FullSimplify[SSmirror[ArcSin[1],\[Theta]pmax,\[Pi]/2,0]];
Print["Minimum imprecision (natural units): ", minS];
Print["Minimum imprecision (natural units): ", maxS];
Print[""]
Print["differential imprecision"]
Print["smirror = ", smirror[\[Theta],\[Phi],\[Pi]/2,0,\[Theta]0,\[Phi]0]]
Print["smirrorcirc = ", smirrorcirc[\[Theta],\[Phi],\[Theta]0,\[Phi]0]]
Print["smirroratom = ", smirroratom[\[Theta],\[Phi],\[Theta]0,\[Phi]0]]
Print[""]
Print["Total imprecision \[Theta]0=0"]
Print["Smirror(\[Theta]0=0)= ", FullSimplify[Smirror[\[Theta]D,\[Pi]/2,0,0,0]]]
Print["Smirrorcirc(\[Theta]0=0) = ", FullSimplify[Smirrorcirc[\[Theta]D,0,0]]]
Print[""]
Print["Total imprecision \[Theta]0=\[Pi]/2"]
Print["Smirror(\[Theta]0=\[Pi]/2)= ", FullSimplify[Smirror[\[Theta]D,\[Pi]/2,0,\[Pi]/2,\[Phi]0]]]
Print["Smirrorcirc(\[Theta]0=\[Pi]/2) = ", FullSimplify[Smirrorcirc[\[Theta]D,\[Pi]/2,\[Phi]0]]]
Print["Smirroratom= ", FullSimplify[Smirroratom[\[Theta]D,\[Theta]0,\[Phi]0]]]
Print[""]
Print["Total imprecision \[Theta]D=\[Pi]/2"]
Print["Smirror(\[Theta]=\[Pi]/2)= ", FullSimplify[Smirror[\[Pi]/2,\[Pi]/2,0,\[Theta]0,\[Phi]0]]]
Print["Smirrorcirc(\[Theta]=\[Pi]/2) = ", FullSimplify[Smirrorcirc[\[Pi]/2,\[Theta]0,\[Phi]0]]]
Print["Smirroratom(\[Theta]=\[Pi]/2)= ", FullSimplify[Smirroratom[\[Pi]/2,\[Theta]0,\[Phi]0]]]
Print[""]
Print["SmirrorQPD =", SmirrorQPD[\[Theta],\[Theta]0,\[Phi]0]]


Needs["PlotLegends`"]
(*----------------------------------------------------------------------------------------------*)
(*----------------------------------------------------------------------------------------------*)
(*PLOTS*)
(*----------------------------------------------------------------------------------------------*)
(*standard font size for all plots*)
lgdfontsize = 16;
(*frame flag that can be used for all plots*)
frameflg = False;
(*----------------------------------------------------------------------------------------------*)
(*3D axes label*)
xyzlabl = {Style["x",Bold,Black,lgdfontsize],Style["y",Bold,Black,lgdfontsize],Style["z",Bold,Black,lgdfontsize]};
(*angular information for linear polarization for the minimal imprecision configuration*)
SphericalPlot3D[1/ssmirror[Ttheta,Pphi,\[Theta]pmax,\[Theta]0max,\[Phi]0max] ,{Ttheta,0,\[Pi]/2},{Pphi,0,2*\[Pi]},
	AxesLabel ->xyzlabl, PlotRange->All]
(*angular information for circular polarization for the minimal imprecision configuration*)
SphericalPlot3D[1/ssmirrorcirc[Ttheta,Pphi,\[Theta]0max,\[Phi]0max] ,{Ttheta,0,\[Pi]/2},{Pphi,0,2*\[Pi]},
	AxesLabel ->xyzlabl, PlotRange->All]
(*Normalization by the highest imprecision*)
norm = SSmirror[ArcSin[1],\[Theta]pmax,\[Pi]/2,0];
(*2D axes label*)
xylabl = {Style["NA",Bold,Black,lgdfontsize],Style["1/S",Bold,Black,lgdfontsize]};
(*plot legend spec*)
pltlgdxyz100={Style["x \[Pi]",Black,lgdfontsize],Style["y \[Pi]",Black,lgdfontsize],Style["z \[Pi]",Black,lgdfontsize],
	Style["x \[Sigma]",Black,lgdfontsize],Style["y \[Sigma]",Black,lgdfontsize],Style["z \[Sigma]",Black,lgdfontsize],
	Style["x, y unpol",Black,lgdfontsize],Style["z unpol",Black,lgdfontsize],Style["y unpol",Black,lgdfontsize]};
(*line style spec*)
xyzlinestyle = {Directive[Red, Thick, Dashing[None]], Directive[Green, Thick, Dashing[None]],Directive[Blue, Thick, Dashing[None]],
Directive[Red, Thick, Dashed], Directive[Green, Thick, Dashed],Directive[Blue, Thick, Dashed],Directive[Orange, Thick, Dashing[None]], Directive[Purple, Thick, Dashing[None]]};
(*plotting*)
Plot[{maxS/SSmirror[ArcSin[x],\[Theta]pmax,\[Pi]/2,0], maxS/SSmirror[ArcSin[x],\[Theta]pmax,\[Pi]/2,\[Pi]/2], maxS/SSmirror[ArcSin[x],\[Theta]pmax,0,0], 
	maxS/SSmirrorcirc[ArcSin[x],\[Pi]/2,0], maxS/SSmirrorcirc[ArcSin[x],\[Pi]/2,\[Pi]/2], maxS/SSmirrorcirc[ArcSin[x],0,0],
	maxS/SSmirroratom[ArcSin[x],\[Pi]/2,0], maxS/SSmirroratom[ArcSin[x],0,0]},{x,0,1},
	PlotLegend->pltlgdxyz100,
	LegendPosition->{-0.4,-0.2},LegendShadow->None,LegendSize->1,LegendBorder->Black,
	ImageSize->Large, AxesLabel -> xylabl , Axes -> True, Frame->frameflg, PlotStyle->xyzlinestyle, BaseStyle->{FontSize->16},
	 PlotTheme->"Monochrome", AspectRatio -> 4/4, Background->White]
Plot[{maxS/SSmirror[ArcSin[x],\[Theta]pmax,\[Pi]/2,0], maxS/SSmirror[ArcSin[x],\[Theta]pmax,\[Pi]/2,\[Pi]/2], maxS/SSmirror[ArcSin[x],\[Theta]pmax,0,0], 
	maxS/SSmirrorcirc[ArcSin[x],\[Pi]/2,0], maxS/SSmirrorcirc[ArcSin[x],\[Pi]/2,\[Pi]/2], maxS/SSmirrorcirc[ArcSin[x],0,0],
	maxS/SSmirroratom[ArcSin[x],\[Pi]/2,0], maxS/SSmirroratom[ArcSin[x],0,0]},{x,0.75,0.95},
	PlotLegend->pltlgdxyz100,
	LegendPosition->{0.4,-0.7},LegendShadow->None,LegendSize->0.8,LegendBorder->Black,
	ImageSize->Large, AxesLabel -> xylabl , Axes -> True, Frame->frameflg, PlotStyle->xyzlinestyle, BaseStyle->{FontSize->16},
	PlotTheme->"Monochrome", AspectRatio -> 4/4, PlotRange->{0.1,0.6}, Background->White]


(*-------------------------------------------------------------------------------------------------*)
(*-------------------------------------------------------------------------------------------------*)
(*FORWARD AND BACKWARD DETECTION*)
(*Forward and backward detection scheme from article:

Felix Tebbenjohanns, Martin Frimmer, and Lukas Novotny
"Optimal position detection of a dipolar scatterer in a focused field"
Phys. Rev. A 100, 043821 \[Dash] Published 14 October 2019
url: https://link.aps.org/doi/10.1103/PhysRevA.100.043821
-------------------------------------------------------------
The label "fw" is included in the variable's names for the forward detection scheme.
The label "bw" is included in the variable's names for the bacward detection scheme.
The authors of this article normalized the final result for the Heisenberg limit. 
This must be corrected here for comparison.
*)
(*Parameter of Gaussian focused beam*)
Agaussian[\[Theta]_]:=(12/35-Cos[\[Theta]]^(5/2)/5-Cos[\[Theta]]^(7/2)/7)/(8/15-Cos[\[Theta]]^(3/2)/3-Cos[\[Theta]]^(5/2)/5);
(*These parameter are only temporary. They are defined in Phys. Rev. A 100, 043821*)
Bfwtmp[s_,c_,A_]:=s*Sqrt[c]{s*(1+2*c)/3, s*(2+c)/3, \[Pi]*(c-A)*(1+c)/4};
Bbwtmp[s_,c_]:=s*Sqrt[c]{s*(1+2*c)/3, s*(2+c)/3, \[Pi]*c*(1+c)/2};
Bfw[\[Theta]_, AA_]:=Evaluate[Integrate[Bfwtmp[Sin[\[Theta]\[Theta]],Cos[\[Theta]\[Theta]],A],{\[Theta]\[Theta],0,\[Theta]}, Assumptions->{Cos[\[Theta]]>=0, 0<A<1, 0<\[Theta]<\[Pi]/2, 0<\[Phi]<2*\[Pi]}]/.{A->AA}];
Bbw[\[Theta]_]:=Evaluate[Integrate[Bbwtmp[Sin[\[Theta]\[Theta]],Cos[\[Theta]\[Theta]]],{\[Theta]\[Theta],0,\[Theta]}, Assumptions->{Cos[\[Theta]]>=0, 0<A<1, 0<\[Theta]<\[Pi]/2, 0<\[Phi]<2*\[Pi]}]/.{A->AA}]
(*Efficiency \[Eta] as defined in Phys. Rev. A 100, 043821. Note that this is not the same
normalization convention that we use. In Phys. Rev. A 100, 043821 the Heisenberg limit
corresponds to \[Eta]=1 for all directions.*)
\[Eta]fw[\[Theta]_,\[Theta]f_]:=Evaluate[{30, 15, 15/(1+(5/2)*Agaussian[\[Theta]f]^2)}*Bfw[\[Theta], Agaussian[\[Theta]]]^2/(\[Pi]^2*Sin[\[Theta]]^2)/.{A->Agaussian[\[Theta]f]}];
\[Eta]bw[\[Theta]_]:=Evaluate[{30, 15, 15/(1+(5/2)*Agaussian[\[Theta]]^2)}*Bbw[\[Theta]]^2/(\[Pi]^2*Sin[\[Theta]]^2)/.{A->Agaussian[\[Theta]]}];
(*-------------------------------------------------------------------------------------------------*)
(*This renormalization and redistribution of axes translates the results reported in
 Phys. Rev. A 100, 043821 into our reference system in our natural units.
 The normalization takes advantage of the fact that the self-homodyne detection is equivalent to the ideal detection*)
(*Total imprecision for forward detection assuming an asymmetric NA lens setup*)
SSfwasym[\[Theta]_,\[Theta]f_]:=Evaluate[FullSimplify[{\[Eta]fw[\[Theta],\[Theta]f][[1]]/SSmirror[\[Pi]/2,\[Theta]pmax,\[Pi]/2,0], \[Eta]fw[\[Theta],\[Theta]f][[3]]/SSmirror[\[Pi]/2,\[Theta]pmax,\[Theta]0max,\[Phi]0max], \[Eta]fw[\[Theta],\[Theta]f][[2]]/SSmirror[\[Pi]/2,\[Theta]pmax,\[Pi]/2,\[Pi]/2]}^(-1)]];
(*Total imprecision for forward detection assuming an ymmetric NA lens setup*)
SSfw[\[Theta]_]:=Evaluate[FullSimplify[SSfwasym[\[Theta],\[Theta]]]];
(*Total imprecision for backward detection assuming an ymmetric NA lens setup*)
SSbw[\[Theta]_]:=Evaluate[FullSimplify[{\[Eta]bw[\[Theta]][[1]]/SSmirror[\[Pi]/2,\[Theta]pmax,\[Pi]/2,0], \[Eta]bw[\[Theta]][[3]]/SSmirror[\[Pi]/2,\[Theta]pmax,\[Theta]0max,\[Phi]0max], \[Eta]bw[\[Theta]][[2]]/SSmirror[\[Pi]/2,\[Theta]pmax,\[Pi]/2,\[Pi]/2]}^(-1)]];
Print["Forward imprecision (natural units) = ", SSfw[\[Theta]]]
Print["Backward imprecision (natural units) = ", SSbw[\[Theta]]]


(*-------------------------------------------------------------------------------------------------*)
(*-------------------------------------------------------------------------------------------------*)
(*PLOT SELF-HOMODYNE AND FORWARD BACKWARD*)
(*line spec*)
xyzlinestyle = {Directive[Red, Thick, Dashing[None]], Directive[Green, Thick, Dashing[None]],Directive[Blue, Thick, Dashing[None]],
Directive[Red, Thick, Dashed], Directive[Green, Thick, Dashed],Directive[Blue, Thick, Dashed],Directive[Purple, Thick, Dashed], Directive[Black, Thin, Dashed], Directive[Black, Thin, Dashed]};
(*legend spec*)
pltlgdxyz={Style["x self",Black,lgdfontsize],Style["y self",Black,lgdfontsize],Style["z self",Black,lgdfontsize],
Style["x fw and bw",Black,lgdfontsize],Style["y forward",Black,lgdfontsize],Style["z fw and bw",Black,lgdfontsize],
Style["y backward",Black,lgdfontsize], Style["H-limit-yz",Black,lgdfontsize], Style["H-limit-x",Black,lgdfontsize]};
(*plot*)
Plot[{minS/SSmirror[ArcSin[x],\[Theta]pmax,\[Pi]/2,0], minS/SSmirror[ArcSin[x],\[Theta]pmax,\[Pi]/2,\[Pi]/2],minS/SSmirror[ArcSin[x],\[Theta]pmax,\[Theta]0max,\[Phi]0max],
	 minS/SSfw[ArcSin[x]][[1]],0,minS/SSfw[ArcSin[x]][[3]], minS/SSbw[ArcSin[x]][[2]],minS/minS, minS/maxS},
	{x,0.01,1}, PlotRange->{0,1.1}, AspectRatio -> 3/4, Axes -> True, Frame->frameflg, AxesLabel -> xylabl, ImageSize->Large,
	PlotLegend->pltlgdxyz,LegendShadow->None,LegendPosition->{1,-0.5},LegendSize->1,PlotStyle->xyzlinestyle]


(*-------------------------------------------------------------------------------------------------*)
(*-------------------------------------------------------------------------------------------------*)
(*PLOT SELF-HOMODYNE WITH QUADRANT PHOTO DETECTOR*)
(*line spec*)
xyzlinestyle = {Directive[Red, Thick, Dashing[None]], Directive[Green, Thick, Dashing[None]],Directive[Blue, Thick, Dashing[None]],
Directive[Red, Thick, Dashed], Directive[Green, Thick, Dashed],Directive[Blue, Thick, Dashed], Directive[Black, Thin, Dashed], Directive[Black, Thin, Dashed]};
(*legend spec*)
pltlgdxyz={Style["x self qpd",Black,lgdfontsize],Style["y self qpd",Black,lgdfontsize],Style["z self qpd",Black,lgdfontsize],
	Style["x self ideal",Black,lgdfontsize],Style["y self ideal",Black,lgdfontsize],Style["z self ideal",Black,lgdfontsize],
	Style["H-limit-yz",Black,lgdfontsize], Style["H-limit-x",Black,lgdfontsize]};
(*plot*)
Plot[{minS/SSmirrorQPD[ArcSin[x],Pi/2,0], minS/SSmirrorQPD[ArcSin[x],\[Pi]/2,\[Pi]/2], minS/SSmirrorQPD[ArcSin[x],0,0],
	minS/SSmirror[ArcSin[x],\[Theta]pmax,\[Pi]/2,0], minS/SSmirror[ArcSin[x],\[Theta]pmax,\[Pi]/2,\[Pi]/2],minS/SSmirror[ArcSin[x],\[Theta]pmax,\[Theta]0max,\[Phi]0max],
	minS/minS, minS/maxS},
	{x,0.01,1}, PlotRange->{0,1.1}, AspectRatio -> 3/4, Axes -> True, Frame->frameflg, AxesLabel -> xylabl, ImageSize->Large,
	PlotLegend->pltlgdxyz,LegendShadow->None,LegendPosition->{1,-0.5},LegendSize->1,PlotStyle->xyzlinestyle]
