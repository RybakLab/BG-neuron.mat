global Cm Ena Ek Eleak_msn Eleak_cortex
global Gnaf_msn Gkdr_msn Gkir_msn Gkaf_msn Gleak_msn
global Gnaf_cortex Gkdr_cortex Gkir_cortex Gkaf_cortex Gleak_cortex
global iNaF_M iNaF_H iKdr_M iKaf_M iKaf_H iKir_M MAX_MH
global MH_inf TAU_mh MH

iNaF_M = 1;         %index for Naf activation
iNaF_H = 2;         %index for Naf inactivation
iKdr_M = 3;         %index for Kdr activation
iKaf_M = 4;         %index for Kaf activation
iKaf_H = 5;         %index for Kir inactivation
iKir_M = 6;         %index for Kir activation
MAX_MH = iKir_M;    %total number of activation&inactivation variables

Cm = 1;             %membane capacitance (uF/cm^2)
%Resting potentials (mV)
Ena = 50;           %Na resting potential 
Ek = -90;           %K resting potential
Eleak_msn = -60;    %Leak resting potential (for MSN neurons)
Eleak_cortex = -65; %Leak resting potential (cortical cells)
%Maximal conductances for MSN neuron type1 (values for soma, s/cm^2)
Gnaf_msn = 100;     %Naf maximal conductance
Gkdr_msn = 10;      %Kdr maximal conductance
Gkir_msn = 0.15;    %Kir maximal conductance
Gkaf_msn = 2;       %Kaf maximal conductance
Gleak_msn = 0.1;    %Leak maximal conductance
%Maximal conductances for pyramidal neuron type2 (values for soma, s/cm^2)
%(Ena = 55mV for pyramidal neuron type1)
Gnaf_cortex = 25;   %Naf maximal conductance
Gkdr_cortex = 7;   %Kdr maximal conductance
Gkir_cortex = 0;    %Kir maximal conductance
Gkaf_cortex = 0;    %Kaf maximal conductance
Gleak_cortex = 0.1; %Leak maximal conductance

MH_inf = zeros( MAX_MH, 1 ); % array for stady-states for activation&inactivation variables
TAU_mh = zeros( MAX_MH, 1 ); % array for time constans for activation&inactivation variables
MH = zeros( MAX_MH, 1 );     % array for activation&inactivation variables

Gnaf_msn = Gnaf_msn/Cm;
Gkdr_msn = Gkdr_msn/Cm;
Gkir_msn = Gkir_msn/Cm;
Gkaf_msn = Gkaf_msn/Cm;
Gleak_msn = Gleak_msn/Cm;

Gnaf_cortex = Gnaf_cortex/Cm;
Gkdr_cortex = Gkdr_cortex/Cm;
Gkir_cortex = Gkir_cortex/Cm;
Gkaf_cortex = Gkaf_cortex/Cm;
Gleak_cortex = Gleak_cortex/Cm;
