global Cm Ena Ek Eleak1 Eleak2 Eleak3
global Gnaf1 Gkdr1 Gkir1 Gkaf1 Gleak1
global Gnaf2 Gkdr2 Gkir2 Gkaf2 Gleak2
global Gnaf3 Gkdr3 Gkir3 Gkaf3 Gleak3
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
Eleak1 = -60;       %Leak resting potential (MSN type I)
Eleak2 = -60;       %Leak resting potential (MSN type II)
Eleak3 = -60;       %Leak resting potential (Pyramidal type II)

%Maximal conductances for MSN neuron type1 (values for soma, s/cm^2)
Gnaf1 = 100;        %Naf maximal conductance (original 1.875)
Gkdr1 = 10;         %Kdr maximal conductance (original 0.0005)
Gkir1 = 0.15;       %Kir maximal conductance (original 0.00015)
Gkaf1 = 2;          %Kaf maximal conductance (original 0.21) (increasing Gkaf1 can cause generation of complex spike)
Gleak1 = 0.1;       %Leak maximal conductance (original 0.0000115)

%Maximal conductances for MSN neuron type2 (values for soma, s/cm^2)
%original values are in S/M^2
%naf = 90000
%kdr = 6.04
%kir = 8
%kaf = 765.24
%leak (0.11500000005175000002328750001048 = 1/(RM=8.69565217))
Gnaf2 = 200;        %Naf maximal conductance (original 9) 120
Gkdr2 = 20;         %Kdr maximal conductance (original 0.000604) 50
Gkir2 = 1.0;        %Kir maximal conductance (original 0.0008) 0.03
Gkaf2 = 2;          %Kaf maximal conductance (original 0.76524) 1.5
Gleak2 = 0.1;       %Leak maximal conductance (original 0.0000115)

%Maximal conductances for pyramidal neuron type1 (values for soma, s/cm^2)
%(Ena = 55mV for pyramidal neuron type1)
Gnaf3 = 32;         %Naf maximal conductance (original 1.875)
Gkdr3 = 10;         %Kdr maximal conductance (original 0.0005)
Gkir3 = 0;          %Kir maximal conductance (original 0.0008) 0.03
Gkaf3 = 2;          %Kaf maximal conductance (original 0.21) (increasing Gkaf1 can cause generation of complex spike)
Gleak3 = 0.1;       %Leak maximal conductance (original 0.0000115)

MH_inf = zeros( MAX_MH, 1 ); % array for stady-states for activation&inactivation variables
TAU_mh = zeros( MAX_MH, 1 ); % array for time constans for activation&inactivation variables
MH = zeros( MAX_MH, 1 );     % array for activation&inactivation variables

Gnaf1 = Gnaf1/Cm;
Gkdr1 = Gkdr1/Cm;
Gkir1 = Gkir1/Cm;
Gkaf1 = Gkaf1/Cm;
Gleak1 = Gleak1/Cm;

Gnaf2 = Gnaf2/Cm;
Gkdr2 = Gkdr2/Cm;
Gkir2 = Gkir2/Cm;
Gkaf2 = Gkaf2/Cm;
Gleak2 = Gleak2/Cm;

Gnaf3 = Gnaf3/Cm;
Gkdr3 = Gkdr3/Cm;
Gkir3 = Gkir3/Cm;
Gkaf3 = Gkaf3/Cm;
Gleak3 = Gleak3/Cm;
