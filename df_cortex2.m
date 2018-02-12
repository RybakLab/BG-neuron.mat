% Type 2. Xiao-Jing Wang Calcium Coding and Adaptive Temporal Computation 
% in Cortical Pyramidal Neurons, J Neurophysiol. 1998; 79(3):1549-66.
function dx = df_cortex2(t,x)
% global variables 
    global Iinj2 T0_inj 
    global Ena Ek Eleak4
    global Gnaf4 Gkdr4 Gkir4 Gkaf4 Gleak4
    global iNaF_M iNaF_H iKdr_M iKir_M iKaf_M iKaf_H MAX_MH
    global MH_inf TAU_mh MH
% initialization
    MH = x(2:1+MAX_MH);       % array for activation&inactivation variables
    v = x(1);                 % membrane potentials
    i_inj = 0; 
    if( t > T0_inj )
        i_inj = Iinj2-(t-T0_inj)*0.025;
    end
        
% NaF activation
    alpha = -0.1*(v+33)/(exp(-(v+33)/10)-1);
    beta = 4/exp((v+58)/12);
    MH_inf(iNaF_M) = alpha/(alpha+beta);
    TAU_mh(iNaF_M) = 1/(alpha+beta);

    alpha = 0.07/(exp((v+50)/10));
    beta = 1/(exp(-(v+20)/10)+1);
    MH_inf(iNaF_H) = alpha/(alpha+beta);
    TAU_mh(iNaF_H) = 1/(alpha+beta);

    TAU_mh(iNaF_M) = TAU_mh(iNaF_M)*1.1;
    TAU_mh(iNaF_H) = TAU_mh(iNaF_H)*2;
%   Inaf = Gnaf4*MH(iNaF_M)^3*MH(iNaF_H)*(v-Ena)

%Kdr channel
    alpha = -0.01*(v+34)/(exp(-(v+34)/10)-1);
    beta = 0.125/exp((v+44)/25);
    MH_inf(iKdr_M) = alpha/(alpha+beta);
    TAU_mh(iKdr_M) = 1/(alpha+beta);

    TAU_mh(iKdr_M) = TAU_mh(iKdr_M)*2.5;
%   Ikdr = Gkdr4*MH(iKdr_M)^4*(v-Ek)

%Kaf channel
    MH_inf(iKaf_M) = 0;
    MH_inf(iKaf_H) = 1;
    TAU_mh(iKaf_M) = 1;
    TAU_mh(iKaf_H) = 1;
%   Ikaf = Gkaf4*MH(iKaf_M)^2*MH(iKaf_H)*(v-Ek)

%Kir channel
    MH_inf(iKir_M) = 1;
    TAU_mh(iKir_M) = 1;
%   Ikir = Gkir4*MH(iKir_M)*(v-Ek)

%Leak channel
%   Ileak = Eleak4*(v-Eleak3)

%DIFFERENTIAL EQUATIONS 
    dx = 0*x;
%   for activation and inactivation
    dx(2:1+MAX_MH) = (-MH+MH_inf)./TAU_mh;
%   for membrane potentials    
    dx(1) = -Gnaf4*MH(iNaF_M)^3*MH(iNaF_H)*(v-Ena)-Gkdr4*MH(iKdr_M)^4*(v-Ek)...
            -Gkir4*MH(iKir_M)*(v-Ek)-Gkaf4*MH(iKaf_M)^2*MH(iKaf_H)*(v-Ek)...
            -Gleak4*(v-Eleak4)-i_inj;
%-- THE END
