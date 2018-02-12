% Type 1. Martin Pospischil, Maria Toledo-Rodriguez, Cyril Monier,
% Zuzanna Piwkowska, Thierry Bal, Yves Frégnac, Henry Markram, Alain Destexhe
% Minimal Hodgkin–Huxley type models for different classes
% of cortical and thalamic neurons, Biol Cybern (2008) 99:427–441
function dx = df_cortex1(t,x)
% global variables 
    global Iinj1 T0_inj 
    global VT1 VT2 Ena Ek Eleak3
    global Gnaf3 Gkdr3 Gkir3 Gkaf3 Gleak3
    global iNaF_M iNaF_H iKdr_M iKir_M iKaf_M iKaf_H MAX_MH
    global MH_inf TAU_mh MH
% initialization
    MH = x(2:1+MAX_MH);       % array for activation&inactivation variables
    v = x(1);                 % membrane potentials
    i_inj = 0; 
    if( t > T0_inj )
        i_inj = Iinj1-(t-T0_inj)*0.015;
    end
        
% NaF activation
    alpha = -0.32*(v-VT1-13)./(exp(-(v-VT1-13)/4)-1);
    beta = 0.28*(v-VT1-40)./(exp((v-VT1-40)/5)-1);
%     alpha = -0.32*(v-13.1)/(exp(-(v-13.1)/4)-1);
%     beta = 0.28*(v-40.1)/(exp((v-40.1)/5)-1);
    MH_inf(iNaF_M) = alpha/(alpha+beta);
    TAU_mh(iNaF_M) = 1/(alpha+beta);
    alpha = 0.128*(exp(-(v-VT1-17)/18));
    beta = 4/(exp(-(v-VT1-40)/5)+1);
%     alpha = 0.128./exp((v-17)/18);
%     beta = 4./(exp(-(v-40)/5)+1);
    
    MH_inf(iNaF_H) = alpha/(alpha+beta);
    TAU_mh(iNaF_H) = 1/(alpha+beta);

    TAU_mh(iNaF_M) = TAU_mh(iNaF_M)*0.6;
    TAU_mh(iNaF_H) = TAU_mh(iNaF_H)*2;
%   Inaf = Gnaf3*MH(iNaF_M)^3*MH(iNaF_H)*(v-Ena)

%Kdr channel
    alpha = -0.32*(v-VT2-13)./(exp(-(v-VT2-13)/4)-1);
    beta = 0.28*(v-VT2-40)./(exp((v-VT2-40)/5)-1);
%     alpha = -0.016*(v-35.1)/(exp(-(v-35.1)/5)-1);
%     beta = 0.25./exp((v-20)/40);
    MH_inf(iKdr_M) = alpha/(alpha+beta);
    TAU_mh(iKdr_M) = 1/(alpha+beta);

    TAU_mh(iKdr_M) = TAU_mh(iKdr_M)*7;
%   Ikdr = Gkdr3*MH(iKdr_M)^4*(v-Ek)

%Kaf channel
    MH_inf(iKaf_M) = 0;
    MH_inf(iKaf_H) = 1;
    TAU_mh(iKaf_M) = 1;
    TAU_mh(iKaf_H) = 1;
%   Ikaf = Gkaf3*MH(iKaf_M)^2*MH(iKaf_H)*(v-Ek)

%Kir channel
    MH_inf(iKir_M) = 1;                         % Kir activation (no Kir current)
    TAU_mh(iKir_M) = 1;                         % Time constant for Kir activation
%   Ikir = Gkir3*MH(iKir_M)*(v-Ek)

%Leak channel
%   Ileak = Eleak3*(v-Eleak3)

%DIFFERENTIAL EQUATIONS 
    dx = 0*x;
%   for activation and inactivation
    dx(2:1+MAX_MH) = (-MH+MH_inf)./TAU_mh;
%   for membrane potentials    
    dx(1) = -Gnaf3*MH(iNaF_M)^3*MH(iNaF_H)*(v-Ena)-Gkdr3*MH(iKdr_M)^4*(v-Ek)...
            -Gkir3*MH(iKir_M)*(v-Ek)-Gkaf3*MH(iKaf_M)^2*MH(iKaf_H)*(v-Ek)...
            -Gleak3*(v-Eleak3)-i_inj;
%-- THE END
