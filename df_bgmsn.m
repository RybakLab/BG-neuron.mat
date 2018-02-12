% The simplified model of medium spiny neuron in ventral striatum of
% the basal ganglia in the rat (slice preparation) based on research
% conducting by Finkels lab.
% 1. NMDA/AMPA Ratio Impacts State Transitions and Entrainment to Oscillations 
% in a Computational Model of the Nucleus Accumbens Medium Spiny Projection 
% Neuron. John A. Wolf, Jason T. Moyer, Maciej T. Lazarewicz, Diego Contreras,
% Marianne Benoit-Marand, Patricio O’Donnell, and Leif H. Finkel, 
% The Journal of Neuroscience, October 5, 2005 25(40):9080 –9095
% 2. Effects of Dopaminergic Modulation on the Integrative Properties
% of the Ventral Striatal Medium Spiny Neuron Jason T. Moyer, John A. Wolf,
% and Leif H. Finkel. Neurophysiol 98: 3731–3748, 2007.
% The description of all ion channelse (including tabulated time constants) 
% is below the Matlab code. The parameters of the model were taken from 
% https://senselab.med.yale.edu/ModelDB/showModel.cshtml?model=112834&file=/nacb_msp/#tabs-2
function dx = df_bgmsn(t,x)
% global variables 
    global Iinj1 T0_inj 
    global Ena Ek Eleak_msn 
    global Gnaf_msn Gkdr_msn Gkir_msn Gkaf_msn Gleak_msn
    global iNaF_M iNaF_H iKdr_M iKir_M iKaf_M iKaf_H MAX_MH
    global MH_inf TAU_mh MH
% initialization
    MH = x(2:1+MAX_MH);       % array for activation&inactivation variables
    v = x(1);                 % membrane potentials
    i_inj = 0; 
    if( t > T0_inj )
        i_inj = Iinj1-(t-T0_inj)*0.019;
    end
    i_inj = Iinj1;
% NaF activation (m - generic; h - modified generic (Slope2 for Tm -> -inf))
    MH_inf(iNaF_M) = 1/(1+exp(-(v+23.9)/11.8)); % NaF activation
    TAU_mh(iNaF_M) = 0.1+0.2/cosh((v+1.8)/29);  % Time constant for NaF activation
    MH_inf(iNaF_H) = 1/(1+exp((v+62.9)/10.7));  % NaF inactivation
    TAU_mh(iNaF_H) = 0.3+1/(1+exp((v-8.3)/18)); % Time constant for NaF inactivation
    TAU_mh(iNaF_M) = TAU_mh(iNaF_M)*1;          % 3 times faster for rat (in vitro) - not working
    TAU_mh(iNaF_H) = TAU_mh(iNaF_H)*1;          % 3 times faster for rat (in vitro) - not working
%   Inaf = Gnaf_msn*MH(iNaF_M)^3*MH(iNaF_H)*(v-Ena)
%Kdr channel (activation - alpha/beta)
    alpha =  1/exp(-(v+13)/12.5);
    beta = 1/exp((v+13)/33.3);
    MH_inf(iKdr_M) = alpha/(alpha+beta);        % Kdr activation
    TAU_mh(iKdr_M) = 50*1/(alpha+beta);         % Time constant for Kdr activation
    TAU_mh(iKdr_M) = TAU_mh(iKdr_M)*2;          % 2 times slower for rat (in vitro)
%   Ikdr = Gkdr_msn*MH(iKdr_M)^4*(v-Ek)
%Kir channel (m - generic)
    MH_inf(iKir_M) = 1/(1+exp((v+82.0)/13));    % Kir activation (52+30)
    TAU_mh(iKir_M) = 3.7+3.96/cosh((v+30)/41.5); % Time constant for Kir activation
    TAU_mh(iKir_M) = TAU_mh(iKir_M)*2;          % 2 times slower for rat (in vitro)
%   Ikir = Gkir_msn*MH(iKir_M)*(v-Ek)
%Kaf channel (m & h - modified generic (Slope2 for Tm -> -0))
    MH_inf(iKaf_M) = 1/(1+exp(-(v+10.0)/17.7)); % Kaf activation
    MH_inf(iKaf_H) = 1/(1+exp( (v+75.6)/10.0)); % Kaf inactivation
    TAU_mh(iKaf_M) = 0.85+0.018/exp(v/51.34);% Time constant for Kaf activation
    TAU_mh(iKaf_H) = 14;                        % Time constant for Kaf inactivation
    TAU_mh(iKaf_M) = TAU_mh(iKaf_M)/3;          % 3 times faster for rat (in vitro)
    TAU_mh(iKaf_H) = TAU_mh(iKaf_H)/3;          % 3 times faster for rat (in vitro)
%   Ikaf = Gkaf_msn*MH(iKaf_M)^2*MH(iKaf_H)*(v-Ek)
%Leak channel
%   Ileak = Gleak_msn*(v-Eleak)
%DIFFERENTIAL EQUATIONS 
    dx = 0*x;
%   for activation and inactivation
    dx(2:1+MAX_MH) = (-MH+MH_inf)./TAU_mh;
%   for membrane potentials    
    dx(1) = -Gnaf_msn*MH(iNaF_M)^3*MH(iNaF_H)*(v-Ena)-Gkdr_msn*MH(iKdr_M)^4*(v-Ek)...
            -Gkir_msn*MH(iKir_M)*(v-Ek)-Gkaf_msn*MH(iKaf_M)^2*MH(iKaf_H)*(v-Ek)...
            -Gleak_msn*(v-Eleak_msn)-i_inj;
%-- THE END
