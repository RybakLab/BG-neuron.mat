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
function dx = df_msn1(t,x)
% global variables 
    global Iinj1 T0_inj 
    global Ena Ek Eleak1
    global Gnaf1 Gkdr1 Gkir1 Gkaf1 Gleak1
    global iNaF_M iNaF_H iKdr_M iKir_M iKaf_M iKaf_H MAX_MH
    global MH_inf TAU_mh MH
% initialization
    MH = x(2:1+MAX_MH);       % array for activation&inactivation variables
    v = x(1);                 % membrane potentials
    i_inj = 0; 
    if( t > T0_inj )
        i_inj = Iinj1-(t-T0_inj)*0.019;
    end
        
% NaF activation (m - generic; h - modified generic (Slope2 for Tm -> -inf))
    MH_inf(iNaF_M) = 1/(1+exp(-(v+23.9)/11.8)); % NaF activation
    TAU_mh(iNaF_M) = 0.1+0.2/cosh((v+1.8)/29);  % Time constant for NaF activation
    MH_inf(iNaF_H) = 1/(1+exp((v+62.9)/10.7));  % NaF inactivation
    TAU_mh(iNaF_H) = 0.3+1/(1+exp((v-8.3)/18)); % Time constant for NaF inactivation
    
    TAU_mh(iNaF_M) = TAU_mh(iNaF_M)*1;          % 3 times faster for rat (in vitro) - not working
    TAU_mh(iNaF_H) = TAU_mh(iNaF_H)*1;          % 3 times faster for rat (in vitro) - not working
%   Inaf = Gnaf1*MH(iNaF_M)^3*MH(iNaF_H)*(v-Ena)

%Kdr channel (m - alpha/beta)
%   kdr_a = -(0.616+0.014*v)/exp(-( v+44 )/2.31);
%   kdr_b = 0.0043/exp((v+44)/34);
%   MH_inf(iKdr_M) = kdr_a/(kdr_a+kdr_b);       % Kdr activation
%   TAU_mh(iKdr_M) = 1/(kdr_a+kdr_b);           % Time constant for Kdr activation
%   The original description is failed. Specifically, the time constant
%   changes its sign on [-150,50]mV interval of membrane ponetnial.
%   The Kdr channel was from Migliore et al (1999) Role of an A-type K+
%   conductance in the back-propagation of action potentials in the dendrites 
%   of hippocampal pyramidal neurons. J Comput Neurosci 7: 5–15.
%    kdr_a = exp(-(v+13)/9.09); kdr_b = exp(-(v+13)/12.5);
%    alpha =  1/kdr_b;
%    beta = kdr_a/kdr_b;
    alpha =  1/exp(-(v+13)/12.5);
    beta = 1/exp((v+13)/33.3);
    MH_inf(iKdr_M) = alpha/(alpha+beta);        % Kdr activation
    TAU_mh(iKdr_M) = 50*1/(alpha+beta);         % Time constant for Kdr activation
    TAU_mh(iKdr_M) = TAU_mh(iKdr_M)*2;          % 2 times slower for rat (in vitro)
%   Ikdr = Gkdr1*MH(iKdr_M)^4*(v-Ek)

%Kir channel (m - generic)
    MH_inf(iKir_M) = 1/(1+exp((v+82.0)/13));    % Kir activation (52+30)
    TAU_mh(iKir_M) = 3.7+3.96/cosh((v+30)/41.5); % Time constant for Kir activation
    TAU_mh(iKir_M) = TAU_mh(iKir_M)*2;          % 2 times slower for rat (in vitro)
%   Ikir = Gkir1*MH(iKir_M)*(v-Ek)

%Kaf channel (m & h - modified generic (Slope2 for Tm -> -0))
    MH_inf(iKaf_M) = 1/(1+exp(-(v+10.0)/17.7)); % Kaf activation
    MH_inf(iKaf_H) = 1/(1+exp( (v+75.6)/10.0)); % Kaf inactivation
    TAU_mh(iKaf_M) = 0.85+0.018/exp(v/51.34);% Time constant for Kaf activation
    TAU_mh(iKaf_H) = 14;                        % Time constant for Kaf inactivation
    TAU_mh(iKaf_M) = TAU_mh(iKaf_M)/3;          % 3 times faster for rat (in vitro)
    TAU_mh(iKaf_H) = TAU_mh(iKaf_H)/3;          % 3 times faster for rat (in vitro)
%   Ikaf = Gkaf1*MH(iKaf_M)^2*MH(iKaf_H)*(v-Ek)

%Leak channel
%   Ileak = Gleak1*(v-Eleak)

%DIFFERENTIAL EQUATIONS 
    dx = 0*x;
%   for activation and inactivation
    dx(2:1+MAX_MH) = (-MH+MH_inf)./TAU_mh;
%   for membrane potentials    
    dx(1) = -Gnaf1*MH(iNaF_M)^3*MH(iNaF_H)*(v-Ena)-Gkdr1*MH(iKdr_M)^4*(v-Ek)...
            -Gkir1*MH(iKir_M)*(v-Ek)-Gkaf1*MH(iKaf_M)^2*MH(iKaf_H)*(v-Ek)...
            -Gleak1*(v-Eleak1)-i_inj;
%-- THE END

% % % % %NaF activation                                                          ON
% % % % %    I = Gmax*m^3*h*(V-Ena)
% % % % %    naf_m = 1/(1+exp(-(v+25.9)/11.8); 
% % % % %    naf_h = 1/(1+exp((v+64.9)/10.7);
% % % % %    naf_tm = 0.1+0.2/cosh((v+1.8)/29);
% % % % %    naf_th = 0.3+1/(1+exp((v-8.3)/18));
% % % % 
% % % % %NaP activation                                                         OFF
% % % % %    I = Gmax*m*h*(V-Ena)
% % % % %    nap_m = 1/(1+exp(-(v+52.6)/4.6);
% % % % %    nap_h = 1/(1+exp( (v+48.8)/10.0);
% % % % %    if( v < 40 )
% % % % %       nap_tm = 0.025+0.14*exp((v+40)/10)
% % % % %    else
% % % % %       nap_tm = 0.02+0.145*exp(-(v+40)/10)
% % % % %    endif
% % % % %    nap_th = 1810+4200/cosh(-(V+101)/70.2);
% % % % 
% % % % %Kdr channel                                                             ON
% % % % %    I = Gmax*m^4*(V-Ek);
% % % % %    A = -(0.616+0.014*v)/exp(-( v+44 )/2.31) 
% % % % %    B = 0.0043/exp((v+44)/34)
% % % % %    kdr_m = A(/(A+B)
% % % % %    kdr_tm = 1/(A+B)
% % % % 
% % % % %KIR channel                                                             ON
% % % % %    I = Gmax*m*(V-Ek); 
% % % % %    kir_m = 1/(1+exp((v+52.0)/17.5)
% % % % %    kir_tm = 3.7+3.96/cosh((V+30)/41.5)
% % % % 
% % % % %KAf channel                                                             ON
% % % % %    I = Gmax*m^2*h*(V-Ek);
% % % % %    Gmax = 0.36( 0.0033 dendr)
% % % % %    kaf_m = 1/(1+exp(-(v+10.0)/17.7)
% % % % %    kaf_h = 1/(1+exp( (v+75.6)/10.0)
% % % % %    kaf_tm = 0.85+0.018*exp(-(v+0)/51.34)
% % % % %    kaf_th = 4.67
% % % % 
% % % % %KAs channel                                                            OFF
% % % % %    I = Gmax*m^2*(a*h-(1-a))*(V-Ek); a = 0.996;
% % % % %    Gmax = 0.0104 (9.51*10^-4 dendr)
% % % % %    kas_m = 1/(1+exp(-(v+27.0)/16.0)
% % % % %    kas_h = 1/(1+exp( (v+33.5)/21.5)
% % % % %    kas_tm = 0.378+9.91*exp((-(v+34.3)/30.1)?2)
% % % % %    A = exp(-(v+90.96)/29.01)
% % % % %    B =  exp((v+90.96)/100)
% % % % %    kas_th = 1097.4/(A+B)
% % % % 
% % % % %KRP channel                                                            OFF
% % % % %    I = Gmax*m*(a*h-(1-a))*(V-Ek); a = 0.7; 
% % % % %    Gmax = 1.5*10^-4 (soma only)
% % % % %    krp_m = 1/(1+exp(-(v+13.5)/11.8)
% % % % %    krp_h = 1/(1+exp( (v+54.7)/18.6)
% % % % %    krp_tm = 103/cosh(-(v+76)/67.3)
% % % % %    krp_th = 2085+4935/((1+exp((v-96.7)/13.3)))
% % % % 
% % % % %Leak channel                                                            ON
% % % % %    I = Gmax*(V-Eleak); 
% % % % %    Gmax = 11.5*10^-6
% % % % %    Eleak = -70 mV
% % % % %All time constants are tabulated in [-200,200] mv range
% % % % %1st stage
% % % % %Na channels; Ena = 50 mV
% % % % %NaF channel I = Gmax*m^3*h*(V-Ena)
% % % % %   Gmax= 1.875 (0.0244 dendr)
% % % % %   Minf = 1/(1+exp(-(v+25.9)/11.8)
% % % % %   Hinf = 1/(1+exp( (v+64.9)/10.7)
% % % % %   TauM 
% % % % %       Tm = 0.1+0.2/cosh(0.0614926202219816+0.0345877123808851*V)
% % % % %       Tm = 0.1+0.2/cosh((V+1.8)/29) !!!!!!!!!!!! final 
% % % % %    0.06
% % % % %    0.06
% % % % %    0.06
% % % % %    0.07
% % % % %    0.09
% % % % %    0.11
% % % % %    0.13
% % % % %    0.20
% % % % %    0.32
% % % % %    0.16
% % % % %    0.15
% % % % %    0.12
% % % % %    0.08
% % % % %    0.06
% % % % %    0.06
% % % % %    0.06
% % % % %    0.06
% % % % %   Tauh
% % % % %    Th = 0.286821700366512+1.62748629811592/(1.58575145540076+exp(0.0554564690226297*V))
% % % % %    Th = 0.3+1/(1+exp((V-8.3)/18)) !!!!!!!!!!!! final
% % % % %(mod gen t0+2*t/(exp((v-v12)/slp)+exp(-(v-v12_2)/slp_2)) slp_2->inf)
% % % % %    1.3
% % % % %    1.3
% % % % %    1.3
% % % % %    1.3
% % % % %    1.3
% % % % %    1.3
% % % % %    1.3
% % % % %    1.3
% % % % %    0.85
% % % % %    0.5
% % % % %    0.45
% % % % %    0.32
% % % % %    0.30
% % % % %    0.28
% % % % %    0.28
% % % % %    0.28
% % % % %NaP channel I = Gmax*m*h*(V-Ena); 
% % % % %   Gmax = 5*10^-5 (1.73*10^-7 dendr)
% % % % %   Minf = 1/(1+exp(-(v+52.6)/4.6)
% % % % %   Hinf = 1/(1+exp( (v+48.8)/10.0)
% % % % %   TauM
% % % % %    Tm = 0.025+0.14*exp( (v+40)/10), Vm <40
% % % % %    Tm = 0.02+0.145*exp(-(v+40)/10), Vm >=40
% % % % %   TauH
% % % % %    Th = 1809.94567769994+4197.55537625749/cosh(-1.43737211730468-0.0142508947256976*V)
% % % % %    Th = 1810+4200/cosh(-(V+101)/70.2) !!!!!!!!!!!! final 
% % % % %generic( t0+t/cosh(( v-v12 )/slp )))
% % % % %    4500
% % % % %    4750
% % % % %    5200
% % % % %    6100
% % % % %    6300
% % % % %    5000
% % % % %    4250
% % % % %    3500
% % % % %    3000
% % % % %    2700
% % % % %    2500
% % % % %    2100
% % % % %    2100
% % % % %    2100
% % % % %    2100
% % % % %K channels; Ek = -90 mV
% % % % %Kdr channel I = Gmax*m^4*(V-Ek);
% % % % %   Gmax = 0.0005
% % % % %   Minf = A(/(A+B)
% % % % %   TauM = 1/(A+B)
% % % % %   A = -(0.616+0.014*Vm)/exp(-( Vm+44 )/2.3?1)
% % % % %   B = 0.0043/exp((Vm+44)/34)
% % % % %KIR channel I = Gmax*m*(V-Ek); 
% % % % %   Gmax = 1.4*10^-4
% % % % %   Minf = 1/(1+exp((v+52.0)/17.5)
% % % % %   TauM
% % % % %    Tm = 3.71204354980169+3.95813721509394/cosh(0.720864641782704+0.0241025961231886*V)
% % % % %    Tm = 3.7+3.96/cosh((V+30)/41.5) !!!!!!!!!!!! final
% % % % %generic( t0+t/cosh(( v-v12 )/slp )))
% % % % %    3.7313
% % % % %    4.0000
% % % % %    4.7170
% % % % %    5.3763
% % % % %    6.0606
% % % % %    6.8966
% % % % %    7.6923
% % % % %    7.1429
% % % % %    5.8824
% % % % %    4.4444
% % % % %    4.0000
% % % % %    4.0000
% % % % %    4.0000
% % % % %    4.0000
% % % % %    4.0000
% % % % %    4.0000
% % % % %KAf channel I = Gmax*m^2*h*(V-Ek);
% % % % %   Gmax = 0.36( 0.0033 dendr)
% % % % %   Minf = 1/(1+exp(-(v+10.0)/17.7)
% % % % %   Hinf = 1/(1+exp( (v+75.6)/10.0)
% % % % %   TauM
% % % % %    Tm = 0.845566331075776+0.0176626101976999*exp(-0.0194746130312227*V)
% % % % %    Tm = 0.85+0.018*exp(-(V+0)/51.34) !!!!!!!!!!!! final
% % % % %(mod gen t0+2*t/(exp((v-v12)/slp)+exp(-(v-v12_2)/slp_2)) slp_2->inf)
% % % % %    1.8
% % % % %    1.1
% % % % %    1.0
% % % % %    1.0
% % % % %    0.9
% % % % %    0.8
% % % % %    0.9
% % % % %    0.9
% % % % %    0.9
% % % % %    0.8
% % % % %    0.8
% % % % %   TauH
% % % % %    4.67
% % % % %KAs channel I = Gmax*m^2*(a*h-(1-a))*(V-Ek); a = 0.996;
% % % % %   Gmax = 0.0104 (9.51*10^-4 dendr)
% % % % %   Minf = 1/(1+exp(-(v+27.0)/16.0)
% % % % %   Hinf = 1/(1+exp( (v+33.5)/21.5)
% % % % %   TauM
% % % % %    Tm = 0.378+9.91*exp((-(v+34.3)/30.1)?2)
% % % % %   TauH
% % % % %    A = exp(-(v+90.96)/29.01)
% % % % %    B =  exp((v+90.96)/100)
% % % % %    Th = 1097.4/(A+B)
% % % % %KRP channel I = Gmax*m*(a*h-(1-a))*(V-Ek); a = 0.7; 
% % % % %   Gmax = 1.5*10^-4 (soma only)
% % % % %   Minf = 1/(1+exp(-(v+13.5)/11.8)
% % % % %   Hinf = 1/(1+exp( (v+54.7)/18.6)
% % % % %   TauM
% % % % %    Tm = 102.981176901981/cosh(-1.12973817613827-0.0148684700109162*V)
% % % % %    Tm = 103/cosh(-(V+76)/67.3) !!!!!!!!!!!! final
% % % % %generic( t0+t/cosh(( v-v12 )/slp )))
% % % % %     40
% % % % %     45
% % % % %     48.8
% % % % %     55
% % % % %     64.4
% % % % %     75
% % % % %     83.9
% % % % %     90
% % % % %     93.5
% % % % %     95
% % % % %     95.4
% % % % %     97
% % % % %     99.2
% % % % %     95
% % % % %     79.7
% % % % %     60
% % % % %     44.5
% % % % %     35
% % % % %     29.3
% % % % %     25
% % % % %     20
% % % % %     15
% % % % %     11.6
% % % % %     10
% % % % %     9.6
% % % % %     10
% % % % %     10.5
% % % % %     10
% % % % %     8
% % % % %     5
% % % % %     5
% % % % %   TauH
% % % % %    Th = 2084.72602496305+7077270.02009008/(1434.2075451696+exp(0.0751944587921081*V))
% % % % %    Th = 2085+4935/((1+exp((V-96.7)/13.3))) !!!!!!!!!!!! final
% % % % %(mod gen t0+2*t/(exp((v-v12)/slp)+exp(-(v-v12_2)/slp_2)) slp_2->inf)
% % % % %     7000.0
% % % % %     7000.0
% % % % %     7000.0
% % % % %     7000.0
% % % % %     7000.0
% % % % %     7000.0
% % % % %     7000.0
% % % % %     7000.0
% % % % %     7000.0
% % % % %     7000.0
% % % % %     7000.0
% % % % %     7000.0
% % % % %     7000.0
% % % % %     7000.0
% % % % %     7000.0
% % % % %     7000.0
% % % % %     7000.0
% % % % %     7000.0
% % % % %     7000.0
% % % % %     7000.0
% % % % %     6742.5
% % % % %     6000.0
% % % % %     4740.2
% % % % %     3500.0
% % % % %     2783.3
% % % % %     2500.0
% % % % %     2336.3
% % % % %     2200.0
% % % % %     2083.5
% % % % %     2000.0
% % % % %     2000.0
% % % % %Leak channel I = Gmax*(V-Eleak); 
% % % % %   Gmax = 11.5*10^-6
% % % % %   Eleak = -70 mV
% % % % 
% % % % %2nd stage
% % % % %Kca channels; Ek = -90 mV
% % % % %BK_KCa ;Gmax = 0.12
% % % % %SK_KCa ;Gmax = 0.1885 (dendrite only)
% % % % 
% % % % %Ca channels; Eca = 140 mV (in average)
% % % % %CaL1.2 channel I = gmax*m^2*(a*h-(1-a))*(V-Eca); a= 0.17; Gmax = 6.7*10^-5
% % % % %CaL1.3 channel I = gmax**m^2*h*(V-Eca); Gmax = 3.19*10^-5; (4.25*10^-6 dendr)
% % % % %CaN channel I = gmax*m^2*(a*h-(1-a))*(V-Eca); a= 0.21; Gmax = 1*10^-4
% % % % %CaQ channel I = gmax*m^2*(V-Eca); Gmax = 6*10^-5
% % % % %CaR channel I = gmax*m^3*h*(V-Eca); Gmax = 2.6*10^-4
% % % % %CaT channel I = gmax*m^3*h*(V-Eca); Gmax = 4*10^-4
