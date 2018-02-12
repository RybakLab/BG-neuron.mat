% The simplified model of medium spiny neuron in dorsal striatum of
% the basal ganglia in the mouse (slice preparation) based on research
% conducting by Blackwell1 lab.
% The Effects of NMDA Subunit Composition on Calcium Influx and Spike 
% Timing-Dependent Plasticity in Striatal Medium Spiny Neurons 
% Rebekah C. Evans, Teresa Morera-Herreras, Yihui Cui, Kai Du, Tom Sheehan,
% Jeanette Hellgren Kotaleski, Laurent Venance, Kim T. Blackwell1
% PLoS Comput Biol 8(4): e1002493. doi:10.1371/journal.pcbi.1002493 2012
% The description of all ion channelse (including tabulated time constants) 
% is below the Matlab code. The parameters of the model were taken from 
% https://senselab.med.yale.edu/ModelDB/ShowModel.cshtml?model=145917&file=/MSPNv7.6/#tabs-2
function dx = df_msn2(t,x)
% global variables
    global Iinj2 T0_inj 
    global Ena Ek Eleak2
    global Gnaf2 Gkdr2 Gkir2 Gkaf2 Gleak2
    global iNaF_M iNaF_H iKdr_M iKir_M iKaf_M iKaf_H MAX_MH
    global MH_inf TAU_mh MH
% initialization
    MH = x(2:1+MAX_MH);       % array for activation&inactivation variables
    v = x(1);                 % membrane potentials
    i_inj = 0; 
    if( t > T0_inj )
        i_inj = Iinj2-(t-T0_inj)*0.04;
    end

% NaF activation (m - generic; h - modified generic (Slope2 for Tm -> -inf))
    MH_inf(iNaF_M) = 1/(1+exp(-(v+25)/9.2));    % NaF activation
    TAU_mh(iNaF_M) = 0.1+0.5/cosh((v+73.5)/18.8);% Time constant for NaF ictivation
% original parameters for Th are not working
%     TAU_mh(iNaF_H) = (0.3+1.25/(1+exp((v+30)/8.5))); % Time constant for NaF inactivation
    MH_inf(iNaF_H) = 1/(1+exp((v+62)/6));       % NaF inactivation
    TAU_mh(iNaF_H) = 0.3+1.25/(1+exp((v+15)/8.5)); % Time constant for NaF inactivation
    TAU_mh(iNaF_M) = TAU_mh(iNaF_M)/1.2;        % 1.2 times faster for mouse (in vitro)
    TAU_mh(iNaF_H) = TAU_mh(iNaF_H)*2;          % 1.2 times faster for mouse (in vitro) - not working
%   Inaf = Gnaf2*MH(iNaF_M)^3*MH(iNaF_H)*(v-Ena)

%Kdr channel (m - alpha/beta)
%     kdr_a = exp(-(v+13)/9.09); kdr_b = exp(-(v+13)/12.5);
%     MH_inf(iKdr_M) = 1/(1+kdr_a);               % Kdr activation
%     TAU_mh(iKdr_M) = 50*kdr_b/(1+kdr_a);        % Time constant for Kdr activation
    alpha =  1/exp(-(v+13)/12.5);
    beta = 1/exp((v+13)/33.3);
    MH_inf(iKdr_M) = alpha/(alpha+beta);        % Kdr activation
    TAU_mh(iKdr_M) = 50*1/(alpha+beta);         % Time constant for Kdr activation
    TAU_mh(iKdr_M) = TAU_mh(iKdr_M)*2;          % 2 times slower for mouse (in vitro)
%   Ikdr = Gkdr2*MH(iKdr_M)^4*(v-Ek)

%Kir channel (m - generic)
    MH_inf(iKir_M) = 1/(1+exp((v+102)/13));     % Kir activation
    TAU_mh(iKir_M) = 6.2/cosh(( v+24)/49);      % Time constant for Kir activation (in vivo) 
    TAU_mh(iKir_M) = TAU_mh(iKir_M)*2;          % 2 times slower for mouse (in vitro)
%   Ikir = Gkir2*MH(iKir_M)*(v-Ek)

%Kaf channel (m & h - alpha/beta)
    kafm_a = 1.5/(1+exp(-(v-4)/17));
    kafm_b = 0.6/(1+exp( (v-10)/9));
    kafh_a = 105/(1+exp( (v+121)/22));
    kafh_b = 65/(1+exp(-(v+55)/11.0));
    MH_inf(iKaf_M) = kafm_a/(kafm_a+kafm_b);    % Kaf activation
    MH_inf(iKaf_H) = kafh_a/(kafh_a+kafh_b);    % Kaf inactivation
    TAU_mh(iKaf_M) = 1/(kafm_a+kafm_b);         % Time constant for Kaf activation (in vivo)
    TAU_mh(iKaf_H) = 1/(kafh_a+kafh_b);         % Time constant for Kaf inactivation (in vivo)
    TAU_mh(iKaf_M) = TAU_mh(iKaf_M)/1.5;        % 1.5 times faster for mouse (in vitro)
    TAU_mh(iKaf_H) = TAU_mh(iKaf_H)*1.5;        % 1.5 times slower for mouse (in vitro)
%   Ikaf = Gkaf2*MH(iKaf_M)^2*MH(iKaf_H)*(v-Ek)

%Leak channel
%   Ileak = Gleak2*(v-Eleak)

%DIFFERENTIAL EQUATIONS 
    dx = 0*x;
%   for activation and inactivation
    dx(2:1+MAX_MH) = (-MH+MH_inf)./TAU_mh;
%   for membrane potentials    
    dx(1) = -Gnaf2*MH(iNaF_M)^3*MH(iNaF_H)*(v-Ena)-Gkdr2*MH(iKdr_M)^4*(v-Ek)...
            -Gkir2*MH(iKir_M)*(v-Ek)-Gkaf2*MH(iKaf_M)^2*MH(iKaf_H)*(v-Ek)...
            -Gleak2*(v-Eleak2)-i_inj;
    
%-- THE END

% % % % NaF current
% % % % %    I = Gmax*m^3*h*(V-Ena)
% % % % %    Naf_m = 1/(1+(exp(-(V+25)/9.2)))
% % % % %    naf_h = 1/(1+(exp((V+62)/6)))
% % % % %    Tm = 1.2*(0.1+0.2/cosh((1.8+V)/29))
% % % % %    Th = 1.2*(0.3+1.25/(1+exp((30+V)/8.5)))
% % % % Time constants
% % % % Tm = 0.1+0.52509376022924/cosh(3.91591626832967 + 0.0533116767718528*V)
% % % %      Tm = 0.1+0.5/cosh((73.5+v)/18.8) !!! final
% % % % Time constant Tm(ms) defined as table
% % % % The table corresponds to -100 mV to 50 mV
% % % % 0.3162
% % % % 0.3162
% % % % 0.3162
% % % % 0.4074
% % % % 0.6166
% % % % 0.3548
% % % % 0.2399
% % % % 0.1585
% % % % 0.1047
% % % % 0.0871
% % % % 0.0851
% % % % 0.0813
% % % % 0.0832
% % % % 0.0832
% % % % 0.0832
% % % % 0.0832
% % % % Inactivation 
% % % % Th = 0.290914865865996+0.0372481099913059/(0.0300679383762799+exp(0.11811537092687*V))
% % % %      Th = 0.3+1.25/(1+exp((30+V)/8.5)) !!! final
% % % % Time constant Th(ms) defined as table
% % % % The table corresponds to -100 mV to 50 mV
% % % % 1.5
% % % % 1.5
% % % % 1.5
% % % % 1.5
% % % % 1.5
% % % % 1.5
% % % % 1.5136
% % % % 0.6761
% % % % 0.5129
% % % % 0.4365
% % % % 0.3715
% % % % 0.3388
% % % % 0.2951
% % % % 0.2884
% % % % 0.2754
% % % % 0.2754

% % % % %Kdr channel                                                             ON
% % % % %    I = Gmax*m^4*(V-Ek);
% % % % %    A = exp((V+13)/9.09)
% % % % %    B = exp((V+13)/12.5)
% % % % %    kdr_m =  1/(1+A) 
% % % % %    kdr_tm = 0.5*0.05*B/(1-a)

% % % % Kir current
% % % % Activation 
% % % %      Minf = 1/(1+exp((v+102.0)/13)
% % % % Tm = 2.90015672950307+0.0164576601103018*V+(15.6873447633184+0.18600937955856*V)/cosh(3.02441191029527+0.0495169223038139*V) !!! the best fit
% % % % Tm = 6.18813242007527/cosh( 0.497049538677394+0.0204878944366467*V)
% % % %      Tm = 6.2/cosh(( V+24)/49) !!! final
% % % % Time constant Tm(ms)
% % % % The table corresponds to -150 mV to 50 mV
% % % % 0.2
% % % % 0.2
% % % % 0.2
% % % % 0.2
% % % % 0.2
% % % % 0.38
% % % % 0.97
% % % % 1.486
% % % % 5.3763
% % % % 6.0606
% % % % 6.8966
% % % % 7.6923
% % % % 7.1429
% % % % 5.8824
% % % % 4.4444
% % % % 4.0
% % % % 4.0
% % % % 4.0
% % % % 4.0
% % % % 4.0
% % % % 4.0

% % % % %KAf channel                                                             ON
% % % % %    I = Gmax*m^2*h*(V-Ek);
% % % % %    Gmax = 
% % % % %    kaf_ma = 1.5/(1+exp(-(v-4)/17))
% % % % %    kaf_mb = 0.6/(1+exp( (v-10)/9))
% % % % %    kaf_ha = 105/(1+exp( (v+121)/22))
% % % % %    kaf_hb = 65/(1+exp(-(v+55)/11.0))
% % % % %    kaf_m = kafm_a/(kafm_a+kafm_b)
% % % % %    kaf_h = kafh_a/(kafh_a+kafh_b)
% % % % %    kaf_tm = 1.5*(1/(kafm_a+kafm_b))
% % % % %    kaf_th = 0.67*(1/(kafh_a+kafh_b))

% % % % %Leak channel I = Gmax*(V-Eleak); 
% % % % %   Gmax = 11.5*10^-6
% % % % %   Eleak = -70 mV

