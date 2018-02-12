clear all

    nf = 1;
    V = -150:0.1:50;
    M = zeros( 2, size( V, 2));
    H = zeros( 2, size( V, 2));
    Tm = zeros( 2, size( V, 2));
    Th = zeros( 2, size( V, 2));

%--- NaF channels (type 1 vs type 2 )
%--- steady-state
    M( 1, : ) = 1./(1+exp(-(V+23.9)/11.8));
    M( 2, : ) = 1./(1+exp(-(V+25)/9.2));
    H( 1, : ) = 1./(1+exp((V+62.9)/10.7));
    H( 2, : ) = 1./(1+exp((V+62)/6));
%--- time constants
    Tm( 1, : ) = 0.1+0.2./cosh((V+1.8)/29);
    Th( 1, : ) = 0.3+1./(1+exp((V-8.3)/18));
    Tm( 2, : ) = 0.1+0.5./cosh((V+73.5)/18.8);
    Th( 2, : ) = 0.3+1.25./(1+exp((V+30)/8.5));
%--- plot
    nf = plot_steadystate( 'NaF channels', nf, V, M, H, Tm, Th );
%--- Kdr channels
%--- steady-state
%    Kdr (type I) is unstable solution. It sould be excluded from consideration
%     kdr_a = -(0.616+0.014*V)./exp(-( V+44 )/2.31-1); 
%     kdr_b = 0.0043./exp((V+44)/34); 
%     M( 1, : ) = kdr_a./(0.1+kdr_a+kdr_b);
%     Tm( 1, : ) = 1./(0.1+kdr_a+kdr_b);
    kdr_a = exp(-(V+13)./9.09); kdr_b = exp(-(V+13)./12.5);
    M( 1, : ) =  1./(1+kdr_a);
    Tm( 1, : ) = 50*kdr_b./(1+kdr_a);
    alpha =  1./exp(-(V+13)./12.5);
    beta = 1./exp((V+13)./33.3);
    M( 2, : ) = alpha./(alpha+beta);
    Tm( 2, : ) = 50*1./(alpha+beta);
    
    H( 1, : ) = nan;
    Th( 1, : ) = nan;
    H( 2, : ) = nan;
    Th( 2, : ) = nan;
%--- plot
    nf = plot_steadystate( 'Kdr channels', nf, V, M, H, Tm, Th );
%--- Kir channels
    H( 1, : ) = 1./(1+exp((V+82.0)/13));
    H( 2, : ) = 1./(1+exp((V+102)/13));
    M( 1, : ) = nan;
    M( 2, : ) = nan;
%--- time constants
    Th( 1, : ) = (3.7+3.96./cosh((V+30)/41.5)); 
    Th( 2, : ) = (6.2./cosh(( V+24)/49)); 
    Tm( 1, : ) = nan;
    Tm( 2, : ) = nan;
%--- plot
    nf = plot_steadystate( 'Kir channels', nf, V, M, H, Tm, Th );
%--- Kaf channels
    M( 1, : ) = 1./(1+exp(-(V+10.0)/17.7));  % Kaf activation
    H( 1, : ) = 1./(1+exp( (V+75.6)/10.0));  % Kaf inactivation
    Tm( 1, : ) = (0.85+0.018./exp(V/51.34)); % Time constant for Kaf activation
    Th( 1, : ) = 14;                         % Time constant for Kaf inactivation

    kafm_a = 1.5./(1+exp(-(V-4)/17));
    kafm_b = 0.6./(1+exp( (V-10)/9));
    kafh_a = 0.105./(1+exp( (V+121)/22));
    kafh_b = 0.065./(1+exp(-(V+55)/11.0));
    M( 2, : ) = kafm_a./(kafm_a+kafm_b);  % Kaf activation
    H( 2, : ) = kafh_a./(kafh_a+kafh_b);  % Kaf inactivation
    Tm( 2, : ) = 1./(kafm_a+kafm_b);      % Time constant for Kaf activation
    Th( 2, : ) = 1./(kafh_a+kafh_b);      % Time constant for Kaf inactivation
%--- plot
    nf = plot_steadystate( 'Kaf channels', nf, V, M, H, Tm, Th );


