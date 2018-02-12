clear all

    nf = 1;
    V = -150:0.1:50;
    M = zeros( 2, size( V, 2));
    H = zeros( 2, size( V, 2));
    Tm = zeros( 2, size( V, 2));
    Th = zeros( 2, size( V, 2));
    VT_na = -65;
    VT = -45;
%--- NaF channels (type 1 vs type 2 )
%--- steady-state type 1
    alpha_m = -0.32*(V-VT_na-13)./(exp(-(V-VT_na-13)/4)-1);
    beta_m = 0.28*(V-VT_na-40)./(exp((V-VT_na-40)/5)-1);
    alpha_h = 0.128.*(exp(-(V-VT_na-17)/18));
    beta_h = 4./(exp(-(V-VT_na-40)/5)+1);
%     alpha_m = -0.32*(V-13.1)./(exp(-(V-13.1)/4)-1);
%     beta_m = 0.28*(V-40.1)./(exp((V-40.1)/5)-1);
%     alpha_h = 0.128./exp((V-17)/18);
%     beta_h = 4./(exp(-(V-40)/5)+1);
    
    M( 1, : ) = alpha_m./(alpha_m+beta_m);
    Tm( 1, : ) = 1./(alpha_m+beta_m);
    H( 1, : ) = alpha_h./(alpha_h+beta_h);
    Th( 1, : ) = 1./(alpha_h+beta_h);
%--- steady-state type 2
    alpha_m = -0.1*(V+33)./(exp(-(V+33)/10)-1);
    beta_m = 4./exp((V+58)/12);
    alpha_h = 0.07./(exp((V+50)/10));
    beta_h = 1./(exp(-(V+20)/10)+1);
    M( 2, : ) = alpha_m./(alpha_m+beta_m);
    Tm( 2, : ) = 1./(alpha_m+beta_m);
    H( 2, : ) = alpha_h./(alpha_h+beta_h);
    Th( 2, : ) = 1./(alpha_h+beta_h);
%--- plot
    nf = plot_steadystate( 'NaF channels', nf, V, M, H, Tm, Th );
%--- Kdr channels
%--- steady-state
%--- steady-state type 1
    alpha_m = -0.32*(V-VT-15)./(exp(-(V-VT-15)/5)-1);
    beta_m = 0.5./exp((V-VT-10)/40);
%     alpha_m = -0.016*(V-35.1)./(exp(-(V-35.1)/5)-1);
%     beta_m = 0.25./exp((V-20)/40);
    M( 1, : ) = alpha_m./(alpha_m+beta_m);
    Tm( 1, : ) = 1./(alpha_m+beta_m);
%--- steady-state type 2
    alpha_m = -0.01*(V+34)./(exp(-(V+34)/10)-1);
    beta_m = 0.125./exp((V+44)/25);
    M( 2, : ) = alpha_m./(alpha_m+beta_m);
    Tm( 2, : ) = 1./(alpha_m+beta_m);
   
    H( 1, : ) = nan;
    Th( 1, : ) = nan;
    H( 2, : ) = nan;
    Th( 2, : ) = nan;
%--- plot
    nf = plot_steadystate( 'Kdr channels', nf, V, M, H, Tm, Th );
%--- Kir channels
    H( 1, : ) = nan;
    H( 2, : ) = nan;
    M( 1, : ) = nan;
    M( 2, : ) = nan;
%--- time constants
    Th( 1, : ) = nan;
    Th( 2, : ) = nan;
    Tm( 1, : ) = nan;
    Tm( 2, : ) = nan;
%--- plot
%   nf = plot_steadystate( 'Kir channels', nf, V, M, H, Tm, Th );
%--- Kaf channels
    M( 1, : ) = nan;
    H( 1, : ) = nan;
    Tm( 1, : ) = nan;
    Th( 1, : ) = nan;

    M( 2, : ) = nan;
    H( 2, : ) = nan;
    Tm( 2, : ) = nan;
    Th( 2, : ) = nan;
%--- plot
%    nf = plot_steadystate( 'Kaf channels', nf, V, M, H, Tm, Th );
