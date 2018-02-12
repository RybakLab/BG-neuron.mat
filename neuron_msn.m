clear all
global Iinj1 Iinj2 T0_inj Tmax
    msn_init;

%setup initial conditions for neurons
    T0_inj = 200;   %start to apply injected current
    Tmax = 2200;
    time = 0:0.1:Tmax;

    Iinj1 = -0.31; %min = -0.31 max -32; Fmin = 2; Fmax = 120
    Iinj2 = -1.245; %min = -1.25 max -45; Fmin = 2; Fmax = 100

    x0 = zeros(1, 1+MAX_MH);
    x0(1) = -68; 
    [t,n1] = ode15s( 'df_msn1',time, x0 );
    x0 = zeros(1, 1+MAX_MH);
    x0(1) = -78; 
    [t,n2] = ode15s( 'df_msn2',time, x0 );
    
    figure(5);
    subplot( 2, 1, 1);
    hold on
    plot( time, n1(:,1), 'r-', 'LineWidth', 2 );
    hold off
    subplot( 2, 1, 2);
    hold on
    plot( time, n2(:,1), 'b-', 'LineWidth', 2 );
    hold off
