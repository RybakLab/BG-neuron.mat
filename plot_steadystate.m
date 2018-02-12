%--- NaF channels
function NF = plot_steadystate( name, nf, V, m, h, tm, th )
    NF = nf;
    figure( NF ); 
    subplot( 2, 1, 1);
    hold on
    plot( V, m( 1, : ), 'r-', 'LineWidth', 2 );
    plot( V, m( 2, : ), 'b-', 'LineWidth', 1 );
    plot( V, h( 1, : ), 'r--', 'LineWidth', 2 );
    plot( V, h( 2, : ), 'b--', 'LineWidth', 1 );
    title( strcat( name,': activation & inactivation') );
    hold off
    subplot( 2, 1, 2);
    hold on
    plot( V, tm( 1, : ), 'r-', 'LineWidth', 2 );
    plot( V, tm( 2, : ), 'b-', 'LineWidth', 1 );
    plot( V, th( 1, : ), 'r--', 'LineWidth', 2 );
    plot( V, th( 2, : ), 'b--', 'LineWidth', 1 );
    title( strcat( name,': time constants for activation & inactivation') );
    hold off
    NF = NF+1;
