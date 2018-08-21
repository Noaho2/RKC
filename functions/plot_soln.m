function [] = plot_soln( parms, soln)

if strcmp( parms.timestep, 'no' )
    soln.umat = reshape( soln.u, [parms.nx, parms.ny])';


    contourf( parms.xx, parms.yy, soln.umat, 'edgecolor', 'none' )
    shading flat

    xlabel( '$x$', 'fontsize', 16, 'interpreter', 'latex' )
    ylabel( '$y$', 'fontsize', 16, 'interpreter', 'latex'  )

    set(gca, 'fontsize', 14, 'ticklabelinterpreter', 'latex')
    
else

    tvect = parms.t_save : parms.t_save : parms.T;
    for j = 1 : length(tvect)
        soln.umat = reshape( soln.u(:,j), [parms.nx, parms.ny])';


        figure( j )
        contourf( parms.xx, parms.yy, soln.umat, 'edgecolor', 'none' )
        shading flat

        xlabel( '$x$', 'fontsize', 16, 'interpreter', 'latex' )
        ylabel( '$y$', 'fontsize', 16, 'interpreter', 'latex'  )

        set(gca, 'fontsize', 14, 'ticklabelinterpreter', 'latex')
    end
end