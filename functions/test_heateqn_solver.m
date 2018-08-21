clear all, close all, clc

%Test order of convergence of Poisson solver
%Should be second order

%---User specified parameters
    %# of grid points in x and y dirns
    parms.dxvect = 0.0025;%[0.02, 0.01, 0.005];

    %specify boundaries of domain:
    parms.xbds = [-0.4 1.2];
    parms.ybds = [-0.2 0.3];
    
    %**are we stepping in time?
        %   options are:
        %           'no' : no timestepping
        %           'FE' : Forward Euler
        %           'BE' : Backward Euler
        parms.timestep = 'BE';

            %only if parms.timestep is something other than 'no':
            parms.T = 1; %How long to run
            parms.dtvect = [4 2 1 1/2]* parms.dxvect; %time step
            parms.t_save = parms.T; %save every t_save interval
    %**
    
    %specify exact soln
    parms.exact_soln = @(x,y,t) x.^2 + 4*x.^3.*y.^4 + y.^2 + t.^4.*x.*y;
    
    %Is there a source term? (if not set to 0)
    parms.gexact = @(x,y,t) -4 - 24.*x.*y.^4 - 16 * 3 * y.^2 .* x.^3 + ...
        4*t.^3.*x.*y ;
    
    %Is there an initial condition? (only relevant when parms.timestep is
    %set to something other than 'no').
    parms.u0 = @(x,y) parms.exact_soln( x, y, 0 );
%---


err = zeros( size(parms.dtvect) );
for j = 1 : length( parms.dtvect )
    
    parms.dt = parms.dtvect( j );
    parms.dx = parms.dxvect;%( j );
    
    [parms, soln] = run_solver( parms );
    
    tfinal = soln.t(end);
    
    exact_soln = parms.exact_soln( parms.xx, parms.yy, tfinal );
    exact_soln = reshape( exact_soln', [parms.ntot, 1] );
    
    err(j) = max( abs( soln.u(:, end) - exact_soln ) );
    
    errs = err( j )

end



loglog( parms.dtvect, err, 'ko' ),  hold on
loglog( parms.dtvect, err(end)/(parms.dtvect(end)) * parms.dtvect, 'k--' )

set( gca, 'fontsize', 16, 'TickLabelInterpreter', 'latex' )

xlabel( '$\Delta t$', 'interpreter', 'latex', 'fontsize', 18 )
ylabel( '$||u - u_{exact}||_\infty$', 'interpreter', 'latex', 'fontsize', 18 )


