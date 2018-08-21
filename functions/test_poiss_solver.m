clear all, close all, clc

%Test order of convergence of Poisson solver
%Should be second order

%---User specified parameters
    %# of grid points in x and y dirns
    parms.dxvect = [0.02, 0.01, 0.005, 0.0025];

    %specify boundaries of domain:
    parms.xbds = [-0.4 1.2];
    parms.ybds = [-0.2 0.3];
    
    %are we stepping in time?
    parms.timestep = 'no';
    
    %specify exact soln
    parms.exact_soln = @(x,y) x.^2 + 4*x.^3.*y.^4 + y.^2;
    
    %Is there a source term? (if not set to 0)
    parms.gexact = @(x,y) 4 + 24.*x.*y.^4 + 16 * 3 * y.^2 .* x.^3;
    
    %Is there an initial condition? (only relevant when parms.timestep is
    %set to something other than 'no').
    parms.u0 = @(x,y) 0;
%---

err = zeros( size(parms.dxvect) );
for j = 1 : length( parms.dxvect )
    
    parms.dx = parms.dxvect( j );
    
    [parms, soln] = run_solver( parms );
        
    exact_soln = parms.exact_soln( parms.xx, parms.yy);
    exact_soln = reshape( exact_soln', [parms.ntot, 1] );
    
    err(j) = max(max(abs( soln.u - exact_soln )));
    
    errs = err( j )

end

loglog( parms.dxvect, err, 'ko' ),  hold on
loglog( parms.dxvect, err(end)/(parms.dxvect(end)^2) * parms.dxvect.^2, 'k--' )

set( gca, 'fontsize', 16, 'TickLabelInterpreter', 'latex' )

xlabel( '$\Delta x$', 'interpreter', 'latex', 'fontsize', 18 )
ylabel( '$||u - u_{exact}||_\infty$', 'interpreter', 'latex', 'fontsize', 18 )