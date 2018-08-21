function [parms, soln] = run_solver( parms )

%build grid from user info
parms = build_grid( parms );

%build laplacian
parms.L = build_lap( parms );

if strcmp( parms.timestep, 'no' )
    
    parms = add_bc( parms );
    soln.u = parms.L \ parms.f;
    
else
    
    [parms, soln] = advance( parms );
    
end
