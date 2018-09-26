function parms = add_bc( parms , c)

%incorporate BCs into rhs term

% No time stepping
if strcmp( parms.timestep, 'no' )
    if isfield( parms, 'exact_soln' )
        f_mat = parms.gexact( parms.xx, parms.yy );
        bcs = parms.exact_soln;
    else
        f_mat = parms.g( parms.xx, parms.yy );
        bcs = @(x,y) 0 .* (x .* y );
    end
    scl = -1;
    
% Forward Euler
elseif strcmp( parms.timestep, 'FE' )
    if isfield( parms, 'exact_soln' )
        f_mat = parms.gexact( parms.xx, parms.yy, parms.t );
        bcs = @(x,y) parms.exact_soln( x,y, parms.t - parms.dt );
    else
        f_mat = parms.g( parms.xx, parms.yy, parms.t - parms.dt );%Added '- parms.dt'
        bcs = @(x,y) 0 .* (x .* y );
    end
    scl = 1;
    
% Backward Euler
elseif strcmp( parms.timestep, 'BE' )
    if isfield( parms, 'exact_soln' )
        f_mat = parms.gexact( parms.xx, parms.yy, parms.t );
        bcs = @(x,y) parms.exact_soln( x,y, parms.t );
    else
        f_mat = parms.g( parms.xx, parms.yy, parms.t );
        bcs = @(x,y) 0 .* (x .* y );
    end
    scl = 1;

% Runge-Kutta Fehlberg (does this need to be changed?
%   YES -- I need boundary conditions for 5 different step sizes!
elseif strcmp( parms.timestep, 'RKF45' )
    if isfield( parms, 'exact_soln' )                                   %--------------
        f_mat = parms.gexact( parms.xx, parms.yy, parms.t );            % not updated  |
        bcs = @(x,y) parms.exact_soln( x,y, parms.t);                   %--------------
    else
        f_mat = parms.g( parms.xx, parms.yy, parms.t - (1.0-c) * parms.dt);
        bcs = @(x,y) 0 .* (x .* y );
    end
    scl = 1;
    
end

parms.f = reshape( f_mat', [parms.ntot, 1] );

%bottom
parms.f( inds_mat2vect( 1 : parms.nx, 1, parms) ) =  ...
    parms.f( inds_mat2vect( 1 : parms.nx, 1, parms) ) + scl * ...
    1/(parms.dx^2) * bcs( ...
    parms.x, parms.y(1) - parms.dx )';

%top
parms.f( inds_mat2vect( 1 : parms.nx, parms.ny, parms) ) =  ...
    parms.f( inds_mat2vect( 1 : parms.nx, parms.ny, parms) ) + scl * ...
    1/(parms.dx^2) * bcs( ...
    parms.x, parms.y(end) + parms.dx )';

%left
parms.f( inds_mat2vect( 1 , 1 : parms.ny, parms) ) =  ...
    parms.f( inds_mat2vect( 1, 1 : parms.ny, parms) ) + scl * ...
    1/(parms.dx^2) * bcs(  ...
    parms.x(1) - parms.dx, parms.y )';

%right
parms.f( inds_mat2vect( parms.nx, 1 : parms.ny, parms) ) =  ...
    parms.f( inds_mat2vect( parms.nx, 1 : parms.ny, parms) ) + scl * ...
    1/(parms.dx^2) * bcs( ...
    parms.x(end) + parms.dx, parms.y )';


end




