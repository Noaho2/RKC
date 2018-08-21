function parms = build_grid( parms )

%x and y points
parms.x = parms.xbds(1) + parms.dx : parms.dx : parms.xbds(2) - parms.dx;
parms.y = parms.ybds(1) + parms.dx : parms.dx : parms.ybds(2) - parms.dx;

%# of x and y points
parms.nx = length(parms.x);
parms.ny = length(parms.y);
parms.ntot = parms.nx * parms.ny;

%meshgrid of soln
[parms.xx, parms.yy] = meshgrid( parms.x, parms.y );
