function parms = get_IC( parms )


uold = parms.u0( parms.xx , parms.yy );
parms.uold = reshape( uold', [parms.ntot, 1] );



