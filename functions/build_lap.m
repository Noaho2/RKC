function L = build_lap( parms )


L = sparse( parms.ntot, parms.ntot );


%diagonal terms
indr = inds_mat2vect( 1: parms.nx, 1 : parms.ny, parms );
indc = indr;
L = L + sparse( indr, indc, -4* ones(size(indc)), parms.ntot, parms.ntot );

%above
indr = inds_mat2vect( 1 : parms.nx, 1 : (parms.ny - 1), parms );
indc = inds_mat2vect( 1 : parms.nx, 2 : parms.ny, parms );
L = L + sparse( indr, indc, 1 * ones(size(indc)), parms.ntot, parms.ntot );

%below
indr = inds_mat2vect( 1 : parms.nx, 2 : parms.ny, parms );
indc = inds_mat2vect( 1 : parms.nx, 1 : (parms.ny - 1), parms );
L = L + sparse( indr, indc, 1 * ones(size(indc)), parms.ntot, parms.ntot );

%left
indr = inds_mat2vect( 2 : parms.nx, 1 : parms.ny, parms );
indc = inds_mat2vect( 1 : (parms.nx - 1), 1 : parms.ny, parms );
L = L + sparse( indr, indc, 1 * ones(size(indc)), parms.ntot, parms.ntot );

%right
indr = inds_mat2vect( 1 : (parms.nx - 1), 1 : parms.ny, parms );
indc = inds_mat2vect( 2 : parms.nx, 1 : parms.ny, parms );
L = L + sparse( indr, indc, 1 * ones(size(indc)), parms.ntot, parms.ntot );


L = L / (parms.dx^2);
