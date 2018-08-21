function ivect = inds_mat2vect( imat, jmat, parms )

%take indices on x,y grid and convert them to the appropriate vector index


indx = repmat(imat,[1,length(jmat)]);
indy = repelem(jmat, length(imat));

ivect = (indy-1) * parms.nx + indx;

