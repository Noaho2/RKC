function [ parms, soln ] = advance( parms )

parms = get_IC( parms );
count = 1;

tvect =  parms.dt : parms.dt : parms.T;

for t = tvect
    
    parms.t = t;
    parms = add_bc( parms ); %incorporate BCs
    
    if strcmp( parms.timestep, 'FE' ) %forward euler
        
        u = parms.uold + parms.dt * ( parms.L*parms.uold + parms.f );
        
    elseif strcmp( parms.timestep, 'BE' ) %backward euler
        
        %precompute and store LU-factorization of L
        if t == parms.dt
            [LL,UU,pp,qq,rr] = lu( speye(parms.ntot) - parms.dt * parms.L);
            Linvfun = @(vect) qq*(UU\(LL\(pp*(rr\ vect ) ) ) );
        end
        
        u = Linvfun( parms.uold + parms.dt * parms.f );

        
    end
    
    parms.uold = u;
    
    if mod( t, parms.t_save ) == 0 | t == tvect(end)
        
        soln.uold = parms.uold;
        soln.t(count) = t;
        soln.u(:,count) = u;
        count = count + 1;
        
    end
    
end
    

