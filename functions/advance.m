function [ parms, soln ] = advance( parms )

parms = get_IC( parms );
count = 1;

t = 0;
fprintf('\nFraction complete:\n0.000\n');
t_save_it = parms.t_save;
while t < parms.T
    
    parms.t = t;
    
    if strcmp( parms.timestep, 'FE' ) %forward euler
        
        parms = add_bc( parms ); %incorporate BCs
        u = parms.uold + parms.dt * ( parms.L*parms.uold + parms.f );
        
    elseif strcmp( parms.timestep, 'BE' ) %backward euler
        
        parms = add_bc( parms ); %incorporate BCs
        %precompute and store LU-factorization of L
        if t == 0
            [LL,UU,pp,qq,rr] = lu( speye(parms.ntot) - parms.dt * parms.L);
            Linvfun = @(vect) qq*(UU\(LL\(pp*(rr\ vect ) ) ) );
        end
        
        u = Linvfun( parms.uold + parms.dt * parms.f );
        
        
        
    elseif strcmp( parms.timestep, 'RKF45' ) %Runge-Kutta Fehlberg 
        
        % store Butcher tableau for RKF
        if t == 0
            a = [0         0           0           0          0      ;...
                1/4        0           0           0          0      ;...
                3/32       9/32        0           0          0      ;...
                1932/2197  -7200/2197  7296/2197   0          0      ;...
                439/216    -8          3680/513    -845/4104  0      ;...
                -8/27      2           -3544/2565  1859/4104  -11/40];
            
            c = [0       1/4  3/8         12/13        1      1/2];
            b = [25/216  0    1408/2565   2197/4104    -1/5   0; ...
                 16/135  0    6656/12825  28561/56430  -9/50  2/55];
        end
        
        while true
            
            %Compute the 6 K_i,j values for each element of the soln u_j
            for i = 1:6
                % aK_sum = a_1*K_1,j + a_2*K_2,j + ... + a_i*K_i,j
                %    initialized to a 1D array of zeros for K_1,j
                if i == 1
                   aK_sum = 0.0*ones(parms.ntot, 1); 
                else
                   for j = 1:parms.ntot
                       aK_sum(j) = sum( a(i, 1:i-1) .* K(j, 1:i-1) );
                   end
                end
                
                parms = add_bc( parms , c(i)); %incorporate BCs
                K(:, i) = parms.dt .* ( parms.L*(parms.uold + aK_sum) + parms.f );
              
            end

            % Analyze error between solutions and adjust dt accordingly by a scalar factor s.
            %   Note: ideally we would have seperate dt for each solution, but for
            %   now we'll just leave it as one dt for all elements of u.
            % Find smallest value of s for all u_i :
            s_old = 1;
            for i = 1:parms.ntot
                u_RK4(i) = parms.uold(i) + sum( b(1, :) .* K(i,:) );
                u_RK5(i) = parms.uold(i) + sum( b(2, :) .* K(i,:) );
        
                s = (parms.tol/(.5*abs(u_RK4(i) - u_RK5(i))))^0.25;
                if s < s_old
                    s_old = s;
                end
            end
            
            parms.dt = s_old*parms.dt;
            
            if s > .99 % if dt was adequate, keep results, exit loop.
                u = u_RK5';
                break
            else       % else if dt was too large, scrap results, try again.
            end
        end
        
    end
    
    parms.uold = u;
    
    if t > t_save_it || t_save_it >= parms.T
        
        t_save_it = t_save_it + parms.t_save;
        soln.uold = parms.uold;
        soln.t(count) = t;
        soln.u(:,count) = u;
        count = count + 1;
        
    end
    t = t + parms.dt;
    
    %print progress
    fprintf('\b\b\b\b\b\b');
    fprintf('%.3f\n',t/parms.T);
end
    

