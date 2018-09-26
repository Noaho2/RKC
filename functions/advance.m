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
        
        if t == 0
        % store Butcher tableau for RKF
            a = [0         0           0           0          0      ;...
                1/4        0           0           0          0      ;...
                3/32       9/32        0           0          0      ;...
                1932/2197  -7200/2197  7296/2197   0          0      ;...
                439/216    -8          3680/513    -845/4104  0      ;...
                -8/27      2           -3544/2565  1859/4104  -11/40];
            
            c = [0       1/4  3/8         12/13        1      1/2];
            b = [25/216  0    1408/2565   2197/4104    -1/5   0; ...
                 16/135  0    6656/12825  28561/56430  -9/50  2/55];
             
         % initialize table of K_i values for each element of u    
            K = 0.*ones(parms.ntot,6);
        end
        
        it = 0;
        while true

        %Compute K_1, K_2, K_3, etc.
            for i = 1:6
                if i == 1
                %initialized to an array of zeros for K_1
                   aK_sum = 0.*ones(parms.ntot, 1); 
                else
                    for j = 1:i-1   
                        aK_sum = aK_sum + a(i, j) .* K(:,j);
                    end
                end
               
                parms = add_bc( parms , c(i)); %incorporate BCs
                K(:,i) = parms.dt .* ( parms.L * ( parms.uold + aK_sum ) + parms.f );
            end

            bK_sum1 = 0; %RK4 b*K sum
            bK_sum2 = 0; %RK5 b*K sum
            for i = 1:6
                bK_sum1 = bK_sum1 + ( b(1,i) .* K(:,i) ); 
                bK_sum2 = bK_sum2 + ( b(2,i) .* K(:,i) ); 
            end 
            
            u_RK4 = parms.uold + bK_sum1;
            u_RK5 = parms.uold + bK_sum2;
            
            s_array = parms.tol .* (0.5 .* abs(u_RK4 - u_RK5)).^(-0.25);
            s_min = min(s_array);
            parms.dt = s_min * parms.dt;
            %Note: We have the info to give unique dt values for each 
            % element of u, but trying to keep everything working together
            % while stepping each element a different amount in time would 
            % be tricky... I would guess do-able though.
            
            if s_min > .99 % if our step size was adequate, keep results, exit loop.
                u = u_RK5;
                break
            else % else if dt was too large, scrap results, try again.
                %Just a catch to make sure we don't get stuck in inf loop
                it = it + 1;
                if it > 100
                    disp("Houston, we have a problem.");
                    exit();
                end
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
    
    %print progress:
    fprintf('\b\b\b\b\b\b');
    fprintf('%.3f\n',t/parms.T);
end
    

