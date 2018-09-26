clear all, clc%, close all, clc

addpath('functions/')

%Solve Poisson problem or heat equation on rectangular domain
%   -uses 2nd order FD stencil for Laplacian
%   -assumes Dirichlet boundary conditions
%   -options for time stepping: 
%           'FE' : Forward Euler
%           'BE' : Backward Euler


%---User specified parameters
    %# of grid points in x and y dirns
    parms.dx = 0.005;

    %specify boundaries of domain:
    parms.xbds = [-0.4 1.2];
    parms.ybds = [-0.2 0.3];
    
    %**are we stepping in time?
        %   options are:
        %           'no'    : no timestepping (i.e. solve Poisson problem)
        %           'FE'    : Forward Euler
        %           'BE'    : Backward Euler
        %           'RKF45' : Runge-Kutta Fehlberg  (under construction)
        parms.timestep = 'RKF45';
        
            % if using RKf45, specify the error control tolerance
            parms.tol = 1.0e-3;

            %only if parms.timestep is something other than 'no':
            parms.T = .01; %How long to run
            parms.dt = 0.001; %time step
            parms.t_save = 0.001; %save every t_save interval
            parms.u0 = @(x,y) (x .* y ); %initial condition
    %**
    
    %Is there a source term? (If not, set to 0. If parms.timestep is not 
    %                         'no', f can be a function of time)
    parms.g = @(x,y,t) sin(6.*x).^2 + y.^2 + cos(14.*t);%0.*x;%sin(6.*x).^2 + y.^2 + cos(14.*t);
    
%---

[parms, soln] = run_solver( parms );


plot_soln( parms, soln );

disp("Done");
