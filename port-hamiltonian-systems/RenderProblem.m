function [R, J, Q, W, Winv, B, H, grad_H, t0, tf, x_0, u, par] = ...
                                                 RenderProblem(name, solve)
% Port-hamiltonian-Systems
% Alexandro Jedaidi
% 30.09.2020
%
% This function is based on the system (z)=(J-R)*gradient(H(x))+B*u 
%
% Input:
% name - String, decides which equation has to be created
% solve - boolean, decides if selected equation has to be solved
%
% Output:
% J, skew-symmetric matrix of dimension R^(2n)x(2n)
% R, positive semi-definite, symmetric matrix of dimension R^(2n)x(2n)
% Q, positive semi-definite, symmetric matrix of dimension R^(2n)x(2n)
% B, port-matrix, of dimension R^(2n)x(m)
% H, function handle, R^(2n) -> R and a matrix B
% grad_H, function handle, R^(2n) -> R^(2n)
% W, function handle, R^(2n) -> R^(2n)
% Winv, function handle, R^(2n) -> R^(2n)
% u, function handle, [t0,tf] -> R^(m)
% t0, scalar
% tf, scalar
% x_0, column-vector of dimension R^(2n)
% par, physical parameter, see individual equations for mor info

switch name %switches between different algorithms 
    case 'DWE'          % 1-D wave equation; source see own function

        N = 500;
        % creates equation with all its outputs       
        [R, J, Q, W, Winv, B, H, grad_H, t0, tf, x_0, u, par, f, g] ...
                                                          = DWE(N);
        
        % checks definition and dimension of given matrices
        check_definition(R, J, Q, B, N);                                              
        if solve == 1                          
            h = par.ht;
            % solves given problem with semi-implicit euler       
            symplectic_euler(H, t0, tf, x_0, N, f, g, h);      
        end
        disp('done')

    case 'NWE_linear'                         % non-linear network-equation

        N = 500;
        % creates equation with all its outputs
        [R, J, Q, W, Winv, B, H, grad_H, t0, tf, x_0, u, par, f, g]=...
                                               NWE_linear(N);

        %check_definition(R, J, Q, B, N); 
        
        if solve == 1
            h = 0;
            % solves given problem with semi-implicit euler
            symplectic_euler(H, t0, tf, x_0, N, f, g, h);
        end
        
        disp('done')   
        
    case 'NWE'                                    % linear network-equation

        N = 500;
        % creates equation with all its outputs
        [R, J, Q, W, Winv, B, H, grad_H, t0, tf, x_0, u, par, f, g] ...
                                                 = NWE(N);
        
        check_definition(R, J, Q, B, N);
        
        if solve == 1
            h = 0;
            % solves given problem with semi-implicit euler        
            symplectic_euler(H, t0, tf, x_0, N, f, g, h);
        end
        
        disp('done') 
    case 'Schroedinger'                             % schroedinger-equation
        
        N = 256;
        % creates equation with all its outputs
        [R, J, Q,  W, Winv, B, H, grad_H, t0, tf, x_0, u, par, f, g]=...
                                                  Schroedinger(N);

        check_definition(R, J, Q, B, N);
        
        if solve == 1
            h = 0;
            % solves given problem with semi-implicit euler        
            %Solve_schroedinger(J, H, grad_H, t0, tf, x_0, par, N);
            symplectic_euler(H, t0, tf, x_0, N, f, g, h);
        end
        
        disp('done')
        
    case 'SGE'                                       % sine-gordon-equation

        N = 500;
         % creates equation with all its outputs
        [R, J, Q, W, Winv, B, H, grad_H, t0, tf, x_0, u, par, f, g]...
                                                                = SGE(N);
        
        check_definition(R, J, Q, B, N);
    
        if solve == 1
            h = par.ht;
            % solves given problem with semi-implicit euler 
            symplectic_euler(H, t0, tf, x_0, N, f, g, h);
        end
        
        disp('done')
        
    case 'TLS'                                        % Toda-Lattice-System
    
        % creates list to select between two options for input function u
        list = {'0.1','0.1*sin(t)'};
        [indx,~] = listdlg('SelectionMode','single','ListString',list);

        sel_u = indx;  % indx==1 means u=0.1 and indx==2 means u=0.1*sin(t)
        
        
        N = 500;
        % creates equation with all its outputs
        [R, J, Q, W, Winv, B, H, grad_H, t0, tf, x_0, u, par, f, g] ...
                                                           = TLS(N, sel_u);
        
        check_definition(R, J, Q, B, N);
             
        if solve == 1
            h = 0;
            % solves given problem with semi-implicit euler  
            symplectic_euler(H, t0, tf, x_0, N, f, g, h);
        end
        
        disp('done')
        
    otherwise %occurs, if no valid algorithm was selected
        disp('no valid method')
   
end
