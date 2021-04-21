function [R, J, Q, W, Winv, B, H, grad_H, t0, tf, x_0, u, par, f, g]=NWE(N)
% This function creates the constants for the non-linear-network-equation.
% R, J, Q, B are matrices and t0, tf are scalar and x_0 is a vector. H,
% grad_H and u are function-handles depending on z=[q p].
% This function also creates f and g, which are needed to solve this
% equation with the semi-implicit euler method.
% source: [scan from Prof. Dr. Tatjana Stykel &
% Chaturantabut/Beattie/Gugercin section 2.4]

n = 2*N;                                                        % constants

x_0 = sparse(n,1);                                      % initial condition

t0 = 0;                                    % time intervall in 1^-5 seconds
tf = 10;                                   

par.hx = 0.01;                                                  % increment

par.c = 5;                                                 % in 10^-9 Farad
par.v = 1;                                                        % in Volt
par.l = 2;                                                 % in 10^-9 Henry
par.g_con = 10;                                          % in 10^-9 Siemens
par.r = 1;                                                         % in Ohm

% F is used to create parts of J
F = sparse(diag(ones(N,1))) + sparse(diag(-ones(N-1,1),1));

% skew-symmetric Matrix J
J = [sparse(N,N) F; -F' sparse(N,N)];

% symmetric, positive semidefinite Matrix R
R = [par.g_con.*speye(N) sparse(N,N); sparse(N,N) par.r.*speye(N)];

% Port-Matrix B with dimensions (n x m) and m = 2
B = [sparse(1,N) 1 sparse(1,N-1); sparse(1,N-1) 1 sparse(1,N)]';

% symmetric, positive semidefinite Matrix Q
Q = [sparse(N,N) sparse(N,N); sparse(N,N) speye(N)./par.l];



% This function calculates the product of W*z as a function of z
% without saving the matrix W.
% Input: z as a column-vector
% Output: W(z) as a column-vector
W = @(z,~) [sparse(tril(ones(N,N))) sparse(N,N); sparse(N,N) speye(N)]*z;


% This function calculates the product of inv(W)*z as a function of z
% without saving the matrix inv(W).
% Input: z as a column-vector
% Output: Winv(z) as a column-vector
Winv = @(z,~) [sparse(eye(N)-diag(ones(1,N-1),-1)) sparse(N,N); sparse(N,N) speye(N)]*z;

% function handles of u and the Hamiltonian and its Gradient

% entry-function u: [t0,tf] -> R^(m) with m = 2
u = @(t) [sin(t);cos(t)];


% this function calculates the Hamiltonian as a function of z = [q p].
% Input: q and p as column-vectors;
% Output: H(q,p) as a scalar
H = @(q,p) [q; p]'*(Q*[q; p]) + ...  
                       sum(par.c*par.v^2*(exp(q/(par.c*par.v)))-1-par.v*q);


% this function calculates the Gradient of the Hamiltonian as a
% function of z = [q p].
% Input: q and p as column-vectors;
% Output: grad_Ham(q,p) as a column-vector
grad_H = @(q,p) [(par.v*exp(q/(par.c*par.v))-1); p./par.l];


function func = f_function(q,p,t)
    % This function calculates the top half of z=(J-R)*grad_H(z)+B*u and is
    % used in syplectic_euler to solve the equation with the symplectic
    % euler method.
    % Input: q and p as column-vectors
    % Output: f(q,p,t) as a column-vector
    
    z = (J-R)*grad_H(q,p) + B *u(t);
    func = z(1:N);
end
f = @f_function;

function func = g_function(q,p,t)
    % This function calculates the bottom half of z=(J-R)*grad_H(z)+B*u and
    % is used in syplectic_euler to solve the equation with the symplectic
    % euler method.
    % Input: q and p as column-vectors
    % Output: g(q,p,t) as a column-vector
    z = (J-R)*grad_H(q,p) + B*u(t);
    func = z(N+1:end);
end
g = @g_function;

end

