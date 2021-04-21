function [R_m, J, Q, W, Winv, B, H, grad_H, t0, tf, x_0, u, par, f, g] ...
                                              = NWE_linear(N)
% This function creates the constants for the linear-network-equation.
% R, J, Q, B are matrices and t0, tf are scalar and x_0 is a vector.
% H, grad_H and u are function-handles depending on z=[q p].
% This function also creates f and g, which are needed to solve this
% equation with the semi-implicit euler method.
% source: [e-mail from Prof. Dr. Tatjana Stykel &
% Afkham/Hesthaven17 section 4.3]

n = 2*N;                                                        % constants

t0 = 0;
tf = 10;                                   % time intervall in 1^-5 seconds

par.hx = 0.01;                                                  % increment

x_0 = sparse(n,1);                                      % initial condition

par.R = repmat(0.2,N,1);                                           % in Ohm
par.C = ones(N,1);                                               % in Farad
par.l = ones(N,1);                                               % in Henry
par.R_l = 0.4;                                                     % in Ohm

% R_a is used to create R_m
R_a = [par.R(1:N-1); (par.R(N) + par.R_l)];

% symmetric, positive semidefinite Matrix R with dimension (n x n)
R_m = [sparse(N,N) sparse(N,N); sparse(N,N) sparse(diag(R_a))];

% Port-Matrix B with dimensions (n x m) and m = 1
B = [1; sparse(n-1,1)];

% S is used to create J
S = sparse(diag(-ones(N,1)) + diag(ones(N-1,1),-1));

% skew-symmetric Matrix J
J = [sparse(N,N) S; -S' sparse(N,N)];

% symmetric, positive semidefinite Matrix Q
Q = [sparse(diag(1./par.C)) sparse(N,N); sparse(N,N) diag(1./par.l)];


% This function calculates the product of W*z as a function of z
% without saving the matrix W.
% Input: z as a column-vector
% Output: W(z) as a column-vector
W = @(z,~)[sparse(tril(ones(N,N))) sparse(N,N); sparse(N,N) speye(N)]*z;


% This function calculates the product of inv(W)*z as a function of z
% without saving the matrix inv(W).
% Input: z as a column-vector
% Output: Winv(z) as a column-vector
Winv = @(z,~) [sparse(eye(N)-diag(ones(1,N-1),-1)) sparse(N,N); sparse(N,N) speye(N)]*z;

% function handles of u and the Hamiltonian and its Gradient

% entry-function u: [t0,tf] -> R^(m) with m = 1
u = @(t) sin(t);


% this function calculates the Hamiltonian as a function of z = [q p].
% Input: q and p as column-vectors;
% Output: H(q,p) as a scalar
H = @(q,p)  [q;p]'*(Q*[q;p]);


% this function calculates the Gradient of the Hamiltonian as a
% function of z = [q p]. Input: q and p as column-vectors;
% Output: grad_Ham(q,p) as a column-vector
grad_H = @(q,p) Q*[q;p];


% function handles of f and g for symplectic euler
function func = f_function(q,p,t)
    % This function calculates the top half of z=(J-R)*grad_H(z)+B*u and is
    % used in syplectic_euler to solve the equation with the symplectic
    % euler method.
    % Input: q and p as column-vectors
    % Output: f(q,p,t) as a column-vector
    
    %func = (diag(L)*((diag(R)+diag(R(1:N-1),-1))*q')+[1 zeros(1,N-1)]')';
    z = (J-R_m)*grad_H(q,p) + B*u(t);
    func = z(1:N);
end
f = @f_function;

function func = g_function(q,p,t)
    % This function calculates the bottom half of z=(J-R)*grad_H(z)+B*u and
    % is used in syplectic_euler to solve the equation with the symplectic
    % euler method.
    % Input: q and p as column-vectors
    % Output: g(q,p,t) as a column-vector
    
    %func = zeros(1,N);
    z = (J-R_m)*grad_H(q,p) + B * u(t);
    func = z(N+1:end);
end
g = @g_function;

end

