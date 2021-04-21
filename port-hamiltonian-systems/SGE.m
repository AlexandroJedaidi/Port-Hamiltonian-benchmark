function [R, J, Q, W, Winv, B, H, grad_H, t0, tf, x_0, u, par, f, g]...
                                                                = SGE(N)
% % This function creates the constants for the Sine-Gordon-Equation. R, J,
% B are matrices and t0, tf are scalar and x_0 is a vector. H, grad_H
% and u are function-handles depending on z=[q p]. This function also
% creates f and g, which are needed to solve this equation with the
% semi-implicit euler method.
% source: [e-mail from Prof. Dr. Tatjana Stykel] &
% [Peng/Mohseni'15 section 6.3] 

n = 2*N;                                                        % constants

t0 = 0;                                                    % time intervall
tf = 10;
par.l = 50;                                                 % domain-length
par.hx = par.l/N;                                             % x-step-size
par.ht = 0.0001;                                              % t-step-size

% no R given
R = sparse(n,n); 

% no Port given
B = sparse(n,1);

% no Q given
Q = sparse(n,n);

x_0 = sparse(n,1);                                      % initial condition

% skew-symmetric Matrix J
J = [sparse(N,N) speye(N)./par.hx; -speye(N)./par.hx sparse(N,N)];

% L is used to calculate the hamiltonian and its gradient
e_ = ones(N,1);
L = spdiags([e_ -2*e_ e_], -1:1, N, N);
par.L = L/(par.hx*par.hx);


% This function calculates the product of W*z as a function of z
% without saving the matrix W.
% Input: z as a column-vector
% Output: W(z) as a column-vector
W = @(z,hx) [hx*z(1:length(z)/2);z(length(z)/2+1:end)];


% This function calculates the product of inv(W)*z as a function of z
% without saving the matrix inv(W).
% Input: z as a column-vector
% Output: Winv(z) as a column-vector
Winv = @(z,hx) [z(1:length(z)/2)/hx; z(length(z)/2 +1:end)];

% function handles of u and the Hamiltonian and its Gradient

% entry-function u: [t0,tf] -> R^(m) 
u = @(t) 0;


% this function calculates the Hamiltonian as a function of z = [q p].
% Input: q and p as column-vectors;
% Output: H(q,p) as a scalar
H = @(q,p) par.hx/2*(p'*p - q'*par.L*q) - ...
                                      (2*pi*q(N))/par.hx + (2*pi^2)/par.hx;


% this function calculates the Gradient of the Hamiltonian as a
% function of z = [q p]. Input: q and p as column-vectors;
% Output: grad_Ham(q,p) as a row-vector
grad_H = @(q,p) [-par.hx*par.L*q;par.hx*p] + par.hx*[sin(q(1:N-1)); ...
        (sin(q(N)) - (2*pi)/(par.hx*par.hx)); sparse(N,1)];

% function handles of f and g for symplectic euler
function func = f_function(q,p,~)
    % This function calculates the top half of z=(J-R)*grad_H(z)+B*u and is
    % used in syplectic_euler to solve the equation with the symplectic
    % euler method.
    % Input: q and p as column-vectors
    % Output: f(q,p,t) as a column-vector
    
    z = (J-R)*grad_H(q,p) + B *u();
    func = z(1:N);
end
f = @f_function;

function func = g_function(q,p,~)
    % This function calculates the bottom half of z=(J-R)*grad_H(z)+B*u and
    % is used in syplectic_euler to solve the equation with the symplectic
    % euler method.
    % Input: q and p as column-vectors
    % Output: g(q,p,t) as a column-vector
    
    z = (J-R)*grad_H(q,p) + B*u();
    func = z(N+1:end);
end
g = @g_function;

end

