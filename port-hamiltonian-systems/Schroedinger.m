function [R, J, Q, W, Winv, B, H, grad_H, t0, tf, x_0, u, par, f, g]...
                                                 = Schroedinger(N)
% This function creates the constants for the Schroedinger-equation. R, J,
% Q, B are matrices and t0, tf are scalar and x_0 is a vector. H, grad_H
% and u are function-handles depending on z=[q p]. This function also
% creates f and g, which are needed to solve this equation with the
% semi-implicit euler or with the discrete-gradient-method.
% source: [scan from Prof. Dr. Tatjana Stykel]

n = 2*N;                                                        % constants 

par.l = 0.11;                                           
par.epsilon = 1;
par.c = 1;

ll = (2*pi)/par.l;
par.hx = ll/N;                                                  % step size

% no R given
R = sparse(n,n);

% no Port given
B = sparse(n,1);

% no Q needed
Q = sparse(n,n);

t0 = 0;                                                    % time intervall
tf = 20;

% L is used to create Q
e_ = ones(N,1);
L = spdiags([e_ -2*e_ e_], -1:1, N, N);
L(1,end) = 1;
L(end,1) = 1;
par.L = L/(par.hx*par.hx);

% skew-symmetric Matrix J
J = [sparse(N,N) -speye(N)/par.hx; speye(N)/par.hx sparse(N,N)];

x_0 = zeros(n,1);                                       % initial condition
for o=1:N
    xi = -ll/2 + o*par.hx;
    x_0(o) = (sqrt(2)*sin(par.c*xi)/2)/(cosh(xi));
    x_0(N+o) = (sqrt(2)*cos(par.c*xi)/2)/(cosh(xi));
end


% This function calculates the product of W*z as a function of z
% without saving the matrix W.
% Input: z as a column-vector
% Output: W(z) as a column-vector
W = @(z,hx) [-hx*z(1:length(z)/2); z(length(z)/2+1:end)];


% This function calculates the product of inv(W)*z as a function of z
% without saving the matrix inv(W).
% Input: z as a column-vector
% Output: Winv(z) as a column-vector
Winv = @(z,hx) [-z(1:length(z)/2)/hx; z(length(z)/2+1:end)];

% function handles of u and the Hamiltonian and its Gradient

% entry-function u: [t0,tf] -> R^(m) 
u = @(t) [];


% this function calculates the Hamiltonian as a function of z = [q p].
% Input: q and p as column-vectors;
% Output: H(q,p) as a scalar
H = @(q,p) -par.hx*(q'*(par.L*q))/2 - par.hx*(p'*(par.L*p))/2 - ...
    par.hx*par.epsilon*sum((p.^2+q.^2).^2)/4;

% this function calculates the Gradient of the Hamiltonian as a
% function of z = [q p]. Input: q and p as column-vectors;
% Output: grad_Ham(q,p) as a column-vector
grad_H = @(q,p) (-par.hx)*[par.L*q; par.L*p] - ...
    (par.hx*par.epsilon)*[(p.^2+q.^2).*q; (p.^2+q.^2).*p];

% This function calculates the top half of z=(J-R)*grad_H(z)+B*u and is
% used in symplectic_euler to solve the equation with the symplectic 
% euler method.
% Input: q and p as column-vectors
% Output: f(q,p,t) as a column-vector
f = @(q,p,t) par.L*p+par.epsilon*((p.^2+q.^2).*p);


% This function calculates the bottom half of z=(J-R)*grad_H(z)+B*u and
% is used in symplectic_euler to solve the equation with the symplectic 
% euler method.
% Input: q and p as column-vectors
% Output: g(q,p,t) as a column-vector
g = @(q,p,t) -par.L*q-par.epsilon*((p.^2+q.^2).*q);

end

