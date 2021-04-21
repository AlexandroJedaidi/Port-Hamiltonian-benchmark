function [R, J, Q, W, Winv, B, H, grad_H, t0, tf, x_0, u, par, f, g]...
                                                                   = DWE(N)
% This function creates the constants for the 1-D-wave-equation. R, J,
% Q, B are matrices and t0, tf are scalar and x_0 is a vector. H, grad_H
% and u are function-handles depending on z=[q p]. This function also
% creates f and g, which are needed to solve this equation with the
% semi-implicit euler method.
% source: [Peng/Mohseni'15, sections 6.1, 6.2] &
% [Afkham/Hesthaven'17, section 4.1]
% comment: the wave speed c is given with c_q and is already squared

% constants


par.l = 1;                                                  % domain length
par.c_q = 0.1;                                                 % wave speed  

par.hx = par.l/N;                                             % x-increment
par.ht = 0.0005;                                              % t-increment

n = 2*N;

% no Port given
B = sparse(n,1);

t0 = 0;
tf = 10;                                              % time-intervall in s

x = zeros(N,1);                                                      % Grid
for i=1:N
    x(i) = (i*par.l)/N;
end

% used to create bottom half of R
r = zeros(N,1);
for m=1:N
    r(m) = 0.1 + 0.9*(m/N);
end

q_a = zeros(N,1);                                       % initial condition
for k=1:N                                       
    q_a(k) = cubic_spline(10*abs(x(k)-0.5));    
end                                             
p_a = zeros(N,1);                               
                                                
x_0 = [q_a; p_a];   

% symmetric, positive semidefinite Matrix R
r_ = [zeros(N,1); r];
R = sparse(diag(r_));

% skew-symmetric Matrix J
J = [sparse(N,N) speye(N); -speye(N) sparse(N,N)];

% L is used to create the top part of Q
e_ = ones(N,1);
L = spdiags([e_ -2*e_ e_], -1:1, N, N);
par.L = L/(par.hx*par.hx);

% symmetric, positive semidefinite Matrix Q
Q =[-(par.c_q*par.L) sparse(N,N); sparse(N,N) speye(N)];

% This function calculates the product of W*z as a function of z
% without saving the matrix W.
% Input: z as a column-vector
% Output: W(z) as a column-vector
W = @(z,h_x) [hx*z(1:length(z)/2); z(length(z)/2+1:end)];


% This function calculates the product of inv(W)*z as a function of z
% without saving the matrix inv(W).
% Input: z as a column-vector
% Output: Winv(z) as a column-vector
Winv = @(z,hx) [z(1:length(z)/2)/hx;z(length(z)/2+1:end)];

% function handles of u and the Hamiltonian and its Gradient

% entry-function u: [t0,tf] -> R^(m) 
u = @(t) 0;


% this function calculates the Gradient of the Hamiltonian as a
% function of z = [q p]. Input: q and p as column-vectors;
% Output: grad_Ham(q,p) as a column-vector                                           
grad_H = @(q,p) Q * [q;p];


% this function calculates the Hamiltonian as a function of z = [q p].
% Input: q and p as column-vectors;
% Output: H(q,p) as a scalar
H = @(q,p) ([q;p]'*(Q*[q;p]))./2 ;

function func = f_func(q,p,t)
    % This function calculates the top half of z=(J-R)*grad_H(z)+B*u and is
    % used in syplectic_euler to solve the equation with the symplectic
    % euler method.
    % Input: q and p as column-vectors
    % Output: f(q,p,t) as a column-vector
    %func = p;
    
    z = (J-R)*grad_H(q,p) + B * u(t);
    func = z(1:length(z)/2);
end
f = @f_func;

function func = g_func(q,p,t)
    % This function calculates the bottom half of z=(J-R)*grad_H(z)+B*u and
    % is used in syplectic_euler to solve the equation with the symplectic
    % euler method.
    % Input: q and p as column-vectors
    % Output: g(q,p,t) as a column-vector

    z = (J-R)*grad_H(q,p) + B* u(t);
    func = z(length(z)/2+1:end);

    %func =  par.c_q .* (par.L * q) - sparse(diag((r)))*p;
end
g = @g_func;


end

