function [R, J, Q, W, Winv, B, H, grad_H, t0, tf, x_0, u, par, f, g]...
                                                            = TLS(N, sel_u)
% % This function creates the constants for the Toda-Lattice-System. R, J,
% Q, B are matrices and t0, tf are scalar and x_0 is a vector. H, grad_H
% and u are function-handles depending on z=[q p]. This function also
% creates f and g, which are needed to solve this equation with the
% semi-implicit euler method.
% source: [Chaturantabut/Beattie/Gugercin`16, section 3.4.2]

n = 2*N;                                          % N-particle Toda lattice

% skew-symmetric Matrix J equals canonical form of J (no W needed)
J = [sparse(N,N) speye(N);-speye(N) sparse(N,N)];

gam = repmat(0.1,N,1);                                  % init gamma values

% symmetric, positive semidefinite Matrix R
R = [sparse(N,N) sparse(N,N); sparse(N,N) sparse(diag(gam))];     

% Port-Matrix with dimensions (n x m) and m = 1
B = [sparse(N,1); 1; sparse(N-1,1)];         

x_0 = sparse(n,1);                                      % initial condition

t0 = 0;                                                    % time intervall
tf = 100;                        

par.hx = 0.01;                                                  % step size

% symmetric, positive semidefinite Matrix Q
Q = [sparse(N,N) sparse(N,N); sparse(N,N) speye(N)];

% This function calculates the product of W*z as a function of z
% without saving the matrix W.
% Input: z as a column-vector
% Output: W(z) as a column-vector
% comment: this time no Transformation-matrix is needed because J is in
% canonical form [0 I; -I 0]
W = @(z,hx) [];


% This function calculates the product of inv(W)*z as a function of z
% without saving the matrix inv(W).
% Input: z as a column-vector
% Output: Winv(z) as a column-vector
% comment: this time no Transformation-matrix is needed because J is in
% canonical form [0 I; -I 0]
Winv  = @(z,hx) [];


% entry-function u: [t0,tf] -> R^(m) with m = 1
function func = U(t)                      
    if sel_u == 1
        func = 0.1;
    elseif sel_u == 2
        func = 0.1*sin(t);
    end
end
u = @U;

function g = Hamiltonian(q,p)
    % this function calculates the Hamiltonian as a function of z = [q p].
    % Input: q and p as column-vectors;
    % Output: H(q,p) as a scalar
    
    sum_ = 0;
    for k=1:N-1                    
        sum_ = sum_ + exp(q(k)-q(k+1));
    end
    g = 1/2 * sum(p.^2) + sum_ + exp(q(N)) - q(1) - N;
end
H = @Hamiltonian;


function func = grad_Ham(q,p)
    % this function calculates the Gradient of the Hamiltonian as a
    % function of z = [q p]. Input: q and p as column-vectors;
    % Output: grad_Ham(q,p) as a column-vector
    
    func = zeros(2*N,1);                       
    func(1) = exp(q(1)-q(2)) - 1;
    func(N) = exp(q(N)) - exp(q(N-1) - q(N));

    for j=2:N-1
        func(j) = exp(q(j) - q(j+1)) - exp(q(j-1) -q(j));
    end

    func(N+1:2*N) = p;
end
grad_H = @grad_Ham;

% function handles of f and g for symplectic euler
function func = f_function(q,p,t)
    % This function calculates the top half of z=(J-R)*grad_H(z)+B*u and is
    % used in syplectic_euler to solve the equation with the symplectic
    % euler method.
    % Input: q and p as column-vectors
    % Output: f(q,p,t) as a column-vector
    
    %func = flip(p);
    z = (J-R)*grad_H(q,p) + B*u(t);
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
    %for k1=1:N
    %    func(k1) = -(gam(k1)*q(k1) - q(N-k1+1));
    %end
    %func = func + B*u(t);
    z = (J-R)*grad_H(q,p) + B *u(t);
    func = z(N+1:end);
end
g = @g_function;

end

