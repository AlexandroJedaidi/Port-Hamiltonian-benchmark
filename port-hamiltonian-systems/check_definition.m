function [] = check_definition(R, J, Q, B, N)
% checking properties of R, J, Q and B

n = 2*N;

J_2n = [sparse(N,N) speye(N);-speye(N) sparse(N,N)];

disp("checking if J has canonical form [0 I; -I 0]")
i = isequal(J,J_2n);

disp(" ")

if i == 1
    disp("no transformationmatrix W needed")
    j = 1;
    J_sym = 1;
else
    j = 0;
    disp("checking if J is skew symmetric")
    J_sym = issymmetric(J, 'skew');
end

disp(" ")

disp("checking if R is symmectric and positive semi-definite")
i =  isequal(R,sparse(n,n));
if i == 1
    r_z = 1;
    R_sym = 1;
    issemi_R = 1;
else
    r_z = 0;
    R_sym = issymmetric(R);
    d1 = eig(R);
    issemi_R = all(d1 >= 0);
end

disp(" ")

disp("checking B")
i = isequal(B,sparse(n,1));
if i == 1
    disp("no port given")
    b_z = 1;
else
    b_z = 0;
end
disp(" ")

if b_z == 1 && r_z == 1 && j == 1
    disp("canonical hamiltonian system given")
    disp("the hamiltonian should be constant (check graph while solving!)")
end
if b_z == 1 && r_z == 1
    disp("dissipative hamiltonian system given")
    disp("the hamiltonian should be constant (check graph while solving!)")
end

disp(" ")

disp("checking Q")
i = isequal(Q,sparse(n,n));
if i == 1
    disp("no Q given")
    Q_sym = 1;
    issemi_Q = 1;
else
    Q_sym = issymmetric(R);
    d2 = eig(R);
    issemi_Q = all(d2 >= 0);
end

disp(" ")

error = [J_sym R_sym issemi_R Q_sym issemi_Q];
test = [1 1 1 1 1];
t = isequal(error,test);
if t == 1
    disp("no errors")
end

disp(" ")

% display dimensions of given matrices
disp("dimensions:")
matrix = {'Q'; 'R'; 'J'; 'B'};
dimension = [size(Q); size(R); size(J); size(B)];
T = table(matrix,dimension)
 
end

