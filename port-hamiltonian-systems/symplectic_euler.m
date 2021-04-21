function [] = symplectic_euler(H, t0, tf, x_0, N, f, g, h)
%This function solves the given system with the matrices J,H,B, the 
%function H and starting point x_0,

% dq = f(t,p)
% dp = g(t,q)

nt = 1000;

% initial condition
q0 = x_0(1:length(x_0)/2);  
p0 = x_0(length(x_0)/2+1:end);

% step size
if h == 0
    ht = (tf-t0)/(nt*100);
else
    ht = h;
end

% x-axis as time
t = linspace(t0,tf,nt+1);

% memory-allocation of q and p
q = zeros(N,nt+1); q(:,1) = q0;
p = zeros(N,nt+1); p(:,1) = p0;

for k = 2:nt+1
    
    k = k;
    
    % symplectic euler

    p(:,k) = p(:,k-1) + ht*g(q(:,k-1),p(:,k-1),t(:,k));
    q(:,k) = q(:,k-1) + ht*f(q(:,k-1),p(:,k),t(:,k));
    
end

% y-axis as hamiltonian
for i=1:length(t)
  Ht(i)=H(q(:,i),p(:,i));
end

% plot H(t)
figure(1); clf;
plot(t,Ht);
ylabel('H(t)');
xlabel('t (sec)');

% plot |H(t) - H(0)|
errH = abs(Ht-Ht(1)*ones(1,nt+1));
figure(2); clf;
semilogy(t(2:end),errH(2:end));
ylabel('|H(t)-H(0)|');
xlabel('t (sec)');

end