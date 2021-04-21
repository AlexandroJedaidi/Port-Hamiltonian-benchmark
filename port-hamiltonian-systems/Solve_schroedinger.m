function [] = Solve_schroedinger(J, H, grad_H, t0, tf, x_0, par, N)
%This function solves the given system with the matrices J,H,B, the 
%function H and starting point x_0.

nt = 1000;    

ht = (tf-t0)/(nt);
hx = par.hx; 
epsilon = par.epsilon;

% initial condition
q0 = x_0(1:N);  
p0 = x_0(N+1:end);

t = linspace(t0,tf,nt+1);

q = zeros(N,nt+1); q(:,1) = q0;
p = zeros(N,nt+1); p(:,1) = p0;

hessH = @(q,p) -hx*[par.L sparse(N,N); sparse(N,N) par.L] - ...
     (par.hx*par.epsilon)*[spdiags(3*q.^2+p.^2,0,N,N) spdiags(2*q.*p,0,N,N);
     spdiags(2*q.*p,0,N,N) spdiags(q.^2+3*p.^2,0,N,N)];

for k = 2:nt+1
    
    k = k;
    
    A = speye(2*N)-(ht/2)*J*hessH(q(:,k-1),p(:,k-1));
    x = (A)\(J*grad_H(q(:,k-1),p(:,k-1)));
    q(:,k) = q(:,k-1) + ht*x(1:N);
    p(:,k) = p(:,k-1) + ht*x(N+1:end);
  
end

for i=1:length(t)
  Ht(i)=H(q(:,i),p(:,i));
end

figure(1); clf;
plot(t,Ht);
ylabel('H(t)');
xlabel('t (sec)');

errH = abs(Ht-Ht(1)*ones(1,nt+1));
figure(2); clf;
semilogy(t(2:end),errH(2:end));
ylabel('|H(t)-H(0)|');
xlabel('t (sec)');


end



