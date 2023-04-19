function [u,duu,ymax,x]=blasius_u(x)
N1=0.3;
N2=0.2;
n=300;
u0=[0,0,N1];
u1=[0,0,N2];
nspan=[0 n];
sol=ode45(@blf,nspan,u0);
U=deval(sol,n);
phi0=U(2)-1;
sol=ode45(@blf,nspan,u1);
U=deval(sol,n);
phi1=U(2)-1;
N3=N2-((phi1*(N2-N1))/(phi1-phi0));
u3=[0,0,N3];
while N3<0.3320
   sol=ode45(@blf,nspan,u3);
    U=deval(sol,n);
    phi0=phi1;
    phi1=U(2)-1;
   N1=N2;
   N2=N3;
   N3=N2-((phi1*(N2-N1))/(phi1-phi0));
   u3=[0,0,N3];
 
end

ymax = 15;
y_new = 0.5 * (x + 1) * ymax;
u_temp=deval(sol,1.721*y_new);
display(u_temp)
u = u_temp(2,:);
figure(4);
f = u_temp(1,:);
fd_dash = u_temp(3,:);
duu = -0.5 * f .* fd_dash*(1.721)^2; 
end