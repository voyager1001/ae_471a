zi=sqrt(-1); 

% input data 
iflow=input('Poiseuille Orr-Sommerfeld (1) or Couette Orr-Sommerfeld (2) or Blasius Orr-Sommerfeld (3) or Poiseuille Squire (4) or Couette Squire (5) or Blasius Squire (6) '); 
nosmod=input('Enter N the number of OS modes: '); 
R=input('Enter the Reynolds number: ');
alp=input('Enter alpha: ');

beta=0;
ak2=sqrt(alp^2+beta^2);

%Collocation points
y=zeros(nosmod+1,1);
for i=0:nosmod
    y(i+1)=cos(i*pi/nosmod);
end

% generate Chebyshev differentiation matrices 
[D0,D1,D2,D4]=Dmat(nosmod,y); 

% set up Orr-Sommerfeld matrices A and B
if iflow==1
    [A,B]=pois(nosmod,alp,beta,R,D0,D1,D2,D4); 
elseif iflow==2
    [A,B]=couet(nosmod,alp,beta,R,D0,D1,D2,D4); 
elseif iflow==3
    [A,B,ymax]=blasius(nosmod,alp,beta,R,D0,D1,D2,D4,y);
elseif iflow==4
    [A,B]=pois_squire(nosmod,alp,beta,R,D0,D1,D2,D4); 
elseif iflow==5
    [A,B]=couet_squire(nosmod,alp,beta,R,D0,D1,D2,D4); 
elseif iflow==6
    [A,B,ymax]=blasius_squire(nosmod,alp,beta,R,D0,D1,D2,D4,y);
end

% compute eigenvalues and eigenvectors
% [d,v]=eig(inv(B)*A); 
    Q=B\A;
    %[a,b]=eig(Q);
    [d,v]=eigs(Q,125,'largestimag','Tolerance',1e-12,'Display',true,'MaxIterations',500);

%plotting eigenvalues
eigenvals = diag(v);
imaginary_part = imag(eigenvals);
real_part = real(eigenvals);
figure(1)
plot(real_part,imaginary_part,".k")
xlabel('Real Part') 
ylabel('Imaginary Parts')
if(iflow==2||iflow==5)
    xlim([-1 1])
    ylim([-1 0.1])
else
xlim([0 1])
ylim([-1 0.1])
end

%Plotting eigenvectors
idx = find((real(eigenvals) > 0.2 & real(eigenvals) < 0.5) & (imag(eigenvals) > -0.2 & imag(eigenvals) < 0.1), 1);
desired_eigenvec = d(:,idx);
desired_eigenval = eigenvals(idx);

% Computing velocity or vorticity
product = D0'*desired_eigenvec;

if(iflow==3||iflow==6)
    product1 = D1'*desired_eigenvec/(1i*ak2);
    y_new = 0.5 * (y + 1) * ymax;
    figure(2)
    plot(abs(product),y_new,'k',LineWidth=1)
    hold on
    plot(-real(product),y_new,'c')
    hold on
    plot(imag(product),y_new,'g')
    hold on
    xlabel('v') 
    ylabel('y') 
    figure(3)
    plot(abs(product1),y_new,'k',LineWidth=1)
    hold on
    plot(-real(product1),y_new,'c')
    hold on
    plot(imag(product1),y_new,'g')
    xlabel('u') 
    ylabel('y') 
else 
    figure(2)
    plot(y, abs(product),'k',LineWidth=1)
    hold on
    plot(y, real(product),'c')
    hold on
    plot(y, imag(product),'g')
end
