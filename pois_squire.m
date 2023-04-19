function [A,B]=pois_squire(nosmod,alp,beta,R, T, T1, T2, T4) 

ak2=sqrt(alp^2+beta^2);

vec=(0:1:nosmod)'; 
u=(ones(length(vec),1)-cos(pi*vec/nosmod).^2); 
duu=-2*(ones(length(vec),1));

A=zeros(nosmod+1,nosmod+1);
B=zeros(nosmod+1,nosmod+1);
p=complex(0,-5);
    for i=2:nosmod
        for j=1:nosmod+1
                A(1,j)=p*T(j,1);
                A(nosmod+1,j)=p*T(j,nosmod+1);
                A(i,j)=(complex(0,1)*ak2*u(i)+(ak2^2)/R)*T(j,i)-T2(j,i)/R;
        end
    end
    for i=1:nosmod+1
        for j=1:nosmod+1
            B(i,j)=T(j,i)*complex(0,1)*alp;
        end
    end 
end


