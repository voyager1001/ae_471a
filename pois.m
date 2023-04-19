function [A,B]=pois(nosmod,alp,beta,R, T, T1, T2, T4)
 
ak2=sqrt(alp^2+beta^2);
vec=(0:1:nosmod)';

u=(ones(length(vec),1)-cos(pi*vec/nosmod).^2); 
duu =-2*(ones(length(vec),1));

A=zeros(nosmod+1,nosmod+1);
B=zeros(nosmod+1,nosmod+1);
p=complex(0,-5);
    for i=3:nosmod-1
        for j=1:nosmod+1
                A(1,j)=p*T(j,1);
                A(2,j)=p*T1(j,1);
                A(nosmod,j)=p*T1(j,nosmod+1);
                A(nosmod+1,j)=p*T(j,nosmod+1);
                A(i,j)=L1(i,alp,ak2,R,u,duu)*T(j,i)+L2(i,alp,ak2,R,u)*T2(j,i)-T4(j,i)/(complex(0,1)*alp*R);
        end
    end
    for i=3:nosmod-1
        for j=1:nosmod+1
                B(1,j)=T(j,1);
                B(2,j)=T1(j,1);
                B(nosmod,j)=T1(j,nosmod+1);
                B(nosmod+1,j)=T(j,nosmod+1);
                B(i,j)=T2(j,i)-ak2*ak2*T(j,i);
        end
    end
end

function x=L1(y,a,a1,Re,U,U2)
    x=-U(y)*a1*a1-U2(y)-(a1^4)/(complex(0,1)*a*Re);
end

function x=L2(y,a,a1,Re,U)
    x=U(y)+(2*a1*a1)/(complex(0,1)*a*Re);
end



