function [T,T1,T2,T4]=Dmat(N,y)

T=zeros(N+1,N+1);
T1=zeros(N+1,N+1);
T2=zeros(N+1,N+1);
T3=zeros(N+1,N+1);
T4=zeros(N+1,N+1);

for j=1:N+1
    T(1,j)=1;
    T(2,j)=y(j);
    T1(2,j)=T(1,j);
    T(3,j)=2*y(j)*T(2,j)-T(1,j);
    T1(3,j)=4*T(2,j);
    T2(3,j)=4*T1(2,j);
end

for i=4:N+1
    for j=1:N+1
        T(i,j)=2*y(j)*T(i-1,j)-T(i-2,j);
        T1(i,j)=2*(i-1)*T(i-1,j)+((i-1)/(i-3))*T1(i-2,j);
        T2(i,j)=2*(i-1)*T1(i-1,j)+((i-1)/(i-3))*T2(i-2,j);
        T3(i,j)=2*(i-1)*T2(i-1,j)+((i-1)/(i-3))*T3(i-2,j);
        T4(i,j)=2*(i-1)*T3(i-1,j)+((i-1)/(i-3))*T4(i-2,j);
    end
end
end