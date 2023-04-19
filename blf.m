function dfdn=blf(~,f)
    dfdn=zeros(3,1);
    dfdn(1,1)=f(2);
    dfdn(2,1)=f(3);
    dfdn(3,1)=-0.5*f(1)*f(3);
end
