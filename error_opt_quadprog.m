clear all
close all
a1 = -1; a2 = -2;
Am = [0 1; a1 a2];
Bm = [0; 1];
Cm = [1 0];
Dm = 0;
Ts = 0.1;
Nsim =100;

[Adn,Bdn,Cdn,Ddn] = c2dm(Am,Bm,Cm,Dm,Ts);   
Np=10;
Nc=3;
refy=ones(Np,1)*1;
refu=ones(Nc,1)*1;
ref=[refy;refu];
[Phi, F]= mpcgain1(Adn,Bdn,Cdn,Nc,Np);
Q1=diag(ones(1,Np))*1;
Q2=diag(ones(1,Nc))*0.01;
H=blkdiag(Q1,Q2);
x0=[0;0];
Aeq= [-eye(Np) Phi];

x=x0;
for i=1:Nsim
    beq= Aeq*ref+F*x0;
    [U, fval, exitflag, output] = quadprog(H,[],[],[],Aeq,beq);
    z_opt  = ref - U;
    Y_opt  = z_opt(1:Np);
    Uc_opt = z_opt(Np+1:end);
    if i>30 
        Uc_opt(1)=Uc_opt(1)+1;
    end
    x=Adn*x+Bdn*Uc_opt(1);
    y(i)=Cdn*x;
    x0=x;
end
plot(y)
    
function [Phi,F]=mpcgain1(Ap,Bp,Cp,Nc,Np)
    [m1,n1]=size(Cp);
    [n1,n_in]=size(Bp);
    n=n1+m1;
    h(1,:)=Cp;
    F(1,:)=Cp*Ap;
    for kk=2:Np
        h(kk,:)=h(kk-1,:)*Ap;
        F(kk,:)= F(kk-1,:)*Ap;
    end
    v=h*Bp;
    Phi=zeros(Np,Nc); %declare the dimension of Phi
    Phi(:,1)=v;
    for i=2:Nc
        Phi(:,i)=[zeros(i-1,1);v(1:Np-i+1,1)]; %Toeplitz matrix
    end
end

