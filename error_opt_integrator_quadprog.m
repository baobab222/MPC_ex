clear; clc;

a1 = -1; a2 = -2;
A = [0 1; a1 a2];
B = [0; 1];
C = [1 0];
D = 0;
Ts = 0.1;
Nsim =100;

[Adn,Bdn,Cdn,Ddn] = c2dm(A,B,C,D,Ts);   % вместо c2dm
Ad= [Adn zeros(2,1);Cdn*Adn 1];
Bd= [Bdn;Cdn*Bdn];
Cd= [zeros(1,2) 1];
Dd=Ddn;
Np = 10;
Nc = 3;                 % для абсолютного входа лучше Nc=Np
refy = ones(Np,1);
refu = ones(Nc,1)*0;       % u_ss = 1 для y=1
ref = [refy; refu];
% for i=1:Nsim+Np
%    REFY(i,1)= 0.01*sin(0.1*i); 
% end

[Phi,F] = mpcgain1(Ad,Bd,Cd,Nc,Np);

Q1 = 10*eye(Np);    
Q2 = 1*eye(Nc);
H = blkdiag(Q1,Q2);
f = zeros(Np+Nc,1);

x0 = [0;0];
Aeq = [-eye(Np), Phi];

x = [x0;Cdn*x0];
xd=x0;

y = zeros(1,Nsim);
u = zeros(1,Nsim);
dist= 0;
u=0;
for i = 1:Nsim
%     refy= REFY(i:i+Np-1);
%     ref = [refy; refu];
    beq = Aeq*ref + F*x;

    [U, fval, exitflag] = quadprog(H, f, [], [], Aeq, beq);
    if exitflag ~= 1
        warning('quadprog exitflag = %d', exitflag);
    end

    z_opt = ref - U;
    Uc_opt = z_opt(Np+1:end);
    if i>50
        dist=1;
    end
    du(i) = Uc_opt(1);
    u(i+1)=u(i)+du(i);
    xd_prev =xd;
    xd= Adn*xd+Bdn*(u(i+1)+dist);
%     x = Ad*x + Bd*(u(i+1)+dist);
    y(i) = Cdn*xd;
    x=[xd-xd_prev;y(i)];
%     x = x+[randn(1,1)*0.1;0];
end
% plot((0:Nsim-1)*Ts,REFY(1:Nsim))
% hold on
plot((0:Nsim-1)*Ts, y, 'LineWidth', 1.2);
grid on;
xlabel('t, s');
ylabel('y');

function [Phi,F] = mpcgain1(Ap,Bp,Cp,Nc,Np)
    [~,n1] = size(Cp);

    h = zeros(Np,n1);
    F = zeros(Np,n1);

    h(1,:) = Cp;
    F(1,:) = Cp*Ap;

    for kk = 2:Np
        h(kk,:) = h(kk-1,:)*Ap;
        F(kk,:) = h(kk,:)*Ap;   % C*A^kk
    end

    Phi = zeros(Np,Nc);
    for i = 1:Nc
        for kk = i:Np
            Phi(kk,i) = h(kk-i+1,:)*Bp;  % C*A^(kk-i)*B
        end
    end
end

