function [u1,y1,deltau1,k]= simuob(xm,u,y,sp,Ap,Bp,Cp,A,B,C,N_sim,Omega,Psi,K_ob,Lzerot)
%closed-loop simulation without constraints
% [m1,n1]=size(Cp);
% [n1,n_in]=size(Bp);
Xf=[xm;(y-sp(:,1))];
[ny,n]=size(C);
[n,n_in]=size(B);
X_hat=rand(n,1);
for kk=1:N_sim;
    Xsp=[zeros(n-ny,1);sp(:,kk)];
    %eta=-(Omega\Psi)*Xf;
    eta=-(Omega\Psi)*(X_hat-Xsp);
    deltau=Lzerot*eta;
    u=u+deltau;
    deltau1(:,kk)=deltau;
    u1(1:n_in,kk)=u;
    y1(1:ny,kk)=y;
    %X_hat=A*X_hat+K_ob*(y-C*X_hat)+B*deltau;
    X_hat=A*X_hat+K_ob*(C*Xf-C*X_hat)+B*deltau;
    %%%%
    %plant simulation
    %%%%%%
    xm_old=xm;
    xm=Ap*xm+Bp*u; % calculate xm(k+1)
    y=Cp*xm; %calculate y(k+1)
    %updating feedback state variable Xf
    Xf=[xm-xm_old;(y-sp(:,kk+1))];
end
k=0:(N_sim-1);
end