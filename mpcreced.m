function [y1,u1]=mpcreced(Ap,Bp,Cp,Dp,Nc,Np,rw,r,u,y,xm,N_sim)
    %mpcreced() implements the Receding Horizon Control algorithm
    %
    %Inputs:
    %
    %   Ap : 
    %   Bp : 
    %   Cp : 
    %   Dp : 
    %   Nc : 
    %   Np : 
    %   rw : 
    %   r :  
    %   u : 
    %   y : 
    %   xm : 
    %   N_sim : 
    %
    %Outputs:
    %
    %   y1 : 
    %   u1 : 
    %
    [Phi,F,Phi_Phi,Phi_F,Phi_R,Ae,Be,Ce]=mpcgains(Ap,Bp,Cp,Dp,Nc,Np);
    [n,n_in]=size(Be);
    Xf=zeros(n,1);
    for kk=1:N_sim
            dU=inv(Phi_Phi+rw*eye(Nc,Nc))*(Phi_R*r(kk)-Phi_F*Xf);
            du=dU(1,1);
            u=u+du;
            u1(kk)=u;
            y1(kk)=y;

            xm_old=xm;
            xm=Ap*xm+Bp*u;
            y=Cp*xm;
            Xf=[xm-xm_old;y];
    end