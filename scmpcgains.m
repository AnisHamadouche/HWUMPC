function [Phi,F,Phi_Phi,Phi_F,Phi_R,Ae,Be,Ce]=scmpcgains(Ap,Bp,Cp,Nc,Np)
%scmpcgains() calculates Phi'Phi, Phi'F and Phi'Rs control gains of future
%control signals for 1 optimization window of lenght Np.
%See: (1.17) in https://drive.google.com/file/d/1Pfn9IzV24xJ07utEBiJOAV5xY8Kj_0yv/view?usp=sharing
%
%Inputs:
%   Ap : n1 x n1 matrix;
%   Bp : n1 x n_in matrix;
%   Cp : m1 x n1 matrix;
%   Nc : control horizon;
%   Np : prediction horizon;
%
%Outputs:
%   Ae : n1+m1 x n1+m1 matrix
%   Be : n1+m1 x n_in matrix
%   Ce : m1 x n1+m1 matrix
%
%Usage:
%   [Phi_Phi,Phi_F,Phi_R,Ae,Be,Ce]=mpcgain(Ap,Bp,Cp,Nc,Np);

    %[Ad,Bd,Cd,Dd]=c2dm(Ac,Bc,Cc,Dc,Delta_t);
    [m1,n1]=size(Cp);
    [n1,n_in]=size(Bp);
    Ae=eye(n1+m1,n1+m1);
    Ae(1:n1,1:n1)=Ap;
    Ae(n1+1:n1+m1,1:n1)=Cp*Ap;
    Be=zeros(n1+m1,n_in);
    Be(1:n1,:)=Bp;
    Be(n1+1:n1+m1,:)=Cp*Bp;
    Ce=zeros(m1,n1+m1);
    Ce(:,n1+1:n1+m1)=eye(m1,m1);
    n=n1+m1;
    h(1,:)=Ce;
    F(1,:)=Ce*Ae;
    for kk=2:Np
        h(kk,:)=h(kk-1,:)*Ae;
        F(kk,:)= F(kk-1,:)*Ae;
    end
    v=h*Be;
    Phi=zeros(Np,Nc); %declare the dimension of Phi
    Phi(:,1)=v; % first column of Phi
    for i=2:Nc
        Phi(:,i)=[zeros(i-1,1);v(1:Np-i+1,1)]; %Toeplitz matrix
    end
    BarRs=ones(Np,1);
    Phi_Phi= Phi'*Phi;
    Phi_F= Phi'*F;
    Phi_R=Phi'*BarRs;
end