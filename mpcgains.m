function [Phi,F,Phi_Phi,Phi_F,Phi_R,Ae,Be,Ce]=mpcgains(Ap,Bp,Cp,Dp,Nc,Np)
%mpcgains() calculates Phi'Phi, Phi'F and Phi'Rs control gains of future
%control signals for 1 optimization window of lenght Np.
%See: (1.17) in https://drive.google.com/file/d/1Pfn9IzV24xJ07utEBiJOAV5xY8Kj_0yv/view?usp=sharing
%
%Inputs (from continuous ss):
%   Ap : n1 x n1 matrix ;
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
    [Ae,Be,Ce] = ss_augment(Ap,Bp,Cp,Dp);%augmentation
    n=n1+m1;
    h=zeros(Np*size(Ce,1),size(Ce,2));    
    F=zeros(Np*size(Ce,1),size(Ce,2));
    h(1:m1,:)=Ce;
    F(1:m1,:)=Ce*Ae;
    for kk=2:Np
        h((kk-1)*m1+1:(kk)*m1,:)=h((kk-2)*m1+1:(kk-1)*m1,:)*Ae;
        F((kk-1)*m1+1:(kk)*m1,:)=F((kk-2)*m1+1:(kk-1)*m1,:)*Ae;
    end
    v=zeros(Np*size(Ce*Ae*Be,1),size(Ce*Ae*Be,2));
    v(1:m1,:)=h(1:m1,:)*Be;
    for kk=2:Np
        v((kk-1)*m1:(kk)*m1,:)=h((kk-1)*m1:(kk)*m1,:)*Be;
    end
    %v=h*Be; %first column of the Toeplitz matrix Phi
    Phi=zeros(Np*size(Ce*Ae*Be,1),Nc*size(Ce*Ae*Be,2)); %declare the dimension of Phi
    Phi(:,1:size(Ce*Ae*Be,2))=v;
    for i=2:Nc
        Phi(:,(i-1)*size(Ce*Ae*Be,2)+1:i*size(Ce*Ae*Be,2))=[zeros((i-1)*size(Ce*Ae*Be,1),size(Ce*Ae*Be,2))', v(1:(Np-i+1)*size(Ce*Ae*Be,1),:)']'; %Toeplitz matrix
    end
    BarRs=ones(Np*m1,1); %fixed set-point during Np samples
    Phi_Phi=Phi'*Phi;
    Phi_F=Phi'*F;
    Phi_R=Phi'*BarRs;
end