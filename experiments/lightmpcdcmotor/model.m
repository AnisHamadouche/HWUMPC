clear all
num=[0.1];
denum=[1 0 0 0 0 0 0 -0.8];
[Ap,Bp,Cp,Dp]=tf2ss(num,denum);
[n1,n_in]=size(Bp);
[m1,n1]=size(Cp);
A=eye(n1+m1,n1+m1);
A(1:n1,1:n1)=Ap;
A(n1+1:n1+m1,1:n1)=Cp*Ap;
B=zeros(n1+m1,n_in);
B(1:n1,:)=Bp;
B(n1+1:n1+m1,:)=Cp*Bp;
C=zeros(m1,n1+m1);
C(:,n1+1:n1+m1)=eye(m1,m1);