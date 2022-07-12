%1 - Define the state space system parameters
%2 - Set up the input signal, initial conditions, etc.
x0 = [0.1 0.2]';
n1=size(x,1);
Ns = Nc;
dUreg=inv(Phi_Phi+10*eye(size(Phi_Phi,1)))*(Phi_R*0-Phi_F*x0);
dU=inv(Phi_Phi)*(Phi_R*0-Phi_F*x0);
%y=zeros(Ns,1);
%3 - Perform the system simulation
Yreg=F*x+Phi*dUreg;%Output
Y=F*x+Phi*dU;%Output
%stairs(Yreg)
%hold
%stairs(Y)
u=[0,cumsum(dU)']';
x=x0;
X(:,1)=x;
for i=2:Ns
    X(:,i)=Ac*X(:,i-1)+Bc*u(i);%State trajectory
end
plot(X')
hold
u=[0,cumsum(dUreg)']';
x=x0;
X(:,1)=x;
for i=2:Ns
    X(:,i)=Ac*X(:,i-1)+Bc*u(i);%State trajectory
end
plot(X')
%Y=reshape(Y,Np,size(Ac,1))
%dU=reshape(dU,[Np,size(Bc,2)])
%sys = ss(Ac,Bc,Cc,Dc,Ts);
%t=1:1:Ns
%lsim(sys,cumsum(dU),t*Ts,x0)