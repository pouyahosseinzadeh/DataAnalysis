clear all
close all
% Pouya Hosseinzadeh

% Load the data
load dati_silver

% Input data
sig2 = 1;     % Variance of measurement noise
q = 0.00001;    % Variance of White process w for state X3
%P0 = 0.5;     % Initial condition for covariance matrix of estimation error
%x0 =  [3 0.5 0.5]; 
%x0 =  [-1 2 0.5]; % initial condition for the states
x0 =  [-4 5 0.2]; % Initial condition for the states
%x0 =  [-1.2 0.1 0.2]; 
x = x0;

% Matrix G of the system linearized df/dw
G =[0 0 0;0 0 (x(2))^2; 0 0 1];

% Covariance matrix of the disturbance w
Q =zeros(3);
Q(3,3)=q;

% Covariance matrix of v
R = sig2*eye(1);

% Initialization
N=size(U,1);
xp=zeros(3,1);
%P = P0*eye(3);
P=diag([1.5 1 0.7]);

Xest(:,1)=x;
PDiag(:,1)=diag(P);

% Recursion
for k=1:N-1
    
    % Prediction step
    xp(1) = x(1) + x(2);
    xp(2) = -x(1) + 0.5*x(2) - (x(3))*(x(2))^3 + U(k);
    xp(3) = x(3);
    
    F = [1 1 0; -1  x(2)-3*x(2)*x(3) -(x(2))^3; 0 0 1]; 
    Pp = F*P*F'+ G*Q*G';
    
    % Correction step
    H = [0 1 0];
    K = Pp*H'*inv(H*Pp*H'+R);
    xc = xp + K*(Y(k+1,:)'-xp(2));
    Pc = Pp*(eye(3)-H'*K');
    x=xc;
    P=Pc;
    
    Xest(:,k+1)=x;
    PDiag(:,k+1)=diag(P);
    
end

t = 1:150;

% Comparison between state X and the estimation Xest
 
figure, plot(t,Xest(1,:),t,X(:,1),'g'),legend('Estimation for X1','True value of X1')
ylabel('Value of state'), xlabel('Time')

figure, plot(t,Xest(2,:),t,X(:,2),'g'),legend('Estimation for X2','True value of X2')
ylabel('Value of state'), xlabel('Time')



%Comparing the estimation error with the confidence intervals +/- 3*sqrt(P(K|K))

figure, plot(t,X(:,1)-Xest(1,:)','r',t,3*sqrt(PDiag(1,:)),'b--',t,-3*sqrt(PDiag(1,:)),'b--'), grid on
legend('Estimation error of state variable X1','Confidence interval')
ylabel('Estimation error'), xlabel('Time')

figure, plot(t,X(:,2)-Xest(2,:)','r',t,3*sqrt(PDiag(2,:)),'b--',t,-3*sqrt(PDiag(2,:)),'b--'), grid on
legend('Estimation error of state variable X2','Confidence interval')
ylabel('Estimation error'), xlabel('Time')


figure, plot(t,Xest(3,:),'g'),legend('Estimation for X3')
ylabel('Value of state'), xlabel('Time')


figure, plot(t,0.1*ones(size(t))-Xest(3,:)','r',t,3*sqrt(PDiag(3,:)),'b--',t,-3*sqrt(PDiag(3,:)),'b--'), grid on
legend('Estimation error of state variable X2','Confidence interval')
ylabel('Estimation error'), xlabel('Time')
