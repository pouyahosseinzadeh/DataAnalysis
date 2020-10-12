clear all
close all
% Pouya Hosseinzadeh

% Load the data
load dati_silver

% Input data
L = 0.1;
sig2 = 1;   % Variance of measurement noise
x0 = [-1 2];    % initial condition -6,4
x = x0;

% Matrix G of the system linearized df/dw
G =zeros(2); % w = 0

% Covariance matrix of w
Q =zeros(2); % No disturbance, w = 0

% Covariance matrix of v
R = sig2*eye(1);

% Initialization
N=size(U,1);
xp=zeros(2,1);
P=diag([10 4]);         % 0.5, 0.5


Xest(:,1)=x;
PDiag(:,1)=diag(P);

% Recursion
for k=1:N-1
    
    % Prediction step
    xp(1) = x(1) + x(2);
    xp(2) = -x(1) + 0.5*x(2) - L*(x(2))^3 + U(k);
    F = [1 1; -1 0.5-3*L*(x(2))^2];    
    Pp = F*P*F'+ G*Q*G';
   
    % Correction step
    H = [0 1];
    K = Pp*H'*inv(H*Pp*H'+R);
    xc = xp + K*(Y(k+1,:)'-[xp(2)]);
    Pc = Pp*(eye(2)-H'*K');
    x=xc;
    P=Pc;
    
    Xest(:,k+1)=x;
    PDiag(:,k+1)=diag(P);
    
end


t = 1:150;

% Comparing the true state X and the estimation Xest

figure, plot(t,Xest(1,:),t,X(:,1),'g'),legend('Estimation for X1','True value of X1')
ylabel('Value of state'), xlabel('Time')

figure, plot(t,Xest(2,:),t,X(:,2),'g'),legend('Estimation for X2','True value of X2')
ylabel('Value of state'), xlabel('Time')



% Comparing the estimation error with the confidence intervals +/- 3*sqrt(P(K|K))

figure, plot(t,X(:,1)-Xest(1,:)','r',t,3*sqrt(PDiag(1,:)),'b--',t,-3*sqrt(PDiag(1,:)),'b--'), grid on
legend('Estimation error of state variable X1','Confidence interval')
ylabel('Estimation error'), xlabel('Time')

figure, plot(t,X(:,2)-Xest(2,:)','r',t,3*sqrt(PDiag(2,:)),'b--',t,-3*sqrt(PDiag(2,:)),'b--'), grid on
legend('Estimation error of state variable X2','Confidence interval')
ylabel('Estimation error'), xlabel('Time')
