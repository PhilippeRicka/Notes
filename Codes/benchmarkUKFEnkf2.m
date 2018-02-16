clear

format long

T = 100;
dt = pi/10;
time = (0:T-1)*dt;

R1 = 1600;
R2 = 13000;
C = 2.5e-5;

p = 5 + rand*ones(1,T);
q=0*p;
for t=2:T
    p(t) = p(t-1) + (1+randn(1)*sqrt(p(1))/10)*cos(t*dt)*dt;    % nodal pressure is a slightly noised sinusoid
end
for t=1:T-1
    q(t+1) = (p(t+1)/R2+C*(p(t+1)-p(t))/dt); % corresponding RCR input current
end
q(1)=q(T);

q = q + randn(1,T)*sqrt(mean(q))/500;

        % UKF

mu = ones(1,T)*q(1)*(C+R2*dt)/(C*R2+dt*R2^2);       % first state estimate
Sigma = 1;                  % first covariance estimate is arbitrary or fed by extra knowledge

n = 1;                      % UKF parameter
kappa = 2;                  %      ''

W = [kappa/(n+kappa) , 1/(2*(n+kappa)) , 1/(2*(n+kappa))];  % constant weights


for t = 1:T-1

    % sigma-points around state
    X = mu(t) + [0 ; sqrt((n+kappa)*Sigma) ; -sqrt((n+kappa)*Sigma)];

    % unscented transform to get predicted state (time update)
    Xbar = R2*dt*(q(t+1)+C*X/dt)/(C*R2+dt);

    % mean and covariance predictions
    mubar = W*Xbar;
    Sigmabar = W*((Xbar-mubar).^2);

    % measurement predictions
    Yf = (1/R2+C/dt)*Xbar - mu(t)*C/dt;
    yhat = W*Yf;

    % innovation covariance
    Cpy = W*((Xbar - mubar).*(Yf - yhat));
    Cyy = W*((Yf - yhat).*(Yf - yhat)) + Sigma;
    K = Cpy/Cyy;

    mu(t+1) = mubar + K*(q(t+1) - yhat);
    %Sigma = Sigmabar - K*Cyy*K';

end

subplot(121)
plot(time,p,time,mu)


    % EnKF

qens = 20;

muen = ones(1,T)*q(1)*(C+R2*dt)/(C*R2+dt*R2^2);       % first state estimate
Sigma = 1;                  % first covariance estimate is arbitrary or fed by extra knowledge

W = ones(1,qens)/qens;

Nloop=1;
meanenerror = 0;
for N = 1:Nloop
    
    Xf = randn(qens,1)*sqrt(abs(q(1))) + q(1);
    Xa = Xf;
    Xabar = W*Xa;
    
    y_kal = zeros(qens,1);
    
for t = 2:T

    % propagation
    Xf = R2*dt*(q(t)+C*Xa/dt)/(C*R2+dt);

    % mean predictions
    muf = W*Xf;
    
    % measurement predictions
    Yf = (1/R2+C/dt)*Xf - muen(t-1)*C/dt;
    yhat = W*Yf;

    % measurement noising
    E = randn(qens,1)*sqrt(abs(q(t)))/100;
    y_kal = E + q(t)*ones(qens,1);

    % innovation covariance
    Cpy = W*((Xf - muf).*(Yf - yhat));
    Cyy = var(Yf) + var(E);
    K = Cpy/Cyy;

    for i = 1:qens
        Xa = Xf + K*(y_kal(i) - Yf(i));
    end
    
    muen(t) = mean(Xa);

end

meanenerror = meanenerror + sqrt(mean((p-muen).^2));
end
meanenerror = meanenerror/Nloop;

subplot(122)
plot(time,p,time,mu)

disp('Mean L2 norm of estimation deviation :')
disp(['    UKF : ',num2str(sqrt(mean((p-mu).^2)))])
disp(['    EnKF : ',num2str(meanenerror)])