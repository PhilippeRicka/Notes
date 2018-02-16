clear

format long

T = 100;
dt = pi/10;
time = (0:T-1)*dt;

R1 = 1600;
R2 = 13000;
C = 2.5e-5;



% state and measurement simulations

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




    % Freibourg

    tic();
    
for t = 2:T
    
    % sigma-points
    Xa = mu(t-1) + [0 ; sqrt((n+kappa)*Sigma) ; -sqrt((n+kappa)*Sigma)];
    
    % time update
    Xf = R2*dt*(q(t)+C*Xa/dt)/(C*R2+dt);
    
    muXf = W*Xf;
    SigmaXf = W*((Xf-muXf).^2);
    
    % measurement predictions
    Ya = muXf + [0 ; sqrt((n+kappa)*SigmaXf) ; -sqrt((n+kappa)*SigmaXf)];
    Yf = (1/R2+C/dt)*Ya - mu(t-1)*C/dt;
    
    muYf = W*Yf;
    
    % innovation and Kalman gain
    Cyy = W*((Yf - muYf).^2) + Sigma;
    Cxy = W*((Xf - muXf).*(Yf - muYf));
    K = Cxy*inv(Cyy);
    
    % estimation
    mu(t) = muXf + K*(q(t) - muYf);
    % Sigma = SigmaXf - K*Cyy*K';
    
end

disp(['Freibourg : ',num2str(sqrt(mean((p-mu).^2))),'    time : ',num2str(toc())])
subplot(221)
plot(time,q)
subplot(222)
plot(time,p,time,mu)


    % Oregon
    
    tic();
    
for t = 2:T
    
    % sigma-points
    Xa = mu(t-1) + [0 ; sqrt((n+kappa)*Sigma) ; -sqrt((n+kappa)*Sigma)];
    
    % time update
    Xf = R2*dt*(q(t)+C*Xa/dt)/(C*R2+dt);
    
    muXf = W*Xf;
    % SigmaXf = W*((Xf-muXf).^2);
    
    % measurement predictions
    Yf = (1/R2+C/dt)*Xf - mu(t-1)*C/dt;
    
    muYf = W*Yf;
    
    % innovation and Kalman gain
    Cyy = W*((Yf - muYf).^2);
    Cxy = W*((Xf - muXf).*(Yf - muYf));
    K = Cxy*inv(Cyy);
    
    % estimation
    mu(t) = muXf + K*(q(t) - muYf);
    % Sigma = SigmaXf - K*Cyy*K';
    
end

disp(['   Oregon : ',num2str(sqrt(mean((p-mu).^2))),'    time : ',num2str(toc())])
subplot(223)
plot(time,q)
subplot(224)
plot(time,p,time,mu)