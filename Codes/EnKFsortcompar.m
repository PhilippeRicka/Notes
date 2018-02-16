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

            % EnKF unsorted

qens = 100;

mu = ones(1,T)*q(1)*(C+R2*dt)/(C*R2+dt*R2^2);       % first state estimate
Sigma = 1;                  % first covariance estimate is arbitrary or fed by extra knowledge

W = ones(1,qens)/qens;

Nloop=1;
meanenerror = 0;

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
    Yf = (1/R2+C/dt)*Xf - mu(t-1)*C/dt;
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
    
    mu(t) = mean(Xa);

end

unserror = sqrt(mean((p-mu).^2));



            % EnKF sorted

qens = 100;

mu = ones(1,T)*q(1)*(C+R2*dt)/(C*R2+dt*R2^2);       % first state estimate
Sigma = 1;                  % first covariance estimate is arbitrary or fed by extra knowledge

W = ones(1,qens)/qens;

Nloop=1;
meanenerror = 0;

Xf = randn(qens,1)*sqrt(abs(q(1))) + q(1);
Xf = sort(Xf);
Xa = Xf;
Xabar = W*Xa;

y_kal = zeros(qens,1);

for t = 2:T

    % propagation
    Xf = R2*dt*(q(t)+C*Xa/dt)/(C*R2+dt);
    Xf = sort(Xf);

    % mean predictions
    muf = W*Xf;
    
    % measurement predictions
    Yf = (1/R2+C/dt)*Xf - mu(t-1)*C/dt;
    Yf = sort(Yf);
    yhat = W*Yf;

    % measurement noising
    E = randn(qens,1)*sqrt(abs(q(t)))/100;
    y_kal = E + q(t)*ones(qens,1);
    y_kal = sort(y_kal);

    % innovation covariance
    Cpy = W*((Xf - muf).*(Yf - yhat));
    Cyy = var(Yf) + var(E);
    K = Cpy/Cyy;

    for i = 1:qens
        Xa = Xf + K*(y_kal(i) - Yf(i));
    end
    
    Xa = sort(Xa);
    
    mu(t) = mean(Xa);

end

serror = sqrt(mean((p-mu).^2));

disp(['Error unsorted filter : ',num2str(unserror)])
disp(['  Error sorted filter : ',num2str(serror)])
disp(['Relative accuracy gain via sorting ensemble : ',num2str(unserror/serror)])