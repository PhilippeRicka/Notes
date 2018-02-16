clear

T = 100;
dt = pi/10;
time = (0:T-1)*dt;

Nloop = 10;
ukfeff = 0;
ukft = 0;
enkfeff = 0;
enkft = 0;

for L=1:Nloop

    x = 5 + rand*ones(1,T);
    for t=2:T
        x(t) = x(t-1) + (1+randn(1)/10)*cos(t*dt)*dt;    % state is a slightly noised sinusoid
    end

    sd_mes = 0.4;               % measurements are this much accurate
    y = x + randn(1,T)*sd_mes;  % simulated measurements

            % UKF

    mu = ones(1,T)*y(1);        % first state estimate is the measurement
    Sigma = 1;                  % first covariance estimate is arbitrary or fed by extra knowledge

    n = 1;                      % UKF parameter
    kappa = 2;                  %      ''
    
    tic();

    W = [kappa/(n+kappa) , 1/(2*(n+kappa)) , 1/(2*(n+kappa))];  % constant weights

    for t = 1:T-1

        % sigma-points
        X = mu(t) + [0 ; sqrt((n+kappa)*Sigma) ; -sqrt((n+kappa)*Sigma)];

        % unscented transform
        Xf = X + cos(t*dt)*dt;

        % mean and covariance predictions
        muf = W*Xf;
        Sigmabar = W*((Xf-muf).^2);

        % measurement predictions
        Yf = Xf;%mubar + [0 ; sqrt((n+kappa)*Sigmabar) ; -sqrt((n+kappa)*Sigmabar)];
        yhat = W*Yf;

        % innovation covariance
        S = W*((Yf - yhat).^2) + Sigma;
        SigmabarX = W*((Xf-muf).*(Yf-yhat));
        K = SigmabarX/S;

        mu(t+1) = muf + K*(y(t) - yhat);
        %Sigma = Sigmabar - K*S*K;

    end
    
    ukftime = toc();
    ukfefficiency = sqrt(mean((x-y).^2)/mean((x-mu).^2));   % efficiency is the L2 norm ratio of measurement deviation by estimation deviation
    subplot(121)
    plot(time,x,time,y,time,mu)

            % EnKF

    Ntest = 10;
    enkftime = 0;
    enkfefficiency = 0;
    qens = 30;

    mu = ones(1,T)*y(1);        % first state estimate is the measurement
    W = ones(1,qens)/qens;
    
    for I=1:Ntest
        
        tic();

        Xf = randn(qens,1)*sqrt(abs(y(1)))/100 + y(1);
        Xa = Xf;
        Xa_bar = mean(Xa);

        y_kal = zeros(qens,1); % initialisation

        for t=2:T
 
            % propagation
            Xf = Xa + cos(t*dt)*dt;

            % mean predictions
            muf = W*Xf;

            % measurement predictions
            Yf = Xf;
            yhat = W*Yf;
            
            % measurement noising
            E = randn(qens,1)*sqrt(abs(y(t)))/100;
            y_kal = E + y(t)*ones(qens,1);

            % innovation covariance
            S = var(Yf) + var(E);
            SigmabarX = W*((Xf-muf).*(Yf-yhat));
            K = SigmabarX/S;
            
            for q = 1:qens
                Xa(q) = Xf(q) + K*(y_kal(q) - Yf(q));
            end

            mu(t) = mean(Xa);

        end

            enkftime = enkftime + toc();

            eff = sqrt(mean((x-y).^2)/mean((x-mu).^2));   % efficiency is the L2 norm ratio of measurement deviation by estimation deviation
            enkfefficiency = enkfefficiency + eff;
        
    end

    enkftime = enkftime/Ntest;
    enkfefficiency = enkfefficiency/Ntest;
    subplot(122)
    plot(time,x,time,y,time,mu)
    
    ukfeff = ukfeff + ukfefficiency;
    enkfeff = enkfeff + enkfefficiency;
    ukft = ukft + ukftime;
    enkft = enkft + enkftime;
    
end

ukfeff = ukfeff/Nloop;
enkfeff = enkfeff/Nloop;
ukft = ukft/Nloop;
enkft = enkft/Nloop;

disp(' ')
disp('Ensemble Kalman Filter')
disp(' ')
disp(['    mean time = ',num2str(enkft),' s'])
disp(['    mean efficiency = ',num2str(enkfeff)])
disp(' ')
disp('Unscented Kalman Filter')
disp(' ')
disp(['    mean time = ',num2str(ukft),' s'])
disp(['    mean efficiency = ',num2str(ukfeff)])
disp(' ')
plot(time,p,time,mu)