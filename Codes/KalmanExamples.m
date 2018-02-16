clear
clf
format long

nloop = 5;
T = 0.8;
dt = 0.01;
C = 0.000025;
R1 = 1600;
R2 = 13000;
alpha = 1/(C/(dt) + 1/R2);
beta = 1 + R1/R2 + R1*C/(dt);
gamma = R1*C/(dt);

time = 1:T/dt;
q = (10*exp(-100*(time*dt-0.38).^2)+5)/100000;
p = 0*time;

T = length(q)*dt;
time = 1:T/dt;




for i=1:nloop
    for t=1:T/dt-1
        p(t+1) = alpha*(beta*q(t+1)-gamma*q(t)+C*p(t)/(dt));
    end
p(1)=p(T/dt);
end



    % Premier exemple : on bruite p (sd=4mmHg) puis on le filtre
    
    y = p + randn(1,T/dt)*4/133.322;
    est_p = ones(1,T/dt)*y(1);
    
    qens = 10;
    
    pf = randn(1,qens)*sqrt(y(1))/100 + y(1);
    pa = pf;
    pa_bar = mean(pa);
    
    y_kal = zeros(1,qens); % initialisation
    
    for t = 1:T/dt
        Pxy = 0;
        Pyy = 0;
        
        pf = pa + randn(1,qens).*sqrt(pf)/100;
        
        E = randn(1,qens)*sqrt(y(t))/100;
        y_kal = E + y(t)*ones(1,qens); % bruitage artificiel
        
        Pxy = var(pf);
        Pyy = Pxy + var(E);
        K = Pxy/Pyy;
        
        for i = 1:qens
            pa(i) = pf(i) + K*(y_kal(i) - pf(i));
        end
        
        est_p(t) = mean(pa);
    end
    
%     subplot(211)
%     plot(time,p,time,y,time,estim_p)
%     eff = mean((p-y).^2)/mean((p-estim_p).^2);
%     aff = ['efficacite pour p : ',num2str(eff)];
%     disp(aff)
%     
    
    
    % Deuxieme exemple : p est une mesure, on cherche q connaissant les parametres
    
    y = p;
    estim_q = ones(1,T/dt)*5/100000;
    
    
    qf = 5*ones(1,qens)/100000;
    qa = qf;
    qa_bar = mean(qa);
    
    y_kal = zeros(1,qens);

    for t = 2:T/dt
        Pxy = 0;
        Pyy = 0;
        
        qf = estim_q(t-1) + randn(1,qens).*sqrt(abs(qa))/100;
        
        E = randn(1,qens)*sqrt(y(t))/100;
        y_kal = E + y(t)*ones(1,qens); % bruitage artificiel
        
        Pxy = cov(qf);
        Pyy = Pxy + cov(E);
        K = Pxy/Pyy;
        
        for i = 1:qens
            qa(i) = qf(i) + K*(y_kal(i) - alpha*(beta*qf(i)-gamma*qa(i)+C*p(t-1)/dt));
        end
        
        test(t) = Pxy;
        estim_q(t) = mean(qa);
    end
    
%     qudev = mean((q-estim_q).^2);
%     aff = ['ecart quadratique moyen avec q : ',num2str(qudev)];
%     disp(aff)
%     subplot(212)
%     plot(time,q,time,estim_q)
%    
    
    
    
    % Troisieme exemple : UKF

    clear
    
    T = 100;
    dt = 10*pi/T;
    time = (0:T-1)*dt;
    
    trajectories = zeros(3,T);
    
    sd = 1;
    sd_mes = 0.4;
    n = 1;
    kappa = 2;
    
    x = zeros(1,T);
    for t=2:T
        x(t) = x(t-1) + (1+randn(1)/4)*cos(t*dt)*dt;
    end
    
    y = x + randn(1,T)*sd_mes;
    mu = ones(1,T)*y(1);
    
    Sigma = sd^2;
    W = [kappa/(n+kappa) , 1/(2*(n+kappa)) , 1/(2*(n+kappa))];
    tic();
    for t = 1:T-1
        
        % sigma-points
        X = mu(t) + [0 ; sqrt((n+kappa)*Sigma) ; -sqrt((n+kappa)*Sigma)];
   
        % unscented transform
        Xbar = X + cos(t*dt)*dt;
        
        % mean and covariance predictions
        mubar = W*Xbar;
        Sigmabar = W*((Xbar-mubar).^2);
        
        % measurement predictions
        Ybar = mubar + [0 ; sqrt((n+kappa)*Sigmabar) ; -sqrt((n+kappa)*Sigmabar)];
        yhat = W*Ybar;

        % innovation covariance
        S = W*((Ybar - yhat).^2) + sd_mes^2;
        SigmabarX = W*((Xbar-mubar).*(Ybar-yhat));
        K = SigmabarX/S;
        
        mu(t+1) = mubar + K*(y(t) - yhat);
        Sigma = Sigmabar - K*S*K;
        
    end
    toc();
    
    
    clf
    plot(time,x,time,y,time,mu)%,time,trajectories(1,:),time,trajectories(2,:),time,trajectories(3,:))
    legend('x','y','estimation x')%,'\sigma point 1','\sigma point 2','\sigma point 3')
    eff = mean((x-y).^2)/mean((x-mu).^2);
    aff = ['efficacite : ',num2str(eff)];
    disp(aff)

