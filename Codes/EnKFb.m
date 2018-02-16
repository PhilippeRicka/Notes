clear
clf()

g=-9.81;
qens=25;
I = [1 0 0; 0 1 0; 0 0 1];

% Vecteur temps
dt=0.2;
T=5;
t =[0:dt:T];

% Initialisation vraies valeurs x
x = 0 * [t;t;t];
% Matrice de transition relative ? la dynamique r?elle
A = [ 1 dt dt*dt/2; 0 1 dt; 0 0 1];

% Valeurs mesur?es
y = 0*x;

sd_bruit_x = 1;
sd_bruit_v = 0.5;
sd_bruit_a = 0.25;
Q = [ sd_bruit_x*sd_bruit_x 0 0; 0 sd_bruit_v*sd_bruit_v 0; 0 0 sd_bruit_a*sd_bruit_a ];

% Etat du processus
etat_vrai = [ x(1,1);  x(2,1); x(3,1) ];

% Etat r?el
x(3,1)=g;
for i = 2:length(t)
    x_bruit = grand(1, 'mn', [0;0;0], Q)
    etat_vrai = [ x(1,i-1);  x(2,i-1); x(3,i-1) ];
    etat_vrai = A * etat_vrai + x_bruit;
    x(:,i) = etat_vrai(:);
end

% Mesures
sd_mes_x = 5;
sd_mes_v = 5;
sd_mes_a = 10;
R = [ sd_mes_x*sd_mes_x 0 0; 0 sd_mes_v*sd_mes_v 0; 0 0 sd_mes_a*sd_mes_a ];

for k=1:length(t)
    y_bruit = grand(1, 'mn', [0;0;0], R);
    y_vrai(:,k) = x(:,k) + y_bruit;
end
for i=1:qens
    e=grand(1,'mn', [0;0;0], R.*1/qens);
    y(3*i-2:3*i,1) = y_vrai(:,1)+e;
    E(:,i) = e;
end

% Algo
xf = 0*y ;
xf(:,1) = y(:,1);

% Calcul de la moyenne
xf_bar = [0;0;0];
for i = 1:qens
    xf_bar = xf_bar + xf(3*i-2:3*i,1);
end
xf_bar = xf_bar/qens;

xa = xf;
for k = 2:length(t)
    Dk = [(y(2,1)*dt-g*(2*k-1)*dt)^2 0 0 ; 0 (-g*dt)^2 0 ; 0 0 0.1];
    xf_bar = [0;0;0];
    for i = 1:qens
        delta = grand(1, 'mn', [0;0;0], Dk);
        xf(3*i-2:3*i,k) = xa(3*i-2:3*i,k) + delta;
        xf_bar = xf_bar + xf(3*i-2:3*i,k);
    end
    xf_bar = xf_bar/qens;
     
    Pxyk = 0*I;
    for i = 1:qens
        e=grand(1,'mn', [0;0;0], R/qens);
        y(3*i-2:3*i,k) = y_vrai(:,k)+e;
        E(1:3,i) = e;
        Pxyk = Pxyk + (xf(3*i-2:3*i,k) - xf_bar)*(I*xf(3*i-2:3*i,k) - I*xf_bar)';
    end
    Pxyk = Pxyk/(qens-1);
    % Dans notre cas, Pxyk = Pyyk - Rk
    Rk = E*E'/(qens-1);
    Rk = [Rk(1,1) 0 0 ; 0 Rk(2,2) 0 ; 0 0 Rk(3,3)];
    Pyyk = Pxyk + Rk;
    
    K = Pxyk*inv(Pyyk);

    for i = 1:qens
        if k<length(t)
            xa(3*i-2:3*i,k+1) = xf(3*i-2:3*i,k) + K * (y(3*i-2:3*i,k) - xf(3*i-2:3*i,k));
        else 
            xa(3*i-2:3*i,k+1) = xf(3*i-2:3*i,k);
        end
    end

  % Pr?diction
  pre(1:3,k) = xf_bar;%xa(1:3,k);

  % Erreur de pr?diction
%  err(i) = norm(x(:,i)-xf(:,i),2)/norm(x(:,i)-an(:,i),2);
end

% Bleu  : r?el
% Vert  : mesur?
% Rouge : pr?dit
%subplot(313)
%plot(t, x(3,:), t, y(3,:), t, xa(3,:));  
subplot(212)
plot(t, x(2,:), t, y_vrai(2,:), t, pre(2,:));  
subplot(211)
plot(t, x(1,:), t, y_vrai(1,:), t, pre(1,:));
%subplot(223)
%plot(t, pre(1,:))
%subplot(224)
%plot(t, pre(2,:))
%subplot(223)
%plot(t,err)
