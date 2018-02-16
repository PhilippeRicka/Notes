function P = step_p(p,q,r1,r2,c,t,dt)
    alpha = 1/(c/(dt) + 1/r2);
    beta = 1 + r1/r2 + r1*c/(dt);
    gamma = r1*c/(dt);
    
    P = alpha*(beta*q(t)-gamma*q(t-1)+c*p/(dt));
end