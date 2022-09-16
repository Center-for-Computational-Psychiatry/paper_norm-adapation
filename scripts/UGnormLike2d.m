function [ Like , museries] = UGnormLike2d(offers, rejects, param, x )
%Norm adaptation model in Ultimatum Game - Bayesian Observer
%variable initial Norm
%Andreas Hula, 11.December.2014
n = length(offers);
Like =0;
k = param(3);
mu= 20*exp(x(3))/(1+exp(x(3)));
envy = exp(x(1))/(1+exp(x(1)));
temp = exp(x(2))/(1+exp(x(2)));
museries = zeros(1,n+1);
museries(1) = mu;

    for i = 1:n  
        k = k+1;
        mu =(k-1)/k*mu +1/k*offers(i);
        museries(i+1)=mu;
        U = offers(i)-envy*max(mu-offers(i),0);
        Like = Like -log(exp(temp*U*(rejects(i)))/(1+exp(temp*U)));        
    end

end