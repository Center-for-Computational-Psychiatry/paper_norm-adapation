function [ Like, museries ] = UGnormLike2g(offers, rejects, param, x )
%Norm adaptation model in Ultimatum Game - Rescorla-Wagner based
%Norm adaptation model, variable initial norm
%Andreas Hula, 11.December.2014
n = length(offers);
Like =0;
mu= 20*exp(x(4))/(1+exp(x(4)));
envy = exp(x(1))/(1+exp(x(1)));
adaptation = exp(x(2))/(1+exp(x(2)));
temp = exp(x(3))/(1+exp(x(3)));
museries = zeros(1,n+1);
museries(1) = mu;
    for i = 1:n
        mu = mu +adaptation*(offers(i)-mu); 
        museries(i+1)=mu;
        U = offers(i)-envy*max(mu-offers(i),0);
        Like = Like -log(exp(temp*U*(rejects(i)))/(1+exp(temp*U))); 
    end

end