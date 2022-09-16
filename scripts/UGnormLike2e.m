function [ Like, museries ] = UGnormLike2e(offers, rejects, param, x )
%Norm adaptation model in Ultimatum Game - fixed norm Fehr-Schmidt
%non adapting reference
%Andreas Hula, 11.December.2014
n = length(offers);
Like =0;
mu= 10;
envy = exp(x(1))/(1+exp(x(1)));
temp = exp(x(2))/(1+exp(x(2)));
museries = zeros(1,n+1);
museries(1) = mu;
    for i = 1:n
        U = offers(i)-envy*max(mu-offers(i),0);
        museries(i+1)=mu;        
        Like = Like -log(exp(temp*U*(rejects(i)))/(1+exp(temp*U)));        
    end

end