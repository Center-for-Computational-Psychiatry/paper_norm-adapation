load('UG_mindful_new_59ss_2014.mat')
load('NewInd.mat')
%Implementation of data fit for 
%lesion data set. Ultimatum Game Data
%All Models
%Andreas Hula, 11.December.2014 
norounds = 46; %Number of Rounds
nosubjects=59; %Number of Subjects
options =optimset('TolFun',10^(-14),'TolX',10^(-14),'Largescale','off','DerivativeCheck','off','GradObj','off','MaxFunEvals',2000,'MaxIter',400);

    offers= zeros(nosubjects, norounds); %store offer data
    rejects=zeros(nosubjects, norounds); %store rejection data
    val = zeros(nosubjects,3, 5); %parameter values returned by the optimization
    result = zeros(nosubjects,3, 5); %The transformed parameter values
    Like =1000+zeros(nosubjects,5); %Calculate Likelihoods
    BIC = zeros(nosubjects,5);  %Calculate BIC
    NormCurve = zeros(nosubjects, norounds+1, 5); %Evolution of the norm curves
       for i = 1:nosubjects
            for j=1:length(ug2{i+1,3})%j = find(trialtype(i,:)==1)
                offers(i,j) =ug2{i+1,3}(j); % read out offers 
                rejects(i,j)=ug2{i+1,4}(j); % raw data: reject=1 accept =0; now: reject=0 accept=1
            end
       end
 
    param =[10,2,4,10]; %Starting values for a Bayesian Observer       
       
       
    v = [0 0 0 0 0 0]; %starting point for optimization
for i =1:nosubjects 
    for s =1:18 %try multiple starting points - Determined empiricially
            if s==1
               v =[0 0 0 0 1 1];
            end
              if s==2
              v =[1 0 0 0 1 0];
              end
            if s==3
              v =[1 -1 0 -1 0 0];
            end
             if s==4
              v =[1 1 0 1 0 0];
             end  
             if s==5
                v =[1 1 0 0 1 0];
             end
                if s==6
                v =[-1 -1 0 -2 -2 -2];
                end
                if s==7
                v =[1 1 0 -1 0 0];
                end            
             if s==8
                v =[1 -1 0 -1 1 1];
             end
             if s==9
                v =[1 0 0 0 0 -1];
             end   
             if s==10
                v =[1 0 1 0 0 0];
             end  
           if s==11
                v =[1 0 1 1 1 0];
            end
              if s==12
                v =[-1 -1 0 1 -1 0];
              end
            if s==13
                v =[1 -1 1 1 -1 0];
            end
             if s==14
                v =[1 1 0 1 0 0];
             end  
             if s==15
                v =[1 1 0 -1 1 0];
             end
                if s==16
                v =[1 1 -1 -1 0 0];
                end
             if s==17
                v =[1 1 -1 0 1 0];
             end
             if s==18
                v =[-1 -1 0 1 -1 0];         
             end             
            y = v(1:2);
        [y,Likeold,~,~,~,~] = fminunc(@(y) ... %minimisation of negative log likelihoods 
            UGnormLike2e(ug2{i+1,3},ug2{i+1,4},param, y),y,options); %Fehr-Schmidt
                                                %non adopting model
        if (Like(i,1)>Likeold) %only keep the best likelihood and parameters
            Like(i,1)=Likeold; %replace negative loglikelihood with better 
          val(i,1:length(y),1)=y; %keep better parameters 
        end
              y = v(1:3);
        [y,Likeold,~,~,~,~] = fminunc(@(y) ... %minimisation of negative log likelihoods 
            UGnormLike2a(ug2{i+1,3},ug2{i+1,4},param, y),y,options); %Rescorla-Wagner
                                                %fixed initial norm model
        if (Like(i,2)>Likeold) %only keep the best likelihood and parameters
            Like(i,2)=Likeold; %replace negative loglikelihood with better 
          val(i,1:length(y),2)=y; %keep better parameters 
        end  
            y = v(1:4);
         [y,Likeold,~,~,~,~] = fminunc(@(y) ... %minimisation of negative log likelihoods 
            UGnormLike2g(ug2{i+1,3},ug2{i+1,4},param, y),y,options); %Rescorla-Wagner
                                                %variable initial norm model
        if (Like(i,3)>Likeold) %only keep the best likelihood and parameters
            Like(i,3)=Likeold; %replace negative loglikelihood with better 
           val(i,1:length(y),3)=y; %keep better parameters 
        end    
            y = v(1:2);
        [y,Likeold,~,~,~,~] = fminunc(@(y) ... %minimisation of negative log likelihoods 
            UGnormLike2b(ug2{i+1,3},ug2{i+1,4},param, y),y,options); %Bayes-Observer
                                                %fixed initial norm model
        if (Like(i,4)>Likeold) %only keep the best likelihood and parameters
            Like(i,4)=Likeold; %replace negative loglikelihood with better 
          val(i,1:length(y),4)=y; %keep better parameters 
        end  
            y = v(1:3);
        [y,Likeold,~,~,~,~] = fminunc(@(y) ... %minimisation of negative log likelihoods 
            UGnormLike2d(ug2{i+1,3},ug2{i+1,4},param, y),y,options); %Bayes-Observer
                                                %variable initial norm model
        if (Like(i,5)>Likeold) %only keep the best likelihood and parameters
            Like(i,5)=Likeold; %replace negative loglikelihood with better 
          val(i,1:length(y),5)=y; %keep better parameters 
        end         
    end    

    [~,NormCurve(i,1:(length(ug2{i+1,3})+1),1)] = UGnormLike2e(ug2{i+1,3},ug2{i+1,4},param, val(i,1:2,1));
    [~,NormCurve(i,1:(length(ug2{i+1,3})+1),2)] = UGnormLike2a(ug2{i+1,3},ug2{i+1,4},param, val(i,1:3,2));
    [~,NormCurve(i,1:(length(ug2{i+1,3})+1),3)] = UGnormLike2g(ug2{i+1,3},ug2{i+1,4},param, val(i,1:4,3));
    [~,NormCurve(i,1:(length(ug2{i+1,3})+1),4)] = UGnormLike2b(ug2{i+1,3},ug2{i+1,4},param, val(i,1:2,4));
    [~,NormCurve(i,1:(length(ug2{i+1,3})+1),5)] = UGnormLike2d(ug2{i+1,3},ug2{i+1,4},param, val(i,1:3,5));    
    
    BIC(i,1) = Like(i,1) + 2/2*(log(length(ug2{i+1,3}))-log(2*pi)); %
    BIC(i,2) = Like(i,2) + 3/2*(log(length(ug2{i+1,3}))-log(2*pi)); %      
    BIC(i,3) = Like(i,3) + 4/2*(log(length(ug2{i+1,3}))-log(2*pi)); %      
    BIC(i,4) = Like(i,4) + 2/2*(log(length(ug2{i+1,3}))-log(2*pi)); %      
    BIC(i,5) = Like(i,5) + 3/2*(log(length(ug2{i+1,3}))-log(2*pi)); %    
end
    
    for i =1:4
       BICnorm(1,i) = mean(BIC(find(Lesion_Ind==(i-1)),1));
       BICstd(1,i) = std(BIC(find(Lesion_Ind==(i-1)),1));
       BICnorm(2,i) = mean(BIC(find(Lesion_Ind==(i-1)),2));
       BICstd(2,i) = std(BIC(find(Lesion_Ind==(i-1)),2));
       BICnorm(3,i) = mean(BIC(find(Lesion_Ind==(i-1)),3));
       BICstd(3,i) = std(BIC(find(Lesion_Ind==(i-1)),3));
       BICnorm(4,i) = mean(BIC(find(Lesion_Ind==(i-1)),4));
       BICstd(4,i) = std(BIC(find(Lesion_Ind==(i-1)),4));
       BICnorm(5,i) = mean(BIC(find(Lesion_Ind==(i-1)),5));
       BICstd(5,i) = std(BIC(find(Lesion_Ind==(i-1)),5));      
    end
    
     result(:,1,1) =  exp(val(:,1,1))./(1+exp(val(:,1,1)));
     result(:,2,1) =  exp(val(:,2,1))./(1+exp(val(:,2,1)));
     result(:,1,2) =  exp(val(:,1,2))./(1+exp(val(:,1,2)));
     result(:,2,2) =  exp(val(:,2,2))./(1+exp(val(:,2,2)));
     result(:,3,2) =  exp(val(:,3,2))./(1+exp(val(:,3,2)));
     result(:,1,3) =  exp(val(:,1,3))./(1+exp(val(:,1,3)));
     result(:,2,3) =  exp(val(:,2,3))./(1+exp(val(:,2,3)));
     result(:,3,3) =  exp(val(:,3,3))./(1+exp(val(:,3,3)));    
     result(:,4,3) =  20*exp(val(:,4,3))./(1+exp(val(:,4,3)));   
     result(:,1,4) =  exp(val(:,1,4))./(1+exp(val(:,1,4)));
     result(:,2,4) =  exp(val(:,2,4))./(1+exp(val(:,2,4)));     
     result(:,1,5) =  exp(val(:,1,5))./(1+exp(val(:,1,5)));
     result(:,2,5) =  exp(val(:,2,5))./(1+exp(val(:,2,5)));
     result(:,3,5) =  20*exp(val(:,3,5))./(1+exp(val(:,3,5)));