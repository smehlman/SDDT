% ====================== SDDT_pred_muN_pD.m ========================
% This code accompanies Ehlman et al. (Am Nat) to generate Figures 2,3, and
% 4. Small alterations to this code (by altering sigmas) can be used to
% generate Figure 5. 

% State-dependent detection theory is used to calculate decision 
% thresholds before the introduction of exotic predators. These thresholds
% are then used post-introduction and fitness change is assessed.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% -------------------------------------------------------------------
% ****************************  PARAMETERS  *************************
% -------------------------------------------------------------------
 
clear

muN = 0:0.1:2;        % mean of exotic cue distribution (also must change later)
pD = 0.01:0.01:0.3;   % probability of danger
pD_old = pD;          % used after exotic intro

pN = 0;               % probability of exotic pred (before intro, obvs 0)
pS = 1 - pN - pD;     % probability of safety

muS = 0;       % mean of safe cue distribution
muD = 2;       % mean of dangerous cue distribution 
sigmaS = 1;    % variance of safe cue distribution
sigmaD = 1;    % variance of dangerous cue distribution
sigmaN = 1;    % variance of exotic cue distribution

e = 0.9;        % prob of escape given that RUN when DANGER
m = 0.1;        % prob that pred misses when FORAGE when DANGER
eN = 0.9;       % prob of escape given that RUN when NOVEL
mN = 0.1;       % prob that pred misses when FORAGE when NOVEL


horizon = 100;  % Number of time steps from start to end

reserves = 30;  % Number of reserve levels 

k = 15;         % terminal fitness function - x-value of inflection point
phi = 0.2;      % terminal fitness function - slope at inflection point

pre_rv_pD = [];
post_rv_pD = [];
probmat_pre = [];
new_prob = [];
original_prob = [];   
    
bestsofar = 0;

lambda = 2;           % mean food obtained when foraging
max_food_gain = 10;   % max food gained (negligible prob with small lambda)

food_vector = food_from_foraging_horizon(lambda, max_food_gain);  
% returns vector of probabilities.  
% Note that food_vector(1) is prob of obtaining 0 food  

range = -10:0.001:10;  % range of values over which threshold is optimized

threshold = [];        % originating threshold matrix
threshold(1:(reserves), horizon,numel(pD)) = 0;

% ---Originating 'error'/'misses' matrices--- 
error_s1 = zeros(reserves,horizon,numel(pD));
error_d1 = zeros(reserves,horizon,numel(pD));
error_s2 = zeros(reserves,horizon,numel(pD));
error_d2 = zeros(reserves,horizon,numel(pD));
error_n = zeros(reserves,horizon,numel(pD),numel(muN));

% ---Originating RV matrices---
for u = 1:numel(pD)
  for f = 1:numel(muN)  
    for i = 1:reserves
      expected_rv_pre(i, horizon + 1,u) = 1/( 1 + exp(-phi*(i - k)) );
      expected_rv_post(i, horizon + 1,u,f) = 1/( 1 + exp(-phi*(i - k)) );
      % 1/( 1 + exp(-phi*(i - k)) )
    end
    % initialise all the other payoff values to zero.
    for j = 1:horizon
      for i = 1:reserves
        expected_rv_pre(i,j,u) = 0;
        expected_rv_post(i,j,u,f) = 0;
      end
    end 
  end
end

for u = 1:numel(pD) % Loop through all pD values
    
muN = 2;            
    
pN = 0;
pS = 1-pN-pD(u);

lhs = [];

    % ---Left hand side of signal detection theory equality---
    for j = 1:numel(range)
       lhs(j) = ((pN * normpdf(range(j),muN,sigmaN) + ...
                pD(u) * normpdf(range(j),muD,sigmaD) ) / (pN + pD(u)))...
                / normpdf(range(j),muS,sigmaS);
    end

rhs = [];

    for j = 1:horizon
       timestamp = horizon + 1 - j;  %The time step that we're dealing with

       %  Loop for updating thresholds
       for i = 1:reserves
        
          if i == 1
             threshold(i,timestamp,u) = max(range);    % This applies due to losing one unit of energy each time step.  

          else
              total_val_foraging_horizon_u_pre; % Calculates v_forage_survive: expected RS (if foraging and survive), as probability of each amount of food multiplied by value of resulting state (given starting in state i).  
    
              v_forage_survive

              % ---right hand side of signal detection theory equality---
              rhs(i) = pS/(pD(u) + pN) * (v_forage_survive - expected_rv_pre(i-1, timestamp + 1,u))/...
                            ((pD(u) * (e*expected_rv_pre(i-1, timestamp + 1,u) - m* v_forage_survive )...
                            + pN * (eN*expected_rv_pre(i-1, timestamp + 1,u) - mN* v_forage_survive ))/...
                            (pD(u) + pN));

              tmp = abs(lhs(:)-rhs(i));        
              [qirrelevant n] = min(tmp);          
              threshold(i, timestamp,u) = range(n);            
              tmp = [];
          end % if i case

          error_d1(i,timestamp,u) = pD(u) * normcdf(threshold(i,timestamp,u),muD,sigmaD);
          error_s1(i,timestamp,u) = pS * 1-normcdf(threshold(i,timestamp,u),muS,sigmaS);
          
       end  % for i = 1: reserves
       
       
       % ---Calculating pre-exotic reproductive value---
       for i = 1:reserves

          total_val_foraging_horizon_u_pre; % calculates v_forage_survive: expected RS (if foraging and survive), as probability of each amount of food multiplied by value of resulting state (given starting in state i).  

         if i == 1

            v_run_D = 0;
            v_cont_D = pD(u) * m * v_forage_survive;
            v_run_S = 0;
            v_cont_S = pS * v_forage_survive;
            expected_rv_pre(i, timestamp,u) = v_run_D + v_cont_D + v_run_S + v_cont_S;

         else  % qqp does this deal ok with reserves already being at max?

            v_run_D = pD(u) * (1-normcdf(threshold(i, timestamp,u),muD,sigmaD)) * expected_rv_pre(i-1, timestamp+1,u) * e;
            v_cont_D = pD(u) * normcdf(threshold(i, timestamp,u),muD,sigmaD) * v_forage_survive * m;
            v_run_S = pS * (1-normcdf(threshold(i, timestamp,u),muS,sigmaS)) * expected_rv_pre(i-1, timestamp+1,u);
            v_cont_S = pS * normcdf(threshold(i, timestamp,u),muS,sigmaS) * v_forage_survive;
            
            expected_rv_pre(i, timestamp,u) = v_run_D + v_cont_D + v_run_S + v_cont_S;

         end  % if - else
       end  % for i = 1:reserves
    end   % for each time
    

%%%%%%%%%%%%%%%% AFTER EXOTIC INTRODUCTION %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pD_safe = 0;
pN = pD(u);         % Default is full replacement of native with exotic
pS = 1-pD_safe-pN;
    
        % ***********************************************************    
        %         CALC EXPECTED REPRO SUCCESS
        % *********************************************************** 
        
  muN = 0:0.1:2;      % mean of exotic cue distribution
    
    for f = 1:numel(muN)  
       timestamp = []; 
        
       for j = 1:horizon
            timestamp = horizon + 1 - j;  
            
        for i = 1:reserves
           
            total_val_foraging_horizon_u_post_f; % calculates v_forage_survive: expected RS (if foraging and survive), as probability of each amount of food multiplied by value of resulting state (given starting in state i).  
       
            error_d2(i,timestamp,u) = pD_safe * normcdf(threshold(i, timestamp,u),muD,sigmaD);
            error_s2(i,timestamp,u) = pS * 1-normcdf(threshold(i, timestamp,u),muS,sigmaS);
            error_n(i,timestamp,u,f) = normcdf(threshold(i, timestamp,u),muN(f),sigmaN); % 'Misses' conditional on exotic being present here
            
         if i == 1

            v_run_D = 0;
            v_cont_D = pD_safe * m * v_forage_survive;
            v_run_S = 0;
            v_cont_S = pS * v_forage_survive;
            v_run_N = 0;
            v_cont_N = pN * mN * v_forage_survive;
            expected_rv_post(i,timestamp,u,f) = v_run_D + v_cont_D + v_run_S + v_cont_S + ...
                            v_run_N + v_cont_N;

         else

            v_run_D = pD_safe * (1-normcdf(threshold(i, timestamp,u),muD,sigmaD))...
                            * expected_rv_post(i-1, timestamp+1,u,f) * e;
            v_cont_D = pD_safe * normcdf(threshold(i, timestamp,u),muD,sigmaD) * ...
                            v_forage_survive * m;
            v_run_S = pS * (1-normcdf(threshold(i, timestamp,u),muS,...
                            sigmaS)) * expected_rv_post(i-1, timestamp+1,u,f);
            v_cont_S = pS * normcdf(threshold(i, timestamp,u),muS,...
                            sigmaS) * v_forage_survive;
            v_run_N = pN * (1-normcdf(threshold(i, timestamp,u),muN(f),sigmaN))...
                            * expected_rv_post(i-1, timestamp+1,u,f) * eN;
            v_cont_N = pN * normcdf(threshold(i, timestamp,u),muN(f),sigmaN)...
                            * v_forage_survive * mN;

            expected_rv_post(i,timestamp,u,f) = v_run_D + v_cont_D + v_run_S + v_cont_S + v_run_N + v_cont_N;               
 
         end  % if - else
        end  % for i = 1:reserves
      end % for horizon
    end % for muN
end  % for pD  


%%%%%%%% PLOTTING RESULTS %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% PLOTTING ERROR RATE (Fig 4 in MS) %%%%

    figure(10)
    hold on
    title('reserves low, time early')
    xlabel('pD')
    ylabel('error rate')

    ylim([0 1])
    plot(pD,squeeze(error_n(2,2,:,1)),'--','LineWidth',2)
    plot(pD,squeeze(error_n(2,2,:,11)),'.-','LineWidth',2)
    plot(pD,squeeze(error_n(2,2,:,21)),'LineWidth',2)
    hold off

    figure(11)
    hold on
    title('reserves low, time mid')
    xlabel('pD')
    ylabel('error rate')

    ylim([0 1])
    plot(pD,squeeze(error_n(2,50,:,1)),'--','LineWidth',2)
    plot(pD,squeeze(error_n(2,50,:,11)),'.-','LineWidth',2)
    plot(pD,squeeze(error_n(2,50,:,21)),'LineWidth',2)
    hold off

    figure(12)
    hold on
    title('reserves low, time late')
    xlabel('pD')
    ylabel('error rate')

    ylim([0 1])
    plot(pD,squeeze(error_n(2,98,:,1)),'--','LineWidth',2)
    plot(pD,squeeze(error_n(2,98,:,11)),'.-','LineWidth',2)
    plot(pD,squeeze(error_n(2,98,:,21)),'LineWidth',2)
    hold off

    figure(13)
    hold on
    title('reserves mid, time early')
    xlabel('pD')
    ylabel('error rate')

    ylim([0 1])
    plot(pD,squeeze(error_n(15,2,:,1)),'--','LineWidth',2)
    plot(pD,squeeze(error_n(15,2,:,11)),'.-','LineWidth',2)
    plot(pD,squeeze(error_n(15,2,:,21)),'LineWidth',2)
    hold off

    figure(14)
    hold on
    title('reserves mid, time mid')
    xlabel('pD')
    ylabel('error rate')

    ylim([0 1])
    plot(pD,squeeze(error_n(15,50,:,1)),'--','LineWidth',2)
    plot(pD,squeeze(error_n(15,50,:,11)),'.-','LineWidth',2)
    plot(pD,squeeze(error_n(15,50,:,21)),'LineWidth',2)
    hold off

    figure(15)
    hold on
    title('reserves mid, time late')
    xlabel('pD')
    ylabel('error rate')

    ylim([0 1])
    plot(pD,squeeze(error_n(15,98,:,1)),'--','LineWidth',2)
    plot(pD,squeeze(error_n(15,98,:,11)),'.-','LineWidth',2)
    plot(pD,squeeze(error_n(15,98,:,21)),'LineWidth',2)
    hold off

    figure(16)
    hold on
    title('reserves high, time early')
    xlabel('pD')
    ylabel('error rate')

    ylim([0 1])
    plot(pD,squeeze(error_n(29,2,:,1)),'--','LineWidth',2)
    plot(pD,squeeze(error_n(29,2,:,11)),'.-','LineWidth',2)
    plot(pD,squeeze(error_n(29,2,:,21)),'LineWidth',2)
    hold off

    figure(17)
    hold on
    title('reserves high, time mid')
    xlabel('pD')
    ylabel('error rate')

    ylim([0 1])
    plot(pD,squeeze(error_n(29,50,:,1)),'--','LineWidth',2)
    plot(pD,squeeze(error_n(29,50,:,11)),'.-','LineWidth',2)
    plot(pD,squeeze(error_n(29,50,:,21)),'LineWidth',2)
    hold off

    figure(18)
    hold on
    title('reserves high, time late')
    xlabel('pD')
    ylabel('error rate')

    ylim([0 1])
    plot(pD,squeeze(error_n(29,98,:,1)),'--','LineWidth',2)
    plot(pD,squeeze(error_n(29,98,:,11)),'.-','LineWidth',2)
    plot(pD,squeeze(error_n(29,98,:,21)),'LineWidth',2)
    hold off


%     %%%%%% Plotting THRESHOLDS (Figure 3 in MS) %%%%%%
% 

    g = find(squeeze(threshold(2,1,:)) == min(squeeze(threshold(2,1,:))));
    gg = find(squeeze(threshold(5,1,:)) == min(squeeze(threshold(5,1,:))));
    ggg = find(squeeze(threshold(10,1,:)) == min(squeeze(threshold(10,1,:))));
    gggg = find(squeeze(threshold(15,1,:)) == min(squeeze(threshold(15,1,:))));
    ggggg = find(squeeze(threshold(20,1,:)) == min(squeeze(threshold(20,1,:))));
    gggggg = find(squeeze(threshold(29,1,:)) == min(squeeze(threshold(29,1,:))));

    figure(20)
    hold on
    title('time early')
    xlabel('pD')
    ylabel('threshold')
    %xlim([0 0.5])
    ylim([-1 3.5])
    plot(pD,squeeze(threshold(2,1,:)),'m','LineWidth',2)
    plot(pD(g(1)),min(squeeze(threshold(2,1,:))),'ko')
    plot(pD,squeeze(threshold(5,1,:)),'r','LineWidth',2)
    plot(pD(gg(1)),min(squeeze(threshold(5,1,:))),'ko')
    plot(pD,squeeze(threshold(10,1,:)),'y','LineWidth',2)
    plot(pD(ggg(1)),min(squeeze(threshold(10,1,:))),'ko')
    plot(pD,squeeze(threshold(15,1,:)),'g','LineWidth',2)
    plot(pD(gggg(1)),min(squeeze(threshold(15,1,:))),'ko')
    plot(pD,squeeze(threshold(20,1,:)),'c','LineWidth',2)
    plot(pD(ggggg(1)),min(squeeze(threshold(20,1,:))),'ko')
    plot(pD,squeeze(threshold(29,1,:)),'b','LineWidth',2)
    plot(pD(gggggg(1)),min(squeeze(threshold(29,1,:))),'ko')
    hold off

    g = find(squeeze(threshold(2,50,:)) == min(squeeze(threshold(2,50,:))));
    gg = find(squeeze(threshold(5,50,:)) == min(squeeze(threshold(5,50,:))));
    ggg = find(squeeze(threshold(10,50,:)) == min(squeeze(threshold(10,50,:))));
    gggg = find(squeeze(threshold(15,50,:)) == min(squeeze(threshold(15,50,:))));
    ggggg = find(squeeze(threshold(20,50,:)) == min(squeeze(threshold(20,50,:))));
    gggggg = find(squeeze(threshold(29,50,:)) == min(squeeze(threshold(29,50,:))));
    
    figure(21)
    hold on
    title('time mid')
    xlabel('pD')
    ylabel('threshold')
    %xlim([0 0.5])
    ylim([-1 3.5])
    plot(pD,squeeze(threshold(2,50,:)),'m','LineWidth',2)
    plot(pD(g(1)),min(squeeze(threshold(2,50,:))),'ko')
    plot(pD,squeeze(threshold(5,50,:)),'r','LineWidth',2)
    plot(pD(gg(1)),min(squeeze(threshold(5,50,:))),'ko')
    plot(pD,squeeze(threshold(10,50,:)),'y','LineWidth',2)
    plot(pD(ggg(1)),min(squeeze(threshold(10,50,:))),'ko')
    plot(pD,squeeze(threshold(15,50,:)),'g','LineWidth',2)
    plot(pD(gggg(1)),min(squeeze(threshold(15,50,:))),'ko')
    plot(pD,squeeze(threshold(20,50,:)),'c','LineWidth',2)
    plot(pD(ggggg(1)),min(squeeze(threshold(20,50,:))),'ko')
    plot(pD,squeeze(threshold(29,50,:)),'b','LineWidth',2)
    plot(pD(gggggg(1)),min(squeeze(threshold(29,50,:))),'ko')
    hold off

    g = find(squeeze(threshold(2,99,:)) == min(squeeze(threshold(2,99,:))));
    gg = find(squeeze(threshold(5,99,:)) == min(squeeze(threshold(5,99,:))));
    ggg = find(squeeze(threshold(10,99,:)) == min(squeeze(threshold(10,99,:))));
    gggg = find(squeeze(threshold(15,99,:)) == min(squeeze(threshold(15,99,:))));
    ggggg = find(squeeze(threshold(20,99,:)) == min(squeeze(threshold(20,99,:))));
    gggggg = find(squeeze(threshold(29,99,:)) == min(squeeze(threshold(29,99,:))));
    
    figure(22)
    hold on
    title('time late')
    xlabel('pD')
    ylabel('threshold')
    %xlim([0 0.5])
    ylim([-1 3.5])
    plot(pD,squeeze(threshold(2,99,:)),'m','LineWidth',2)
    plot(pD(g(1)),min(squeeze(threshold(2,99,:))),'ko')
    plot(pD,squeeze(threshold(5,99,:)),'r','LineWidth',2)
    plot(pD(gg(1)),min(squeeze(threshold(5,99,:))),'ko')
    plot(pD,squeeze(threshold(10,99,:)),'y','LineWidth',2)
    plot(pD(ggg(1)),min(squeeze(threshold(10,99,:))),'ko')
    plot(pD,squeeze(threshold(15,99,:)),'g','LineWidth',2)
    plot(pD(gggg(1)),min(squeeze(threshold(15,99,:))),'ko')
    plot(pD,squeeze(threshold(20,99,:)),'c','LineWidth',2)
    plot(pD(ggggg(1)),min(squeeze(threshold(20,99,:))),'ko')
    plot(pD,squeeze(threshold(29,99,:)),'b','LineWidth',2)
    plot(pD(gggggg(1)),min(squeeze(threshold(29,99,:))),'ko')
    hold off

    
    
%%%%%%%%%%%%%%%%%% ADDITIONAL CODE %%%%%%%%%%%%%%%%%    
    
% In order to run code, must establish the following functions separately
% in MATLAB working directory:

% % ****************************************************************
% % total_val_foraging_horizon_u_pre 
% % Sums up the vold values corresponding to reserves which would be reached if foraging without being killed.
% % ****************************************************************
% v_forage_survive = 0;
% for food_gain = 0:max_food_gain
%    resulting_state = food_gain + i - 1;
%    if (resulting_state > reserves) % Would exceed L 
%       resulting_state = reserves;
%    end
%    if (resulting_state > 0)  % if statement checks that we?re not looking up values for a resultant state of zero.
%       v_forage_survive = v_forage_survive + food_vector(food_gain+1)*expected_rv_pre(resulting_state, timestamp + 1,u);
%    end % check we?re not looking up values for a resultant state of zero.
% end % for food_gain

% % ****************************************************************
% % total_val_foraging_u_post_f 
% % ****************************************************************
% v_forage_survive = 0;
% for food_gain = 0:max_food_gain
%    resulting_state = food_gain + i - 1;
%    if (resulting_state > reserves) % Would exceed L 
%       resulting_state = reserves;
%    end
%    if (resulting_state > 0)  % if statement checks that we?re not looking up values for a resultant state of zero.
%       v_forage_survive = v_forage_survive + food_vector(food_gain+1)*expected_rv_post(resulting_state, timestamp + 1,u,f);
%    end % check we?re not looking up values for a resultant state of zero.
% end % for food_gain

% % ****************************************************************
% function food_vector = food_from_foraging(lambda, food_max)
% % ****************************************************************
% total = 0;
% for i = 0:food_max
%   food_vector(i+1) = poisspdf(i,lambda);
%   total = total + food_vector(i+1);
% end
% for i = 0:food_max
%   food_vector(i+1) = food_vector(i+1)/total; 
% end