function periods = backward_recursion_markov(variables,observed_price,last_price,start_inventory,dist,priceDist,solution_method, suboptimal_decision_CEC, suboptimal_decision_R1, suboptimal_decision_R2)

%------------------------<allocate values>---------------------------------
% allocate variabels from struct "variable"
[ch, cp, d_max, I_max, k, m, pi, P_pr_reg, t_max, prices, normal_mu, normal_sigma,poisson_lambda,nbin_r,nbin_p] = allocate(variables);
demand_dist=get_demand_dist(dist, d_max, variables); 

PI_markov=cell(1,t_max); % Variable that will contain the PI values for t=1, 2,...t_max
nP_markov=cell(1,t_max); % Variable that will contain the next prices for t=1, 2,...t_max
markovmodel=getMarkovModel(priceDist);
for t=1:t_max+1
    if t==1
        PI_markov{t}=update_pi_markov(markovmodel,prices,t,pi,k); % Estimated values for regime 1 or regime 2 in period 1 (given); here [0.5 0.5]
        nP_markov{t}=next_price_markov(markovmodel,prices,t,pi);
    else
       PI_markov{t}=update_pi_markov(markovmodel,prices,t,PI_markov,k,solution_method);
       nP_markov{t}=next_price_markov(markovmodel,prices,t,PI_markov);
    end
end
PI_markov=PI_markov(2:end); %PI_Markov{1} is Period 0 (p_t-1).
nP_markov=nP_markov(1:end-1);
%PI=reshapePIMarkov(PI_Markov, t_max, b_t);
PI_needed_markov=needed_pi_markov(PI_markov,observed_price,last_price,prices);
if isequal(PI_needed_markov{1}, [0 0])
    disp('No reachable state')
    periods={};
    return
end
for t=1:t_max
    PI_needed{t}=PI_needed_markov{t}(:,1);
end

inv=0:I_max;

for t=t_max:-1:1  %starting in the last period
   if t==t_max % Prices for last period: No stages of the future must be considered.
        all_comb=combvec(prices,inv,PI_needed{t}')'; % only states that can be reached are under examination
        period_matrix=zeros(length(all_comb),I_max+1); % Matrix containing all costs for all possible stages and every demand ordering
        period_min_cost=zeros(length(all_comb),1); % vector containgin the minimum costs of every stage
        period_min_order=zeros(length(all_comb),1); % vector containing the order decission with minimum costs
        
        for s=1:length(all_comb) %loop over all reachable states s
            for y=0:I_max %loop over all order decision 
                period_matrix(s,y+1)=single_period_cost(y,cp,ch,all_comb(s,1),d_max,all_comb(s,2),demand_dist,dist); %allComb(s,1): actual price p_t, allComb(s,2): actual inventory I_t
            end
            period_min_cost(s)=min(period_matrix(s,:)); % minimal costs of state s
            period_min_order(s)=find(period_matrix(s,:)==min(period_matrix(s,:)))-1; %optimal order decision in state s
        end
   else % Not in last period
       if t~=1
           all_comb=combvec(prices,inv,PI_needed{t}')'; %only states that can be reached are under examination -> PI_needed
       else
           all_comb=combvec(observed_price,start_inventory,PI_needed{t}')'; %in period 1 only the case for p_1=observed price and I_1=startInventory is under examination
       end
       all_comb_size=size(all_comb);
       period_matrix=zeros(all_comb_size(1),I_max+1);
       period_min_cost=zeros(all_comb_size(1),1); 
       period_min_order=zeros(all_comb_size(1),1); % vector containing the order decision with minimum costs
       next_prices=get_next_prices_markov(nP_markov{t+1}, all_comb(:,1),observed_price,last_price,prices, t); 
       if t~=1 
           next_prices=reorderMarkov(next_prices,prices,I_max,t); %adjust the price vector to each possible state in t+1
       end
       for s=1:all_comb_size(1)
           mins=get_mins(PI_needed,prices, t, all_comb(s,3), all_comb(s,1), periods{t+1}.period_min_cost, inv); %get right minimal costs from t+1
            for y=0:I_max
                period_matrix(s,y+1)=backward_cost_function(y,cp,ch,all_comb(s,1),d_max,I_max,all_comb(s,2),next_prices(s,:),mins,demand_dist,dist);
            end
            period_min_cost(s)=min(period_matrix(s,:)); % minimal costs of state s
            period_min_order(s)=find(period_matrix(s,:)==min(period_matrix(s,:)))-1; %optimal order decision in state s
       end
   end
   
   % Assign calculated values to the variable periods
   periods{t}.PI=PI_markov{t};
   periods{t}.PI_needed=PI_needed{t};
   periods{t}.period_matrix=period_matrix;
   
   periods{t}.suboptimal_cost_CEC=period_matrix(s,suboptimal_decision_CEC+1);
   periods{t}.suboptimal_cost_R1=period_matrix(s,suboptimal_decision_R1+1);
   periods{t}.suboptimal_cost_R2=period_matrix(s,suboptimal_decision_R2+1);
   
   periods{t}.period_min_cost=period_min_cost;
   periods{t}.period_min_order=period_min_order;
   periods{t}.nextPrice=nP_markov{t};
end
