function periods = backward_recursion(variables, observed_price, start_inventory,dist,solution_method, suboptimal_decision_CEC, suboptimal_decision_R1, suboptimal_decision_R2)

[ch, cp, d_max, I_max, k, m, pi, P_pr_reg, t_max, prices, normal_mu, normal_sigma,poisson_lambda,nbin_r,nbin_p] = allocate(variables);
demand_dist=get_demand_dist(dist, d_max, variables); 

PI=cell(1,t_max); %variable that will contain the belief values for t=1, 2,...t_max
nP=cell(1,t_max); %variable that will contain the next prices for t=1, 2,...t_max
for t=1:t_max
    if t==1
        PI{t}=pi; %initial beliefs for regime 1 and 2 in period 1 (given) 
    else
        if strcmp(solution_method, 'OLFC')
        PI{t}=update_pi(k,m,pi,P_pr_reg,t,solution_method);
        else 
        PI{t}=update_pi(k,m,PI{t-1},P_pr_reg,t,solution_method);
        end 
    end
    nP{t}=next_price(PI{t}, P_pr_reg); %calculate probabilities for all possible prices in period t
end

PI_needed=needed_pi(PI, observed_price, prices); %only consider beliefs, that can be reached, given the observed price in period 1

inv=0:I_max; 

for t=t_max:-1:1 %starting in the last period
   if t==t_max 
        all_comb=combvec(prices,inv,PI_needed{t}')'; %only states that can be reached are under examination
        period_matrix=zeros(length(all_comb),I_max+1); %matrix containing all costs for all possible periods and every demand ordering
        period_min_cost=zeros(length(all_comb),1); %vector containing the minimum costs of every period
        period_min_order=zeros(length(all_comb),1); %vector containing the order decision of every period with minimum costs
        for s=1:length(all_comb) %loop over all reachable states s
            for y=0:I_max %loop over all order decision
                period_matrix(s,y+1)=single_period_cost(y,cp,ch,all_comb(s,1),d_max,all_comb(s,2),demand_dist,dist); %allComb(s,1): actual price p_t, allComb(s,2): actual inventory I_t
            end
            period_min_cost(s)=min(period_matrix(s,:)); %minimal costs of state s
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
       period_min_order=zeros(all_comb_size(1),1); 
       next_prices=get_next_prices(nP{t+1}, all_comb(:,1), observed_price,prices, t);  %price probabilities where p_1 = observedPrice
       if t~=1 
           next_prices=aux_matrix(next_prices,prices,I_max,t); %adjust the price vector to each possible state in t+1
       end
       for s=1:all_comb_size(1)
           mins=get_mins(PI_needed,prices, t, all_comb(s,3), all_comb(s,1), periods{t+1}.period_min_cost, inv); %get right minimal costs from t+1
            for y=0:I_max
                period_matrix(s,y+1)=backward_cost_function(y,cp,ch,all_comb(s,1),d_max,I_max,all_comb(s,2),next_prices(s,:),mins,demand_dist,dist);
            end
            period_min_cost(s)=min(period_matrix(s,:)); %minimal costs of state s
            period_min_order(s)=find(period_matrix(s,:)==min(period_matrix(s,:)))-1; %optimal order decision in state s
       end
   end
   periods{t}.PI=PI{t};
   periods{t}.PI_needed=PI_needed{t};
   periods{t}.period_matrix=period_matrix;
   periods{t}.suboptimal_cost_CEC=period_matrix(s,suboptimal_decision_CEC+1);
   periods{t}.suboptimal_cost_R1=period_matrix(s,suboptimal_decision_R1+1);
   periods{t}.suboptimal_cost_R2=period_matrix(s,suboptimal_decision_R2+1);
   periods{t}.period_min_cost=period_min_cost;
   periods{t}.period_min_order=period_min_order;
   periods{t}.next_price=nP{t};
end