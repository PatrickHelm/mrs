
clear all
close all hidden
format long
clc
addpath('src');

%_ _SOLUTION METHOD _ _ _ _ _ _ _  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
solution_method='R2';              %exact, CEC, R1, R2
%_ _SUBOPTIMAL DECISIONS _ _ _ _ _ _ _  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
suboptimal_decision_vector_R1=xlsread('Stochastic_Regimes_MR-MO.xlsx',1,'B4:B11');
suboptimal_decision_vector_R2=xlsread('Stochastic_Regimes_MR-MO.xlsx',1,'C4:C11');
suboptimal_decision_matrix_CEC=xlsread('Stochastic_Regimes_MR-MO.xlsx',1,'F4:M11');
%_ _INVENTORY COST _ _ _ _ _ _ _  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
variables.t_max=4;                       
variables.ch=17.5;                            
variables.cp=6000;  
%_ _DEMAND _ _ _ _ _ _ _  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
variables.d_max=1;                        
variables.I_max=4;                       
demand_dist='manual';                 %uniform, normal, poisson, nbin, manual
    variables.normal.mu=15;
    variables.normal.sigma=3;
    variables.poisson.lambda=1;
    variables.nbin.r=22.5;
    variables.nbin.p=0.6;
start_inventory=0;    
%_ _PRICE PROCESS _ _ _ _ _ _ _  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
variables.k=[0.99 0.01; 0.01 0.99];         
variables.m=2;                              
beliefs=[0.5];
variables.prices=[1000 1500 2000 2500 3000 3500 4000 4500]; 
price_dist={'RW','MR'};
variables.P_pr_reg=0; %muss nur definiert sein, da unter variables, aber für Markov-Fall nicht von Bedeutung

%_ _CALCULATIONS START _ _ _ _ _ _ _  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _    

diary myFirstDiaryFile.txt

for i=1:numel(variables.prices)
    observed_price=variables.prices(i);
    for h=1:1%numel(variables.prices)
        last_price=variables.prices(h);
        for j=1:1%numel(beliefs)
            pi1=beliefs(j);
            variables.pi=[pi1 1-pi1];
            
            suboptimal_decision_CEC=suboptimal_decision_matrix_CEC(j,i);
            suboptimal_decision_R1=suboptimal_decision_vector_R1(i);
            suboptimal_decision_R2=suboptimal_decision_vector_R2(i);


            allSave=cell(variables.t_max, 1);           



            for T=variables.t_max 
                variables.t_max=T;
            tic
                periods=backward_recursion_markov(variables,observed_price,last_price,start_inventory,demand_dist,price_dist,solution_method,suboptimal_decision_CEC,suboptimal_decision_R1,suboptimal_decision_R2);
                if ~isempty(periods)
                disp(['<------------Planning horizon : ' num2str(T) ' periods --------------------->']);
                disp(['- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -']);
                disp(['For planning horizon of ' num2str(T) ' periods, and observed price sequence: ' num2str(last_price) '-->' num2str(observed_price) ' and initial belief ' num2str(variables.pi) ':']);
                disp(['Optimal Order Decision: ' num2str(periods{1}.period_min_order)]);
                disp(['Minimum Total Expected Cost: ' num2str(periods{1}.period_min_cost)]);
                disp(['- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -']);
                disp(['Suboptimal Cost (order decision (CEC): ' num2str(suboptimal_decision_CEC) '): ' num2str(periods{1}.suboptimal_cost_CEC)]);
                disp(['Optimality gap in % total cost (CEC): ' num2str((periods{1}.suboptimal_cost_CEC/periods{1}.period_min_cost)*100-100)]);
                disp(['- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -']);
                disp(['Suboptimal Cost (order decision (R1): ' num2str(suboptimal_decision_R1) '): ' num2str(periods{1}.suboptimal_cost_R1)]);
                disp(['Optimality gap in % total cost (R1): ' num2str((periods{1}.suboptimal_cost_R1/periods{1}.period_min_cost)*100-100)]);
                disp(['- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -']);
                disp(['Suboptimal Cost (order decision (R2): ' num2str(suboptimal_decision_R2) '): ' num2str(periods{1}.suboptimal_cost_R2)]);
                disp(['Optimality gap in % total cost (R2): ' num2str((periods{1}.suboptimal_cost_R2/periods{1}.period_min_cost)*100-100)]);
                end 
            toc
                allSave{T}={periods};
            end
        end 
        j+1;
    end 
    h+1;
end 
i+1;
diary off
