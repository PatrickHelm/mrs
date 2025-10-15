clear all
close all hidden
format long
clc
addpath('src');

%_ _SOLUTION METHOD _ _ _ _ _ _ _  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
solution_method='R1';              %exact, CEC, R1, R2
%_ _SUBOPTIMAL DECISIONS _ _ _ _ _ _ _  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
suboptimal_decision_vector_R1=xlsread('Stochastic_Regimes_HL-LL.xlsx',1,'B4:B12');
suboptimal_decision_vector_R2=xlsread('Stochastic_Regimes_HL-LL.xlsx',1,'C4:C12');
suboptimal_decision_matrix_CEC=xlsread('Stochastic_Regimes_HL-LL.xlsx',1,'F4:T12');
%_ _INVENTORY COST _ _ _ _ _ _ _  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
variables.t_max=4;                         
variables.ch=15;                            
variables.cp=40;  
%_ _DEMAND _ _ _ _ _ _ _  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
variables.d_max=30;                        
variables.I_max=120;                        
demand_dist='normal';                 %uniform, normal, poisson, nbin, manual
    variables.normal.mu=15;
    variables.normal.sigma=3;
    variables.poisson.lambda=1;
    variables.nbin.r=22.5;
    variables.nbin.p=0.6;
start_inventory=0;    
%_ _PRICE PROCESS _ _ _ _ _ _ _  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ 
variables.k=[0.9 0.1; 0.1 0.9];         
variables.m=2;  
beliefs=[0.1 0.3 0.5 0.7 0.9];
%variables.prices=[100 150 200 250 300 350 400 450 500 550 600 650 700 750 800];
variables.prices=[10 12.5 15 17.5 20 22.5 25 27.5 30];
        %Real  World:
    %variables.P_pr_reg=[0.00000000000642 0.00000001522976; 0.00000000043814 0.00000486876449; 0.00000002003070 0.00047977596157; 0.00000061330835 0.01457311158330; 0.00001257653612 0.13644580403084; 0.00017272002688 0.39378771955745; 0.00158863395670 0.35031460248286; 0.00978598341421 0.09606128727711; 0.04037240294844 0.00811956293533; 0.11154862280620 0.00021154900723; 0.20641598015534 0.00000169896104; 0.25581258672549 0.00000000420581; 0.21232444231738 0.00000000000321; 0.11802596626598 0.000000000000001; 0.04393945106367 0.00000000000000000005];
        %HL-LL:
    variables.P_pr_reg=[0.0000012590141 0.0842412539718; 0.0000573844503 0.2387343482365; 0.0013060640686 0.3378404161123; 0.0148436719380 0.2387343482365; 0.0842412539718 0.0842412539718; 0.2387343482365 0.0148436719380; 0.3378404161123 0.0013060640686; 0.2387343482365 0.0000573844503; 0.0842412539718 0.0000012590141]; 
        %LV-HV:
    %variables.P_pr_reg=[0.0000018583873 0.1111111111111; 0.0004407417288 0.1111111111111; 0.0219102327357 0.1111111111111; 0.2283098678790 0.1111111111111; 0.4986745985385 0.1111111111111; 0.2283098678790 0.1111111111111; 0.0219102327357 0.1111111111111; 0.0004407417288 0.1111111111111; 0.0000018583873 0.1111111111111];
        %HL-LL (larger discretization):
    %variables.P_pr_reg=[0.000002480179 0.165950028484; 0.002572864946 0.665524597907; 0.165950028484 0.165950028484; 0.665524597907 0.002572864946; 0.165950028484 0.000002480179];

%_ _CALCULATIONS START _ _ _ _ _ _ _  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _    

diary myFirstDiaryFile.txt

for i=1:numel(variables.prices)
    observed_price=variables.prices(i);
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
                periods=backward_recursion(variables, observed_price, start_inventory, demand_dist,solution_method, suboptimal_decision_CEC, suboptimal_decision_R1, suboptimal_decision_R2);
                disp(['<------------Planning horizon : ' num2str(T) ' periods --------------------->']);
                disp(['- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -']);
                disp(['For planning horizon of ' num2str(T) ' periods and observed price: ' num2str(observed_price) ' and initial belief ' num2str(variables.pi) ' in period 1:']);
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
               
            toc
                allSave{T}={periods};
            end
    end 
    j+1;
end 
i+1;
diary off