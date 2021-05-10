
% cooperativity parameters
a_mon_coop = 100;
a_dim_coop = 20;
b_mon_coop = 100;
b_dim_coop = 40;

kp = 0.2;
konc = 0.1;
TimeMaxBind = 10;
tss_time =4;

nr_sim=40;
tfinal = 1600000;

WUS=(0.008:0.01:0.4080);

%mutant specific parameters
if strcmp(str_mutant,'WT')
    
    Kdarray=[0.9571;0.1855;0.3663;0.5652;1.249];
    Kd2array=[0.5;1;0.5;1;0.5].*(Kdarray);
    
    nr_sites=length(Kdarray);
    
    k_on=konc*ones(nr_sites,1);
    
    distance_matrix
    
elseif strcmp(str_mutant,'950M')
    
    Kdarray=[0.1855;0.3663;0.5652;1.249];
    Kd2array=[1;0.5;1;0.5].*(Kdarray);
    
    nr_sites=length(Kdarray);
    
    k_on=konc*ones(nr_sites,1);
    
    distance_matrix
    
elseif strcmp(str_mutant,'970M')
    Kdarray=[0.9571;0.3663;0.5652;1.249];
    Kd2array=[0.5;0.5;1;0.5].*(Kdarray);
    
    nr_sites=length(Kdarray);
    
    k_on=konc*ones(nr_sites,1);
    
    distance_matrix
    
elseif strcmp(str_mutant,'DM')
    Kdarray=[0.9571;0.5652;1.249];
    Kd2array=[0.5;1;0.5].*(Kdarray);
    
    nr_sites=length(Kdarray);
    
    k_on=konc*ones(nr_sites,1);
    
    distance_matrix
    
elseif strcmp(str_mutant,'1060i')
    
    Kdarray=[1.249];
    Kd2array=[0.5].*(Kdarray);
    
    nr_sites=length(Kdarray);
    
    k_on=konc*ones(nr_sites,1);
    
    distance_matrix
    
elseif strcmp(str_mutant,'970M4i')
    Kdarray=[0.05830];
    Kd2array=[1].*(Kdarray);
    
    nr_sites=length(Kdarray);
    
    k_on=konc*ones(nr_sites,1);
    
    distance_matrix
    
end





