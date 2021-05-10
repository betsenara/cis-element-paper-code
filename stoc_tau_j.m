function [j,tau,EVENTS,SITES] = stoc_tau_j(STATE,Kdarray,Kd2array,k_on,wus,Int_Mat_Dimer,Int_Mat_Mon,a_mon_coop,a_dim_coop,b_mon_coop,b_dim_coop,kp,IN_SITE)

 
        [SITES,PROBS,EVENTS] = stochastic_events_probs(STATE,Kdarray,Kd2array,k_on,wus,Int_Mat_Dimer,Int_Mat_Mon,a_mon_coop,a_dim_coop,b_mon_coop,b_dim_coop,kp,IN_SITE);
       
        %% find tau, stochastic time step:
        
        r = rand(2,1);
        r1 = r(1);
        r2 = r(2);
        
        a0 = sum(PROBS);
        tau = (1/a0)*log(1/r1);
        
        %Find j, the event that happens
        
        j = 1;
        
        while sum(PROBS(1:j))<=r2*a0 
            j = j+1;
        end
end

