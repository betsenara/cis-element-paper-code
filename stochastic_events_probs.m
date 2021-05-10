function [SITES,PROBS,EVENTS] = stochastic_events_probs(STATE,Kdarray,Kd2array,k_on,wus,Int_Mat_Dimer,Int_Mat_Mon,a_mon_coop,a_dim_coop,b_mon_coop,b_dim_coop,kp,IN_SITE)


SITES = [];
PROBS =[];
EVENTS = [];
nrSITE=length(STATE);

atanc = 0.3;

f_mon_coop = @(x) -a_mon_coop.*atan(atanc*(x-b_mon_coop)) + a_mon_coop.*atan(10^5) + 1 ;

f_dim_coop = @(x)-a_dim_coop.*atan(atanc*(x-b_dim_coop)) + a_dim_coop.*atan(10^5) + 1 ;
  
    num_mon=0;
    num_dim=0;
    dim_state_vec=zeros(length(STATE),1);
     mon_state_vec=zeros(length(STATE),1);
    for si=1:nrSITE
        if STATE(si)=='M'
            num_mon=num_mon+1;
             mon_state_vec(si)=1;
        end
    end
    
    for si=1:length(STATE)
        if STATE(si)=='D'
            num_dim=num_dim+1;
            dim_state_vec(si)=1;
        end
    end

    
    for i=1:length(STATE)
     
        
        if STATE(i)=='N'
            Site1 = i;
            Site2 = [];
          
             Prob1 = k_on(i)*wus;
            Prob2 =[];
            
            Event1 = 1;
            Event2 =[];
        end
        if STATE(i)=='M'
            Site1 = i;
            Site2 = i;
       
            Prob1 = k_on(i)*wus;

           Coop_i = Int_Mat_Mon(:,i);
           ind_mon_coop = find((Coop_i.*mon_state_vec)>0);
           if isempty(ind_mon_coop)
             coop_prod=1;
           else
    
              coop_prod=prod(f_mon_coop(Coop_i(ind_mon_coop)));
           end

           Prob2 =k_on(i)*Kdarray(i)./coop_prod;
 
            
            Event1 = 2;
            Event2 = 3;
        end
      
        if STATE(i)=='D'
            Site1 = i;
            Site2 = [];

          Coop_i = Int_Mat_Dimer(:,i);
           ind_dim_coop = find((Coop_i.*dim_state_vec)>0);
           if isempty(ind_dim_coop)
             coop_prod=1;
           else
              coop_prod=prod(f_dim_coop(Coop_i(ind_dim_coop)));
           end

           Prob1 = k_on(i)*Kd2array(i)./coop_prod;
           Prob2 =[];
            
            Event1 = 4;
            Event2 =[];
        end
        SITES = [SITES,Site1,Site2];
        PROBS =[PROBS,Prob1,Prob2];
        EVENTS =[EVENTS,Event1,Event2];
        
     
    end

    SITES = [SITES,[1:length(STATE)]];
    PROBS =[PROBS,kp*IN_SITE*mon_state_vec'];
    EVENTS =[EVENTS,zeros(1,length(STATE))];
    
  
end
