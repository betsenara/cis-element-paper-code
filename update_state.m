function [STATE,num_mon,num_dim,mon_state_vec,dim_state_vec] = update_state(STATE,EVENTS,SITES,result,j)

 if EVENTS(j)~=0
     
     STATE(SITES(j)) = result(EVENTS(j));
            
 end
        % number of monomers, number of dimers
           num_mon=0;
           num_dim=0;
           nrSITE=length(STATE);
            mon_state_vec=zeros(size(STATE));
            dim_state_vec=zeros(size(STATE));
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
        
end

