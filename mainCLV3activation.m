%clear all

str_mutant = 'WT';
%str_mutant = '950M';
%str_mutant = '970M';
%str_mutant = 'DM';
%str_mutant = '970Mi';
%str_mutant = '1060i';

parameters


mRNA_nr_sim = cell(length(nr_sim),1);


for ii=1:nr_sim   
    mRNAsave = zeros(1,length(WUS));
    
    
    for wi=1:length(WUS)
        
        %initialize
        time_curr= 0;
        copyfcn = @(N,nr) [repmat(N,[1 nr]) ];
        STATE = copyfcn('N',nr_sites);
        result = STATE;
        EVENTS = [1];
        SITES=[1];
        j=1;
        tau=0;
        mon_state_vec=zeros(size(STATE));
        dim_state_vec=zeros(size(STATE));
        
        % Transcription Start Site (TSS), empty (IN_SITE =1) or full (IN_SITE =0)
        IN_SITE =1;
        mRNA = 0;
        % find stochastic time step and next event
        [j,tau,EVENTS,SITES] = stoc_tau_j(STATE,Kdarray,Kd2array,k_on,WUS(wi),Int_Mat_Dimer,Int_Mat_Mon,a_mon_coop,a_dim_coop,b_mon_coop,b_dim_coop,kp,IN_SITE);
        
        time_remain_mon_vec= zeros(length(STATE),1);
        % the vector of total cooperativity for each site
        coop_number = (Int_Mat_Mon*mon_state_vec').*mon_state_vec';
        % the indices of monomers with no cooperativity
        ind_mon= find(((coop_number==0)&(mon_state_vec'>0))>0);
        % finding the minimum of the remaining times of monomers with no cooperativity
        [time_remain_mon,ind_min] = min(time_remain_mon_vec(ind_mon));
        % find the index of the monomer corresponding to minimum time
        ind_min_mon=ind_mon(ind_min);
        % if there are no monomers with no coop, set time_remain_mon to 0.
        if isempty(time_remain_mon)
            time_remain_mon = ~isempty(time_remain_mon);
        end
        
        while time_curr+tau<tfinal
            
            if ((sum((coop_number==0).*mon_state_vec')>0)&&(time_remain_mon<tau))&&(IN_SITE==1 ||(IN_SITE==0 && time_remain_mon<time_remain_pol))
                time_curr = time_curr+time_remain_mon;
                STATE(ind_min_mon)='N';
                
                
                time_remain_mon_vec= time_remain_mon_vec-time_remain_mon*mon_state_vec'.*(coop_number==0);
                mon_state_vec(ind_min_mon)=0;
                coop_number = (Int_Mat_Mon*mon_state_vec').*mon_state_vec';
                ind_mon= find(((coop_number==0)&(mon_state_vec'>0))>0);
                [time_remain_mon,ind_min] = min(time_remain_mon_vec(ind_mon));
                ind_min_mon=ind_mon(ind_min);
                if isempty(time_remain_mon)
                    time_remain_mon = ~isempty(time_remain_mon);
                end
                
                if IN_SITE==0
                    time_remain_pol = time_remain_pol-time_remain_mon;
                end
                
                
                
            elseif (IN_SITE==0 && time_remain_pol<tau)&&((sum((coop_number==0).*mon_state_vec')==0)||((sum((coop_number==0).*mon_state_vec')>0)&&(time_remain_pol<time_remain_mon)))
                
                time_curr = time_curr+time_remain_pol;
                mRNA = mRNA+1;
                IN_SITE=1;
                
                
                if isempty(time_remain_mon)
                    time_remain_mon = ~isempty(time_remain_mon);
                end
                
                if sum((coop_number==0).*mon_state_vec)>0
                    time_remain_mon_vec= time_remain_mon_vec-time_remain_pol*mon_state_vec'.*(coop_number==0);
                    coop_number = (Int_Mat_Mon*mon_state_vec').*mon_state_vec';
                    ind_mon= find(((coop_number==0)&(mon_state_vec'>0))>0);
                    [time_remain_mon,ind_min] = min(time_remain_mon_vec(ind_mon));
                    ind_min_mon=ind_mon(ind_min);
                end
                
            else
                time_curr = time_curr + tau;
                if sum(mon_state_vec)>0
                    time_remain_mon_vec= time_remain_mon_vec-tau*mon_state_vec'.*(coop_number==0);
                    coop_number = (Int_Mat_Mon*mon_state_vec').*mon_state_vec';
                    ind_mon= find(((coop_number==0)&(mon_state_vec'>0))>0);
                    [time_remain_mon,ind_min] = min(time_remain_mon_vec(ind_mon));
                    ind_min_mon=ind_mon(ind_min);
                    
                end
                
                if IN_SITE==0
                    time_remain_pol = time_remain_pol-tau;
                end
                
                result = 'MDNM';
                [STATE,num_mon,num_dim,mon_state_vec,dim_state_vec] = update_state(STATE,EVENTS,SITES,result,j);
                if (EVENTS(j)==1)||(EVENTS(j)==4)
                    % if a new monomer binds
                    time_remain_mon_vec(SITES(j))=TimeMaxBind;
                    coop_number = (Int_Mat_Mon*mon_state_vec').*mon_state_vec';
                    ind_mon= find(((coop_number==0)&(mon_state_vec'>0))>0);
                    [time_remain_mon,ind_min] = min(time_remain_mon_vec(ind_mon));
                    ind_min_mon=ind_mon(ind_min);
                end
                if EVENTS(j)==0
                    IN_SITE=0;
                    time_remain_pol=tss_time;
                end

            end
             % find stochastic time step and next event
            [j,tau,EVENTS,SITES] = stoc_tau_j(STATE,Kdarray,Kd2array,k_on,WUS(wi),Int_Mat_Dimer,Int_Mat_Mon,a_mon_coop,a_dim_coop,b_mon_coop,b_dim_coop,kp,IN_SITE);
            
        end
        %%%%%%%%%%% end of while
        mRNAsave(1,wi) = mRNA;
        
    end
    
    
    mRNA_nr_sim{ii} =  mRNAsave;
    
end
summRNA = zeros(1,length(WUS));

for ii=1:nr_sim
    mRNA_save = mRNA_nr_sim{ii};
    summRNA = summRNA + mRNA_save(1,:);
    
end



str = [str_mutant,'_pm',num2str(TimeMaxBind),'kp',num2str(kp),'tr_t',num2str(tss_time),'amc',num2str(a_mon_coop),'adc',num2str(a_dim_coop),'bmc',num2str(b_mon_coop),'bdc',num2str(b_dim_coop),'.mat'];

save(str)


