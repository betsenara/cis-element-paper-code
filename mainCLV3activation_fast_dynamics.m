%clear all

str_mutant = '970M4i';
%str_mutant = '1060i';

parameters_fast_dynamics

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
        
        
       while time_curr+tau<tfinal
                    
                    if IN_SITE==0 && time_remain_pol<tau
                        
                        time_curr = time_curr+time_remain_pol;
                        mRNA = mRNA+1;
                        IN_SITE=1;
                    
                    else
                        time_curr = time_curr + tau;
                        
                        if IN_SITE==0
                            time_remain_pol = time_remain_pol-tau;
                        end
                        
                        result = 'MDNM';
                        [STATE,num_mon,num_dim,mon_state_vec,dim_state_vec] = update_state(STATE,EVENTS,SITES,result,j);

                        if EVENTS(j)==0
                            IN_SITE=0;
                            time_remain_pol=tss_time;
                        end       
                    end
                    
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



str = [str_mutant,'fast_pm',num2str(TimeMaxBind),'kp',num2str(kp),'tr_t',num2str(tss_time),'amc',num2str(a_mon_coop),'adc',num2str(a_dim_coop),'bmc',num2str(b_mon_coop),'bdc',num2str(b_dim_coop),'.mat'];

save(str)


