
% 970M4i
str_mutant = '970M4i';

parameters
strparam = ['_pm',num2str(TimeMaxBind),'kp',num2str(kp),'tr_t',num2str(tss_time),'amc',num2str(a_mon_coop),'adc',num2str(a_dim_coop),'bmc',num2str(b_mon_coop),'bdc',num2str(b_dim_coop),'.mat'];
% fast dynamics
%strparam = ['fast_pm',num2str(TimeMaxBind),'kp',num2str(kp),'tr_t',num2str(tss_time),'amc',num2str(a_mon_coop),'adc',num2str(a_dim_coop),'bmc',num2str(b_mon_coop),'bdc',num2str(b_dim_coop),'.mat'];

str = [str_mutant,strparam];
load(str)
M970M4imRNA = (1/nr_sim)*summRNA;
plot(WUS,M970M4imRNA,'Color',[1 0.5 0.5],'LineWidth',3)
hold on

% 1060i
str_mutant = '1060i';
str = [str_mutant,strparam];
load(str)
M1060imRNA = (1/nr_sim)*summRNA;
plot(WUS,M970mRNA,'Color',[1 1 0.5],'LineWidth',3)
hold on




legend({'970M4i','1060i'},'FontSize',10)

ylabel('mRNA')
xlabel('WUS')
set(gca,'FontSize',20)
xlim([0 0.2])
