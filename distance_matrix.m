
if strcmp(str_mutant,'WT')
    
    Int_Mat_Dimer=zeros(length(Kdarray));
    Int_Mat_Dimer(1,2)=50;
    Int_Mat_Dimer(2,1)=50;
    Int_Mat_Dimer(3,2)=32;
    Int_Mat_Dimer(2,3)=32;
    Int_Mat_Dimer(3,4)=18;
    Int_Mat_Dimer(4,3)=18;
    Int_Mat_Dimer(5,4)=84;
    Int_Mat_Dimer(4,5)=84;
    Int_Mat_Dimer(1,3)=82;
    Int_Mat_Dimer(3,1)=82;
    Int_Mat_Dimer(1,4)=100;
    Int_Mat_Dimer(4,1)=100;
    Int_Mat_Dimer(1,5)=184;
    Int_Mat_Dimer(5,1)=184;
    Int_Mat_Dimer(4,2)=50;
    Int_Mat_Dimer(2,4)=50;
    Int_Mat_Dimer(5,2)=134;
    Int_Mat_Dimer(2,5)=134;
    Int_Mat_Dimer(5,3)=102;
    Int_Mat_Dimer(3,5)=102;
    
    Int_Mat_Mon=zeros(length(Kdarray));
    Int_Mat_Mon(2,1)=50;
    Int_Mat_Mon(1,2)=50;
    Int_Mat_Mon(2,3)=32;
    Int_Mat_Mon(3,2)=32;
    Int_Mat_Mon(4,2)=50;
    Int_Mat_Mon(2,4)=50;
    Int_Mat_Mon(5,2)=134;
    Int_Mat_Mon(2,5)=134;
    
elseif strcmp(str_mutant,'950M')
    
    Int_Mat_Dimer=zeros(length(Kdarray));
    Int_Mat_Dimer(1,2)=32;
    Int_Mat_Dimer(2,1)=32;
    Int_Mat_Dimer(3,2)=18;
    Int_Mat_Dimer(2,3)=18;
    Int_Mat_Dimer(3,4)=84;
    Int_Mat_Dimer(4,3)=84;
    
    Int_Mat_Dimer(1,3)=50;
    Int_Mat_Dimer(3,1)=50;
    Int_Mat_Dimer(1,4)=134;
    Int_Mat_Dimer(4,1)=134;
    
    Int_Mat_Dimer(4,2)=102;
    Int_Mat_Dimer(2,4)=102;
    
    Int_Mat_Mon=zeros(length(Kdarray));
    
    Int_Mat_Mon(2,1)=32;
    Int_Mat_Mon(1,2)=32;
    Int_Mat_Mon(1,3)=50;
    Int_Mat_Mon(3,1)=50;
    Int_Mat_Mon(1,4)=134;
    Int_Mat_Mon(4,1)=134;
    
elseif strcmp(str_mutant,'970M')
    Int_Mat_Dimer=zeros(length(Kdarray));
    Int_Mat_Dimer(1,2)=82;
    Int_Mat_Dimer(2,1)=82;
    Int_Mat_Dimer(3,2)=18;
    Int_Mat_Dimer(2,3)=18;
    Int_Mat_Dimer(3,4)=84;
    Int_Mat_Dimer(4,3)=84;
    
    Int_Mat_Dimer(1,3)=100;
    Int_Mat_Dimer(3,1)=100;
    Int_Mat_Dimer(1,4)=184;
    Int_Mat_Dimer(4,1)=184;
    
    Int_Mat_Dimer(4,2)=102;
    Int_Mat_Dimer(2,4)=102;
    
    Int_Mat_Mon=zeros(length(Kdarray));
    
elseif strcmp(str_mutant,'DM')
    Int_Mat_Dimer=zeros(length(Kdarray));
    Int_Mat_Dimer(1,2)=100;
    Int_Mat_Dimer(2,1)=100;
    Int_Mat_Dimer(3,2)=84;
    Int_Mat_Dimer(2,3)=84;   
    Int_Mat_Dimer(1,3)=184;
    Int_Mat_Dimer(3,1)=184;
 
    Int_Mat_Mon=zeros(length(Kdarray));
    
elseif strcmp(str_mutant,'1060i')
    Int_Mat_Dimer=zeros(length(Kdarray));
    Int_Mat_Mon=zeros(length(Kdarray));
elseif strcmp(str_mutant,'970M4i')
    Int_Mat_Dimer=zeros(length(Kdarray));
    Int_Mat_Mon=zeros(length(Kdarray));
    
end