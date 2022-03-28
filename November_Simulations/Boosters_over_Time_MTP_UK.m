function [V3,Phase]=Boosters_over_Time_MTP_UK(UpTake, V2, Date)

%%

Separation=170;   % six month post dose 2.

Delay=7;
Inf_to_Symp=5;

Rollout=1.0; SLOWER=1.0;

load Regional_PP.mat

% Amount available in the England each week for boosting.
Vwf=[2.3 1.9 1.9 1.9 1.8 1.4 1.0 1.0]*Rollout*1e6;
Vwf=Vwf*sum(Region_PP(2:11,:),'all')/sum(Region_PP(2:8,:),'all'); % rescale to UK.
Vwf_date=datenum(2021,11,1)+1-datenum(2020,1,1);


load Vacc_Data.mat

propBUT=UpTake;

if nargin>=2
    Vacc2=V2;
end

if ~exist('Vacc3')
    Vacc3=0*Vacc1;
end

Vacc3(:,Date:end,:)=[];  %start Boosters on 20th Sept.


% Vacc1=Vacc1*SLOWER;
% Vacc2=Vacc2*SLOWER;
% Vacc_Type1=Vacc_Type1*SLOWER;
% Vacc_Type2=Vacc_Type2*SLOWER;
% UpTake=UpTake*SLOWER;
% Rollout=Rollout*SLOWER;

if length(UpTake)==1
    UpTake=ones(1,21)*UpTake;
end

% currently the same in all regions and all ages.


a80=[17:21]; p80=Region_PP(1:11,a80)/sum(Region_PP(2:11,a80),'all');  %Note average for UK data !!
aW=[4:13]; pW=Region_PP(1:11,aW)/sum(Region_PP(2:11,aW),'all');  %Note average for UK workers data !!

T=[datenum(2021,1,4):datenum(2022,11,1)] +1-datenum(2020,1,1);
V3=zeros(11,max(T),21);

% Assume last value of Vwf doses a week if we don't have data
T=[datenum(2021,1,25):datenum(2023,7,1)] +1-datenum(2020,1,1);
v(T)=Vwf(end)/7;

% Use SPI-M's forecasts....
T=Vwf_date + ([1:length(Vwf)]-1)*7;
for t=0:6
    v(T+t)=Vwf/7;
end
%v(1:T)=Vwf(1)/7;

% Using TRUE early data to date
T=300:size(Vacc3,2);
V3(2:8,T,:)=Vacc3(2:8,T,:);
V2(2:8,T,:)=Vacc2(2:8,T,:);
Vacc2(1,:,:)=sum(Vacc2(2:8,:,:),1);
Vacc3(1,:,:)=sum(Vacc3(2:8,:,:),1);
% and guess for W, S, NI
for R=9:11
    V3(R,T,:)=squeeze(Vacc3(1,T,:)).*(ones(length(T),1)*(Region_PP(R,:)./Region_PP(1,:)));
    V2(R,T,:)=squeeze(Vacc2(1,T,:)).*(ones(length(T),1)*(Region_PP(R,:)./Region_PP(1,:)));
end

%FORECASTING FORWARDS
T=(size(Vacc3,2)+1):(datenum(2023,7,1) +1-datenum(2020,1,1));

%v(v>2.5e6/7)=2.5e6/7;
RPPall=sum(Region_PP(1:11,:),2)/sum(Region_PP(2:11,:),'all');   % RPPall = proportion of population in UK
for R=2:11
    RPP(R,1:21)=Region_PP(R,:)./sum(Region_PP(2:11,:),1);
end

for R=2:11
    
    UpTake(1:21)=min(propBUT(1:21).*squeeze(sum(V2(R,:,1:21),2))'./Region_PP(R,1:21),0.99*ones(1,21));
    
    for t=T
        
        V3(R,t,:)=0;  V2(R,t,:)=0; Phase(t)=0;
        Eligible=sum(squeeze(V2(R,1:(t-Separation),:)),1)-sum(squeeze(V3(R,1:(t-1),:)),1);

        Available=RPPall(R)*v(t);
        
        
        TbV=aW;
        Already=sum(V3(R,1:t,TbV),'all');  Required=mean(0.9*(0.7e6+1.6e6)*sum(Region_PP(R,aW),2)/sum(Region_PP(2:11,aW),'all')); % 90% in Health Care Workers
        if Already<Required & Available>0
            Given=min(Required-Already,Available/2);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=1; 
        end

        TbV=a80;
        Already=sum(V3(R,1:t,TbV),'all');  Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % over 80s & CareHomes
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=2; 
        end
        
        TbV=(75/5)+1;
        Already=sum(V3(R,1:t,TbV),'all');  Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 75-79
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=3; 
        end
        
        TbV=(70/5)+1;
        Already=sum(V3(R,1:t,TbV),'all');  Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 70-74
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=4; 
        end
        
        TbV=[5:14];
        Already=sum(V3(R,1:t,TbV),'all');  Required=mean(UpTake(TbV))*(2.3e6+2.2e6)*sum(Region_PP(R,TbV),2)/sum(Region_PP(2:11,TbV),'all'); % CEV
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=4.5; 
        end
 
        TbV=(65/5);
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 65-69
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=5; 
        end
        
             
       TbV=[5:12];
       Already=sum(V3(R,1:t,TbV),'all');   Required=mean(UpTake(TbV))*(2.3e6+2.2e6+8.5e6)*sum(Region_PP(R,TbV),2)/sum(Region_PP(2:11,TbV),'all'); % Health Conditions <65
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=6; 
        end            
        
        TbV=(60/5)+1;
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 60-64
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=7; 
        end                     
        
        TbV=(55/5)+1;
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 55-59
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=8; 
        end 
        
        TbV=(50/5)+1;
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 50-54
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=9; 
        end 
        
        TbV=[9:10];
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 40-49
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=10; 
        end 
        
        TbV=[7:8];
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 30-39
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=11; 
        end 
        
        TbV=[5:6];
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 20-29
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=12;
        end
                                                 
        TbV=[4];
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 15-19
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=13; 
        end 
                                                        
        TbV=3;
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 10-14
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=14; 
        end 
                                                            
        TbV=2;
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 5-9
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=15; 
        end 
                                                                
        TbV=1;
        Already=sum(V3(R,1:t,TbV),'all');   Required=sum(UpTake(TbV).*Region_PP(R,TbV),'all'); % 0-4
        if Already<Required & Available>0
            Given=min(Required-Already,Available);
            Given=min(Eligible(TbV),Region_PP(R,TbV)*Given./sum(Region_PP(R,TbV)));
            V3(R,t,TbV)=Given;
            Available=Available-sum(Given);
            Phase(t)=16; 
        end        
    end
end

% Remove any impossible values
for R=2:11
    for A=1:21
        if squeeze(sum(V3(R,:,A),2))>Region_PP(R,A)
            V3(R,:,A)=V3(R,:,A)*0.99*Region_PP(R,A)/squeeze(sum(V3(R,:,A),2));
        end   
    end
end

% Put in the delay to effect
V3(:,[1:end]+Delay,:)=V3;  V3(:,1:Delay,:)=0;

