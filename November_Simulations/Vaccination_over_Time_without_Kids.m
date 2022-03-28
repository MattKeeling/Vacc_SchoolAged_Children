function [V1,V2,Phase,RatioPfMd]=Vaccination_over_Time_without_Kids(UpTake, Date)

Separation=11*7;
Delay=14;
Inf_to_Symp=5;

load Surrogate_Vacc_Data.mat
load Regional_PP.mat

% Remove the later data and use approximation instead
Vacc1(:,Date:end,:)=[];
Vacc2(:,Date:end,:)=[];

a80=[17:21]; p80=Region_PP(1:11,a80)/sum(Region_PP(2:11,a80),'all');  %Note average for UK data !!
aW=[4:13]; pW=Region_PP(1:11,aW)/sum(Region_PP(2:11,aW),'all');  %Note average for UK data !!

T=[datenum(2021,1,4):datenum(2022,11,1)] +1-datenum(2020,1,1);
V1=zeros(11,2000,21);  V2=zeros(11,2000,21);

% Approximations to the amount of vaccine going forwards if we don't have
% the data (or if we didn't have the data at the time).
v(350:370)=5e5/7;
v(370:440)=2.5e6/7;
v(440:540)=3.2e6/7;
v(540:600)=1.3e6/7;
v(600:2000)=2.5e5/7;

% Using TRUE early data to date
T=300:size(Vacc1,2);
V1(2:8,T,:)=Vacc1(2:8,T,:);
V2(2:8,T,:)=Vacc2(2:8,T,:);
Vacc1(1,:,:)=sum(Vacc1(2:8,:,:),1);
Vacc2(1,:,:)=sum(Vacc2(2:8,:,:),1);

%FORECASTING FORWARDS
T=(size(Vacc1,2)+1):(datenum(2023,7,1) +1-datenum(2020,1,1));

RPPall=sum(Region_PP(1:11,:),2)/sum(Region_PP(2:11,:),'all');   % RPPall = proportion of population in UK
for R=2:11
    RPP(R,1:21)=Region_PP(R,:)./sum(Region_PP(2:11,:),1);
end

Keep_UpTake=UpTake;
UpTake(1:3)=0; % remove 0-14 year olds.

for R=2:11
    
    for t=T
        
        if t>500
            Separation=Separation-1;   %reduce until at 8 weeks.
            Separation(Separation<8*7)=8*7;
        end
        
        V1(R,t,:)=0;
        V2(R,t,:)=sum(V1(R,1:(t-Separation),:),2)-sum(V2(R,1:(t-1),:),2);
        V2(R,t,1:3)=0; 
        
        V2(R,t,(V2(R,t,:)<0))=0;
        if sum(V2(R,t,:),3)>RPPall(R)*v(t)
            V2(R,t,:)=V2(R,t,:)*RPPall(R)*v(t)/sum(V2(R,t,:),3);
        end
        
        v1=v(t)*RPPall(R)-sum(V2(R,t,:),3); % Amount of vaccine left to distribute.
        %v1=v(t)*RPPall(R);
        v1(v1<0)=0;
        
        Phase(t)=2; TbV=a80;
        
        if sum(V1(R,1:t,aW),'all')>UpTake*(0.7e6+1.6e6)*sum(Region_PP(R,aW),2)/sum(Region_PP(2:11,aW),'all') % Done HCW
            VHCW=0;
        else
            VHCW=v1/2;
        end
        
        Vold=v1-VHCW;
        
        
        if sum(V1(R,1:t,TbV),'all')>sum(UpTake(TbV).*Region_PP(R,TbV),'all')  % over 80s & CareHomes
            Phase(t)=3; TbV=(75/5)+1;
            
            if sum(V1(R,1:t,TbV),'all')>sum(UpTake(TbV).*Region_PP(R,TbV),'all') % 75-79
                Phase(t)=4; TbV=(70/5)+1;
                
                if sum(V1(R,1:t,TbV),'all')>sum(UpTake(TbV).*Region_PP(R,TbV),'all') % 70-74
                    Phase(t)=4.5; TbV=[5:14];
                    
                    if sum(V1(R,1:t,TbV),'all')>mean(UpTake(TbV))*(2.3e6+2.2e6)*sum(Region_PP(R,TbV),2)/sum(Region_PP(2:11,TbV),'all') % Extremely vulnerable
                        Phase(t)=5; TbV=(65/5)+1;
                        
                        if sum(V1(R,1:t,TbV),'all')>sum(UpTake(TbV).*Region_PP(R,TbV),'all') % 65-69
                            Phase(t)=6; TbV=[5:12];
                            
                            if sum(V1(R,1:t,TbV),'all')>mean(UpTake(TbV))*(2.3e6+2.2e6+8.5e6)*sum(Region_PP(R,TbV),2)/sum(Region_PP(2:11,TbV),'all') % Health Conditions <65
                                Phase(t)=7; TbV=(60/5)+1;
                                
                                if sum(V1(R,1:t,TbV),'all')>sum(UpTake(TbV).*Region_PP(R,TbV),'all') % 60-64
                                    Phase(t)=8; TbV=(55/5)+1;
                                    
                                    if sum(V1(R,1:t,TbV),'all')>sum(UpTake(TbV).*Region_PP(R,TbV),'all') % 55-59
                                        Phase(t)=9; TbV=(50/5)+1;
                                        
                                        if sum(V1(R,1:t,TbV),'all')>sum(UpTake(TbV).*Region_PP(R,TbV),'all') % 50-54
                                            Phase(t)=10; TbV=[9:10];
                                            
                                            if sum(V1(R,1:t,TbV),'all')>sum(UpTake(TbV).*Region_PP(R,TbV),'all') % 40-49
                                                Phase(t)=10; TbV=[7:8];
                                                
                                                if sum(V1(R,1:t,TbV),'all')>sum(UpTake(TbV).*Region_PP(R,TbV),'all') % 30-39
                                                    Phase(t)=10; TbV=[5:6];
                                                    
                                                    if sum(V1(R,1:t,TbV),'all')>sum(UpTake(TbV).*Region_PP(R,TbV),'all') % 20-29
                                                        Phase(t)=10; TbV=[4];
                                                        
                                                        if sum(V1(R,1:t,TbV),'all')>sum(UpTake(TbV).*Region_PP(R,TbV),'all') % 15-19
                                                            Phase(t)=10; TbV=3;
                                                            
                                                            if sum(V1(R,1:t,TbV),'all')>sum(UpTake(TbV).*Region_PP(R,TbV),'all') % 10-14 not used in this code
                                                                Phase(t)=10; TbV=2;
                                                                
                                                                if sum(V1(R,1:t,TbV),'all')>sum(UpTake(TbV).*Region_PP(R,TbV),'all') % 5-9 not used in this code
                                                                    Phase(t)=10; TbV=1;
                                                                    
                                                                    if sum(V1(R,1:t,TbV),'all')>sum(UpTake(TbV).*Region_PP(R,TbV),'all') % 0-5 not used in this code
                                                                        Phase(t)=0; TbV=aW; Vold=0;
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        V1(R,t,aW)=pW(R,:)*(VHCW);
        
        V1(R,t,TbV)=Vold*Region_PP(R,TbV)/sum(Region_PP(R,TbV),'all');
                
    end
end

% Remove any impossible values
for R=2:11
    for A=1:21
        if squeeze(sum(V1(R,:,A),2))>Region_PP(R,A)
            V1(R,:,A)=V1(R,:,A)*0.99*Region_PP(R,A)/squeeze(sum(V1(R,:,A),2));
        end
        if squeeze(sum(V2(R,:,A),2))>Region_PP(R,A)
            V2(R,:,A)=V2(R,:,A)*0.99*Region_PP(R,A)/squeeze(sum(V2(R,:,A),2));
        end
    end
end


% Put in the delay to effect
V1(:,[1:end]+Delay,:)=V1;  V1(:,1:Delay,:)=0;
V2(:,[1:end]+Delay,:)=V2;  V2(:,1:Delay,:)=0;


