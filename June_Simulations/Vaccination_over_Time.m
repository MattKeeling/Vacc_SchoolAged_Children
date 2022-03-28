function [V1,V2,VET1,VET2,RC,RCH,RCD,Phase,RatioPfMd]=Vaccination_over_Time(UpTake,nVEI,nVES, nVEH, nVED, Rollout, SLOWER, DATE)

%%
%
%%

Separation=11*7;
Delay=14;
Inf_to_Symp=5;

load Surrogate_Vacc_Data.mat

Vacc1(:,DATE:end,:)=[];
Vacc2(:,DATE:end,:)=[];
Vacc_Type1(:,DATE:end,:,:)=[];
Vacc_Type2(:,DATE:end,:,:)=[];

load Other_Data.mat

Vacc1=Vacc1*SLOWER;
Vacc2=Vacc2*SLOWER;
Vacc_Type1=Vacc_Type1*SLOWER;
Vacc_Type2=Vacc_Type2*SLOWER;
UpTake=UpTake*SLOWER;
Rollout=Rollout*SLOWER;

if length(nVEI)==2   % Assuming Pfizer=AZ
    VET1=ones(11,21)*nVEI(1);
    VE1=ones(11,21)*nVES(1);
    VEH1=ones(11,21)*nVEH(1);
    VED1=ones(11,21)*nVED(1);
    VET2=ones(11,21)*nVEI(2);
    VE2=ones(11,21)*nVES(2);
    VEH2=ones(11,21)*nVEH(2);
    VED2=ones(11,21)*nVED(2);
else
    % Vacc_Type1=Pf, AZ, Md, NA  % Vaccination already done
    NumPfMd=squeeze(sum(Vacc_Type1(2:8,:,:,[1 3]),[1 2 4]));
    NumAZ=squeeze(sum(Vacc_Type1(2:8,:,:,2),[1 2 4]));
    NumNA=squeeze(sum(Vacc_Type1(2:8,:,:,4),[1 2 4]));
    NumNA_to_PfMd=NumPfMd.*NumNA./(NumPfMd+NumAZ);
    NumNA_to_AZ=NumAZ.*NumNA./(NumPfMd+NumAZ);
    NumPfMd=NumPfMd + NumNA_to_PfMd;
    NumAZ=NumAZ + NumNA_to_AZ;
    
    %Vaccinations to come
    Add_Vacc = sum(Region_PP(2:8,:),1).*UpTake - (NumPfMd+NumAZ)';
    Add_Vacc(Add_Vacc<0)=0;
    
    AddVaccPfMd=0.4*Add_Vacc;
    AddVaccAZ=0.6*Add_Vacc;
    AddVaccPfMd(1:8)=Add_Vacc(1:8);
    AddVaccAZ(1:8)=0*Add_Vacc(1:8);
    
    RatioPfMd=(NumPfMd'+AddVaccPfMd)./(NumPfMd'+NumAZ'+AddVaccPfMd+AddVaccAZ);
    RatioAZ=1-RatioPfMd;
    
    % VE input = Pf1, Pf2, AZ1, AZ2
    VET1=ones(11,1).*(nVEI(1)*RatioPfMd + nVEI(3)*RatioAZ);
    VE1=ones(11,1).*(nVES(1)*RatioPfMd + nVES(3)*RatioAZ);
    VEH1=ones(11,1).*(nVEH(1)*RatioPfMd + nVEH(3)*RatioAZ);
    VED1=ones(11,1).*(nVED(1)*RatioPfMd + nVED(3)*RatioAZ);
    VET2=ones(11,1).*(nVEI(2)*RatioPfMd + nVEI(4)*RatioAZ);
    VE2=ones(11,1).*(nVES(2)*RatioPfMd + nVES(4)*RatioAZ);
    VEH2=ones(11,1).*(nVEH(2)*RatioPfMd + nVEH(4)*RatioAZ);
    VED2=ones(11,1).*(nVED(2)*RatioPfMd + nVED(4)*RatioAZ);
    
end
    
if length(UpTake)==1
    UpTake=ones(1,21)*UpTake;
end

% currently the same in all regions and all ages.

% Estimated cumulative vaccinations each week
if Rollout==0  % default
    Rollout=1;
end
 
%Estimated weekly doses in England
Vwf=[2.5 2.0 1.6 1.6 1.5 1.4 2.0]*Rollout*1e6;
Vwf_date=datenum(2021,6,14)+1-datenum(2020,1,1);

Vwf=Vwf*sum(Region_PP(2:11,:),'all','omitnan')/sum(Region_PP(2:8,:),'all','omitnan');

a80=[17:21]; p80=Region_PP(1:11,a80)/sum(Region_PP(2:11,a80),'all');  %Note average for UK data !!
aW=[4:13]; pW=Region_PP(1:11,aW)/sum(Region_PP(2:11,aW),'all');  %Note average for UK data !!

T=[datenum(2021,1,4):datenum(2022,11,1)] +1-datenum(2020,1,1);
V1=zeros(11,max(T),21);  V2=zeros(11,max(T),21);

% Assume last value of Vwf doses a week if we don't have data
T=[datenum(2021,1,25):datenum(2023,7,1)] +1-datenum(2020,1,1);
v(T)=Vwf(end)/7;

% Use estimation
T=Vwf_date + ([1:length(Vwf)]-1)*7;
for t=0:6
    v(T+t)=Vwf/7;
end


% Using TRUE early data to date
T=300:size(Vacc1,2);
V1(2:11,T,:)=Vacc1(2:11,T,:);
V2(2:11,T,:)=Vacc2(2:11,T,:);
Vacc1(1,:,:)=sum(Vacc1(2:8,:,:),1);
Vacc2(1,:,:)=sum(Vacc2(2:8,:,:),1);


%FORECASTING FORWARDS
T=(size(Vacc1,2)+1):(datenum(2023,7,1) +1-datenum(2020,1,1));

RPPall=sum(Region_PP(1:11,:),2)/sum(Region_PP(2:11,:),'all');   % RPPall = proportion of population in UK
for R=2:11
    RPP(R,1:21)=Region_PP(R,:)./sum(Region_PP(2:11,:),1);
end

for R=2:11
    
    for t=T
        
        if t>500
            Separation=Separation-1;   %reduce until at 8 weeks.
            Separation(Separation<8*7)=8*7;
        end
        
        
        V1(R,t,:)=0;
        V2(R,t,:)=sum(V1(R,1:(t-Separation),:),2)-sum(V2(R,1:(t-1),:),2);
        V2(R,t,(V2(R,t,:)<0))=0;
        if sum(V2(R,t,:),3)>RPPall(R)*v(t)
            V2(R,t,:)=V2(R,t,:)*RPPall(R)*v(t)/sum(V2(R,t,:),3);
        end
        
        v1=v(t)*RPPall(R)-sum(V2(R,t,:),3);
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
                                                            
                                                            if sum(V1(R,1:t,TbV),'all')>sum(UpTake(TbV).*Region_PP(R,TbV),'all') % 10-14
                                                                Phase(t)=10; TbV=2;
                                                                
                                                                if sum(V1(R,1:t,TbV),'all')>sum(UpTake(TbV).*Region_PP(R,TbV),'all') % 5-9
                                                                    Phase(t)=10; TbV=1;
                                                                    
                                                                    if sum(V1(R,1:t,TbV),'all')>sum(UpTake(TbV).*Region_PP(R,TbV),'all') % 0-5
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

Region_PP(Region_PP==0)=1;
v1=cumsum(V1,2);
v2=cumsum(V2,2);
v1o=v1-v2;
clear W1 W2 RC

for t=(2+Delay+Separation):size(V1,2)
    RC(:,t+Inf_to_Symp,:)=(squeeze(v1o(:,t-Delay,:)).*(VE1-VET1) + squeeze(v2(:,t-Delay,:)).*(VE2-VET2)) ./ (Region_PP-squeeze(v1o(:,t-Delay,:)).*VET1-squeeze(v2(:,t-Delay,:)).*VET2);
    RCD(:,t+Inf_to_Symp,:)=(squeeze(v1o(:,t-Delay,:)).*(VED1-VET1) + squeeze(v2(:,t-Delay,:)).*(VED2-VET2)) ./ (Region_PP-squeeze(v1o(:,t-Delay,:)).*VET1-squeeze(v2(:,t-Delay,:)).*VET2);
    RCH(:,t+Inf_to_Symp,:)=(squeeze(v1o(:,t-Delay,:)).*(VEH1-VET1) + squeeze(v2(:,t-Delay,:)).*(VEH2-VET2)) ./ (Region_PP-squeeze(v1o(:,t-Delay,:)).*VET1-squeeze(v2(:,t-Delay,:)).*VET2);
    
end
RC=1-RC;
RCH=1-RCH;
RCD=1-RCD;

% Put in the delay to effect
V1(:,[1:end]+Delay,:)=V1;  V1(:,1:Delay,:)=0;
V2(:,[1:end]+Delay,:)=V2;  V2(:,1:Delay,:)=0;


