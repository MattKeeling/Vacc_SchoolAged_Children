function []=Generate_Output(DATE_STR)

%%

Matlab_Save=1;

Rmax=11;
Rtot=[2:11];

RUN_UNTIL=datenum(2023,2,1);

START_VACCINATION_12_17=datenum(2021,9,1);
START_VACCINATION_5_11=datenum(2022,3,1);

MaxType=4; % maximum number of variants.

for LOOP=0:2 % LOOP determines the vaccination strategy 
    
    load(['Data_File_' DATE_STR '.mat']);
    
    Date=datenum(DATE_STR,'dd_mm_yy')+1-datenum(2020,1,1); % days, with 01/01/2020 being day 1.
    
    [V1,V2,Phase, RatioPf]=Vaccination_over_Time_without_Kids(UT,Date); V3=0*V2;
    [V3,Phase]=Boosters_over_Time_MTP_UK(BUT,V2,Date);   
    
    load Regional_PP.mat
    
    if LOOP==0 % No Vaccination of Kids
    end
    
    if LOOP==1 || LOOP==2 % Vaccinate 12-17 year olds, with 2 doses
        clear v1 v2
        Kids_Uptake=0.7;  %70%
        Vacc_Speed=150e3;  % 150,000 per week
        
        t=(START_VACCINATION_12_17+1-datenum(2020,1,1)):size(V1,2);
        NumKids=sum(Region_PP(2:11,4)*0.6 + Region_PP(2:11,4)*0.6);
        v1(2:11,t,4)=(Vacc_Speed/7)*Region_PP(2:11,4)*0.6*ones(1,length(t))/NumKids;  % 60% of 15-19 age group
        v1(2:11,t,3)=(Vacc_Speed/7)*Region_PP(2:11,3)*0.6*ones(1,length(t))/NumKids;  % 60% of 10-14 age group
        m=find(cumsum(sum(v1(2:11,:,4)))>Kids_Uptake*sum(Region_PP(2:11,4)*0.6));  % find when > Kids_Uptake have been vaccinated
        v1(2:11,m,4)=0;
        m=find(cumsum(sum(v1(2:11,:,3)))>Kids_Uptake*sum(Region_PP(2:11,4)*0.6));  % find when > Kids_Uptake have been vaccinated
        v1(2:11,m,3)=0;
        v2(:,t,:)=v1(:,t-12*7,:);  % second dose after 12 weeks.    
        V1(:,:,1:4)=V1(:,:,1:4)+v1; V2(:,:,1:4)=V2(:,:,1:4)+v2;
    end
    
    if LOOP==2 % Vaccinate 5-11 year olds, with 2 doses
        clear v1 v2
        Kids_Uptake=0.7;  %70%
        Vacc_Speed=150e3;  % 150,000 per week
        
        t=(START_VACCINATION_5_11+1-datenum(2020,1,1)):size(V1,2);
        NumKids=sum(Region_PP(2:11,2)*1.0 + Region_PP(2:11,3)*0.4);
        v1(2:11,t,2)=(Vacc_Speed/7)*Region_PP(2:11,2)*1.0*ones(1,length(t))/NumKids;  % 100% of 5-9 age group
        v1(2:11,t,3)=(Vacc_Speed/7)*Region_PP(2:11,3)*0.4*ones(1,length(t))/NumKids;  % 40% of 10-14 age group
        m=find(cumsum(sum(v1(2:11,:,2)))>Kids_Uptake*sum(Region_PP(2:11,2)*1.0));  % find when > Kids_Uptake have been vaccinated
        v1(2:11,m,2)=0;
        m=find(cumsum(sum(v1(2:11,:,3)))>Kids_Uptake*sum(Region_PP(2:11,4)*0.6));  % find when > Kids_Uptake have been vaccinated
        v1(2:11,m,3)=0;
        v2(:,t,:)=v1(:,t-12*7,:);  % second dose after 12 weeks.    
        V1(:,:,1:3)=V1(:,:,1:3)+v1; V2(:,:,1:3)=V2(:,:,1:3)+v2;  
    end
    
    
    %%
    clear nALL nDEATHS nHOSP_AD nHOSP_OCC nICU_OCC nICU_AD
    
    Names={'England','East of England','London','Midlands','North East and Yorkshire','North West','South East','South West','Wales','Scotland','Northern Ireland','United Kingdom'};
   
    %%
    % generate run_stops by region
    Comps = zeros(4,Num_Comp+1,11);
    
    tmpRS=RUN_STOP;  Z=(RUN_STOP(end-1)+7):7:(ceil(RUN_UNTIL+1-datenum(2020,1,1)));
    tmpRS(Num_Comp+[1:length(Z)]-1)=Z;
    
    RUN_STOPs = zeros(11,Num_Comp+1);
    for Region = 2:11
        RUN_STOPs(Region,1:length(tmpRS)) = tmpRS;
    end
    
    
    Keep_RUN_STOPs=RUN_STOPs(2,:);
    
    if size(V1,2)<max(RUN_STOPs,[],'all')
        mx=max(RUN_STOPs,[],'all')+100;
        for r=1:11  for a=1:21
                V1(r,end:mx,a)=V1(r,end,a);
                V2(r,end:mx,a)=V2(r,end,a);              
            end
        end
    end
    
    %%
    
    maxtime = max(RUN_STOPs(end,:))+30; % needs to be big enough to take account of any lags
    nALL = zeros(11,maxtime,84);
    nDEATHS = zeros(11,maxtime,84);
    nHOSP_AD = zeros(11,maxtime,84);
    nHOSP_OCC = zeros(11,maxtime,84);
    nICU_AD = zeros(11,maxtime,84);
    nICU_OCC = zeros(11,maxtime,84); 
    
    % generate Comp by region
    clear Comps
    
    RUN_STARTs=[82*ones(11,1) RUN_STOPs(:,1:(end-1))];
     
    % NEED TO ADD IN DECLINE IN PRECAUTIONARY BEHAVIOUR
    NC=Num_Comp;
    for Region = Rtot
        Comps(Region,1:NC)=COMPLIANCE(1:NC,Region);
        Comps(Region,NC:size(RUN_STOPs,2))=Comps(Region,NC);
        CompsO(Region,1:NC)=COMPLIANCEO(1:NC,Region);
        CompsO(Region,NC:size(RUN_STOPs,2))=CompsO(Region,NC);
        m=find(RUN_STOPs(Region,:)>Date);  Decline=0.025; 
        for i=1:length(m)
            Comps(Region,m(i))=max(Comps(Region,m(i)-1)-Decline,0.001);  %Drop compliance slowly each week until zero
        end
        CompsO(Region,m)=CompsO(Region,m(1)-1)*Comps(Region,m)/Comps(Region,m(1)-1);
    end
    
    MAX_REGION=Rmax;
    
    if Date<700 % Don't include Omicron
        IMPORT_DAY(23:33)=1e5; IMPORT_LEVEL(23:33)=0;
    end
    if Date<450 % Don't include Delta
        IMPORT_DAY(12:22)=1e5; IMPORT_LEVEL(12:22)=0;
    end
    
    TXT_STR=['Warwick_Output_Loop' num2str(LOOP)];
                 
    %parfor Region=[2:8]
       for Region=2:8
           
           fprintf(1,'Simulating %s, Vaccine Loop %d\n',Names{Region},LOOP);
           
            [nDC, ~, ~, nHospital, inHospital, nICU, inICU, nDeaths, ~, All, ~]=Simulate_One_Region(Region, TAU, ALPHA, INC_P, SCALING(Region), FACTOR(Region), H_FACTOR(Region:11:end), I_FACTOR(Region:11:end), D_FACTOR(Region:11:end),  ...
            H_STRETCH(Region:11:end), I_STRETCH(Region:11:end), LAG(Region:11:end), START_DATE(Region)+1, 1, Comps(Region,:), CompsO(Region,:), RUN_STOPs(Region,:), NV_BETA(Region:11:end)', NV_SPEED(Region:11:end)', IMPORT_DAY(Region:11:end),IMPORT_LEVEL(Region:11:end),...
            squeeze(V1(Region,:,:)), squeeze(V2(Region,:,:)), squeeze(V3(Region,:,:)), Transmission_Reduction, VEffI, VEffS, VEffH, VEffD, RatioPf,...
            Detection, Susceptibility, GAMMA, WaningSpeed, WaningEfficacy, [0.1 0], MaxType, [], 0);
        
        nDeaths=nDeaths.* ((1+DH_SCALING(Region)*sum(inHospital,2)./sum(Region_PP(Region,:),2))*ones(1,size(nDeaths,2)));
        
        T=1:size(All,1);

        padding = maxtime-length(T);
        nALL(Region,:,:)=[All; zeros(padding,84)];
        nDEATHS(Region,:,:)=[nDeaths; zeros(padding,84)];
        nHOSP_OCC(Region,:,:)=[inHospital; zeros(padding,84)];
        nHOSP_AD(Region,:,:) = [nHospital; zeros(padding,84)];
        nICU_OCC(Region,:,:) = [inICU; zeros(padding,84)];
        nICU_AD(Region,:,:) = [nICU; zeros(padding,84)];
    end
    
    
    
    %% NEED TO COMBINE THE OLD & NEW VARIANT !
    
    A=1:size(nALL,3); a=1:21; b=21+a; c=42+a; d=63+a;
    nALL=nALL(:,:,a)+nALL(:,:,b)+nALL(:,:,c)+nALL(:,:,d);
    nDEATHS=nDEATHS(:,:,a)+nDEATHS(:,:,b)+nDEATHS(:,:,c)+nDEATHS(:,:,d);
    nHOSP_OCC=nHOSP_OCC(:,:,a)+nHOSP_OCC(:,:,b)+nHOSP_OCC(:,:,c)+nHOSP_OCC(:,:,d);
    nHOSP_AD=nHOSP_AD(:,:,a)+nHOSP_AD(:,:,b)+nHOSP_AD(:,:,c)+nHOSP_AD(:,:,d);
    nICU_OCC=nICU_OCC(:,:,a)+nICU_OCC(:,:,b)+nICU_OCC(:,:,c)+nICU_OCC(:,:,d);
    nICU_AD=nICU_AD(:,:,a)+nICU_AD(:,:,b)+nICU_AD(:,:,c)+nICU_AD(:,:,d);

    if Matlab_Save        
        save([TXT_STR '_' DATE_STR '.mat'],'nALL','nDEATHS','nHOSP_OCC','nHOSP_AD','nICU_OCC','nICU_AD');
    end
    
    
end

end