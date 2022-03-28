%% Run a single sample of the July code

clear
Generate_Output;

%% Plot Total Infections, Hospital Admissions and Deaths for 4 different age groupings and 7 different behaviours
% Without childhood vaccination is in grey, with 12-17 year old vaccination is in pale red, and with 5-17 year old vaccination is in pale green
% The outline colours correspond to the assumed behaviour going forwards.

clear

Color=[0.8 0.8 0.8; 1 0.8 0.8; 0.7 0.9 0.7];
ECol=[0.6 0 0.8; 1 0 0; 0.3 0.3 1; 0 0.85 0; 1 0 1; 1 0.7 0; 0 0.9 0.8];


figure(2); clf; set(gcf,'position',[1 541 1174 804]);

for Vacc_Loop=0:2
    
    load(['Warwick_Output_Loop' num2str(Vacc_Loop) '_26_06_21.mat']);
    
    VarNames={'ALL','HOSP_AD','DEATHS'};
    TrueNames={'Infections','Hospital Admissions','Deaths'};
    
    T=1:size(ALL,3); Weighting=1+0*T; Weighting(1:(datenum(2021,7,18)+1-datenum(2020,1,1)))=0;
    Weighting((datenum(2022,1,1)+1-datenum(2020,1,1)):(datenum(2022,12,31)+1-datenum(2020,1,1)))=1-[1:365]/365;
    Weighting((datenum(2023,1,1)+1-datenum(2020,1,1)):end)=0;
    
    X=1/4; x=1/8;
    
    for Names=1:3
        
        for PreCBeh=1:7
            
            eval(['Y=squeeze(' VarNames{Names} '(PreCBeh,:,:,:));']);
            subplot(4,3,Names); y=Y(:,:,2)+0.4*Y(:,:,3);
            fill(PreCBeh+Vacc_Loop*X+x*[-1 1 1 -1]-X,[0 0 1 1]*(sum(y(2:8,:),1)*Weighting'),'r','FaceColor',Color(Vacc_Loop+1,:),'EdgeColor',ECol(PreCBeh,:),'LineWidth',2); hold on
            ylabel('Aged 5-11'); title(TrueNames{Names}); set(gca,'XLim',[0.2 7.8]);
            subplot(4,3,Names+3); y=0.6*Y(:,:,3)+0.6*Y(:,:,4);
            fill(PreCBeh+Vacc_Loop*X+x*[-1 1 1 -1]-X,[0 0 1 1]*(sum(y(2:8,:),1)*Weighting'),'r','FaceColor',Color(Vacc_Loop+1,:),'EdgeColor',ECol(PreCBeh,:),'LineWidth',2); hold on
            ylabel('Aged 12-17'); set(gca,'XLim',[0.2 7.8]);
            subplot(4,3,Names+6); y=Y(:,:,2)+Y(:,:,3)+0.6*Y(:,:,4);
            fill(PreCBeh+Vacc_Loop*X+x*[-1 1 1 -1]-X,[0 0 1 1]*(sum(y(2:8,:),1)*Weighting'),'r','FaceColor',Color(Vacc_Loop+1,:),'EdgeColor',ECol(PreCBeh,:),'LineWidth',2); hold on
            ylabel('Aged 5-17'); set(gca,'XLim',[0.2 7.8]);
            subplot(4,3,Names+9); y=sum(Y(:,:,:),3);
            fill(PreCBeh+Vacc_Loop*X+x*[-1 1 1 -1]-X,[0 0 1 1]*(sum(y(2:8,:),1)*Weighting'),'r','FaceColor',Color(Vacc_Loop+1,:),'EdgeColor',ECol(PreCBeh,:),'LineWidth',2); hold on
            ylabel('All Ages'); xlabel('Vaccinated ages'); set(gca,'XLim',[0.2 7.8]);
        end
    end
    drawnow;
end

a=axes;
for VL=0:2
    h(VL+1)=fill([-1 1 1 -1],[0 0 1 1],'r','FaceColor',Color(VL+1,:),'EdgeColor','k','LineWidth',0.5); hold on
end
for PreCBeh=1:7
    h(PreCBeh+3)=fill([-1 1 1 -1],[0 0 1 1],'r','FaceColor','w','EdgeColor',ECol(PreCBeh,:),'LineWidth',2); hold on
end
axis([10 11 10 11]); axis off
legend(h,'Vaccine >18','Vaccine 12-17','Vaccine 5-17','Behaviour 1','Behaviour 2','Behaviour 3','Behaviour 4','Behaviour 5','Behaviour 6','Behaviour 7','Location','NorthOutside','NumColumns',10);
a.Position=[0.1300 0.1044 0.7750 0.8337];


