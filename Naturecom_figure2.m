%% LOAD BEHAVIOR FILE
loadbehavior =1
if loadbehavior ==1
    server = 'D:\Dropbox'
    cd ([  server '\Work_related\Duration opto_2' ])
    clear all
    load Data_figure2 
end

thismarkerdize =30
ErrobarlineWidth= 1
%%
% Select Animal actions (left 0 vs right 1 turn)
Animal_Action_Left_Right = DataB.Response.Action;
% [6 8 2 3 5 9 16 13]these rats have swithced ruls that is T2>T1 go left
% and the rest T2>T1 go right. 
Animal_Action_Left_Right(ismember(DataB.Info.Rat,[6 8 2 3 5 9 16 13]))= 1-Animal_Action_Left_Right(ismember(DataB.Info.Rat,[6 8 2 3 5 9 16 13]));
thisactionnew=Animal_Action_Left_Right;
NTD= (DataB.Stimulus.Duration{2} -DataB.Stimulus.Duration{1} )./(DataB.Stimulus.Duration{2} +DataB.Stimulus.Duration{1} );
NTD = (round(NTD*100)/100);
opto = DataB.Stimulus.opto;
Pscyometric_dataromance = DataB.Response.Correctness;
Bootstrap_percnet_choice =[];
Fit_Parameteres=[];
theseNTD = [-0.3485   -0.2400      -0.1170       0  0.1170     0.2400    0.3496 ];
Pscyometric_data=[];
Correct=[];




%% creat figure  2C psychomteric curves vs1 photoexcitation duration rats

%%

NTDcounter=0
figure
subplot 221

Rats =unique(DataB.Info.Rat(DataB.Info.Optoinfo==1));
allRats{3} = [9 10];
allRats{2} = [1:4 5  10  9  ];
allx = linspace(-0.35,0.35,1000);
allRats{1} = [ 1:5 ];
[thesepairs, ccc] = unique([DataB.Stimulus.Duration{2} ;DataB.Stimulus.Duration{1} ]','rows');
theseNTD = [-0.3485   -0.2400      -0.1170       0  0.1170     0.2400    0.3496 ];

for Condition = [1:3]
    
    
    thesedays = unique(DataB.Info.Date );
    thuissum=[];
    perf=[];
    Rats = allRats{Condition};
    for Ratid=1:numel(Rats)
    thesedays=unique(DataB.Info.Date(ismember(DataB.Info.Rat,Rats(Ratid))));

        andthiscg =DataB.Info.TrialNumber>0&DataB.Stimulus.Duration{1}<335&...
        DataB.Stimulus.Duration{1}>333&ismember(DataB.Info.Rat ,Rats(Ratid))...
        &ismember(DataB.Info.Date ,thesedays([1:end]))&DataB.Response.Perfromance>0.5;
        for thiss=1:size(theseNTD,2)
            
            NTDcounter=NTDcounter+1;
            thiscondition = NTD>theseNTD(thiss)-0.05&NTD<theseNTD(thiss)+0.05&andthiscg;
            thuissum{1}(Ratid,thiss) = sum((opto==7)&thiscondition);
            thuissum{1}(Ratid,thiss) = sum((opto==5)&thiscondition);
            thuissum{1}(Ratid,thiss) = sum((opto==3|opto==1)&thiscondition);
            perf{1}(Ratid,thiss)=mean(thisactionnew((opto==7)&thiscondition));
            perf{2}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1}<90 ...
                &DataB.Stimulus.Nominal{2}<90&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1}...
                &(opto==3|opto==1)&thiscondition));
            perf{3}(Ratid,thiss)=mean(thisactionnew((opto==5)&thiscondition));
            perf{4}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{2} <DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));
            perf{5}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{2} <90&DataB.Stimulus.Nominal{1} <90 &(opto==3|opto==1)&thiscondition));
            perf{6}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1} >90&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));
            perf{6}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1} >90&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));
            
            
        end
        
    end
    
    
    
    
    
    
    bin=theseNTD;

    thiscolor=[111 203 212; 8 8 8;238 59 35]./255;
      

    F=@(g,I,u,v,x) (g+(1-g-I)./(1+exp(-(x-u)/v)))*100;
    
    UL=[0.5, 0.5, 10, 5];
    SP=[ 0.4, 0.4, 5, 0.2];
    LM=[0, 0, -10, 0.1];
    
    
    UL=[0.5, 0.5, 10, 5];
    SP=[ 0.4, 0.4, 5, 0.3];
    LM=[0, 0, -10, 0.02];
    plot([0 0], [ 0 100],'-','color',[0.5 0.5 0.5],'Linewidth',2)
    hold on
    plot([-0.36 0.36], [ 50 50],'-','color',[0.5 0.5 0.5],'Linewidth',2)
    
    
    psycho=thesepairs(:,2)==334;
    if Condition==2
     i=[ 5]
    elseif Condition==3
    i=[1]
    elseif  Condition==1
     i=[1]
    end
        
        bootfit= (fit(bin',nanmean(perf{i})'*100,F,'StartPoint',SP,'Upper',UL,'Lower',LM));
        hold on
        hL1 = plot(bootfit,'-');
        hold on
        set([hL1 ],'color',thiscolor(Condition,:),'Linewidth',2)
        plot(bin, nanmean(perf{i})*100,'.-','color',thiscolor(Condition,:),'Markersize',thismarkerdize)
         h=errorbar(bin, nanmean(perf{i})*100,std(perf{i}*100)./sqrt(size(perf{i},1)),...
             'color',thiscolor(Condition,:),'Markersize',30);
        h.LineWidth = ErrobarlineWidth;
        h.LineStyle = 'none';
        
    
        
        fitparms(Condition,:) = coeffvalues(bootfit)
        allptoints (Condition,:)= F(  fitparms(Condition,1), fitparms(Condition,2), fitparms(Condition,3), fitparms(Condition,4),allx)
end


title('All Conditions: duration rats')
set(gca,'xtick',[-0.3:0.1:0.3])
set(gca,'ytick',[0:20:100])
set(gca,'xtick',[-0.2 0 0.2])
set(gca,'xtickLabel',[-0.2 0 0.2])
ylabel('Choice T2>T1 (%)')
xlabel('\DeltaT')
set(gca,'FontSize',15)
hold on
axis square
set(gca,'TickDir', 'out' )
box off

xlim([-0.36 0.36])
ylim([0 100])
legend off
figurecorrectionOPTO(gca,[0.02 0.02],12)


%%

subplot 222

Rats = [1:5];
thuissum=[];
perf=[];
NTDcounter = 0;
thisactionnew =Animal_Action_Left_Right;
alldates=[];
for Ratid=1:numel(Rats)
    thesedays=unique(DataB.Info.Date(ismember(DataB.Info.Rat,Rats(Ratid))));
 
        andthiscg =DataB.Info.TrialNumber>0&DataB.Stimulus.Duration{1}<335&...
        DataB.Stimulus.Duration{1}>333&ismember(DataB.Info.Rat ,Rats(Ratid))...
        &DataB.Stimulus.Nominal{2}>78&...
        DataB.Stimulus.Nominal{1}>78&ismember(DataB.Info.Date ,thesedays([ 1:end]))&DataB.Response.Perfromance>0.5;
        alldates = [alldates thesedays];

    for thiss=1:size(theseNTD,2)

        NTDcounter=NTDcounter+1;
        thiscondition = NTD>theseNTD(thiss)-0.05&NTD<theseNTD(thiss)+0.05&andthiscg;
        thuissum{1}(Ratid,thiss) = sum((opto==7)&thiscondition);
        thuissum{3}(Ratid,thiss) = sum((opto==5)&thiscondition);
        thuissum{2}(Ratid,thiss) = sum((opto==3|opto==1)&thiscondition);

        thuissumaction{1}(Ratid,thiss) = sum(thisactionnew((opto==7)&thiscondition));
        thuissumaction{3}(Ratid,thiss) = sum(thisactionnew((opto==5)&thiscondition));
        thuissumaction{2}(Ratid,thiss) = sum(thisactionnew((opto==3|opto==1)&thiscondition));


        perf{1}(Ratid,thiss)=mean(thisactionnew((opto==7)&thiscondition));
        perf{2}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1}<90 ...
            &DataB.Stimulus.Nominal{2}<90&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1}...
            &(opto==3|opto==1)&thiscondition));
        perf{3}(Ratid,thiss)=mean(thisactionnew((opto==5)&thiscondition));
        perf{4}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{2} <DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));
        perf{5}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1} <90&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));
        perf{6}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1} >90&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));


    end

end






bin=theseNTD;
thiscolor=[60 88 167;0 0 0; 93 186 70]./255;

F=@(g,I,u,v,x) (g+(1-g-I)./(1+exp(-(x-u)/v)))*100;

UL=[0.5, 0.5, 10, 5];
SP=[ 0.4, 0.4, 5, 0.2];
LM=[0, 0, -10, 0.1];


UL=[0.5, 0.5, 10, 6];
SP=[ 0.4, 0.4, 5, 0.3];
LM=[0, 0, -10, 0.01];
plot([0 0], [ 0 100],'-','color',[0.5 0.5 0.5],'Linewidth',2)
hold on
plot([-0.36 0.36], [ 50 50],'-','color',[0.5 0.5 0.5],'Linewidth',2)


for i=[1  3  ]

    bootfit= (fit(bin',mean(perf{i})'*100,F,'StartPoint',SP,'Upper',UL,'Lower',LM));
    hold on
    hL1 = plot(bootfit,'-');
    hold on
    set([hL1 ],'color',thiscolor(i,:),'Linewidth',2)
    plot(bin, mean(perf{i})*100,'.','color',thiscolor(i,:),'Markersize',thismarkerdize)
    h=errorbar(bin, mean(perf{i})*100,std(perf{i}*100)./sqrt(size(perf{i},1)),...
        'color',thiscolor(i,:),'Markersize',30);
    h.LineWidth = ErrobarlineWidth;
    h.LineStyle = 'none';

end


title('vS1 photoexcitation: duration rats')
set(gca,'xtick',[-0.3:0.1:0.3])
set(gca,'ytick',[0:20:100])
set(gca,'xtick',[-0.2 0 0.2])
set(gca,'xtickLabel',[-0.2 0 0.2])
ylabel('Choice T2>T1 (%)')
xlabel('\DeltaT')


set(gca,'FontSize',15)
hold on
axis square
set(gca,'TickDir', 'out' )
box off

xlim([-0.36 0.36])
ylim([0 100])
legend off
figurecorrectionOPTO(gca,[0.05 0.05],12)

%% creat figure  2d External light duration rats (Control)
subplot 223
Rats = [1:5];
Rats = [9 10];
Rats =[10    11    12    13    14    15]

thesedays = unique(DataB.Info.Date );
thefirst_7= find(opto==7,1,'first');
thuissum=[]
perf=[];
alldates=[];
% andthiscg =(DataB.Stimulus.Duration{1} <400&DataB.Stimulus.Duration{1} >200&ismember(DataB.Info.Rat ,Rat))&ismember(DataB.Info.Date ,[thesedays])&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1} ;
for Ratid=1:numel(Rats)
    thesedays=unique(DataB.Info.Date(ismember(DataB.Info.Rat,Rats(Ratid))));
 

            alldates = [alldates thesedays];

    andthiscg =DataB.Info.TrialNumber>0&DataB.Stimulus.Duration{1}<335&...
        DataB.Stimulus.Duration{1}>333&ismember(DataB.Info.Rat ,Rats(Ratid))...
        &DataB.Stimulus.Nominal{2}>78&...
        DataB.Stimulus.Nominal{1}>78&ismember(DataB.Info.Date ,thesedays([ 1:end]))&DataB.Response.Perfromance>0.5;



    for thiss=1:size(theseNTD,2)

        NTDcounter=NTDcounter+1;
        thiscondition = NTD>theseNTD(thiss)-0.05&NTD<theseNTD(thiss)+0.05&andthiscg;
        thuissum{1}(Ratid,thiss) = sum((opto==7)&thiscondition);
        thuissum{3}(Ratid,thiss) = sum((opto==5)&thiscondition);
        thuissum{2}(Ratid,thiss) = sum((opto==3|opto==1)&thiscondition);

        thuissumaction{1}(Ratid,thiss) = sum(thisactionnew((opto==7)&thiscondition));
        thuissumaction{3}(Ratid,thiss) = sum(thisactionnew((opto==5)&thiscondition));
        thuissumaction{2}(Ratid,thiss) = sum(thisactionnew((opto==3|opto==1)&thiscondition));


        perf{1}(Ratid,thiss)=mean(thisactionnew((opto==7)&thiscondition));
        %     perf(2,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{2} >DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));
        perf{2}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1}<90 ...
            &DataB.Stimulus.Nominal{2}<90&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1}...
            &(opto==3|opto==1)&thiscondition));
        %      perf(2,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1} <90&DataB.Stimulus.Nominal{2} <90&(opto==3|opto==1)&thiscondition));
        perf{3}(Ratid,thiss)=mean(thisactionnew((opto==5)&thiscondition));
        perf{4}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{2} <DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));
        perf{5}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1} <90&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));
        perf{6}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1} >90&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));
        perf{6}(Ratid,thiss)=mean(thisactionnew(DataB.Stimulus.Nominal{1} >90&DataB.Stimulus.Nominal{2} ==DataB.Stimulus.Nominal{1} &(opto==3|opto==1)&thiscondition));


    end

end




F=@(g,I,u,v,x) (g+(1-g-I)./(1+exp(-(x-u)/v)))*100;

UL=[0.5, 0.5, 10, 5];
SP=[ 0.4, 0.4, 5, 0.2];
LM=[0, 0, -10, 0.1];


UL=[0.5, 0.5, 10, 6];
SP=[ 0.4, 0.4, 5, 0.3];
LM=[0, 0, -10, 0.01];
plot([0 0], [ 0 100],'-','color',[0.5 0.5 0.5],'Linewidth',2)
hold on
plot([-0.36 0.36], [ 50 50],'-','color',[0.5 0.5 0.5],'Linewidth',2)


for i=[1  3  ]

    bootfit= (fit(bin',mean(perf{i})'*100,F,'StartPoint',SP,'Upper',UL,'Lower',LM));
    hold on
    hL1 = plot(bootfit,'-');
    hold on
    set([hL1 ],'color',thiscolor(i,:),'Linewidth',2)
    plot(bin, mean(perf{i})*100,'.','color',thiscolor(i,:),'Markersize',thismarkerdize)
    h=errorbar(bin, mean(perf{i})*100,std(perf{i}*100)./sqrt(size(perf{i},1)),...
        'color',thiscolor(i,:),'Markersize',30);
    h.LineWidth = ErrobarlineWidth;
    h.LineStyle = 'none';

end


title('External light duration rats (Control)')
set(gca,'xtick',[-0.3:0.1:0.3])
set(gca,'ytick',[0:20:100])
set(gca,'xtick',[-0.2 0 0.2])
set(gca,'xtickLabel',[-0.2 0 0.2])
ylabel('Choice T2>T1 (%)')
xlabel('\DeltaT')


set(gca,'FontSize',15)
hold on
axis square
set(gca,'TickDir', 'out' )
box off

xlim([-0.36 0.36])
ylim([0 100])
legend off
figurecorrectionOPTO(gca,[0.05 0.05],12)

%%
%close all
load(fullfile('D:\Dropbox\OptoDurationPaper\A_NATCOMSUBMISSION','WM_Speed_OptoBC_3rats.mat'))
data=TrialSeb;
data.action = ~(xor(data.cat,data.hitmiss))

% rats = [5 7 8 9];%unique(TrialSeb.rat);
% PsychoFun = 'CumGauss'; param = [100 10 0 0.001];
% p0 = [100 10 0 0]; %[mu sigma gamma lambda]
% lb = [90 1 0 0]; %LOWER BOUND VECTOR: [mu sigma gamma lambda]
% ub = [105 30 0 0]; %UPPER BOUND VECTOR :[mu sigma gamma lambda]
% % super_params = zeros(2,numel(rats),4);
% % CONDI=[120 140];
% threshold=[120 140];
%tasks=unique(TrialSeb.task(TrialSeb.condition==CONDI));
%rorder=randperm(12);

rats = [20,31:32];
conditions = [6,4];% 6 opto at stim 1, 4 opto at stim2

%dats=unique(data.date);

% ints2 = unique(data.std_comp);
% ints1 = unique(data.std_base);
% durs2 = unique(data.dur_comp);
% durs1 = unique(data.dur_base);

% cols=brewermap(12,'Spectral');
% colvals=[11,9]
%cols=cols(2:end,:)
pars=nan(2,4);
%cols=[ 0.2 0.65 1; 0 0.65 0.4]
cols=[96/255 195/255 50/255; 69/255 98/255 183/255]

index = data.std_base==80%&data.date==dats(5);%&(data.date>=dats(11)|data.date==dats(9)
ints2 = unique(data.std_comp(data.std_base==80));
NSD=(ints2-80)./(ints2+80)


%% Create psychometric curve out of the mean data between the 3 rats and add SEM.
clear res
subplot 224
hold off
bin=NSD;
for co = 1:numel(conditions)
 clear res

           for iint2 = 1:numel(ints2)
               for irats=1:numel(rats)
            res.act(irats,iint2) = mean(data.action(index & data.std_comp==ints2(iint2)&data.optopos==conditions(co)&data.rat==rats(irats)));%&data.trialID<=150)
            res.var(irats,iint2)= std(data.action(index & data.std_comp==ints2(iint2)&data.optopos==conditions(co)&data.rat==rats(irats)));
            res.count(irats,iint2) = sum(index & data.std_comp==ints2(iint2)&data.optopos==conditions(co)&data.rat==rats(irats));%&data.trialID<=150
               end
               end
           
           res.act(:)=flipud(res.act(:));
           res.var(:)=flipud(res.var(:));
           res.count(:)=flipud(res.count(:));
           
      p0 = [0 0.1 0 0];               % p0 = [mean(durs2) 0.001 0.1 0.1];
        xaxis = linspace(min(NSD),max(NSD),201);
      lb = [min(NSD) -1 0 0];
                      % lb = [min(durs2) 0.0001 0.1 0.1];
      ub = [max(NSD) 1 0.15 0.15];
                       %ub= [max(durs2) 0.1 0.1 0.1];
        pars(co,:)= psychoMLE(p0,NSD,mean(res.act),lb,ub);%... &res(ir).count(iint2,:)
        curve = cumulativegaussianLapse(pars(co,:),xaxis);
        

       % SEMrats(co,:)=std(res.act)/sqrt(3);

%         subplot(2,2,ir)
        % lh(co)=plot(NSD,mean(res.act)*100,'o','color',cols(co,:),'MarkerFaceColor',cols(co,:),'MarkerSize',8);


LM=[0, 0, -10, 0.1];


UL=[0.5, 0.5, 10, 6];
SP=[ 0.4, 0.4, 5, 0.3];
LM=[0, 0, -10, 0.01];
plot([0 0], [ 0 100],'-','color',[0.5 0.5 0.5],'Linewidth',2)
hold on
    bootfit= (fit(bin',mean(res.act)'*100,F,'StartPoint',SP,'Upper',UL,'Lower',LM));
    hold on
    hL1 = plot(bootfit,'-');
    hold on
    set([hL1 ],'color',cols(co,:),'Linewidth',2)
    plot(bin, mean(res.act)'*100,'.','color',cols(co,:),'Markersize',thismarkerdize)
             lh(co)=errorbar(NSD,mean(res.act)*100,std(res.act)*100/sqrt(3),'o','color',cols(co,:),'MarkerFaceColor',cols(co,:),'MarkerSize',8);

        %title(['Rat ',num2str(rats(ir))])
%        colormap(cols)
%   c= colorbar
%         set(c,'ticks',[0.05 0.5 0.95],'ticklabels',{45,101,156},'direction','normal','box','off','TickLength',0.001)
%         ylabel(c,'Speed [mm/s]')
%         xlim([200 800])
  end

        xlabel('NSD')
        ylabel('p (S2>S1)')
        legend(lh,{'light On Stim1','light On Stim2'})
        %xcale('log')
        
        %datashift=datapoints(1,:)-datapoints(2,:)
    figurecorrectionOPTO(gca,[0.05 0.05],12)
xlim([-0.3 0.3])

%%



