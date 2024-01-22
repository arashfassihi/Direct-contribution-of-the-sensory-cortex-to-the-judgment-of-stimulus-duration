clear all
figure
loadbehavior = 1;
if loadbehavior ==1
thesesession = [];

%%
% cd 'D:\Dropbox\Work_related\Duration opto\'
    server = 'D:\Dropbox';
    cd (fullfile( server,'\OptoDurationPaper\A_NATCOMSUBMISSION'))
    
    clear all
    load Data_figure2 


    %%
    Animal_Action_Left_Right = DataB.Response.Action;
    TDI=(DataB.Stimulus.Duration{2} -DataB.Stimulus.Duration{1} )./(DataB.Stimulus.Duration{2} +DataB.Stimulus.Duration{1} );
    TDI = (round(TDI*100)/100);
    opto = DataB.Stimulus.opto;
    Pscyometric_dataromance = DataB.Response.Correctness;
        Animal_Action_Left_Right = DataB.Response.Action;

    Animal_Action_Left_Right(ismember(DataB.Info.Rat,[6 8 2 3 5 9 16 13]))= 1-Animal_Action_Left_Right(ismember(DataB.Info.Rat,[6 8 2 3 5 9 16 13]));
end
%% do bootstrap  (fast version) only mean choice bias
TDIcounter = 0;
minperf = 0.5

%  Rat =[ 11:14 10   16   ];
Rat =[1:5];
% Rat =[ 11 12 13 14   16]
%   Rat =[]

Bootstrap_percnet_choice =[];
Fit_Parameteres=[];
theseTDI = [         0  0.1170 ];
theseTDI = [-0.3485   -0.2400      -0.1170       0  0.1170     0.2400    0.3496 ];
theseTDI = [   -0.2400      -0.1170       0  0.1170       0.2400  ];
theseTDI = [   -0.2400      -0.1170       0  0.1170       0.2400  ];

Pscyometric_data=[];
Correct=[];
thesedatesrat9 =['03-Dec-2018';'06-Dec-2018';'11-Dec-2018';'01-Mar-2019';'04-Mar-2019';'05-Mar-2019';'06-Mar-2019';'11-Mar-2019';'12-Mar-2019';'13-Mar-2019';'14-Mar-2019';'26-Mar-2019';'29-Mar-2019';'01-Apr-2019';'04-Apr-2019']
desireddates9 =datenum(thesedatesrat9);
theseundesiredtrials = ismember(DataB.Info.Rat ,9)&~ismember(DataB.Info.Date,desireddates9);
theseundesiredtrials = ismember(DataB.Info.Rat ,9);

thesedatesrat10 =[ '10-May-2018';'14-May-2018';'16-May-2018';'17-May-2018';'21-May-2018';'22-May-2018';'24-May-2018';'28-May-2018';'29-May-2018']

desireddates10 =datenum(thesedatesrat10);
theseundesiredtrials =theseundesiredtrials|(ismember(DataB.Info.Rat ,10)&~ismember(DataB.Info.Date,desireddates10));
theseundesiredtrials =theseundesiredtrials|ismember(DataB.Info.Rat ,10);
theseundesiredtrials =~true(size(DataB.Info.Rat));

%  theseundesiredtrials = logical(zeros(size(theseundesiredtrials)));
clear Op
clear noOp


Desired_condition =DataB.Info.TrialNumber>0&DataB.Stimulus.Duration{1}<335&...
    DataB.Stimulus.Duration{1}>332&ismember(DataB.Info.Rat ,Rat)&DataB.Stimulus.Nominal{2}>78&...
    DataB.Stimulus.Nominal{1}>78&DataB.Stimulus.Nominal{1}<82&DataB.Stimulus.Nominal{2}<82&~theseundesiredtrials&DataB.Response.Perfromance>minperf;

Desired_condition_Nolightgrouped =DataB.Info.TrialNumber>0&DataB.Stimulus.Duration{1}<335&...
    DataB.Stimulus.Duration{1}>332&ismember(DataB.Info.Rat ,[1:5 10 9 ])&DataB.Stimulus.Nominal{2}>78&...
    DataB.Stimulus.Nominal{1}>78&DataB.Stimulus.Nominal{1}<90&DataB.Stimulus.Nominal{2}<90&~theseundesiredtrials&DataB.Response.Perfromance>minperf;

Desired_condition_Halo =DataB.Info.TrialNumber>0&DataB.Stimulus.Duration{1}<335&...
    DataB.Stimulus.Duration{1}>332&ismember(DataB.Info.Rat ,[10 9   ])&DataB.Stimulus.Nominal{2}>78&...
    DataB.Stimulus.Nominal{1}>78&DataB.Stimulus.Nominal{1}<90&~theseundesiredtrials&DataB.Stimulus.Nominal{2}<90&DataB.Response.Perfromance>minperf;


Desired_condition_ExternalLight =DataB.Info.TrialNumber>0&DataB.Stimulus.Duration{1}<335&...
    DataB.Stimulus.Duration{1}>332&ismember(DataB.Info.Rat ,[ 11:16 ])&DataB.Stimulus.Nominal{2}>78&...
    DataB.Stimulus.Nominal{1}>78&DataB.Stimulus.Nominal{1}<90&DataB.Stimulus.Nominal{2}<90&~theseundesiredtrials&DataB.Response.Perfromance>minperf;

for j=1:numel(theseTDI)

    tic

    TDIcounter=TDIcounter+1;
    updated_desired_condition = TDI>(theseTDI(j))-0.03&TDI<(theseTDI(j))+0.03&Desired_condition;
    updated_desired_condition_externallight = TDI>(theseTDI(j))-0.03&TDI<(theseTDI(j))+0.03&Desired_condition_ExternalLight;

    updated_desired_condition_nolight_grouped = TDI>(theseTDI(j))-0.03&TDI<(theseTDI(j))+0.03&Desired_condition_Nolightgrouped;
    updated_desired_condition_Halo = TDI>(theseTDI(j))-0.03&TDI<(theseTDI(j))+0.03&Desired_condition_Halo;



    DesireTrials.Halo.stim2=find(((opto==7)&updated_desired_condition_Halo));
    DesireTrials.Halo.stim1=find(((opto==5)&updated_desired_condition_Halo));
    DesireTrials.chr2.stim2 = find(((opto==7)&updated_desired_condition));
    DesireTrials.chr2.stim1 = find(((opto==5)&updated_desired_condition));
    DesireTrials.Nolight.currentcond = find(((opto<5)&updated_desired_condition));
    DesireTrials.Nolight.Chr2andHalorats = find(((opto<5)&updated_desired_condition_nolight_grouped));

    DesireTrials.externallight.stim2 = find(((opto==7)&updated_desired_condition_externallight));
    DesireTrials.externallight.stim1 = find(((opto==5)&updated_desired_condition_externallight));




    thisboots =DesireTrials.chr2.stim2(randi(numel(DesireTrials.chr2.stim2),1000,round(numel(DesireTrials.chr2.stim2)./2)));
    Res.OptoStim2.Chr2(j,:)=mean(Animal_Action_Left_Right(thisboots),2);

    thisboots =DesireTrials.chr2.stim1(randi(numel(DesireTrials.chr2.stim1),1000,round(numel(DesireTrials.chr2.stim1)./2)));
    Res.OptoStim1.Chr2(j,:)=mean(Animal_Action_Left_Right(thisboots),2);


    thisboots =DesireTrials.Halo.stim2(randi(numel(DesireTrials.Halo.stim2),1000,round(numel(DesireTrials.Halo.stim2)./2)));
    Res.OptoStim2.Halo(j,:)=mean(Animal_Action_Left_Right(thisboots),2);



    thisboots =DesireTrials.Halo.stim1(randi(numel(DesireTrials.Halo.stim1),1000,round(numel(DesireTrials.Halo.stim1)./2)));
    Res.OptoStim1.Halo(j,:)=mean(Animal_Action_Left_Right(thisboots),2);

    thisboots =  DesireTrials.Nolight.Chr2andHalorats(randi(numel(  DesireTrials.Nolight.Chr2andHalorats ),1000,round(numel(  DesireTrials.Nolight.Chr2andHalorats )./2)));
    Res.Nolight.Chr2andHalorats(j,:)=mean(Animal_Action_Left_Right(thisboots),2);


 


   
    
    
    thisboots =DesireTrials.externallight.stim1(randi(numel(DesireTrials.externallight.stim1),1000,round(numel(DesireTrials.externallight.stim1)./2)));
    Res.External.stim1(j,:)=mean(Animal_Action_Left_Right(thisboots),2);


    thisboots =DesireTrials.externallight.stim2(randi(numel(DesireTrials.externallight.stim2),1000,round(numel(DesireTrials.externallight.stim2)./2)));
    Res.External.stim2(j,:)=mean(Animal_Action_Left_Right(thisboots),2);



    %     Op{j}.Correct_stim2=mean(Pscyometric_dataromance(thisboots),2);
    %     Op{j}.ActualPercentT2_T1_stim2=mean(Animal_Action_Left_Right(thistrials),2);
    %     Op{j}.ActualPercentHaloT2_T1_stim2=mean(Animal_Action_Left_Right(thistrialshalo),2);
    %
    %     % Opto on stimulus1
    %     thistrials = find(((opto==5)&updated_desired_condition));
    %     thisboots = thistrials(randi(numel(thistrials),1000,round(numel(thistrials))));
    %     Op{j}.PercentT2_T1_stim1 = mean(Animal_Action_Left_Right(thisboots),2);
    %     Op{j}.Correct_stim1 = mean(Pscyometric_dataromance(thisboots),2);
    %     Op{j}.ActualPercentT2_T1_stim1=mean(Animal_Action_Left_Right(thistrials),2);
    %
    %     % no opto
    %     thistrials = find(((opto<5)&updated_desired_condition));
    %     thisboots =thistrials(randi(numel(thistrials),1000,round(numel(thistrials)./2)));
    %     noOp{j}.PercentT2_T1=mean(Animal_Action_Left_Right(thisboots),2);
    %     noOp{j}.Correct=mean(Pscyometric_dataromance(thisboots),2);
    %     noOp{j}.ActualPercentT2_T1=mean(Animal_Action_Left_Right(thistrials),2);
    %     % no opto all rats figure 3E
    %     thistrials = find(((opto<5)&updated_desired_condition_nolight_grouped));
    %     thisboots =thistrials(randi(numel(thistrials),1000,round(numel(thistrials)./2)));
    %     noOp{j}.PercentT2_T1_figure3E=mean(Animal_Action_Left_Right(thisboots),2);
    %     noOp{j}.Correct_figure3E=mean(Pscyometric_dataromance(thisboots),2);
    %     noOp{j}.ActualPercentT2_T1_figure3E=mean(Animal_Action_Left_Right(thistrials),2);
    %     noOp{j}.figure3Enumberoftrials=numel(thistrials);
    %     noOp{j}.ActualPercentT2_T1_figure3E=mean(Animal_Action_Left_Right(thistrials),2);
end


%%
pavlues.HalovsContorl=sum ((sum(Res.OptoStim2.Halo)-sum( Res.Nolight.Chr2andHalorats))>0)./1000;
pavlues.Chr2vsContorl=1-sum ((sum(Res.OptoStim2.Chr2)-sum(Res.Nolight.Chr2andHalorats))>0)./1000;
pavlues.Chr2stim2vsstim1=1-sum (sum(Res.OptoStim2.Chr2)-sum( Res.OptoStim1.Chr2)>0)./1000;
pavlues.stim2vsstim1External=1-sum (sum(Res.External.stim2)-sum(Res.External.stim1)>0)./1000;
pavlues.Halostim2vsstim1=sum((sum(Res.OptoStim2.Halo)-sum(Res.OptoStim1.Halo))>0)./1000;
% pavlues.chr2stim1vsContorl=1-sum ((sum(Res.OptoStim1.Chr2)-sum( Res.Nolight.Chr2andHalorats))>0)./1000;
% pavlues.HalovsContorl=sum (sum(Res.OptoStim2.Halo)-sum( Res.Nolight.Chr2andHalorats)>0)./1000;
pavlues.Ch2S2vsS1=1-sum((sum(Res.OptoStim2.Chr2)-sum(Res.OptoStim1.Chr2))>0)./1000;

%%
hold off
 plot(theseTDI,mean(Res.OptoStim1.Chr2'))
h= errorbar(theseTDI,mean(Res.OptoStim1.Halo'),std(Res.OptoStim1.Halo'))
hold on
h.Color= [0 0.8 0.1]
h2= errorbar(theseTDI,mean(Res.OptoStim2.Halo'),std(Res.OptoStim2.Halo'))
figurecorrectionOPTO(gca,[0.05 0.05],12)
h2.Color= [0 0.2 0.8]

