istim = 2
counter=0
clear a
clear c
clear spn
cd 'D:\Dropbox\OptoDurationPaper\A_NATCOMSUBMISSION\related to figure3\'
load STIM1
load STIM2
load STIM2high.mat
load STIM2low
figure
for j=1:size(STIM1,2)
    
     i=6
        %              if istim
        %                  i =4;
        %              end
        if ~isempty(STIM1{i,j})&&~isempty(STIM2{i,j})
            
            counter = counter+1;
            b = STIM1{i,j};
            a(j )= mean(b.count{istim});
            b = STIM2{i,j};
            c(j) = mean(b.count{istim});
            spn(j) = mean(b.ratesp{istim});
            
            
        else
            spn(i,j) =NaN;
            c(i,j) =NaN;
            a(i,j )= NaN;
        end
end
%     end
% hold off
% 
% plot(a(spn>2),c(spn>2),'.k','Markersize',12)
% hold on
% plot([1 200] ,[1 200],'--','color',[0.4 0.4 0.4])
% 
% axis square
%     set(gca,'FontSize',15)
% set(gca,'TickDir', 'out','TickLength',[0.05, 0.01])
% box off
% hold on
% 
% plot(mean(a(spn>2)),mean(c(spn>2)),'+r','markersize',20,'Linewidth',4)

%%


% X = -40:20:120
% bar(X,hist(c(spn>2)-a(spn>2),X))



%%
figure

thiscolor = jet(7)
for      i=1:7;

for permi=1:200

istim = 2;
counter=0;
clear a
clear c
clear spn
for j=1:size(STIM1,2)
    
        %              if istim
        %                  i =4;
        %              end
        if ~isempty(STIM1{i,j})&&~isempty(STIM2{i,j})
            
            counter = counter+1;
            b = STIM1{i,j};
              thiscount = b.count{istim};
            a(j )= thiscount(randi(numel(thiscount)));
            b = STIM2{i,j};
            thiscount = b.count{istim};
            c(j) = thiscount(randi(numel(thiscount)));
            spn(j) = mean(b.ratesp{istim});
            
            
        else
            spn(i,j) =NaN;
            c(i,j) =NaN;
            a(i,j )= NaN;
        end
    end
thesevalues = find(spn>2);
subplot (2,2,3)
thesevalues = thesevalues(randi(numel(thesevalues),1,numel(thesevalues)));

if permi>100
plot(log(mean(a(thesevalues))),log(mean(c(thesevalues))),'.','color',thiscolor(i,:),'markersize',15)
hold on
end
meanspikecount.control(i,permi)=mean(a(thesevalues));
meanspikecount.opto(i,permi)=mean(c(thesevalues));

end
end
%%
subplot (2,2,1)
    hold off
subplot (2,2,4)
    hold off

for i=1:7
[f,xi] = ksdensity(log(meanspikecount.control(i,:))); 
subplot (2,2,1)
plot(xi,f,'color',thiscolor(i,:),'markersize',10)
hold on
[f,xi] = ksdensity(log(meanspikecount.opto(i,:))); 
subplot (2,2,4)
plot(f,xi,'color',thiscolor(i,:),'markersize',10)
hold on
end
subplot (2,2,3)

plot([1 200] ,[1 200],'--','color',[0.4 0.4 0.4])

axis square
    set(gca,'FontSize',15)
set(gca,'TickDir', 'out','TickLength',[0.05, 0.01])
box off
hold on


limit_axis = [5 50];
subplot (2,2,3)
% set(gca, 'yScale', 'log')
% set(gca, 'xScale', 'log')
xlim(log([ limit_axis]))
ylim(log([ limit_axis]))
set(gca,'xtick',log(10:10:50))
set(gca,'ytick',log(10:10:50))
set(gca,'xticklabel',(10:10:50))
set(gca,'yticklabel',(10:10:50))
subplot (2,2,1)
    hold off
    set(gca , 'Visible','off')
xlim(log([limit_axis]))
axis square

subplot (2,2,4)
    hold off
set(gca , 'Visible','off')
ylim(log([ limit_axis]))
axis square


figure
for i=1:7
[f,xi] = ksdensity((meanspikecount.control(i,:))-(meanspikecount.opto(i,:))); 


plot(xi,f,'color',thiscolor(i,:),'markersize',10)
hold on


end
%%
plotmeanpoints=0;
if plotmeanpoints ==1
for i=1:7
plot([1 mean(meanspikecount.control(i,:))],[mean(meanspikecount.opto(i,:)) mean(meanspikecount.opto(i,:))],'--','color',thiscolor(i,:))
plot([mean(meanspikecount.control(i,:)) mean(meanspikecount.control(i,:))],[1 mean(meanspikecount.opto(i,:))],'--','color',thiscolor(i,:))


end
end
%%
D = [    161   205   264   334   422   545   694  ];

figure(2)
thiscolor = jet(7)
for      i=1:7;

for permi=1:200

istim = 2;
counter=0;
clear a
clear c
clear spn
for j=1:size(STIM1,2)
    
        %              if istim
        %                  i =4;
        %              end
        if ~isempty(STIM1{i,j})&&~isempty(STIM2{i,j})
            
            counter = counter+1;
            b = STIM1{i,j};
              thiscount = sum((b.PSTH{2}(:,1000+D(i)-100:1000+D(i))),2)*10;
            a(j )= thiscount(randi(numel(thiscount)));
            b = STIM2{i,j};
            thiscount =sum((b.PSTH{2}(:,1000+D(i)-100:1000+D(i))),2)*10;
            c(j) = thiscount(randi(numel(thiscount)));
            spn(j) = mean(b.ratesp{istim});
            
            
        else
            spn(i,j) =NaN;
            c(i,j) =NaN;
            a(i,j )= NaN;
        end
    end
thesevalues = find(spn>1);
subplot (2,2,3)

thesevalues = thesevalues(randi(numel(thesevalues),1,numel(thesevalues)));
if permi>100

plot(log(mean(a(thesevalues))),log(mean(c(thesevalues))),'.','color',thiscolor(i,:),'markersize',10)

hold on
end




hold on
meanfiringrate.control(i,permi)=mean(a(thesevalues));
meanfiringrate.opto(i,permi)=mean(c(thesevalues));
% meanfiringrate.allcontrol{i,permi)=mean(a(thesevalues));
% meanfiringrate.allopto{i,permi}=mean(c(thesevalues));
end
end
%%
subplot (2,2,1)
    hold off
subplot (2,2,4)
    hold off

for i=1:7
[f,xi] = ksdensity(log(meanfiringrate.control(i,:))); 
subplot (2,2,1)
plot(xi,f,'color',thiscolor(i,:),'markersize',10)
hold on
[f,xi] = ksdensity(log(meanfiringrate.opto(i,:))); 
subplot (2,2,4)
plot(f,xi,'color',thiscolor(i,:),'markersize',10)
hold on
end
subplot (2,2,3)

plot([1 200] ,[1 200],'--','color',[0.4 0.4 0.4])

axis square
    set(gca,'FontSize',15)
set(gca,'TickDir', 'out','TickLength',[0.05, 0.01])
box off
hold on


limit_axis = [30 80];
subplot (2,2,3)
% set(gca, 'yScale', 'log')
% set(gca, 'xScale', 'log')
xlim(log([ limit_axis]))
ylim(log([ limit_axis]))
set(gca,'xtick',log(10:20:80))
set(gca,'ytick',log(10:20:80))
set(gca,'xticklabel',(10:20:80))
set(gca,'yticklabel',(10:20:80))
subplot (2,2,1)
    hold off
    set(gca , 'Visible','off')
xlim(log([limit_axis]))
axis square

subplot (2,2,4)
    hold off
set(gca , 'Visible','off')
ylim(log([ limit_axis]))
axis square

%%
figure
for i=1:7
plot([0 mean(meanfiringrate.control(i,:))],[mean(meanfiringrate.opto(i,:)) mean(meanfiringrate.opto(i,:))],'--','color',thiscolor(i,:))
plot([mean(meanfiringrate.control(i,:)) mean(meanfiringrate.control(i,:))],[0 mean(meanfiringrate.opto(i,:))],'--','color',thiscolor(i,:))


end
%%
figure
for i=1:7
[f,xi] = ksdensity((meanfiringrate.control(i,:))-(meanfiringrate.opto(i,:))); 


plot(xi,f,'color',thiscolor(i,:),'markersize',10)
hold on


end
%% for differnt speed


%%
figure

thiscolor = jet(7);
for      i=1:7;

for permi=1:100

istim = 2;
counter=0;
clear a
clear c
clear spn
for j=1:size(STIM2high,2)
    
        %              if istim
        %                  i =4;
        %              end
        if ~isempty(STIM2{i,j})&&~isempty(STIM2high{i,j})&&~isempty(STIM2low{i,j})
            
            counter = counter+1;
            b = STIM2low{i,j};
              thiscount = b.count{istim};
            a(j )= thiscount(randi(numel(thiscount)));
            b = STIM2high{i,j};
            thiscount = b.count{istim};
            c(j) = thiscount(randi(numel(thiscount)));
                        b = STIM2{i,j};

            spn(j) = mean(b.ratesp{istim});
            
            
        else
            spn(i,j) =NaN;
            c(i,j) =NaN;
            a(i,j )= NaN;
        end
    end
thesevalues = find(spn>0);
thesevalues = thesevalues(randi(numel(thesevalues),1,numel(thesevalues)));
plot(mean(a(thesevalues)),mean(c(thesevalues)),'.','color',thiscolor(i,:),'markersize',10)
hold on
meanspikecount.control(i,permi)=mean(a(thesevalues));
meanspikecount.opto(i,permi)=mean(c(thesevalues));

end
end
plot([1 200] ,[1 200],'--','color',[0.4 0.4 0.4])

axis square
    set(gca,'FontSize',15)
set(gca,'TickDir', 'out','TickLength',[0.05, 0.01])
box off
hold on
%%

for i=1:7
plot([1 mean(meanspikecount.control(i,:))],[mean(meanspikecount.opto(i,:)) mean(meanspikecount.opto(i,:))],'--','color',thiscolor(i,:))
hold on

plot([mean(meanspikecount.control(i,:)) mean(meanspikecount.control(i,:))],[1 mean(meanspikecount.opto(i,:))],'--','color',thiscolor(i,:))

end

 set(gca, 'xScale', 'log')
set(gca, 'yScale', 'log')
xlim([5 45])
ylim([5 45])