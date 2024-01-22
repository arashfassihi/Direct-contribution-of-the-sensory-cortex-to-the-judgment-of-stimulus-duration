
clear Spikecount
clear inds
istim = 2
clear a
clear c
counter=0
clear spn
for j=1:size(STIM1,2)
    
    for i=7
        %              if istim
        %                  i =4;
        %              end
        if ~isempty(STIM1{i,j})&&~isempty(STIM2{i,j})
            
            counter = counter+1;
            b = STIM1{i,j};
            a(counter,: )= mean(b.PSTH{istim}(:,1000:1000+D(i)));
            b = STIM2{i,j};
            c(counter,:) = mean(b.PSTH{istim}(:,1000:1000+D(i)));
            spn(counter) = mean(b.ratesp{2});
            
            
        end
    end
end


[rr ,cc ] =find(isnan(c));

%


%%

optocumsum =mean(cumsum(c(spn>2,:)'),2);
controlcumsum =mean(cumsum(a(spn>2,:)'),2);


%%
close all

plot(optocumsum,controlcumsum,'.')


[Intersectvalue, indopto,indcontrol]=intersect(round(optocumsum*1),round(controlcumsum*1))
colorcode = cool (numel(Intersectvalue));
for i=1:numel(Intersectvalue)
    
    if Intersectvalue(i)==20
        text ( indcontrol(i)+40,indopto(i)+10,num2str([indcontrol(i);indopto(i)]))
         plot(indcontrol(i),indopto(i),'.','color',[0 1 0],'markersize',25)
        hold on
        plot(indcontrol(i),indopto(i),'o','color',[0 1 0],'markersize',7,'markeredgecolor','k')
     
    else
          plot(indcontrol(i),indopto(i),'.','color',colorcode(i,:),'markersize',25)
        hold on
        plot(indcontrol(i),indopto(i),'o','color',colorcode(i,:),'markersize',7,'markeredgecolor','k')
    end
end

hold on

plot([0 700],[0 700],'k')

set(gca,'FontSize',15)
set(gca,'TickDir', 'out','TickLength',[0.05, 0.01])
box off
%%
figure
for i=1:30
        if Intersectvalue(i)==20
                plot([0 1],[i i]-1,'Linewidth',15,'color',[ 0 1 0])

        else
    plot([0 1],[i i]-1,'Linewidth',15,'color',colorcode(i,:))
    hold on
        end
end
box off
set(gca, 'Xticklabel',[])

