function steady()
%first run plot_data.m in ../final_time to get height and base values
clc
close all
global Height Base
CaseList = [1];

close all
figure(1); clf; hold on
set(gcf,'Position',[0 680 2*560 420])
set(gcf,'FileName',['steady1'])
sigma=1.1:.1:2;

for Case=CaseList
    com=importdata(['segregation.' num2str(Case) '.com']);
    ix=1:size(com,1);
    t=com(ix,1)';
    COMX=(com(ix,2)'-Base(Case))/Height(Case);
    COMY=(com(ix,3)'-Base(Case))/Height(Case);
    COMZ=(com(ix,4)'-Base(Case))/Height(Case);
    SmallCOMX=(com(ix,5)'-Base(Case))/Height(Case);
    SmallCOMY=(com(ix,6)'-Base(Case))/Height(Case);
    SmallCOMZ=(com(ix,7)'-Base(Case))/Height(Case);
    LargeCOMX=(com(ix,8)'-Base(Case))/Height(Case);
    LargeCOMY=(com(ix,9)'-Base(Case))/Height(Case);
    LargeCOMZ=(com(ix,10)'-Base(Case))/Height(Case);
    
    subplot(1,length(CaseList),find(CaseList==Case))
    %Time resets in this code, so you want the second time it starts
    tstart=max(find(t<1));
    plot(t(tstart:end),LargeCOMZ(tstart:end),'b-','DisplayName','COM Large');
    plot(t(tstart:end),COMZ(tstart:end),'k-','DisplayName','COM');
    plot(t(tstart:end),SmallCOMZ(tstart:end),'r-','DisplayName','COM Small');
    
   % plot(t(tstart:end),[LargeCOMZ(tstart:end); COMZ(tstart:end); SmallCOMZ(tstart:end)],'DisplayName',{'COM large','COM total','COM small'});
    xlabel('$t$')
    %ylim([.25,.75])
    %xlim([0,100])
    title(['$\sigma=' num2str(sigma(Case)) '$'])
    if (Case==CaseList(end))
        legend('show')
        set(legend,'Location','East')
    end
    
    %set(gca,'xtick',0:20:100,'ytick',0.25:.05:.75);
    set(gca,'ytickLabel',[]);
    if Case==CaseList(1)
      set(gca,'ytickLabelMode', 'auto')
      ylabel('$\hat{z}_f$','Rotation',90)
    end
     %grid on
     pos = get(gca,'Position');
     set(gca,'Box','on');
     set(gca,'Position',get(gca,'Position')-.05*[+1.75*pos(3) -2*pos(4) -3.5*pos(3) 2*pos(4)]);
end

print_figures()
return
