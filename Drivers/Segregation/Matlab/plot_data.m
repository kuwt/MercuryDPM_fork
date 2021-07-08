function plot_data()
addpath('..')
%load_data('new')
load_data();
clear global Pcolor
close all
print_data();
plot_nu(1);
%plot_nu(1);
% plot_nu(10);

plot_phi();
%plot_P();
plot_phi_topcut();

print_figures();
return


function print_data()
global stat small large Height Base
t=unique(cellfun(@(d) diff(d.time), stat));
disp(['averaging time intervall: ' num2str(t)]);
w=unique(cellfun(@(d) d.w, stat));
disp(['averaging width: ' num2str(w)]);
Height=cellfun(@(d) d.FlowHeight(1), stat);
disp(['flow height: ' num2str(Height)]);
Base=cellfun(@(d) d.Base(1), stat);
disp(['flow height: ' num2str(Base)]);
return

function plot_phi_topcut()
load_data([],'topcut');
plot_nu(1);
%plot_nu(6);
%plot_nu(10);
plot_phi();
set(gcf,'FileName',['phifitall_topcut'])

%plot_P();
load_data();
return

function plot_phi()
figure();
set(gcf,'Position',[0 0 560*2 420*.8])
set(gcf,'FileName',['phifitallshort'])
CaseList=[1];
for i=CaseList
  subplot(1,ceil(length(CaseList)/1),find(CaseList==i,1));
  P(i)=fit_phi(i,'subplot'); 
  set(gca,'xtick',0:.25:1,'ytick',0:.25:1);
  set(gca,'xtickLabel',[],'ytickLabel',[]);
  if i==CaseList(1)
      text(-0.22,0.5,'$\hat{z}_f$','Rotation',90,'VerticalAlignment','Middle','HorizontalAlignment','Center')
      set(gca,'ytickLabelMode', 'auto')
  end
  text(0.5,-0.15,'$\phi$','VerticalAlignment','Middle','HorizontalAlignment','Center')
  set(gca,'xtickLabelMode', 'auto')
  %grid on
  pos = get(gca,'Position');
  set(gca,'Box','on');
  set(gca,'Position',get(gca,'Position')-.05*[+1.75*pos(3) -2*pos(4) -3.5*pos(3) 2*pos(4)]);
end
return

function plot_P()
global Pcolor
if isempty(Pcolor)|Pcolor==[.5 .5 .5];
  Pcolor = [0 0 0];
  marker = 'x';
  DN='fit to full flow height';
else
  Pcolor = [.5 .5 .5];
  marker = 'd';
  DN='fit to dense base layer';
end  
for i=1:10, [P(i) Pinf(i,:)]=fit_phi(i,'nofigure'); end

figure(21);
set(gcf,'Position',[0 0 560 420])
set(gcf,'FileName','P')
%clf
hold on
sigmainv = (0.0033:0.0003:0.0060)/0.003;
%plot(sigmainv,P,'ko');
errorbar(sigmainv,P,P-Pinf(:,1)',P-Pinf(:,2)',marker,'Color',Pcolor,'DisplayName',DN)
ylabel('$P_s$')
xlabel('$\sigma^{-1}$')
if Pcolor==[.5 .5 .5];
  C=get(gca,'Children');
  legend(C([3 1]),'Location','SouthEast'); 
end
% 'coefficients' The coefficient names. Use a cell array if there are multiple names. The following names are not allowed: i, j, pi, inf, nan, eps.
% 'dependent'    The dependent (response) variable name
% 'independent'  The independent (predictor) variable name
% 'options'      The default fit options for the object
% 'problem'      The problem-dependent (fixed) parameter names. Use a cell array if there are multiple names. The default is none.

%MONOD fit: half saturation constant K, saturation value Pmax
% s=sigmainv-1;
% f = fittype('rat01');
% c = fit(s',(P./s)',f, 'StartPoint',[10 1] );
% coeff = coeffvalues(c);
% Pmax = coeff(1); K=coeff(2);
% sinv=1:.05:2;
% plot(sinv,Pmax*(sinv-1)./(sinv-1+K),'Color',Pcolor);
%f = fittype('Pmax*(sigmainv-1)/(K+sigmainv-1)','coefficients',{'Pmax';'K'},'independent','sigmainv');
f = fittype('Pmax*(1-exp(-K*(sigmainv-1)))','coefficients',{'Pmax';'K'},'independent','sigmainv');
c = fit(sigmainv',P',f,'StartPoint',[7 5] );
coeff = coeffvalues(c);
Pmax = coeff(1); 
K=coeff(2);
sinv=1:.02:2;
%plot(sinv,Pmax*(sinv-1)./(sinv-1+K),'Color',Pcolor);
plot(sinv,Pmax*(1-exp(-K*(sinv-1))),'Color',Pcolor);
disp(['coefficients Pmax=' num2str(coeff(1)) ' , K=' num2str(coeff(2))])
% confidence interval
%var = predint(c,sinv);
%plot(sinv,var(:,1),':','Color',Pcolor);
%plot(sinv,var(:,2),':','Color',Pcolor);

% another plotting method
%plot(c,sigmainv,P);


axis tight
return


function [P, Pinf]=fit_phi(Case,opt)
global stat small 
sigma=1.1:.1:2;

phi = small{Case}.VolumeFraction./stat{Case}.VolumeFraction;
z=(stat{Case}.z-stat{Case}.Base)/stat{Case}.FlowHeight;

ind=z>0&z<1; 
z=z(ind); 
phi=phi(ind);

phim = mean(phi);
phifit = @(P) (1-exp(-phim*P))*exp((phim-z)*P)./...
  (1-exp(-(1-phim)*P)+(1-exp(-phim*P))*exp((phim-z)*P));

s = fitoptions('Method','NonlinearLeastSquares',...
  'Lower',0,...
  'Upper',Inf,...
  'Startpoint',1);

f = fittype('(1-exp(-phim*P))*exp((phim-z)*P)./(1-exp(-(1-phim)*P)+(1-exp(-phim*P))*exp((phim-z)*P))','options',s,'coefficients','P','independent','z','problem','phim');
c = fit(z,phi,f,'phim',phim);
P = coeffvalues(c);
Pinf = confint(c);

if ~exist('opt','var')
  figure(10+Case);
  set(gcf,'Position',[0 0 560 420])
  set(gcf,'FileName',['phifit' num2str(Case,'_%d')])

  hold on

  plot(phifit(P),z,'DisplayName',['$\phi_{fit}(' num2str(P,'%.2f') ')$'])
  plot(phifit(Pinf(1)),z,':')
  plot(phifit(Pinf(2)),z,':')

  legend('show')
  axis tight
elseif strcmp(opt,'subplot')
  hold on

  plot(phi,z,'k','DisplayName','$\phi$')
  plot(phifit(P),z,'DisplayName',['$\phi_{fit}(' num2str(P,'%.2f') ')$'])
%   legend('show')
   plot(phifit(Pinf(1)),z,':')
   plot(phifit(Pinf(2)),z,':')

 text(.9,.9,['$\sigma^{-1}=' num2str(sigma(Case)) '$'],'VerticalAlignment','Top','HorizontalAlignment','Right')
 %xlabel('$\phi$','Position',[.5 -.01 1])
 %ylabel('$z$')
  set(gca,'YTick',[0 .5 1])
  set(gca,'XTick',[0 .5 1])
%   set(gca,'YTickLabel',[0 1])
%   set(gca,'XTickLabel',[0 1])
  %axis tight
end

return

% function plot_phi(Case)
% global stat small 
% 
% figure(10+Case);
% set(gcf,'Position',[0 0 560 420])
% set(gcf,'FileName',['phi' num2str(Case,'_%d')])
% clf
% hold on
% 
% phi = small{Case}.VolumeFraction./stat{Case}.VolumeFraction;
% z=(stat{Case}.z-stat{Case}.Base)/stat{Case}.FlowHeight;
% 
% ind=z>0&z<1; 
% z=z(ind); 
% phi=phi(ind);
% 
% plot(phi,z,'k','DisplayName','$\phi$')
% 
% title(stat{Case}.name)
% xlabel('$\phi^s$')
% ylabel('$z$')
% axis([0 1 0 1]);
% 
% return

function plot_nu(Case)
global stat small large

figure(Case);
set(gcf,'Position',[0 0 560 420])
set(gcf,'FileName',['nu' num2str(Case,'_%d')])
clf
hold on

kappa=0.02;
stat{Case}.StressZZ=sum(stat{Case}.VolumeFraction)-cumsum(stat{Case}.VolumeFraction);
stat{Case}.Surface2=max(stat{Case}.z.*(stat{Case}.StressZZ>kappa*stat{Case}.StressZZ(1)),[],1);


plot(stat{Case}.VolumeFraction,stat{Case}.z/stat{Case}.Surface2,'k',...
  'DisplayName','$\Phi$')
plot(small{Case}.VolumeFraction,small{Case}.z/stat{Case}.Surface2,'r',...
  'DisplayName','$\Phi^s$')
plot(large{Case}.VolumeFraction,large{Case}.z/stat{Case}.Surface2,'--b',...
  'DisplayName','$\Phi^l$')


%title(stat{Case}.name)
legend('show')
set(legend,'Location','East')
ylabel('$\hat{z}$')
xlim([0 max(xlim)]);
plot(xlim,[1 1]*stat{Case}.Surface/stat{Case}.Surface2,':','Color',.5*[1 1 1]);
ylim([0 1]);


return

function load_data(opt,opt2)

if exist('opt','var')
  if (strcmp(opt,'new'))
    clear global
    
  end
end
  
% where we store the data
%clear global
global stat small large

% get load_statistics
if ~exist('loadStatistics','file')
  addpath('../../../Matlab/');
  addpath('..');
end




stat=loadStatistics({'Segregation.1.stat'});
small=loadStatistics({'Segregation.1.small.stat'});
large = loadStatistics({'Segregation.1.large.stat'});

%if isempty(stat ), 
%     stat  = loadStatistics({...
%         '../final_time/segregation.1.stat',...
%         '../final_time/segregation.2.stat',...
%         '../final_time/segregation.3.stat',...
%         '../final_time/segregation.4.stat',...
%         '../final_time/segregation.5.stat',...
%         '../final_time/segregation.6.stat',...
%         '../final_time/segregation.7.stat',...
%         '../final_time/segregation.8.stat',...
%         '../final_time/segregation.9.stat',...
%         '../final_time/segregation.10.stat',...
%         }); 
%     small  = loadstatistics({...
%         '../final_time/segregation.1.small.stat',...
%         '../final_time/segregation.2.small.stat',...
%         '../final_time/segregation.3.small.stat',...
%         '../final_time/segregation.4.small.stat',...
%         '../final_time/segregation.5.small.stat',...
%         '../final_time/segregation.6.small.stat',...
%         '../final_time/segregation.7.small.stat',...
%         '../final_time/segregation.8.small.stat',...
%         '../final_time/segregation.9.small.stat',...
%         '../final_time/segregation.10.small.stat',...
%         });
%     large  = loadstatistics({...
%         '../final_time/segregation.1.large.stat',...
%         '../final_time/segregation.2.large.stat',...
%         '../final_time/segregation.3.large.stat',...
%         '../final_time/segregation.4.large.stat',...
%         '../final_time/segregation.5.large.stat',...
%         '../final_time/segregation.6.large.stat',...
%         '../final_time/segregation.7.large.stat',...
%         '../final_time/segregation.8.large.stat',...
%         '../final_time/segregation.9.large.stat',...
%         '../final_time/segregation.10.large.stat',...
%         }); 
%    save stat.mat stat small large
    %load stat.mat
%end


for Case=1:length(stat)
  if (false)
    stat{Case}.Base=min(stat{Case}.z(stat{Case}.VolumeFraction>0.5*max(stat{Case}.VolumeFraction)));
    stat{Case}.Surface=max(stat{Case}.z(stat{Case}.VolumeFraction>0.5*max(stat{Case}.VolumeFraction)));
  else
    kappa=0.02;
    stat{Case}.StressZZ=sum(stat{Case}.VolumeFraction)-cumsum(stat{Case}.VolumeFraction);
    stat{Case}.Base=min(stat{Case}.z(stat{Case}.StressZZ<(1-kappa)*stat{Case}.StressZZ(1)),[],1);
    if (~exist('opt','var'))
        stat{Case}.Surface=max(stat{Case}.z.*(stat{Case}.StressZZ>kappa*stat{Case}.StressZZ(1)),[],1);
    else %cut off pure region at top
%         stat{Case}.Surface=...
%           min(stat{Case}.z(small{Case}.VolumeFraction./stat{Case}.VolumeFraction<min(small{Case}.VolumeFraction./stat{Case}.VolumeFraction)+0.02),[],1);
         stat{Case}.Surface=...
           max(stat{Case}.z(stat{Case}.VolumeFraction>median(stat{Case}.VolumeFraction)));
    end
  end
  stat{Case}.FlowHeight=stat{Case}.Surface-stat{Case}.Base;
end

return


