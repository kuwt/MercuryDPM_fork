clear all; close all;
disp('running chute_periodic_stat.m ...');

theta_maxstop = [];
h_maxstop = [];
% load data from file to variable stat
filename = dir(['*.*.stat']);
num_per_figure = 36;
for fignum =1:ceil(length(filename)/num_per_figure)
  h = figure(fignum);
  clf(h);
  set(h,'Position',[1 1 1679*.7 939*.7]);
  directory = cd;
  theta_stop=[];
  h_stop=[];
  d = 0.72e-3;
  ilist = ((fignum-1)*num_per_figure + 1) : min(length(filename), fignum*num_per_figure)
  for i=ilist
    disp(filename(i).name);
    % gather data
    stat = importdata(filename(i).name,'\t',1);
    if (false)%(size(stat.data,1)>=11)
      stat.data = (stat.data(1:end-10,:)+stat.data(2:end-9,:)+stat.data(3:end-8,:)+stat.data(4:end-7,:)+stat.data(5:end-6,:)+stat.data(6:end-5,:) ...
        +stat.data(7:end-4,:)+stat.data(8:end-3,:)+stat.data(9:end-2,:)+stat.data(10:end-1,:)+stat.data(11:end,:))/11;
    end
    t = stat.data(:,1);
    ekin = stat.data(:,2);
    h = stat.data(:,end);
    theta = stat.data(:,end-1);
    n_max = size(stat.data,1);

    %plot data
    subplot(ceil(sqrt(length(ilist))),ceil(sqrt(length(ilist))),mod(i-1,num_per_figure)+1);
    if (n_max<10) linespec='-x'; else linespec='-'; end
    plot(t,log(ekin)/log(10),linespec);
    ylabel(['log(' stat.colheaders{2} ')'],'Interpreter','none');
    xlabel(stat.colheaders{1},'Interpreter','none');
    title(filename(i).name(1:end-5));
    axis([0 6 -15 0]);
    
    % gather data for figure 2 and title
    stop = min(find(ekin<1e-11));
    if (stop)
      theta_stop = [theta_stop theta(stop)];
      h_stop = [h_stop h(end)-d/2];
      hold on; plot(t(stop),log(ekin(stop))/log(10),'ro'); hold off
      title([filename(i).name(1:end-5) ' th=' num2str(theta_stop(end)) ...
        ',hs/d=' num2str(h_stop(end)/d) ],'Interpreter','none');
    end
  end

  figure(ceil(length(filename)/num_per_figure)+1);
  hold on;
  plot(theta_stop,h_stop/d,'o');
  %plot(theta_init,2 COM_z/d,'-or');
  title('\circ denotes angles and heights for which the flow stops')
  xlabel('\theta');
  ylabel('h_{stop}/d');
  if (theta_stop)
    theta_maxstop(end+1) = theta_stop(end);
    h_maxstop(end+1) = h_stop(end);
  end

end

d=0.72e-3;
figure(fignum+1);
hold on;
plot(theta_maxstop,h_maxstop/d,':');
coeff = findfit(theta_maxstop(3:end), h_maxstop(3:end)/d);
title(['A=' num2str(coeff(1)) ', \theta_1=' num2str(atan(coeff(2))/pi*180) ', \theta_2=' num2str(atan(coeff(3))/pi*180) ])
set(figure(fignum+1),'PaperType','A4');