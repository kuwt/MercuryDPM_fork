% stat file is assumed to have the following structure
% headerline: denotes t, nx, ny, nz, and names of variables
% other lines: lists t, nx, ny, nz, and variables
clear; close all;
disp('running chute_evals_stats.m ...');

% load data from file to variable stat
filename = dir('*.stat')
for k=1:length(filename)
  disp(['... loading data from file' filename(k).name ' (creating filename, stat, header, header_double)']);
  stat = importdata(filename(k).name,' ',5);
  name = stat.textdata{1};
  CG = stat.textdata{3};
  n = cell2mat(textscan(stat.textdata{2},'%f'));
  ntotal = prod(max(1,n));
  xyz = zeros(3,ntotal);
  xyz(:) = cell2mat(textscan(stat.textdata{4},'%f '));
  titles = textscan(stat.textdata{5},'%s ');
  titles = titles{1};
  numvar = round((size(stat.data,2)-1)/ntotal);
  t = stat.data(:,1);

  % average data over N_avg time steps
  disp(['... average data over N_avg time steps (creating N_avg, N_total, prefactor, stat.avg)']);
  N_total = size(stat.data,1);
  N_avg = 1;%ceil(.1*N_total);
  stat.avg = stat.data(:,2:end);
  for i=1:N_avg-1
    stat.avg = stat.avg + [zeros(i,ntotal*numvar); stat.data(1:end-i,2:end)];
  end
  prefactor = [1:N_avg N_avg*ones(1,length(t)-N_avg)]'*ones(1,ntotal*numvar);
  stat.avg = stat.avg ./ prefactor;

  % set subplot arrangment
  disp(['... set subplot arrangment (creating nx, ny)']);
  if (numvar<=2) ny = 1;
  elseif (numvar<=6) ny = 2;
  elseif (numvar<=12) ny = 3;
  else ny = 4;
  end
  nx = ceil(numvar/ny);

  % get dim
  disp(['... get dim (creating dim)']);
  dim = 0;
  if (n(1)>1)	dim = dim+1; end
  if (n(2)>1)	dim = dim+1; end
  if (n(3)>1)	dim = dim+1; end

  % plot statistics
  if (true)
    disp(['... plot statistics(x,y) for last timestep in figure 1']);
    label = 'xyz';
    label = label(find(n));
    plot_i = 1:numvar;
    plot_i = [1 2 3 4 5 9 13 6 7 8 11 12 16];
    h = figure(10*k+1);
    set(h,'Position',[1 1 1679 939],'PaperType','A4','PaperOrientation','landscape','PaperPosition',[0 0 11.69 8.26]);
    for i=1:numvar, handle(i) = subplot(ny,nx,plot_i(i)); end 
    htextbox = annotation('textbox',[0 0 1 1]);
    set(htextbox,'String',[name ', ' CG ', n ' num2str(n')] ,'Interpreter','none')
    
    switch dim
    case 1
      X = xyz(find(n>1),:);
      for i=1:numvar
        subplot(ny,nx,plot_i(i));
        bar(X,stat.avg(end,i:numvar:end),1);
        title(titles{i});
        xlabel(label);
        axis('tight');
      end
      saveas(h ,['figure' num2str(h ) '.pdf']);
    case 2
      %for x-averaging
      hX = figure(10*k+2);
      set(hX,'Position',[1 1 1679 939],'PaperType','A4','PaperOrientation','landscape','PaperPosition',[0 0 11.69 8.26]);
      for i=1:numvar, handleX(i) = subplot(ny,nx,plot_i(i)); end 
      htextboxX = annotation('textbox',[0 0 1 1]);
      set(htextboxX,'String',[get(htextbox,'String') ': This graph is a contourplot of the average along x'],'Interpreter','none')
      %for y-averaging
      hY = figure(10*k+3);
      set(hY,'Position',[1 1 1679 939],'PaperType','A4','PaperOrientation','landscape','PaperPosition',[0 0 11.69 8.26]);
      for i=1:numvar, handleY(i) = subplot(ny,nx,plot_i(i)); end 
      htextboxY = annotation('textbox',[0 0 1 1]);
      set(htextboxY,'String',[get(htextbox,'String') ': This graph is a contourplot of the average along y'],'Interpreter','none')
      figure(10*k+1); 
      
      m = find(n>1);
      X=zeros(n(m)')';
      Y=zeros(n(m)')';
      Z=zeros(n(m)')';
      X(:)=xyz(m(1),:);
      Y(:)=xyz(m(2),:);
      for i=1:numvar
        subplot(ny,nx,plot_i(i));
        Z(:) = stat.avg(end,i:numvar:end);
        contourf(X,Y,Z);
        view(0,90);
        axis('tight');
        title(titles{i});
        colorbar;
        xlabel(label(1));
        ylabel(label(2));
        %x-averaging
        bar(handleX(i),X(1,:),sum(Z,1)/size(X,2),1);
        title(handleX(i),titles{i});
        xlabel(handleX(i),label(1));
        axis(handleX(i),'tight');
        %y-averaging
        bar(handleY(i),Y(:,1),sum(Z,2)/size(X,1),1);
        title(handleY(i),titles{i});
        xlabel(handleY(i),label(2));
        axis(handleY(i),'tight');
      end
      saveas(h ,['figure' num2str(h ) '.pdf']);
      saveas(hX,['figure' num2str(hX) '.pdf']);
      saveas(hY,['figure' num2str(hY) '.pdf']);
    case 3
      %for x-averaging
      hX = figure(10*k+2);
      set(hX,'Position',[1 1 1679 939],'PaperType','A4','PaperOrientation','landscape','PaperPosition',[0 0 11.69 8.26]);
      htextboxX = annotation('textbox',[0 0 1 1]);
      set(htextboxX,'String',[get(htextbox,'String') ': This graph is a contourplot of the average along x'],'Interpreter','none')
      for i=1:numvar, handleX(i) = subplot(ny,nx,plot_i(i)); end 
      %for y-averaging
      hY = figure(10*k+3);
      set(hY,'Position',[1 1 1679 939],'PaperType','A4','PaperOrientation','landscape','PaperPosition',[0 0 11.69 8.26]);
      htextboxY = annotation('textbox',[0 0 1 1]);
      set(htextboxY,'String',[get(htextbox,'String') ': This graph is a contourplot of the average along y'],'Interpreter','none')
      for i=1:numvar, handleY(i) = subplot(ny,nx,plot_i(i)); end 
      %for z-averaging
      hZ = figure(10*k+4);
      set(hZ,'Position',[1 1 1679 939],'PaperType','A4','PaperOrientation','landscape','PaperPosition',[0 0 11.69 8.26]);
      htextboxZ = annotation('textbox',[0 0 1 1]);
      set(htextboxZ,'String',[get(htextbox,'String') ': This graph is a contourplot of the average along z'],'Interpreter','none')
      for i=1:numvar, handleZ(i) = subplot(ny,nx,plot_i(i)); end 
      figure(10*k+1); 

      X=zeros(n(end:-1:1)');
      Y=zeros(n(end:-1:1)');
      Z=zeros(n(end:-1:1)');
      C=zeros(n(end:-1:1)');
      X(:)=xyz(1,:);
      Y(:)=xyz(2,:);
      Z(:)=xyz(3,:);
      yval = floor(n(2)/2);
      set(htextbox,'String',[get(htextbox,'String') ': This graph is a contourplot of a cross-section along y=' num2str(Y(1,yval,1))],'Interpreter','none')
      for i=1:numvar
        subplot(ny,nx,plot_i(i));
        C(:) = stat.avg(end,i:numvar:end)';
        contourf(permute(X(:,yval,:),[1 3 2]),permute(Z(:,yval,:),[1 3 2]),permute(C(:,yval,:),[1 3 2]),permute(C(:,yval,:),[1 3 2]));%,'EdgeAlpha',0);
        view(0,90);
        axis('tight');
        title(titles{i});
        colorbar;
        xlabel(label(1));
        ylabel(label(3));

        %x-averaging
        contourf(handleX(i),Y(:,:,1),Z(:,:,1),sum(C,3)/size(X,3));
        axis(handleX(i),'tight');
        title(handleX(i),titles{i});
        xlabel(handleX(i),label(2));
        ylabel(handleX(i),label(3));
        %y-averaging
        contourf(handleY(i),permute(X(:,1,:),[1 3 2]),permute(Z(:,1,:),[1 3 2]),permute(sum(C,2)/size(X,2),[1 3 2]));
        axis(handleY(i),'tight');
        title(handleY(i),titles{i});
        xlabel(handleY(i),label(1));
        ylabel(handleY(i),label(3));
        %z-averaging
        contourf(handleZ(i),permute(X(1,:,:),[2 3 1]),permute(Y(1,:,:),[2 3 1]),permute(sum(C,1)/size(X,1),[2 3 1]));
        axis(handleZ(i),'tight');
        title(handleZ(i),titles{i});
        xlabel(handleZ(i),label(1));
        ylabel(handleZ(i),label(2));
        
      end
      figure(10*k+2); for i=1:numvar, subplot(ny,nx,plot_i(i)); hold on; colorbar; set(gca,'Position',get(handle(i),'Position')); end 
      figure(10*k+3); for i=1:numvar, subplot(ny,nx,plot_i(i)); hold on; colorbar; set(gca,'Position',get(handle(i),'Position')); end 
      figure(10*k+4); for i=1:numvar, subplot(ny,nx,plot_i(i)); hold on; colorbar; set(gca,'Position',get(handle(i),'Position')); end 
      saveas(h ,['figure' num2str(h ) '.pdf']);
      saveas(hX,['figure' num2str(hX) '.pdf']);
      saveas(hY,['figure' num2str(hY) '.pdf']);
      saveas(hZ,['figure' num2str(hZ) '.pdf']);
    end
  end



  % plot statistics for each bucket
  if (false)
    disp('... plot statistics(t) for each bucket (x,y) to see if statistics have converged in figure 2')
    h = figure(10*k+2);
    for i=1:numvar
      subplot(ny,nx,i);
      plot(t,stat.avg(:,i:numvar:end));
      title(titles{i});
      xlabel('t');
    end
    saveas(gcf,'figure2.eps','psc2')
  end

  % plot statistics over each bucket after relaxation
  if (false)
    disp('... plot statistics(n*x+y) for last timestep in figure 3')
    h = figure(10*k+3);
    for i=1:numvar
      subplot(ny,nx,i);
      k = bar(stat.avg(end,i:numvar:end),1);
      xlabel('n*x+y');
      title(titles{i});
      axis tight;
    end
  end

  %display data of last timestep
  if (false)
    disp('data of last timestep');
    M = zeros(numvar,ntotal);
    M(:) = stat.avg(end,:);
    disp(cell2mat(titles'));
    disp(num2str(M','%11.5g '));
  end

  disp(['time: ' num2str(t(end))]);
  set(figure(10*k+1),'PaperType','A4');

end