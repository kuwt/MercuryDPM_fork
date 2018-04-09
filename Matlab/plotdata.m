function plotdata(data,info)
% plots data of one timestep

if (isfield(data,'variance')), variance = true; else, variance = false; end

%set h
if (isfield(info,'h')), h = info.h; else, h=1; end
disp(['plotting into figure' num2str(h) ' ...']);

if (~isfield(info,'plottype')), info.plottype = '-'; end

% checks which dimensions are really used and stores their coordinates
dim = ~[ ...
  all(data.coordinates(1,1)==data.coordinates(:,1)) ...
  all(data.coordinates(1,2)==data.coordinates(:,2)) ...
  all(data.coordinates(1,3)==data.coordinates(:,3)) ];
dim_label_full = 'xyz';
dim_label = dim_label_full(dim);
coordinates = data.coordinates(:,dim);
disp(['dim=' num2str(dim) ' ...']);

% stores the grid dimensions
for i=1:sum(dim)
  n(i) = size(coordinates,1) / sum(coordinates(1,i)==coordinates(:,i));
end

% stores coordinates in grid
for i=1:sum(dim)
  if (length(n)==1)
    X{i} = coordinates(:,i);
  else
    X{i} = zeros(n(end:-1:1));
    X{i}(:) = coordinates(:,i);
  end
end
Z = zeros(size(X{1}));
if (variance), Zvar = zeros(size(X{1})); end

%create figure
figure(h);
%clear figure if info.hold does not exist
if (~isfield(info,'hold')), clf; end
%resize figure     
set(gcf,'Position',[300 1 1379 939],'PaperType','A4','PaperOrientation','landscape','PaperPosition',[0 0 11.69 8.26]);
%plot info text
%info_text = [info.text ' t=' mat2str(data.time)];
%global htextbox
%if ~isempty('htextbox'), htextbox = annotation('textbox',[0 0 1 1]); end
%set(htextbox,'String',info_text ,'Interpreter','none')

% defines number of variables and plot order
numvar = size(data.variables,2);
plot_i = 1:numvar;
if (numvar==1), 
  ynum = 1;
  xnum = 1;
elseif (numvar==2), 
  ynum = 1;
  xnum = 2;
elseif (numvar<=4), 
  ynum = 2;
  xnum = 2;
% elseif (numvar==6), %plotting a symmetric matrix
%   ynum = 3;
%   xnum = 3;
%   plot_i = [1:3 5 6 9];
elseif (numvar<=9), 
  ynum = 3;
  xnum = 3;
elseif (numvar<=16), 
  ynum = 4;
  xnum = 4;
  plot_i = [1 5 9 13 6 7 8 11 12 16 2 3 4 10 14 15];
  %ynum = 4;
  %xnum = 3;
elseif (numvar<=24)
  ynum = 4;
  xnum = 6;
  plot_i = [1:2 4:numvar+1];
else (numvar<=36)
  ynum = 6;
  xnum = 6;
%  plot_i = [1:2 4:15 19:36 16:18 3];
end

if (sum(dim)==1)
  % for 1D output
  %set(htextbox,'String',['plotting 1-dimensional data' info_text],'Interpreter','none')
  
	for i=1:numvar
    subplot(ynum,xnum,plot_i(i)); 
		hold on
    Z(:) = data.variables(:,i);
    plot(X{1},Z,info.plottype);
    plot(sum(X{1})/length(X{1}),sum(Z)/length(Z),'rx');
    if (variance), 
      Zvar(:) = sqrt(data.variance(:,i)); 
      plot(X{1},Z-Zvar,':');
      plot(X{1},Z+Zvar,':');
    end
    title(info.variable_names{i},'Interpreter','none');
    xlabel(dim_label);
		hold off
	end

elseif (sum(dim)==2)
  % for 2D output
  
	if (~(isfield(info,'dim_avg') && dim(info.dim_avg)) )
		%plotting without averaging
		%set(htextbox,'String',['plotting 2-dimensional data' info_text],'Interpreter','none')
		for i=1:numvar
			subplot(ynum,xnum,plot_i(i)); hold on
			Z(:) = data.variables(:,i);
			binplot(X{1},X{2},Z);
			colorbar;
			shading flat;
			axis equal;
			title(info.variable_names{i},'Interpreter','none');
			xlabel(dim_label(1));
			ylabel(dim_label(2));
		end
	else
		%plotting averaged data
 
		dim_avg = info.dim_avg(1);
    dim_avg_bool = zeros(1,3);
    dim_avg_bool(dim_avg) = true;
		dim_avg_local = find(find(dim)==dim_avg);
 		%set(htextbox,'String',['plotting ' dim_label_full(dim_avg) '-averaged 2-dimensional data' info_text],'Interpreter','none')
		dim_other = find(dim-dim_avg_bool);
		dim_other_local = find(find(dim)==dim_other);
 		Xcut = squeeze(sum(X{dim_other_local},3-dim_avg_local)/n(dim_avg_local));
		for i=1:numvar
			subplot(ynum,xnum,plot_i(i));
			hold on 
			Z(:) = data.variables(:,i);
			Zcut = squeeze(sum(Z,3-dim_avg_local)/n(dim_avg_local));
			plot(Xcut,Zcut,info.plottype);
			plot(sum(Xcut)/length(Xcut),sum(Zcut)/length(Zcut),'rx');
			xlabel(dim_label(dim_other_local));
			title(info.variable_names{i},'Interpreter','none');
			hold off
		end
	end

elseif (sum(dim)==3)
	%requires at least one averaging dimension
	if (~isfield(info,'dim_avg')), info.dim_avg = 2; end
  if (~isfield(info,'avg_not_cut')), info.avg_not_cut = true; end

	if (length(info.dim_avg)==1)
		%plotting 2D averaged data
		dim_avg = info.dim_avg(1);

% 		if (info.avg_not_cut)
% 			set(htextbox,'String',[...
% 				'plotting ' dim_label(dim_avg) '-averaged 3-dimensional data' ...
% 				info_text],'Interpreter','none')
% 		else
% 			set(htextbox,'String',[...
% 				'plotting center ' dim_label(dim_avg) '-cut of 3-dimensional data' ...
% 				info_text],'Interpreter','none')
% 		end

		dim_other = [1:1:dim_avg-1 dim_avg+1:1:3];
		Xcut{1} = squeeze(sum(X{dim_other(1)},4-dim_avg)/n(dim_avg));
		Xcut{2} = squeeze(sum(X{dim_other(2)},4-dim_avg)/n(dim_avg));

		for i=1:numvar
			subplot(ynum,xnum,plot_i(i)); 
			hold on
			Z(:) = data.variables(:,i);
			if (info.avg_not_cut)
				Zcut = squeeze(sum(Z,4-dim_avg)/n(dim_avg));
			else
				Zperm = permute(Z,4-[dim_other(end:-1:1) dim_avg]);
				Zcut = Zperm(:,:,round(end/2));
			end
			binplot(Xcut{1},Xcut{2},Zcut);
			colorbar;
			shading flat;
			axis equal;
			title(info.variable_names{i},'Interpreter','none');
			xlabel(dim_label(dim_other(1)));
			ylabel(dim_label(dim_other(2)));
			hold off
		end
	else
		%plotting 1D averaged data

% 		set(htextbox,'String',[...
% 		'plotting ' dim_label(info.dim_avg(1)) '-' dim_label(info.dim_avg(2)) '-averaged 3-dimensional data' ...
% 		info_text],'Interpreter','none');
% 		
		dim_other =  [1 2 3]~=info.dim_avg(1) & [1 2 3]~=info.dim_avg(2);
		Xcut =  squeeze(sum(sum(X{dim_other},4-info.dim_avg(2)),4-info.dim_avg(1)) ...
			/prod(n(info.dim_avg)));

		for i=1:numvar
			subplot(ynum,xnum,plot_i(i)); 
			hold on
			Z(:) = data.variables(:,i);
			Zcut = squeeze(sum(sum(Z,4-info.dim_avg(2)),4-info.dim_avg(1)) ...
				/prod(n(info.dim_avg)));
			plot(Xcut,Zcut,info.plottype);
			plot(sum(Xcut)/length(Xcut),sum(Zcut)/length(Zcut),'rx');
			%axis tight
			xlabel(dim_label(dim_other));
			title(info.variable_names{i},'Interpreter','none');
			hold off
		end
		
	end
	
end

if (numvar==16)
  %scale
  for i=[1:16]
    subplot(4,4,i);
    axis tight
  end
  % PVF vertical range: [0:0.7]
  subplot(4,4,1);
  v = axis;
  axis([v(1:2) 0 .7]);
  % PVF vertical range: [0:0.7]
  subplot(4,4,4);
  v = axis;
  axis([v(1:3) .66]);
  % vy and vz with the same vertical axis range, e.g. [-0.02:0.02]
  same_scale([9 13],xnum,ynum);
  % fix vertical axis range for all diagonal stress: e.g. [-20:200]
  same_scale([6 11 16],xnum,ynum);
  % use the same range for Sxy and Syz: e.g. [-10:10]
  same_scale([7 12],xnum,ynum);
  % pls. zoom in [40:50] or [35:55] and plot a line at 45 degrees
  subplot(4,4,10);
  v = axis;
  hold on 
  plot(v(1:2),[45 45])
  %axis([v(1:2) 40 50]);
  subplot(4,4,14);
  v = axis;
  hold on 
  plot(v(1:2),[90 90])
  %axis([v(1:2) 80 91]);
  subplot(4,4,15);
  v = axis;
  hold on 
  plot(v(1:2),[45 45])
  %axis([v(1:2) 40 50]);
elseif (numvar==23 || numvar==24)
  same_scale([4:6],xnum,ynum);
  same_scale([7:12],xnum,ynum);
  same_scale([13:18],xnum,ynum);
  same_scale([19:24],xnum,ynum);
end

return