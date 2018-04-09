function print_figures(figlist,opt)
%print all figures by default
if ~exist('figlist','var')||isempty(figlist), figlist=get(0,'Children'); end
if ~exist('opt'); opt=struct([]); end
if ~exist('opt'); opt=struct([]); end

for i=figlist'
  if ~strcmp(get(i,'Name'),'Printed by LaPrint')
    set(get(get(i,'CurrentAxes'),'Title'),'String','')
    %set default file name
    if isempty(get(i,'FileName')), set(i,'FileName',['figure' num2str(i)]); end
    P=get(i,'Position');
    [~,~]=system('mkdir -v Figures');
    saveas(i,['Figures/' get(i,'FileName') '.fig'])
    PrintLaTeX(i,['Figures/' get(i,'FileName')],P(3)*8/560);
    Axes=get(i,'CurrentAxes');
    title(Axes(1),['saved as Figures/' get(i,'FileName')],'Interpreter','none')
  else
    disp('already printed');
  end
  system(['~/Documents/Work/bin/myMatlabepstopdf ' pwd '/Figures/' get(i,'FileName') '.matlab.tex'])
  if isfield(opt,'scp')
    %system(['scp Figures/' get(i,'FileName') '.matlab.*  Figures/' get(i,'FileName') '.fig $CTW:~/Dropbox/ChuteRheology/Figures'])
    system(['scp Figures/' get(i,'FileName') '.matlab.*  Figures/' get(i,'FileName') '.fig /Users/weinhartt/Dropbox/ChuteRheology/Figures'])
  end
end

return
