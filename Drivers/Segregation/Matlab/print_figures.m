function print_figures(figlist)
%mkdir Figures
%print all figures by default
if ~exist('figlist','var'), figlist=get(0,'Children'); end


for i=figlist'
  %set default file name
  if isempty(get(i,'FileName')), set(i,'FileName',['figure' num2str(i)]); end
  P=get(i,'Position');
  PrintLaTeX(get(i,'Number'),['Figures/' get(i,'FileName')],P(3)*8/560);
  Axes=get(i,'CurrentAxes');
  title(Axes(1),['saved as Figures/' get(i,'FileName')],'Interpreter','none')
end
return