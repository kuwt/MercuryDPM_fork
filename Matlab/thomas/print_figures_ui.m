function print_figures_ui(figlist)
%print all figures by default
if ~exist('figlist','var'), figlist=get(0,'Children'); end
figure(9999); clf
order_figures(3)
figure(9999); 
figlist=sort(figlist)
y=length(figlist)*20;
for i=figlist'
  uicontrol('Style', 'pushbutton', 'String', ['Print ' num2str(i) ' (' get(i,'FileName') ')'],'Position', [20 y 150 20],'Callback',['print_figures(' num2str(i) ')']);
  uicontrol('Style', 'pushbutton', 'String', ['Print&scp' num2str(i)],'Position', [160 y 150 20],'Callback',['print_figures(' num2str(i) ',struct(''scp'',true))']);
  y=y-20;
end

return

