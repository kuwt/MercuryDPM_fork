function plotover(ynum,xnum,num1,num2)
figure(1);
h1Axis=subplot(ynum,xnum,num1);
h1Children=get(h1Axis,'Children');
h1Title=get(h1Axis,'Title');

h2Axis=subplot(ynum,xnum,num2);
h2Children=get(h2Axis,'Children');
h2Title=get(h2Axis,'Title');

figure(2);    % Create a new figure
clf
h3 = axes;    % Create an axes object in the figure
hold on
for i = 1:length(h1Children)
  h = plot(get(h1Children(i),'YData'),get(h2Children(i),'YData'),'-o');
  set(h,'Color',get(h1Children(i),'Color'));
  set(h,'LineStyle', get(h1Children(i),'LineStyle'));
end
xlabel(get(h1Title,'String'))
ylabel(get(h2Title,'String'))
set(h3,'XLim', get(h1Axis,'YLim'));
set(h3,'YLim', get(h2Axis,'YLim'));
set(h3,'XScale','log');
set(h3,'YScale','log');
saveas(2 ,['figure2.pdf']);
return