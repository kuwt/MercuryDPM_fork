figchildren = get(figure(1),'Children');
for h_axes = figchildren'
  if (get(get(h_axes,'Title'),'String')=='height')
    break;
  end
end
h_data = get(h_axes,'Children');
delete(h_data(1:end-1));
h_data = get(h_axes,'Children');
x = get(h_data,'XData');
y = get(h_data,'YData');
range = round([0.3 0.7]*length(x));
[m,b] = findlinearfit(x(range(1):range(2)),y(range(1):range(2)));
axes(h_axes);
hold(h_axes,'on');
plot(h_axes,x,m*x+b,'r','LineWidth',1);
text(x(range(1)),y(range(1)),['     \theta=' num2str(atan(m)*180/pi) '{}^\circ']);
saveas(gcf,'figure1.eps','psc2')
saveas(gcf,'figure1.pdf')
