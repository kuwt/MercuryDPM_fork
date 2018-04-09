function enlargesubplot(ynum,xnum,subnum)
figure(1);
k=subplot(ynum,xnum,subnum);
h=get(k,'Children');

figure    % Create a new figure
axes    % Create an axes object in the figure
new_handle = copyobj(h,gca);
return