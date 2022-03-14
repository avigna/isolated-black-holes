function test_plot_KDE(data, label, densityInverval)

% Values for creatings contours
% [bandwith, density, X, Y] = kde2d(data);
% Template
% data = [SNe_mass_SN(index_), SNe_component_speed_SN(index_)];
% [bandwith, density, X, Y] = kde2d(data,n,MIN_XY,MAX_XY);
% v =max(max(density)).*[0.01,1.0];

n=2^10;
MIN_XY = [1, 0];
MAX_XY = [100, 1000];

[bandwith, density, X, Y]   = kde2d(data,n,MIN_XY,MAX_XY);
v                           = max(max(density)).*[1.0-densityInverval,1.0];

% Plot
figure()
fs=18;
lw=1.5;
sz=1.0;
solar=char(9737);
string1=['Mass [ M_',solar,' ]'];  
string2='Speed [ km s^{-1} ]';
alphaNum = 0.2;
Ylim = [0 400];
Xlim = [1 100];

chosenColor = [138 186 237]./255;

hold on
text(1.15,375,strcat('Z=',label),'FontSize',fs,'FontName','Helvetica')
ylim(Ylim)
xlim(Xlim)
set(gca,'XScale','log')
xticks([1 10 100])
xticklabels({'1','10','100'})
xlabel(string1,'FontSize',fs,'FontName','Helvetica')
ylabel(string2,'FontSize',fs,'FontName','Helvetica')
ax=gca;
ax.FontSize=fs;
box on

contour(X,Y,density,v,'Color',chosenColor,'Linewidth',lw)
scatter(data(:,1), data(:,2), sz, chosenColor, 'Filled', 'MarkerFaceAlpha', alphaNum, 'MarkerEdgeAlpha', alphaNum,'HandleVisibility','off')

end