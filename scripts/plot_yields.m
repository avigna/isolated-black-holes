function plot_yields(saveFlag)
tic;

Z_0_0200 = load('../data/Z_0.02.mat');
Z_0_0142 = load('../data/Z_0.0142.mat');
Z_0_0047 = load('../data/Z_0.0047.mat');
Z_0_0021 = load('../data/Z_0.0021.mat');
Z_0_0100 = load('../data/Z_0.01.mat');
Z_0_0010 = load('../data/Z_0.001.mat');
Z_0_0002 = load('../data/Z_0.0002.mat');
Z_0_0001 = load('../data/Z_0.0001.mat');

yield_Z_0_0200 = Z_0_0200.M./Z_0_0200.M(1);
yield_Z_0_0142 = Z_0_0142.M./Z_0_0142.M(1);
yield_Z_0_0100 = Z_0_0100.M./Z_0_0100.M(1);
yield_Z_0_0047 = Z_0_0047.M./Z_0_0047.M(1);
yield_Z_0_0021 = Z_0_0021.M./Z_0_0021.M(1);
yield_Z_0_0010 = Z_0_0010.M./Z_0_0010.M(1);
yield_Z_0_0002 = Z_0_0002.M./Z_0_0002.M(1);
yield_Z_0_0001 = Z_0_0001.M./Z_0_0001.M(1);

metallicity = [0.02 0.0142 0.01 0.0047 0.0021 0.001 0.0002 0.0001];

% M = [   realMassZAMS, number_zero, number_one, number_three, number_five, number_six, number_ten, number_eleven, number_twelve]

% 0. Isolated BH
num = 2;
yield_number_zero               = [yield_Z_0_0200(num) yield_Z_0_0142(num) yield_Z_0_0100(num) yield_Z_0_0047(num) yield_Z_0_0021(num) yield_Z_0_0010(num) yield_Z_0_0002(num) yield_Z_0_0001(num)];

% 1. First born BH (when disrupted)
num = 3;
yield_number_one                = [yield_Z_0_0200(num) yield_Z_0_0142(num) yield_Z_0_0100(num) yield_Z_0_0047(num) yield_Z_0_0021(num) yield_Z_0_0010(num) yield_Z_0_0002(num) yield_Z_0_0001(num)];

% 3. second born BH (when disrupted)
num = 4;
yield_number_three              = [yield_Z_0_0200(num) yield_Z_0_0142(num) yield_Z_0_0100(num) yield_Z_0_0047(num) yield_Z_0_0021(num) yield_Z_0_0010(num) yield_Z_0_0002(num) yield_Z_0_0001(num)];

% 5. first born BH, when not disrupted (X-ray binarieS at birth)
num = 5;
yield_number_five               = [ yield_Z_0_0200(num) yield_Z_0_0142(num) yield_Z_0_0100(num) yield_Z_0_0047(num) yield_Z_0_0021(num) yield_Z_0_0010(num) yield_Z_0_0002(num) yield_Z_0_0001(num)];

% 6. first born NS, when not disrupted (X-ray binarieS at birth)
num = 6;
yield_number_six                = [yield_Z_0_0200(num) yield_Z_0_0142(num) yield_Z_0_0100(num) yield_Z_0_0047(num) yield_Z_0_0021(num) yield_Z_0_0010(num) yield_Z_0_0002(num) yield_Z_0_0001(num)];

% 10. NS-BH, BH-NS merger product
num = 7;
yield_number_ten                = [yield_Z_0_0200(num) yield_Z_0_0142(num) yield_Z_0_0100(num) yield_Z_0_0047(num) yield_Z_0_0021(num) yield_Z_0_0010(num) yield_Z_0_0002(num) yield_Z_0_0001(num)];

% 11. BBH merger product
num = 8;
yield_number_eleven             = [yield_Z_0_0200(num) yield_Z_0_0142(num) yield_Z_0_0100(num) yield_Z_0_0047(num) yield_Z_0_0021(num) yield_Z_0_0010(num) yield_Z_0_0002(num) yield_Z_0_0001(num)];

% 12. NSNS merger product
num = 9;
yield_number_twelve             = [yield_Z_0_0200(num) yield_Z_0_0142(num) yield_Z_0_0100(num) yield_Z_0_0047(num) yield_Z_0_0021(num) yield_Z_0_0010(num) yield_Z_0_0002(num) yield_Z_0_0001(num)];
  
% Plot
clf
fs=18;
lw=2;
solar=char(9737);
string1='Z';  
string2=['dN_{form}/dM_{SFR} [ M_',solar,'^{-1} ]'];
Ylim = [5*10^-7 10^-3];
Xlim = [0.0001 0.02];

color0 = 105.*[1 1 1]./255;

% color1 = [138 186 237]./255;
% color2 = [154 204 108]./255;
% color3 = [161 71 114]./255;    
% color4 = [235 201 164]./255;
% color5 = [213 94 0]./255;
% color6 = [240 228 66]./255;
% color7 = [0 114 178]./255;

color1 = [0    0.4470    0.7410];
color2 = [    0.8500    0.3250    0.0980];
color3 = [    0.9290    0.6940    0.1250];
color4 = [    0.4940    0.1840    0.5560];
color5 = [    0.4660    0.6740    0.1880];
color6 = [    0.3010    0.7450    0.9330];
color7 = [    0.6350    0.0780    0.1840];

hold on
ylim(Ylim)
xlim(Xlim)
set(gca,'YScale','log')
xticks([0.0001 0.005 0.01 0.015 0.02])
xlabel(string1,'FontSize',fs)
ylabel(string2,'FontSize',fs)
ax=gca;
ax.FontSize=fs;
ax.FontName = 'Times New Roman';
box on

plot(metallicity,yield_number_zero,'Color',color0,'LineWidth',lw)
plot(metallicity,yield_number_one,'Color',color1,'LineWidth',lw)
plot(metallicity,yield_number_three,'Color',color2,'LineWidth',lw)
plot(metallicity,yield_number_five,'-','Color',color3,'LineWidth',lw)
plot(metallicity,yield_number_six,'-','Color',color4,'LineWidth',lw)
plot(metallicity,yield_number_ten,':','Color',color5,'LineWidth',lw)
plot(metallicity,yield_number_eleven,':','Color',color6,'LineWidth',lw)
plot(metallicity,yield_number_twelve,':','Color',color7,'LineWidth',lw)

legend( 'BH_{sin}',...
        'BH_1',...
        'BH_2',...
        'BH-MS',...
        'NS-MS',...
        'BH-NS',...
        'BH-BH',...
        'NS-NS',...
        'NumColumns',1,...
        'box','off',...
        'location','eastoutside')

if saveFlag
    print(gcf,strcat('../plots/png/figure1.png'),'-dpng','-r300');
    saveas(gcf,strcat('../plots/pdf/figure1.pdf'))
end

toc;
end