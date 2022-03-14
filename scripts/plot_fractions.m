function plot_fractions(saveFlag)
tic;

Z_0_0200 = load('../data/fraction_Z_0.02.mat')
Z_0_0142 = load('../data/fraction_Z_0.0142.mat');
Z_0_0100 = load('../data/fraction_Z_0.01.mat');
Z_0_0047 = load('../data/fraction_Z_0.0047.mat');
Z_0_0021 = load('../data/fraction_Z_0.0021.mat');
Z_0_0010 = load('../data/fraction_Z_0.001.mat');
Z_0_0002 = load('../data/fraction_Z_0.0002.mat');
Z_0_0001 = load('../data/fraction_Z_0.0001.mat');

metallicity             = [ 0.0200 ...
                            0.0142 ...
                            0.0100 ...
                            0.0047 ...
                            0.0021 ...
                            0.0010 ...
                            0.0002 ...
                            0.0001];

% Speed limit 50 km/s
fraction_mass_5_speed_50   = [  Z_0_0200.fraction_mass_5_speed_50 ...
                                Z_0_0142.fraction_mass_5_speed_50 ...
                                Z_0_0100.fraction_mass_5_speed_50 ...
                                Z_0_0047.fraction_mass_5_speed_50 ...
                                Z_0_0021.fraction_mass_5_speed_50 ...
                                Z_0_0010.fraction_mass_5_speed_50 ...
                                Z_0_0002.fraction_mass_5_speed_50 ...
                                Z_0_0001.fraction_mass_5_speed_50];

fraction_mass_10_speed_50   = [ Z_0_0200.fraction_mass_10_speed_50 ...
                                Z_0_0142.fraction_mass_10_speed_50 ...
                                Z_0_0100.fraction_mass_10_speed_50 ...
                                Z_0_0047.fraction_mass_10_speed_50 ...
                                Z_0_0021.fraction_mass_10_speed_50 ...
                                Z_0_0010.fraction_mass_10_speed_50 ...
                                Z_0_0002.fraction_mass_10_speed_50 ...
                                Z_0_0001.fraction_mass_10_speed_50];

fraction_mass_20_speed_50   = [ Z_0_0200.fraction_mass_20_speed_50 ...
                                Z_0_0142.fraction_mass_20_speed_50 ...
                                Z_0_0100.fraction_mass_20_speed_50 ...
                                Z_0_0047.fraction_mass_20_speed_50 ...
                                Z_0_0021.fraction_mass_20_speed_50 ...
                                Z_0_0010.fraction_mass_20_speed_50 ...
                                Z_0_0002.fraction_mass_20_speed_50 ...
                                Z_0_0001.fraction_mass_20_speed_50];

fraction_mass_100_speed_50   = [ Z_0_0200.fraction_mass_100_speed_50 ...
                                Z_0_0142.fraction_mass_100_speed_50 ...
                                Z_0_0100.fraction_mass_100_speed_50 ...
                                Z_0_0047.fraction_mass_100_speed_50 ...
                                Z_0_0021.fraction_mass_100_speed_50 ...
                                Z_0_0010.fraction_mass_100_speed_50 ...
                                Z_0_0002.fraction_mass_100_speed_50 ...
                                Z_0_0001.fraction_mass_100_speed_50];




clf
fs=18;
lw=2.0;
string1='Z';  

hold on
xticks([0.0001 0.005 0.01 0.015 0.02])
ylim([0.1 1])
xlabel(string1,'FontSize',fs,'FontName','Times New Roman')
ylabel('N_{bin} / N_{tot}','FontSize',fs,'FontName','Times New Roman')
ax=gca;
ax.FontSize=fs;
ax.FontName = 'Times New Roman';
box on

plot(metallicity,fraction_mass_5_speed_50,'k','LineWidth',lw)
plot(metallicity,fraction_mass_10_speed_50,'--k','LineWidth',lw)
plot(metallicity,fraction_mass_20_speed_50,':k','LineWidth',lw)
plot(metallicity,fraction_mass_100_speed_50,'-.k','LineWidth',lw)

legend( '$M \leq 5\ \rm{M_{\odot}}$',...
        '$M \leq 10\ \rm{M_{\odot}}$',...
        '$M \leq 20\ \rm{M_{\odot}}$',...
        '$M \leq 100\ \rm{M_{\odot}}$',...
        'Location','SouthEast',...
        'Interpreter','latex',...
        'box','off')

if saveFlag
    print(gcf,strcat('../plots/png/figure3.png'),'-dpng','-r300');
    saveas(gcf,strcat('../plots/pdf/figure3.pdf'))
end


toc;
end