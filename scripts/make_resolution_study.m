function make_resolution_study(debugFlag, savePlotFlag)
tic;

% Load Data
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

% Some numerical values
metallicity = [0.02 0.0142 0.01 0.0047 0.0021 0.001 0.0002 0.0001];
number_of_binaries = 5*10^6;
number_of_metallicities = length(metallicity);

% 1. First born BH (when disrupted)
num = 3;
yield_number_one                = [yield_Z_0_0200(num) yield_Z_0_0142(num) yield_Z_0_0100(num) yield_Z_0_0047(num) yield_Z_0_0021(num) yield_Z_0_0010(num) yield_Z_0_0002(num) yield_Z_0_0001(num)];
N_number_one                    = [Z_0_0200.M(num) Z_0_0142.M(num) Z_0_0100.M(num) Z_0_0047.M(num) Z_0_0021.M(num) Z_0_0010.M(num) Z_0_0002.M(num) Z_0_0001.M(num)];

matrix_first_BH = false(number_of_binaries,number_of_metallicities);
matrix_first_BH(1:Z_0_0200.M(num),1)=true;
matrix_first_BH(1:Z_0_0142.M(num),2)=true;
matrix_first_BH(1:Z_0_0100.M(num),3)=true;
matrix_first_BH(1:Z_0_0047.M(num),4)=true;
matrix_first_BH(1:Z_0_0021.M(num),5)=true;
matrix_first_BH(1:Z_0_0010.M(num),6)=true;
matrix_first_BH(1:Z_0_0002.M(num),7)=true;
matrix_first_BH(1:Z_0_0001.M(num),8)=true;

% 12. NSNS merger product
num = 9;
yield_number_twelve             = [yield_Z_0_0200(num) yield_Z_0_0142(num) yield_Z_0_0100(num) yield_Z_0_0047(num) yield_Z_0_0021(num) yield_Z_0_0010(num) yield_Z_0_0002(num) yield_Z_0_0001(num)];
N_number_twelve                 = [Z_0_0200.M(num) Z_0_0142.M(num) Z_0_0100.M(num) Z_0_0047.M(num) Z_0_0021.M(num) Z_0_0010.M(num) Z_0_0002.M(num) Z_0_0001.M(num)];

matrix_NSNS = false(number_of_binaries,number_of_metallicities);
matrix_NSNS(1:Z_0_0200.M(num),1)=true;
matrix_NSNS(1:Z_0_0142.M(num),2)=true;
matrix_NSNS(1:Z_0_0100.M(num),3)=true;
matrix_NSNS(1:Z_0_0047.M(num),4)=true;
matrix_NSNS(1:Z_0_0021.M(num),5)=true;
matrix_NSNS(1:Z_0_0010.M(num),6)=true;
matrix_NSNS(1:Z_0_0002.M(num),7)=true;
matrix_NSNS(1:Z_0_0001.M(num),8)=true;

if debugFlag
    N_number_one
    N_number_twelve
end


% Plot
clf
fs=18;
lw=2;
solar=char(9737);
string1='Z';  
string2=['Number of Systems'];
Xlim = [0.0001 0.02];

gray = 0.5.*[1 1 1];

color1 = [0    0.4470    0.7410];
color7 = [    0.6350    0.0780    0.1840];

hold on
xlim(Xlim)
set(gca,'YScale','log')
xticks([0.0001 0.005 0.01 0.015 0.02])
xlabel(string1,'FontSize',fs)
ylabel(string2,'FontSize',fs)
ax=gca;
ax.FontSize=fs;
ax.FontName = 'Times New Roman';
box on

for i=1:100
r_index = randi([1 number_of_binaries],number_of_binaries,1);
bootstrap_number_one    = [sum(matrix_first_BH(r_index,1)) sum(matrix_first_BH(r_index,2)) sum(matrix_first_BH(r_index,3)) sum(matrix_first_BH(r_index,4)) sum(matrix_first_BH(r_index,5)) sum(matrix_first_BH(r_index,6)) sum(matrix_first_BH(r_index,7)) sum(matrix_first_BH(r_index,8))];
bootstrap_number_twelve = [sum(matrix_NSNS(r_index,1)) sum(matrix_NSNS(r_index,2)) sum(matrix_NSNS(r_index,3)) sum(matrix_NSNS(r_index,4)) sum(matrix_NSNS(r_index,5)) sum(matrix_NSNS(r_index,6)) sum(matrix_NSNS(r_index,7)) sum(matrix_NSNS(r_index,8))];

plot(metallicity,bootstrap_number_one,'Color',gray,'HandleVisibility','off')
plot(metallicity,bootstrap_number_twelve,'Color',gray,'HandleVisibility','off')

clear r_index bootstrap_number_one bootstrap_number_twelve
end

plot(metallicity,N_number_one,'Color',color1,'LineWidth',lw)
plot(metallicity,N_number_twelve,':','Color',color7,'LineWidth',lw)

legend( 'BH_1',...
        'NS-NS',...
        'NumColumns',1,...
        'box','off',...
        'location','eastoutside')

if savePlotFlag
    print(gcf,strcat('../plots/png/resolution_study.png'),'-dpng','-r300');
    saveas(gcf,strcat('../plots/pdf/resolution_study.pdf'))
end

toc;
end