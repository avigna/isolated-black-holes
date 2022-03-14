function analyse_COMPAS_population(metallicity, dir, debugFlag, plotFlag, saveFlag, densityInverval)
tic;
 
if exist(dir)~=1
    dir = '../data/COMPAS/';
end

if debugFlag~=1
    debugFlag = 0;
end

if debugFlag
    dir
    debugFlag
end

display('Chosen metallicity:')
num2str(metallicity)

if metallicity==0.03
    metallicityLabel = 'Z_0.03';
elseif metallicity==0.02
    metallicityLabel = 'Z_0.02';
elseif metallicity==0.0142
    metallicityLabel = 'Z_0.0142';
elseif metallicity==0.01
    metallicityLabel = 'Z_0.01';
elseif metallicity==0.0047
    metallicityLabel = 'Z_0.0047';
elseif metallicity==0.0021
    metallicityLabel = 'Z_0.0021';
elseif metallicity==0.001
    metallicityLabel = 'Z_0.001';
elseif metallicity==0.0002
    metallicityLabel = 'Z_0.0002';
elseif metallicity==0.0001
    metallicityLabel = 'Z_0.0001';
else
    warning('There is no available data with that metallicity.')
end

filename = strcat(dir,metallicityLabel,'/COMPAS_Output.h5');

% Choose analysis
calculate_and_plot_contours(filename,metallicity,debugFlag,plotFlag, saveFlag, densityInverval)
calculate_fractions(filename,num2str(metallicity),debugFlag, saveFlag)


toc;
end