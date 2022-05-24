function analyse_COMPAS_population_full(dir, debugFlag, plotFlag, savePlotFlag, savaDataFlag, densityInverval, lowerLimit)
tic;

if exist(dir)~=1
    % You can download it from: https://zenodo.org/record/6346444
    dir = '/set_the_path_to_your_COMPAS_data';

end

if debugFlag~=1
    debugFlag = 0;
end

if debugFlag
    dir
    debugFlag
end

call_functions_for_specific_metallicity('Z_0.03',0.03,dir, debugFlag, plotFlag, savePlotFlag, savaDataFlag, densityInverval, lowerLimit)
call_functions_for_specific_metallicity('Z_0.02',0.02,dir, debugFlag, plotFlag, savePlotFlag, savaDataFlag, densityInverval, lowerLimit)
tcall_functions_for_specific_metallicity('Z_0.0142',0.0142,dir, debugFlag, plotFlag, savePlotFlag, savaDataFlag, densityInverval, lowerLimit)
call_functions_for_specific_metallicity('Z_0.01',0.01,dir, debugFlag, plotFlag, savePlotFlag, savaDataFlag, densityInverval, lowerLimit)
call_functions_for_specific_metallicity('Z_0.0047',0.0047,dir, debugFlag, plotFlag, savePlotFlag, savaDataFlag, densityInverval, lowerLimit)
call_functions_for_specific_metallicity('Z_0.0021',0.0021,dir, debugFlag, plotFlag, savePlotFlag, savaDataFlag, densityInverval, lowerLimit)
call_functions_for_specific_metallicity('Z_0.001',0.001,dir, debugFlag, plotFlag, savePlotFlag, savaDataFlag, densityInverval, lowerLimit)
call_functions_for_specific_metallicity('Z_0.0002',0.0002,dir, debugFlag, plotFlag, savePlotFlag, savaDataFlag, densityInverval, lowerLimit)
call_functions_for_specific_metallicity('Z_0.0001',0.0001,dir, debugFlag, plotFlag, savePlotFlag, savaDataFlag, densityInverval, lowerLimit)

toc;
end

function call_functions_for_specific_metallicity(metallicityLabel,metallicity, dir, debugFlag, plotFlag, savePlotFlag, savaDataFlag, densityInverval, lowerLimit)
filename = strcat(dir,metallicityLabel,'/COMPAS_Output.h5');
calculate_and_plot_contours(filename, metallicity, debugFlag, plotFlag, savePlotFlag, savaDataFlag, densityInverval, lowerLimit)
calculate_fractions(filename,num2str(metallicity),debugFlag, savaDataFlag, lowerLimit)
calculate_mass_distribution(filename, metallicity, debugFlag, plotFlag, savePlotFlag, lowerLimit)
end