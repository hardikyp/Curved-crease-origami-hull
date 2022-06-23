%% OPENS Geom_Prop.mpr, READS GEOMETRIC PROPERTIES AND SAVES IT TO VARIABLES
clc;clear;close all;
fid = fopen('CAD_Files/Geom_Prop.mpr','r');
numLines = lineCount(fid);
fclose(fid);

fid = fopen('CAD_Files/Geom_Prop.mpr','r');
i = 1;

while i <= numLines
    line = fgetl(fid);
    if string(regexp(line,'Volume','match')) == "Volume"
        volume = str2double(cell2mat(regexp(line,'\d*[.]\d*','match')));
    elseif string(regexp(line,'Centroid','match')) == "Centroid"
        centroidX = str2double(cell2mat(regexp(line,'\d*[.]\d*','match')));
        line = fgetl(fid);
        i = i + 1;
        centroidY = str2double(cell2mat(regexp(line,'\d*[.]\d*','match')));
        line = fgetl(fid);
        i = i + 1;
        centroidZ = str2double(cell2mat(regexp(line,'\d*[.]\d*','match')));
    elseif string(regexp(line,'Radii of gyration','match')) == "Radii of gyration"
        radiusGyrationX = str2double(cell2mat(regexp(line,'\d*[.]\d*','match')));
        line = fgetl(fid);
        i = i + 1;
        radiusGyrationY = str2double(cell2mat(regexp(line,'\d*[.]\d*','match')));
        line = fgetl(fid);
        i = i + 1;
        radiusGyrationZ = str2double(cell2mat(regexp(line,'\d*[.]\d*','match')));
    end
    i = i + 1;
end

density = 0.00000271;               % density of Aluminium (kg/mm^3)
mass = density * volume;

clearvars ans fid i line numLines;
%% Function lineCount(fildID)
% Pre-determines the number of lines in a file
function n = lineCount(fileID)
n = 0;
tline = fgetl(fileID);
    while ischar(tline)
    tline = fgetl(fileID);
    n = n+1;
    end
end
