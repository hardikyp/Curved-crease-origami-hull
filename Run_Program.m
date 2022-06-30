%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                Bar and Hinge for Curved-Crease Origami                %%
%         Modification by: Steven R. Woodruff (stevenrw@umich.edu)        %
%         Modification by: Hardik Patil (hardikyp@umich.edu)              %
%                                                                         %
%             Base code written by: Ke Liu (ke.liu@gatech.edu)            %
% Ref: S. R. Woodruff, E. T. Filipov (2020). 'A bar and hinge formulation %
%          for structural analysis of curved-crease origami.' Submitted.  %
%      K. Liu, G. H. Paulino (2017). 'Nonlinear mechanics of non-rigid    %
%          origami - An efficient computational approach.' PRSA.          %
%      K. Liu, G. H. Paulino (2016). 'MERLIN: A MATLAB implementation to  %
%          capture highly nonlinear behavior of non-rigid origami.'       %
%          Proceedings of IASS Annual Symposium 2016.                     %
%      E. T. Filipov, K. Liu, T. Tachi, M. Schenk, G. H. Paulino (2017).  %
%          'Bar and hinge models for scalable analysis of origami.'  IJSS %
%      K. Liu, G. H. Paulino (2018). 'Highly efficient nonlinear          %
%          structural analysis of origami assemblages using the MERLIN2   %
%          software.' Origami^7.                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reset workspace
clearvars; close all; clc;

% Start timing the script
tic 

% Add required paths to analysis functions and fold pattern functions
addpath('Functions');
addpath('Geometries');
addpath('CAD_Files');

% Select Geometry

% testType = 'Planing';
testType = 'PlaningVariable';
% testType = 'QuickBoat';
% testType = 'QuickBoatVariable';


InputData = struct('restAngle1', -90, ...            % Prescribed folding rest dihedral angle (\phi_R) [deg] (Chine fold angle; Planing = -90, Quickboat = 90)
    'restAngle2', 60, ...                            % Prescribed folding rest dihedral angle (\phi_R) [deg] (Keel / Deadrise angle. 180-prescribed = angle between adj faces of hull; Planing = 60, Quickboat = 60))
    'lengthStraight', 500, ...                       % Length of the straight part of the hull
    'lengthCurve', 454, ...                          % Length of the curved part of the hull
    'tipHeight', 71, ...                             % Height of the nose tip of crease pattern
    'panelWidth', 162, ...                           % Height of the 
    'elasticModulus', 4000, ...                      % Elastic modulus of the sheet (E)
    'thickness', 0.127, ...                          % Thickness of the sheet (t)
    'lengthScale1', 50,...                           % folding length scale factor (L*) [Chine, Planing = 50]
    'lengthScale2', 1,...                            % folding length scale factor (L*) [Keel, Planing = 1]
    'loadMagnitude', 1, ...                          % Load magnitude (displacement or force)
    'loadType', 'Displacement', ...                  % Either 'Force' or 'Displacement' controlled test
    'numberIncrements', 20, ...                      % Specifies the number of displacement or force steps (increments)
    'maxIterations', 50, ...                         % Specifies the maximum number of step iterations
    'testType', testType, ...                        % Specifies specific test (see example options above)
    'numberDivisions', 80, ...                       % This controls the discretization size
    'plotNodes', 'no', ...                           % Plots the flat nodes, bars, and boundary conditions
    'plotEnergy', 'no', ...                          % Plots the strain energy in the elements
    'plotDeformedShape', 'yes', ...                  % Plots the deformed shape after analysis
    'plotReferenceShape', 'yes',...                  % Plots the flat, referece configuration under deformed shape
    'plotIncrement','end');                          % Value determines which step is plotted (or end)
clear testType;

%% Test Inputs
[InputData, inputTest] = Test_Input(InputData);

if inputTest == 0
    fprintf('Input variable error(s). See above.\n\n')
    return
elseif inputTest == 2
    fprintf('Possible error. See above. Hit enter to proceed.\n\n')
    pause
end
clear inputTest;

%% Generate mesh
if strcmp(InputData.testType,'Planing')
    InputData = Geometry_Planing(InputData);
    InputData.foldStep = true;
elseif strcmp(InputData.testType,'PlaningVariable')
    InputData = Geometry_Planing_Variable(InputData);
    InputData.foldStep = true;
elseif strcmp(InputData.testType,'QuickBoat')
    InputData = Geometry_QuickBoat(InputData);
    InputData.foldStep = true;
elseif strcmp(InputData.testType,'QuickBoatVariable')
    InputData = Geometry_QuickBoat_Variable(InputData);
    InputData.foldStep = true;
end

%% Perform analysis
[PreprocessData] = Prepare_Data(InputData);

%% Run Analysis
if strcmp(InputData.testType,'Planing') || strcmp(InputData.testType,'PlaningVariable') 
    [PostprocessData] = Path_Analysis_A(InputData, PreprocessData);
elseif strcmp(InputData.testType,'QuickBoat') || strcmp(InputData.testType,'QuickBoatVariable') 
    [PostprocessData] = Path_Analysis_B(InputData, PreprocessData);
end

%% Prepare Data for Post Processing
PostprocessData = Post_Process(InputData, PreprocessData, PostprocessData);

%% Plot results
Plot_Results(InputData, PreprocessData, PostprocessData);

%% Matching CCO surface with CAD surface (works with planing hull geometry only)
% Enrich with points
% Mesh_0 = enrich_nodes(5, InputData, PostprocessData);
% Mesh_0([180 356 531 708 894 174 353 533 711],:) = [];
% x0 = Mesh_0(:,1);
% y0 = Mesh_0(:,2);
% z0 = Mesh_0(:,3);
% Mesh_1 = readmatrix('planing_surface.txt');
% Mesh_1 = unique([Mesh_1;Mesh_1(:,1) -Mesh_1(:,2) Mesh_1(:,3)],'sorted','rows');
% x1 = Mesh_1(:,1);
% y1 = Mesh_1(:,2);
% z1 = Mesh_1(:,3);

% run('Woodruff_OBJMatcher_RODROT.m')

%% Organise points according to the edges of the hull
% curves = cell(1,10);
% curves{3} = PostprocessData.deformedNodes{end}([InputData.numberDivisions+2:2*InputData.numberDivisions+1 3*InputData.numberDivisions+2],:);
% curves{2} = PostprocessData.deformedNodes{end}(2*InputData.numberDivisions+2:3*InputData.numberDivisions+2,:);
% curves{1} = PostprocessData.deformedNodes{end}([3*InputData.numberDivisions+3:4*InputData.numberDivisions+2 3*InputData.numberDivisions+2],:);
% curves{4} = curves{1};
% curves{4}(:,2) = curves{1}(end,2);
% curves{5} = curves{3};
% curves{5}(:,2) = curves{3}(end,2);
% curves{6} = [curves{2}(1,:);curves{1}(1,:)];
% curves{7} = [curves{2}(1,:);curves{3}(1,:)];
% curves{8} = [curves{1}(1,:);curves{4}(1,:)];
% curves{9} = [curves{3}(1,:);curves{5}(1,:)];
% curves{10} = [curves{4}(1,:);curves{5}(1,:)];

%% Write AutoCAD script file
% run('acad_script.m');

%% Open AUTOCAD CoreConsole and run script for creating a 3D solid
% !"C:\Program Files\Autodesk\AutoCAD 2022\accoreconsole.exe" /i "D:\RESEARCH\Navy_Boat_Project\CCO_MATLAB_Code_for_Planing_Hull\CAD_Files\Drawing1.dwg" /s "D:\RESEARCH\Navy_Boat_Project\CCO_MATLAB_Code_for_Planing_Hull\CAD_Files\Planing_Hull_CAD_Script.scr"

%% Add aft perpendicular value to geom_prop file
% filepath = sprintf("CAD_Files/Geom_Prop_hTip_%d.mpr", item);
% fid = fopen(filepath, 'a');
% fprintf(fid, 'aftPerp: %f\n', curves{1}(end, 1) - curves{1}(1, 1));
% fprintf(fid, 'hullLength: %f\n', curves{1}(end, 1));
% fclose(fid);

%% Setup & Run POWERSEA Simulation
% !python pwrs.py

% Finish timing the script
toc
