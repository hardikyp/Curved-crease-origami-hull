%% GENERATES SCRIPT FILE FOR AUTOCAD TO CREATE SOLID PLANING GEOMETRY
fid = fopen('CAD_Files/Planing_Hull_CAD_Script.scr', 'wt');

numSplines = 9;

curvComb = [4, 1;
            1, 2;
            2, 3;
            3, 5;
            7, 6;
            9, 8;
            4, 5];

numSurfaces = size(curvComb, 1);

fprintf(fid, 'ATTDIA\n');
fprintf(fid, '0\n');
fprintf(fid, 'CMDDIA\n');
fprintf(fid, '0\n');
fprintf(fid, 'FILEDIA\n');
fprintf(fid, '0\n');

% SET DRAWING UNITS
fprintf(fid, "-DWGUNITS\n");
fprintf(fid, "3\n");
fprintf(fid, "\n");
fprintf(fid, "\n");
fprintf(fid, "\n");
fprintf(fid, "\n");

% GENERATE ALL SPLINES
for j = 1:numSplines    
    %GENERATE SPLINE    
    fprintf(fid, 'SPLINE\n');
    for i = 1:size(curves{j}, 1)
        fprintf(fid, '%d,%d,%d\n', curves{j}(i, :));
    end
    fprintf(fid, '\n');
    fprintf(fid, '\n');    
    fprintf(fid, '\n');
    
    % DEFINE VARIABLE NAME 'SPLINE#' FOR LAST DRAWN SPLINE
    fprintf(fid, '(setq spline%d (entlast))\n', j);
end

% LOFT SURFACE BETWEEN ADJACENT SPLINES
for i = 1:numSurfaces
    fprintf(fid, 'LOFT\n');
    % CREATE OBJECT SELECTION
    fprintf(fid, '!spline%d\n', curvComb(i, 1));
    fprintf(fid, '!spline%d\n', curvComb(i, 2));
    fprintf(fid, '\n\n');

    % DEFINE VARIABLE NAME 'SURFACE#' FOR LAST LOFTED SURFACE
    fprintf(fid, '(setq surface%d (entlast))\n', i);
end      

% PERFORM SURFACE SCULPT TO GET A 3D SOLID
fprintf(fid, "SURFSCULPT\n");
for i = 1:numSurfaces
    fprintf(fid, '!surface%d\n', i);
end
fprintf(fid, '\n');

% SET VARIABLE NAME OF THE NEWLY CREATED 3D SOLID TO 'solidboat'
fprintf(fid, '(setq solidboat (entlast))\n');

% ROTATE THE 3D SOLID ALONG X AXIS SO THAT THE TOP IS IN +Z DIRECTION
fprintf(fid, 'ROTATE3D\n');
fprintf(fid, 'ALL\n\n');
fprintf(fid, 'Xaxis\n');
fprintf(fid, '0,0,0\n');
fprintf(fid, '90\n');

% ROTATE THE 3D SOLID AND SPLINES ALONG Z AXIS SO THAT IT FACES -VE X DIRECTION
fprintf(fid, 'ROTATE3D\n');
fprintf(fid, 'ALL\n\n');
fprintf(fid, 'Zaxis\n');
fprintf(fid, '0,0,0\n');
fprintf(fid, '180\n');

% MOVE THE BOAT AND SPLINES IN X DIRECTION TO ALIGN IT AS PER POWERSEA REQUIREMENTS
fprintf(fid, "MOVE\n");
fprintf(fid, "ALL\n\n");
fprintf(fid, "0,0,0\n");
fprintf(fid, "%f,0,0\n", curves{4}(end, 1));

% FIND MASS PROPERTIES FOR THE GIVEN THICKNESS OF THE 3D SOLID
fprintf(fid, 'MASSPROP\n');
fprintf(fid, '!solidboat\n\n');
fprintf(fid, 'Yes\n');
% fprintf(fid, '"D:\\RESEARCH\\Navy_Boat_Project\\CCO_MATLAB_Code_for_Planing_Hull\\CAD_Files\\Geom_Prop.mpr"\n');
fprintf(fid, '"D:\\RESEARCH\\Navy_Boat_Project\\CCO_MATLAB_Code_for_Planing_Hull\\CAD_Files\\Geom_Prop_hTip_%d.mpr"\n', item);
% ERASE EVERYTHING
fprintf(fid, 'ERASE\n');
fprintf(fid, 'ALL\n');
fprintf(fid, '\n');

% GENERATE NECESSARY SPLINES
for j = 1:2    
    %GENERATE SPLINE    
    fprintf(fid, 'SPLINE\n');
    for i = 1:size(curves{j}, 1)
        fprintf(fid, '%d,%d,%d\n', curves{j}(i, :));
    end
    fprintf(fid, '\n');
    fprintf(fid, '\n');    
    fprintf(fid, '\n');
    
    % DEFINE VARIABLE NAME 'SPLINE#' FOR LAST DRAWN SPLINE
    fprintf(fid, '(setq spline%d (entlast))\n', j);
end

% ROTATE THE 3D SOLID ALONG X AXIS SO THAT THE TOP IS IN +Z DIRECTION
fprintf(fid, 'ROTATE3D\n');
fprintf(fid, 'ALL\n\n');
fprintf(fid, 'Xaxis\n');
fprintf(fid, '0,0,0\n');
fprintf(fid, '90\n');

% ROTATE THE 3D SOLID AND SPLINES ALONG Z AXIS SO THAT IT FACES -VE X DIRECTION
fprintf(fid, 'ROTATE3D\n');
fprintf(fid, 'ALL\n\n');
fprintf(fid, 'Zaxis\n');
fprintf(fid, '0,0,0\n');
fprintf(fid, '180\n');

% MOVE THE BOAT AND SPLINES IN X DIRECTION TO ALIGN IT AS PER POWERSEA REQUIREMENTS
fprintf(fid, "MOVE\n");
fprintf(fid, "ALL\n\n");
fprintf(fid, "0,0,0\n");
fprintf(fid, "%f,0,0\n", curves{4}(end, 1));

% SAVE AUTOCAD DWG
fprintf(fid, 'SAVEAS\n');
fprintf(fid, '2018\n');
fprintf(fid, '"D:\\RESEARCH\\Navy_Boat_Project\\CCO_MATLAB_Code_for_Planing_Hull\\CAD_Files\\Planing_Hull_CCO_hTip_%d.dwg"\n', item);
% fprintf(fid, '"D:\\RESEARCH\\Navy_Boat_Project\\CCO_MATLAB_Code_for_Planing_Hull\\CAD_Files\\Planing_Hull_CCO.dwg"\n');

% EXPORT THE GENERATED SOLID AS ACIS (.igs) FILE
fprintf(fid, '_EXPORT\n');
fprintf(fid, '"D:\\RESEARCH\\Navy_Boat_Project\\CCO_MATLAB_Code_for_Planing_Hull\\CAD_Files\\Planing_Hull_CCO_hTip_%d.igs"\n', item);
% fprintf(fid, '"D:\\RESEARCH\\Navy_Boat_Project\\CCO_MATLAB_Code_for_Planing_Hull\\CAD_Files\\Planing_Hull_CCO.igs"\n');
fprintf(fid, 'ALL\n\n');

fprintf(fid, 'ATTDIA\n');
fprintf(fid, '1\n');
fprintf(fid, 'CMDDIA\n');
fprintf(fid, '1\n');
fprintf(fid, 'FILEDIA\n');
fprintf(fid, '1\n');
fprintf(fid, 'QUIT\n');
fprintf(fid, 'Y\n');

fclose(fid);
clearvars i j curvComb numSurfaces numSplines fid;