function [InputData, inputTest] = Test_Input(InputData)

    inputTest = 1;
           
    if ~isnumeric(InputData.restAngle1)
        fprintf('ERROR: <restAngle1> must be a number\n\n')
        inputTest = 0;
        restAngleRad1 = nan;
    elseif ~isreal(InputData.restAngle1)
        fprintf('ERROR: <restAngle1> must be a real number\n\n')
        inputTest = 0;
        restAngleRad1 = nan;
    elseif InputData.restAngle1 < -180 || InputData.restAngle1 > 180
            fprintf('ERROR: <restAngle1> must be between -180 and 180 [deg]\n\n')
            inputTest = 0;
            restAngleRad1 = nan;
    elseif InputData.restAngle1 < 0
            restAngleRad1 = -(180 - abs(InputData.restAngle1))*pi/180;
    else
            restAngleRad1 = (180 - InputData.restAngle1)*pi/180;
    end

    if ~isnumeric(InputData.restAngle2)
        fprintf('ERROR: <restAngle2> must be a number\n\n')
        inputTest = 0;
        restAngleRad2 = nan;
    elseif ~isreal(InputData.restAngle2)
        fprintf('ERROR: <restAngle2> must be a real number\n\n')
        inputTest = 0;
        restAngleRad2 = nan;
    elseif InputData.restAngle2 < -180 || InputData.restAngle2 > 180
            fprintf('ERROR: <restAngle2> must be between -180 and 180 [deg]\n\n')
            inputTest = 0;
            restAngleRad2 = nan;
    elseif InputData.restAngle2 < 0
            restAngleRad2 = -(180 - abs(InputData.restAngle2))*pi/180;
    else
            restAngleRad2 = (180 - InputData.restAngle2)*pi/180;
    end
    
    if ~isnumeric(InputData.elasticModulus)
        fprintf('ERROR: <elasticModulus> must be a number\n\n')
        inputTest = 0;
    elseif ~isreal(InputData.elasticModulus)
        fprintf('ERROR: <elasticModulus> must be a real number\n\n')
        inputTest = 0;
    elseif InputData.elasticModulus < 0
        fprintf('ERROR: <elasticModulus> must be greater than zero\n\n')
        inputTest = 0;
    end

    if ~isnumeric(InputData.thickness)
        fprintf('ERROR: <thickness> must be a number\n\n')
        inputTest = 0;
    elseif ~isreal(InputData.thickness)
        fprintf('ERROR: <thickness> must be a real number\n\n')
        inputTest = 0;
    elseif InputData.thickness < 0
        fprintf('ERROR: <thickness> must be greater than zero\n\n')
        inputTest = 0;
    end

    if ~isnumeric(InputData.lengthScale1)
        fprintf('ERROR: <lengthScale1> must be a number\n\n')
        inputTest = 0;
    elseif ~isreal(InputData.lengthScale1) 
        fprintf('ERROR: <lengthScale1> must be a real number\n\n')
        inputTest = 0;
    elseif InputData.lengthScale1 < 0
        fprintf('ERROR: <lengthScale1> must be greater than zero\n\n')
        inputTest = 0;
    end

    if ~isnumeric(InputData.lengthScale2)
        fprintf('ERROR: <lengthScale2> must be a number\n\n')
        inputTest = 0;
    elseif ~isreal(InputData.lengthScale2) 
        fprintf('ERROR: <lengthScale2> must be a real number\n\n')
        inputTest = 0;
    elseif InputData.lengthScale2 < 0
        fprintf('ERROR: <lengthScale2> must be greater than zero\n\n')
        inputTest = 0;
    end    
    
    if ~isnumeric(InputData.loadMagnitude)
        fprintf('ERROR: <loadMagnitude> must be a number\n\n')
        inputTest = 0;
    elseif ~isreal(InputData.loadMagnitude)
        fprintf('ERROR: <loadMagnitude> must be a real number\n\n')
        inputTest = 0;
    end

%     if ~ischar(InputData.testType)
%         fprintf('ERROR: <testType> must be one of the following:\n\n\t\tFullAnnulus\n\t\tSineWave\n\t\tCurvedSquare\n\t\tCanopy\n\t\tCutAnnulus\n\n')
%         inputTest = 0;
%     elseif ~strcmp(InputData.testType,'FullAnnulus') && ~strcmp(InputData.testType,'SineWave') && ~strcmp(InputData.testType,'CurvedSquare') && ~strcmp(InputData.testType,'Canopy') && ~strcmp(InputData.testType,'CutAnnulus')
%         fprintf('ERROR: <testType> must be one of the following:\n\n\t\tFullAnnulus\n\t\tSineWave\n\t\tCurvedSquare\n\t\tCanopy\n\t\tCutAnnulus\n\n')
%         inputTest = 0;
%     end

    if ~isnumeric(InputData.numberDivisions) || ceil(InputData.numberDivisions) ~= InputData.numberDivisions
        fprintf('ERROR: <numberDivisions> must be an integer\n\n')
        inputTest = 0;
    elseif InputData.numberDivisions < 0
        fprintf('ERROR: <numberDivisions> must be a positive integer\n\n')
        inputTest = 0;
    end

    if ~ischar(InputData.plotNodes)
        fprintf('ERROR: <plotNodes> must be a character\n\n')
        inputTest = 0;
    elseif ~strcmpi(InputData.plotNodes,'yes') && ~strcmpi(InputData.plotNodes,'no')
        fprintf('ERROR: <plotNodes> must be either "yes" or "no" (not case-sensitive)\n\n')
        inputTest = 0;
    end
 
    if ~ischar(InputData.plotEnergy)
        fprintf('ERROR: <plotEnergy> must be a character\n\n')
        inputTest = 0;
    elseif ~strcmpi(InputData.plotEnergy,'yes') && ~strcmpi(InputData.plotEnergy,'no')
        fprintf('ERROR: <plotEnergy> must be either "yes" or "no" (not case-sensitive)\n\n')
        inputTest = 0;
    end

    if ~ischar(InputData.plotDeformedShape)
        fprintf('ERROR: <plotDeformedShape> must be a character\n\n')
        inputTest = 0;
    elseif ~strcmpi(InputData.plotDeformedShape,'yes') && ~strcmpi(InputData.plotDeformedShape,'no')
        fprintf('ERROR: <plotDeformedShape> must be either "yes" or "no" (not case-sensitive)\n\n')
        inputTest = 0;
    end

    if ~ischar(InputData.plotReferenceShape)
        fprintf('ERROR: <plotReferenceShape> must be a character\n\n')
        inputTest = 0;
    elseif ~strcmpi(InputData.plotReferenceShape,'yes') && ~strcmpi(InputData.plotReferenceShape,'no')
        fprintf('ERROR: <plotReferenceShape> must be either "yes" or "no" (not case-sensitive)\n\n')
        inputTest = 0;
    end
    
    if ischar(InputData.plotIncrement)
        if ~strcmpi(InputData.plotIncrement,'end')
            fprintf('ERROR: <plotIncrement> must be an integer between 1 and the specified <incrementNumber> or "end"\n\n')
            inputTest = 0;
        end
    elseif InputData.plotIncrement > InputData.numberIncrements
        fprintf('ERROR: <plotIncrement> might exceed available number of increments; proceed with caution or choose "end" instead\n\n')
        if inputTest == 1
            inputTest = 2;
        end
    end
    
    InputData.restAngleRad1 = restAngleRad1;
    InputData.restAngleRad2 = restAngleRad2;
 
    
           
    
    
    
