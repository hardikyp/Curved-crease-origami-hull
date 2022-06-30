function [PostprocessData] = Path_Analysis_Planing(InputData, PreprocessData)
tol = 1e-6; 
truss.U0 = zeros(3*size(InputData.nodes,1),1);
U = truss.U0;

f_inc1 = InputData.restAngleRad1/InputData.numberIncrements; % Fold angle increment size [rad]
f_inc2 = InputData.restAngleRad2/InputData.numberIncrements; % Fold angle increment size [rad]

angles.Kf = PreprocessData.Kf; % Fold stiffness
angles.Kb = PreprocessData.Kb; % Bending stiffness
truss.A = PreprocessData.barAreas; % Bar cross-sectional areas

angles.CMbend = @(he,h0,kb,L0)EnhancedLinear(he,h0,kb,L0,0,360);
angles.CMfold = @(he,h0,kf,L0)EnhancedLinear(he,h0,kf,L0,0,360);
truss.CM = @(Ex)Ogden(Ex, InputData.elasticModulus); % Constitutive model

angles.fold = PreprocessData.foldingHinges;
angles.bend = PreprocessData.bendingHinges;
truss.Bars = PreprocessData.trussBars;

truss.Node = InputData.nodes;

truss.B = PreprocessData.B; 
truss.L = PreprocessData.L;

angles.pf0 = PreprocessData.initialFoldAngles;
angles.pb0 = PreprocessData.initialBendAngles;
angles.Panel = InputData.panels;

truss.FixedDofs = PreprocessData.fixedDOFs; 
truss.Trigl = PreprocessData.Trigl;

if strcmpi(InputData.loadType, 'Force')
    MaxIcr = InputData.numberIncrements;
%     b_lambda = InputData.loadMagnitude/InputData.numberIncrements;
    b_lambda = 1/InputData.numberIncrements;
    Uhis = zeros(3*size(InputData.nodes,1),MaxIcr);
    FreeDofs = setdiff(1:3*size(InputData.nodes,1),truss.FixedDofs);
    lmd = 0; icrm = 0; MUL = [U,U];
    lamHis = zeros(InputData.numberIncrements,1);
    angles.pf0_Orig=angles.pf0;
    F = PreprocessData.load;
    while icrm<MaxIcr && lmd < 1 % Analysis stops when applied load is below assigned load (\lambda < 1)
        icrm = icrm+1;
        %% Equivalent decrement/increment in angles
        if InputData.foldStep
            angles.pf0(PreprocessData.mountainFolds)=angles.pf0_Orig(PreprocessData.mountainFolds) - f_inc1 * icrm;
            angles.pf0(PreprocessData.valleyFolds)=angles.pf0_Orig(PreprocessData.valleyFolds) + f_inc2 * icrm;
        end
        iter = 0; err = 1;
        fprintf('icrm = %d, lambda = %6.4f\n',icrm,lmd);
        while err>tol && iter<InputData.maxIterations
            iter = iter+1; 
            [IF,K] = GlobalK_fast_ver(U,InputData.nodes,truss,angles);
            R = lmd*F-IF;   
            MRS = [F,R];
            MUL(FreeDofs,:) = K(FreeDofs,FreeDofs)\MRS(FreeDofs,:);
            dUp = MUL(:,1); 
            dUr = MUL(:,2);
            if iter==1
                dUr = 0*dUr; 
            end
            % Modified Generalized Displacement Control Method]
            if iter==1
                if icrm==1
                    sinal=sign(dot(dUp,dUp));
                    dlmd=b_lambda;
                    numgsp=dot(dUp,dUp);   
                else
                    sinal=sinal*sign(dot(dupp1,dUp));
                    gsp=numgsp/dot(dUp,dUp);
                    dlmd=sinal*b_lambda*sqrt(gsp);
                end 
                dupp1=dUp;
                dupc1=dUp;
            else
                dlmd=-dot(dupc1,dUr)/dot(dupc1,dUp);
            end
            dUt = dlmd*dUp+dUr;
            U = U+dUt;
            err = norm(dUt(FreeDofs));
            lmd = lmd+dlmd;
            fprintf('    iter = %d, err = %6.4f, dlambda = %6.4f\n',iter,err,dlmd);
            if err > 1e8
                disp('Divergence!')
                break
            end
        end

        if iter>15
            b_lambda = b_lambda/2;
            disp('Reduce constraint radius...')
            icrm = icrm-1;
            U = Uhis(:,max(icrm,1));  % restore displacement
            lmd = lamHis(max(icrm,1));   % restore load
        elseif iter<3
            disp('Increase constraint radius...')
            b_lambda = b_lambda*1.5;
            Uhis(:,icrm) = U;
            lamHis(icrm) = lmd; 
        else
            Uhis(:,icrm) = U;
            lamHis(icrm) = lmd; 
        end
    end

elseif strcmpi(InputData.loadType, 'Displacement')
    Uhis = zeros(3*size(InputData.nodes,1),InputData.numberIncrements*2);
    Fdsp = PreprocessData.load/InputData.numberIncrements;
    ImpDofs = find(Fdsp~=0);
    FreeDofs = setdiff(setdiff(1:3*size(InputData.nodes,1),truss.FixedDofs),ImpDofs);
    icrm = 0;  
    dspmvd = 0;  
    attmpts = 0;
    mvstepsize = 1;  
    damping = 1;
    Fhis = zeros(InputData.numberIncrements,numel(ImpDofs)); 
    angles.pf0_Orig=angles.pf0;
    while dspmvd <= 1 && attmpts <= 20      
        icrm = icrm+1;
        %% Equivalent decrement/increment in angles
        if InputData.foldStep && icrm <= InputData.numberIncrements
            angles.pf0(PreprocessData.mountainFolds)=angles.pf0_Orig(PreprocessData.mountainFolds) - f_inc1 * icrm;
            angles.pf0(PreprocessData.valleyFolds)=angles.pf0_Orig(PreprocessData.valleyFolds) + f_inc2 * icrm;
        end
        iter = 0; err = 1;   
        fprintf('icrm = %d, dspimps = %6.4f\n',icrm,dspmvd);
        U = U+mvstepsize*Fdsp;
        U(truss.FixedDofs)=0;
        while err>tol && iter<InputData.maxIterations
            iter = iter+1;
            [IF,K] = GlobalK_fast_ver(U,InputData.nodes,truss,angles);
            dU = zeros(3*size(InputData.nodes,1),1);
            dU(FreeDofs) = K(FreeDofs,FreeDofs)\(-IF(FreeDofs));
            err = norm(dU(FreeDofs));
            U = U+damping*dU; 
            fprintf('    iter = %d, err = %6.4f\n',iter,err);
        end

        if iter>=((mvstepsize>1)+1)*InputData.maxIterations/(damping+1)  
            % an aggressive step needs more iterations
            attmpts = attmpts+1;
            icrm = icrm-1;
            if attmpts<=10  
                mvstepsize = mvstepsize*0.5; 
                disp('Take a more conservative step...')
            else
                mvstepsize = max(mvstepsize,1)*1.5;  
                damping = damping*0.75;
                disp('Take a more aggressive step...')
            end
            U = Uhis(:,max(icrm,1)); % restore displacement            
        else
            dspmvd = dspmvd+mvstepsize/InputData.numberIncrements;
            attmpts = 0;
            damping = 1;
            if mvstepsize<1
                mvstepsize = min(mvstepsize*1.1,1); % gradually go back to 1
            else
                mvstepsize = max(mvstepsize*0.9,1);
            end
            Uhis(:,icrm) = U;
            [Fend,~] = GlobalK_fast_ver(U,InputData.nodes,truss,angles);
            Fhis(icrm,:) = Fend(ImpDofs)'; 
        end
    end
else
    disp('Unknown load type!!!')
end

if strcmpi(InputData.loadType,'Force') % lamHis gives the proportion of the load, not the actual load
    Fhis = lamHis*nonzeros(F)';
    lamHis(icrm+1:end,:) = [];
    PostprocessData.forceAchieved = [InputData.loads(:,1), lamHis(end)*InputData.loads(:,2:4)];
end

Uhis(:,icrm+1:end) = [];
Fhis(icrm+1:end,:) = [];

PostprocessData.Uhis = real(Uhis);
PostprocessData.Fhis = real(Fhis);
PostprocessData.foldAngles = angles.pf0;
if strcmpi(InputData.loadType,'Force')
    PostprocessData.lamHis = lamHis;
end
end
