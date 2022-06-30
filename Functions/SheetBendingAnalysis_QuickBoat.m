%% Sheet bending analysis
theta = pi - PostprocessData.bendingAngle(:, end);
E = InputData.elasticModulus;
t = InputData.thickness;
L = PreprocessData.L(1 : size(theta, 1));
bendingEnergy = PostprocessData.bendingEnergy(:, end);

E_Al = 69000; %Mpa
sigmaY_Al = 276; %Mpa

%% Side Panel
thetaSide = theta(1 : 2*InputData.numberDivisions - 2);
LSide = L(1 : 2*InputData.numberDivisions - 2);
bendingEnergySide = bendingEnergy(1 : 2*InputData.numberDivisions - 2);
rhoSide = (24 / (E * t^3)) .* (bendingEnergySide ./ (thetaSide .* LSide));
thetaSideDiag = zeros(size(thetaSide,1) / 2, 1);
thetaSideStraight = thetaSideDiag;
LSideDiag = thetaSideDiag;
LSideStraight = LSideDiag;
bendingEnergySideDiag = LSideDiag;
bendingEnergySideStraight = LSideDiag;
j = 1;
k = 1;
for i = 1 : size(thetaSide, 1)
    if i == 2*InputData.numberDivisions - 3
        thetaSideDiag(j) = thetaSide(i);
        LSideDiag(j) = LSide(i);
        bendingEnergySideDiag(j) = bendingEnergySide(i);
        j = j + 1;
    elseif i == 2*InputData.numberDivisions - 2
            thetaSideStraight(k) = thetaSide(i);
            LSideStraight(k) = LSide(i);
            bendingEnergySideStraight(k) = bendingEnergySide(i);
            k = k + 1;
    else
        if rem(ceil(i / 2), 2) == 0
            thetaSideStraight(k) = thetaSide(i);
            LSideStraight(k) = LSide(i);
            bendingEnergySideStraight(k) = bendingEnergySide(i);
            k = k + 1;
        else
            thetaSideDiag(j) = thetaSide(i);
            LSideDiag(j) = LSide(i);
            bendingEnergySideDiag(j) = bendingEnergySide(i);
            j = j + 1;
        end
    end
end

% Figure for variation of bending hinge angle
% figure(1)
% set(1, "WindowStyle", "docked")
% plot(thetaSideDiag, "LineWidth", 1.5); hold on;
% plot(thetaSideStraight, "LineWidth", 1.5);
% xlabel("Bending Element")
% ylabel("Bending Hinge Angle, \theta (rad)");
% title("Bending angle along the length of side panel");
% legend("Diagonal Member", "Straight Member");
% grid on;

% Figure for variation of absolute bending hinge angle
figure(2)
% set(2, "WindowStyle", "docked")
plot(abs(thetaSideDiag), "LineWidth", 1.5); hold on;
plot(abs(thetaSideStraight), "LineWidth", 1.5);
xlabel("Bending Element")
ylabel("Bending Hinge Angle, |\theta| (rad)");
title("Bending angle along the length of side panel");
legend("Diagonal Member", "Straight Member");
grid on;

figure(3)
plot(abs(rhoSide) .* (E_Al * 0.0005 * 1000), "LineWidth", 1.5); hold on;
yline(sigmaY_Al, "LineWidth", 1.5); hold off;
xlabel("Bending elements");
ylabel("Stress in sheet (MPa)");
title("Sheet yielding analysis for side panel (t = 1mm)");
legend("Sheet Stress", "Yield Stress");
grid on;

%% Hull Panel
thetaHull = theta(2*InputData.numberDivisions - 1 : ...
                  4*InputData.numberDivisions - 4);
LHull = L(2*InputData.numberDivisions - 1 : ...
          4*InputData.numberDivisions - 4);
bendingEnergyHull = bendingEnergy(2*InputData.numberDivisions - 1 : ...
                                  4*InputData.numberDivisions - 4);
rhoHull = (24 / (E * t^3)) .* (bendingEnergyHull ./ (thetaHull .* LHull));
thetaHullDiag = zeros(size(thetaHull,1) / 2, 1);
thetaHullStraight = thetaHullDiag;
LHullDiag = thetaHullDiag;
LHullStraight = LHullDiag;
bendingEnergyHullDiag = LHullDiag;
bendingEnergyHullStraight = LHullDiag;
j = 1;
k = 1;
for i = 1 : size(thetaSide, 1)
    if i == 1
        thetaHullDiag(j) = thetaHull(i);
        LHullDiag(j) = LHull(i);
        bendingEnergyHullDiag(j) = bendingEnergyHull(i);
        j = j + 1;
    elseif i == 2    
        thetaHullStraight(k) = thetaHull(i);
        LHullStraight(k) = LHull(i);
        bendingEnergyHullStraight(k) = bendingEnergyHull(i);
        k = k + 1;
    else
        if rem(ceil(i / 2), 2) == 0
            thetaHullDiag(j) = thetaHull(i);
            LHullDiag(j) = LHull(i);
            bendingEnergyHullDiag(j) = bendingEnergyHull(i);
            j = j + 1;
        else
            thetaHullStraight(k) = thetaHull(i);
            LHullStraight(k) = LHull(i);
            bendingEnergyHullStraight(k) = bendingEnergyHull(i);
            k = k + 1;
        end
    end
end

% Figure for variation of bending hinge angle
% figure(3)
% % set(3, "WindowStyle", "docked")
% plot(thetaHullDiag, "LineWidth", 1.5); hold on;
% plot(thetaHullStraight, "LineWidth", 1.5);
% xlabel("Bending Element")
% ylabel("Bending Hinge Angle, \theta (rad)");
% title("Bending angle along the length of hull (bottom panel)");
% legend("Diagonal Member", "Straight Member");
% grid on;

% Figure for variation of absolute bending hinge angle
figure(4)
% set(4, "WindowStyle", "docked")
plot(abs(thetaHullDiag), "LineWidth", 1.5); hold on;
plot(abs(thetaHullStraight), "LineWidth", 1.5);
xlabel("Bending Element")
ylabel("Bending Hinge Angle, |\theta| (rad)");
title("Bending angle along the length of hull (bottom panel)");
legend("Diagonal Member", "Straight Member");
grid on;

figure(3)
plot(abs(rhoHull) .* (E_Al * 0.0005 * 1000), "LineWidth", 1.5); hold on;
yline(sigmaY_Al, "LineWidth", 1.5); hold off;
xlabel("Bending elements");
ylabel("Stress in sheet (MPa)");
title("Sheet yielding analysis for hull panel (t = 1mm)");
legend("Sheet Stress", "Yield Stress");
grid on;
%%

y = 0.0004/2; %m (thickness / 2)
rho = (24 / (E * t^3)) .* (bendingEnergy ./ (theta .* L));
sheetThickness = [0.0004, 0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003];
stressSide = (E_Al * (sheetThickness ./ 2) * 1000) * prctile(rhoSide, 95);
stressHull = (E_Al * (sheetThickness ./ 2) * 1000) * prctile(rhoHull, 95);


% figure(7)
% plot(abs(rho(1:size(rho)/2) * (E_Al * 0.0005 * 1000)), "LineWidth" , 1.5); hold on;
% yline(sigmaY_Al, "LineWidth" , 1.5); hold off;
% xlabel("Bending elements of the boat");
% ylabel("Stress, \sigma (MPa)");
% title("Stress at location of bending element");

figure(5)
% set(5, "WindowStyle", "docked")
plot(abs(stressSide), "LineWidth" , 1.5);hold on;
plot(abs(stressHull), "LineWidth" , 1.5);
yline(sigmaY_Al, "LineWidth" , 2); hold off;
xticks([1 2 3 4 5 6 7])
xticklabels({"0.4", "0.5", "1", "1.5", "2", "2.5", "3"})
xlabel("Sheet thickness (mm)")
ylabel("Stress in Sheet (Mpa)")
title("Sheet yielding for varying sheet thicknesses")
legend("Side", "Hull", "Yield Stress")
grid on;

%%
length  = 0.93;
widthHull = 0.15;
widthSide = 0.08;
constantSide = (mean(rhoSide) * E_Al * 10e9 / 3) * widthSide / length;
constantHull = (mean(rhoHull) * E_Al * 10e9 / 3) * widthHull / length;
forceSide = abs(constantSide' * sheetThickness.^3);
forceHull = abs(constantHull' * sheetThickness.^3);


figure(6)
% set(6, "WindowStyle", "docked")
plot(forceSide, "LineWidth" , 1.5); hold on;
plot(forceHull, "LineWidth", 1.5);
xticks([1 2 3 4 5 6 7])
xticklabels({"0.4", "0.5", "1", "1.5", "2", "2.5", "3"})
xlabel("Sheet thickness (mm)")
ylabel("Approx. Force (N)")
title("Actuation force required vs sheet thickness")
legend("Side", "Hull")

%% max curvature for sheet along diag on the side
[thetaMaxSideDiag, idx] = max(abs(thetaSideDiag));
rhoMaxSideDiag = (24 / (E * t^3)) * (bendingEnergySideDiag(idx) / ...
                 (thetaMaxSideDiag * LSideDiag(idx)));
stressSideDiag = E_Al * y * rhoMaxSideDiag * 1000;

% max curvature for sheet along straight on the side
[thetaMaxSideStraight, idx] = max(abs(thetaSideStraight));
rhoMaxSideStraight = (24 / (E * t^3)) * (bendingEnergySideStraight(idx) / ...
                     (thetaMaxSideStraight * LSideStraight(idx)));
stressSideStraight = E_Al * y * rhoMaxSideStraight * 1000;

% max curvature for sheet along diag on the hull
[thetaMaxHullDiag, idx] = max(abs(thetaHullDiag));
rhoMaxHullDiag = (24 / (E * t^3)) * (bendingEnergyHullDiag(idx) / ...
                 (thetaMaxHullDiag * LHullDiag(idx)));
stressHullDiag = E_Al * y * rhoMaxHullDiag * 1000;

% max curvatue for sheet along straight on the hull
[thetaMaxHullStraight, idx] = max(abs(thetaHullStraight));
rhoMaxHullStraight = (24 / (E * t^3)) * (bendingEnergyHullStraight(idx) / ...
                     (thetaMaxHullStraight * LHullStraight(idx)));
stressHullStraight = E_Al * y * rhoMaxHullStraight * 1000;

function result = secondLarge(vector)
vector = sort(abs(vector), 'descend');
result = vector(2);
end
