%% Sheet bending analysis
% load("curvatureAnalysis.mat");
bendingAngles = PostprocessData.bendingAngle(2 * InputData.numberDivisions: ...
                                             4 * InputData.numberDivisions - 3, ...
                                             end);
theta = pi - bendingAngles;
bendingEnergy = PostprocessData.bendingEnergy(2 * InputData.numberDivisions: ...
                                              4 * InputData.numberDivisions - 3, ...
                                              end);
L = PreprocessData.L(2 * InputData.numberDivisions: ...
                     4 * InputData.numberDivisions - 3);
barAreas = PreprocessData.barAreas(2 * InputData.numberDivisions: ...
                                   4 * InputData.numberDivisions - 3, end);
E = InputData.elasticModulus;
t = InputData.thickness;
E_Al = 69000; %MPa
y = 0.001/2;
sigmaY_Al = 276;

rho = (24 / (E * t^3)) .* (bendingEnergy ./ (theta .* L));
stress = (E_Al * y * 1000) .* rho;

thetaDiag = zeros(size(bendingAngles,1)/2,1);
thetaStraight = thetaDiag;
LDiag = thetaDiag;
LStraight = LDiag;
bendingEnergyDiag = thetaDiag;
bendingEnergyStraight = bendingEnergyDiag;

j = 1;
k = 1;
for i=1:size(theta)
    if i == 1
        thetaDiag(j) = theta(i);
        LDiag(j) = L(i);
        bendingEnergyDiag(j) = bendingEnergy(i);
        j = j+1;
    elseif i == 2
        thetaStraight(k) = theta(i);
        LStraight(k) = L(i);
        bendingEnergyStraight(k) = bendingEnergy(i);
        k = k + 1;
    else
        if rem(ceil(i/2),2) == 0 
            thetaDiag(j) = theta(i);
            LDiag(j) = L(i);
            bendingEnergyDiag(j) = bendingEnergy(i);
            j = j + 1;
        else
            thetaStraight(k) = theta(i);
            LStraight(k) = L(i);
            bendingEnergyStraight(k) = bendingEnergy(i);
            k = k + 1;
        end
    end
end


% figure(1)
% set(1, "WindowStyle", "docked")
% plot(thetaDiag, "LineWidth", 1.5); hold on;
% plot(thetaStraight, "LineWidth", 1.5); hold off;
% xlabel("Bending element");
% ylabel("\theta (rad)");
% title("Bending angle along the length of the hull")
% legend("Diagonal","Straight");
% grid on;

figure(2)
% set(2, "WindowStyle", "docked")
plot(abs(thetaDiag), "LineWidth", 1.5); hold on;
plot(abs(thetaStraight), "LineWidth", 1.5); hold off;
xlabel("Bending element");
ylabel("|\theta| (rad)");
title("Absolute bending angle along the length of the hull")
legend("Diagonal","Straight");
grid on;

figure(3)
% set(3, "WindowStyle", "docked")
plot(abs(stress), "LineWidth", 1.5); hold on;
yline(sigmaY_Al, "LineWidth", 1.5); hold off;
xlabel("Bending Elements");
ylabel("\sigma (MPa)");
legend("Bending Elements", "Yield Stress");
title("Stress in the Aluminium sheet at each bending element, t = 1mm")
grid on;
%%
sheetThickness = [0.0004 0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003]; %mm
stress = (E_Al * (sheetThickness ./ 2) * 1000) * prctile(rho, 95);
figure(4)
% set(4, "WindowStyle", "docked")
plot(abs(stress), "LineWidth" , 1.5); hold on;
yline(sigmaY_Al, "LineWidth", 1.5); hold off;
xticks([1 2 3 4 5 6 7])
xticklabels({"0.4", "0.5", "1", "1.5", "2", "3.5", "3"})
xlabel("Sheet thickness (mm)")
ylabel("Stress in Sheet (Mpa)")
title("Sheet yielding for various sheet thicknesses")
legend("Hull", "Yield Stress")
grid on;

length  = 0.93;
width = 0.15;
constant = (mean(rho, "omitnan") * E_Al * 10e9 / 3) * width / length;
force = constant' * sheetThickness.^3;


figure(6)
% set(6, "WindowStyle", "docked")
plot(force, "LineWidth" , 1.5); hold on;
xticks([1 2 3 4 5 6 7])
xticklabels({"0.4", "0.5", "1", "1.5", "2", "2.5", "3"})
xlabel("Sheet thickness (mm)")
ylabel("Approx. Force (N)")
title("Actuation force required vs sheet thickness");
grid on;