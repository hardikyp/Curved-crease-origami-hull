t = 0.27432;

Pick3 = 'Yes';
Convert2mm_0 = 'No';
Convert2mm_1 = 'No';
AutoPick = 'Yes';
ShowAxes = 'No';
ColorHausdorff = 'Yes';
SampleRate_0 = 1;
SampleRate_1 = 2;
UnitScale = 1;

% Plot Positions
plot_w = 500;
plot_h = 300;
plot_s = 50;
plot_y1 = 485;
plot_y2 = 100;

%% Combine coordinates to make node matrices
if strcmp(Convert2mm_0,'Yes')
    Unit_0 = 25.4;
else
    Unit_0 = 1;
end
if strcmp(Convert2mm_1,'Yes')
    Unit_1 = 25.4;
else
    Unit_1 = 1;
end
M0 = Unit_0*[x0 y0 z0];
M1 = Unit_1*[x1 y1 z1];

%% Plot meshes for point selection
fig = figure('Position',[plot_s plot_y1 plot_w plot_h],'NumberTitle','off');
scatter3(M0(:,1), M0(:,2), M0(:,3), 'MarkerEdgeColor',[0.5 0.5 0.5]);
hold on
scatter3(M1(:,1), M1(:,2), M1(:,3), 'MarkerEdgeColor','k','Marker','*');
dcm_obj = datacursormode(fig);
set(dcm_obj,'UpdateFcn', @ShowIndex);
title('Pick Identical Nodes Between Meshes')
axis equal

%% Input points
if strcmp(AutoPick,'Yes')
    i_Po0 = 881;
    i_Pa0 = 7;
    i_Pb0 = 8;
    i_Po1 = 581;
    i_Pa1 = 33;
    i_Pb1 = 32;
else
    i_Po0 = input('P^o_0 - index of origin of anchor model: ');
    i_Pa0 = input('P^a_0 - index of first match point of anchor model: ');
    i_Pb0 = input('P^b_0 - index of second match point of anchor model: ');
    i_Po1 = input('P^o_1 - index of origin of anchor model: ');
    i_Pa1 = input('P^a_1 - index of first match point of moving model: ');
    i_Pb1 = input('P^b_1 - index of second match point of moving model: ');
end
Po0 = M0(i_Po0,:);
Pa0 = M0(i_Pa0,:);
Pb0 = M0(i_Pb0,:);
Po1 = M1(i_Po1,:);
Pa1 = M1(i_Pa1,:);
Pb1 = M1(i_Pb1,:);

% Check to make sure selected points are not colinear
AO = norm(Pa0-Po0);
AB = norm(Pa0-Pb0);
OB = norm(Po0-Pb0);

theta = acos((-AB^2+AO^2+OB^2)/(2*AO*OB));
if abs(theta - pi) < 0.15
    fprintf('Error: selected points are almost colinear\n\n');
    return
end

% Plot meshes and match points
figure('Position',[plot_s plot_y1 plot_w plot_h],'NumberTitle','off');
scatter3(M0(:,1), M0(:,2), M0(:,3), 'MarkerEdgeColor',[0.5 0.5 0.5]); hold on
scatter3(M1(:,1), M1(:,2), M1(:,3), 'k.'); hold on
scatter3(Po0(1),Po0(2),Po0(3),'md'); hold on
scatter3(Po1(1),Po1(2),Po1(3),'md','MarkerFaceColor','m'); hold on
scatter3(Pa0(1),Pa0(2),Pa0(3),'ys'); hold on
scatter3(Pa1(1),Pa1(2),Pa1(3),'ys','MarkerFaceColor','y'); hold on
scatter3(Pb0(1),Pb0(2),Pb0(3),'c^'); hold on
scatter3(Pb1(1),Pb1(2),Pb1(3),'c^','MarkerFaceColor','c'); hold on
legend('Mesh_0','Mesh_1','O_0','O_1','a_0','a_1','b_0','b_1')
title('Original Meshes')
axis equal

%% Translate both meshes to the origin

% Translate anchor mesh
Tt0 = -Po0;
Tx = Tt0(1);
Ty = Tt0(2);
Tz = Tt0(3);

N0 = size(M0,1);

M0t = [ones(N0,1) zeros(N0,2)]*Tx + [zeros(N0,1) ones(N0,1) zeros(N0,1)]*Ty +...
    [zeros(N0,2) ones(N0,1)]*Tz + M0;

% Translate moving mesh
Tt = - Po1;
Tx = Tt(1);
Ty = Tt(2);
Tz = Tt(3);

N1 = size(M1,1);

M1t = [ones(N1,1) zeros(N1,2)]*Tx + [zeros(N1,1) ones(N1,1) zeros(N1,1)]*Ty +...
    [zeros(N1,2) ones(N1,1)]*Tz + M1;

% Calculate basis of anchor mesh
Pa0 = M0t(i_Pa0,:);
Pb0 = M0t(i_Pb0,:);
ea0 = Pa0/norm(Pa0);
ec0 = cross(Pa0,Pb0);
ec0 = ec0/norm(ec0);
eb0 = cross(ec0,ea0);

% Calculate basis of moving mesh
Pa1 = M1t(i_Pa1,:);
Pb1 = M1t(i_Pb1,:);
ea1 = Pa1/norm(Pa1);
ec1 = cross(Pa1,Pb1);
ec1 = ec1/norm(ec1);
eb1 = cross(ec1,ea1);

% Plot translated meshes
figure('Position',[plot_s+plot_w plot_y1 plot_w plot_h],'NumberTitle','off')
scatter3(M1t(:,1), M1t(:,2), M1t(:,3),'Marker','.','MarkerEdgeColor','k'); hold on
scatter3(M0t(:,1), M0t(:,2), M0t(:,3),'MarkerEdgeColor',[0.5 0.5 0.5]); hold on
line(UnitScale*[0 ea1(1) ],UnitScale*[0 ea1(2) ],UnitScale*[0 ea1(3) ],'Color','r'); hold on
line(UnitScale*[0 eb1(1) ],UnitScale*[0 eb1(2) ],UnitScale*[0 eb1(3) ],'Color','g'); hold on
line(UnitScale*[0 ec1(1) ],UnitScale*[0 ec1(2) ],UnitScale*[0 ec1(3) ],'Color','b'); hold on
line(UnitScale*[0 ea0(1) ],UnitScale*[0 ea0(2) ],UnitScale*[0 ea0(3) ],'Color','r','LineStyle','--'); hold on
line(UnitScale*[0 eb0(1) ],UnitScale*[0 eb0(2) ],UnitScale*[0 eb0(3) ],'Color','g','LineStyle','--'); hold on
line(UnitScale*[0 ec0(1) ],UnitScale*[0 ec0(2) ],UnitScale*[0 ec0(3) ],'Color','b','LineStyle','--'); hold on
axis equal
legend('M_1^t','M_0^t','e^a_1','e^b_1','e^c_1','e^a_0','e^b_0','e^c_0');
title('Translated Meshes')

%% Rotate moveable mesh

% First rotation

% Calculate rotation axis vector
N = cross(ea1,ea0);
n = N/norm(N);

% Calculate rotation magnitude
x = dot(ea1,ea0);
y = dot(n,N);
phi(1) = atan2(y,x);

% Execute first rotation
M1ri = rodrigues_rot(M1t,n,phi(1));

% Second rotation

% Calculate basis of moving mesh
Pa1 = M1ri(i_Pa1,:);
Pb1 = M1ri(i_Pb1,:);
ea1 = Pa1/norm(Pa1);
ec1 = cross(Pa1,Pb1);
ec1 = ec1/norm(ec1);
eb1 = cross(ec1,ea1);

figure('Position',[plot_s+2*plot_w plot_y1 plot_w plot_h],'NumberTitle','off')
scatter3(M0t(:,1), M0t(:,2), M0t(:,3),'MarkerEdgeColor',[0.5 0.5 0.5]); hold on
scatter3(M1ri(:,1), M1ri(:,2), M1ri(:,3),'k.'); hold on
line(UnitScale*[0 ea1(1) ],UnitScale*[0 ea1(2) ],UnitScale*[0 ea1(3) ],'Color','r'); hold on
line(UnitScale*[0 eb1(1) ],UnitScale*[0 eb1(2) ],UnitScale*[0 eb1(3) ],'Color','g'); hold on
line(UnitScale*[0 ec1(1) ],UnitScale*[0 ec1(2) ],UnitScale*[0 ec1(3) ],'Color','b'); hold on
line(UnitScale*[0 ea0(1) ],UnitScale*[0 ea0(2) ],UnitScale*[0 ea0(3) ],'Color','r','LineStyle','--'); hold on
line(UnitScale*[0 eb0(1) ],UnitScale*[0 eb0(2) ],UnitScale*[0 eb0(3) ],'Color','g','LineStyle','--'); hold on
line(UnitScale*[0 ec0(1) ],UnitScale*[0 ec0(2) ],UnitScale*[0 ec0(3) ],'Color','b','LineStyle','--'); hold on
axis equal
legend('M_0^t','M_1^r^i','e^a_1','e^b_1','e^c_1','e^a_0','e^b_0','e^c_0');
xlabel('x');
ylabel('y');
zlabel('z');
title('First Rotation')

% Calculate rotation axis vector
N = cross(ec1,ec0);
n = N/norm(N);

% Calculate rotation magnitude
x = dot(ec1,ec0);
y = dot(n,N);
phi(2) = atan2(y,x);

% Execute second rotation
M1r = rodrigues_rot(M1ri,ea0,phi(2));

% Calculate basis of moving mesh
Pa1 = M1r(i_Pa1,:);
Pb1 = M1r(i_Pb1,:);
ea1 = Pa1/norm(Pa1);
ec1 = cross(Pa1,Pb1);
ec1 = ec1/norm(ec1);
eb1 = cross(ec1,ea1);

if round(norm(eb1-eb0),10) ~= 0 % Don't quite understand why the signs need to be switched, but it works
    phi(2) = -atan2(y,x);
    M1r = rodrigues_rot(M1ri,ea0,phi(2));
    % Calculate basis of moving mesh
    Pa1 = M1r(i_Pa1,:);
    Pb1 = M1r(i_Pb1,:);
    ea1 = Pa1/norm(Pa1);
    ec1 = cross(Pa1,Pb1);
    ec1 = ec1/norm(ec1);
    eb1 = cross(ec1,ea1);
end

figure('Position',[plot_s+3*plot_w plot_y1 plot_w plot_h],'NumberTitle','off')
scatter3(M0t(:,1), M0t(:,2), M0t(:,3),'MarkerEdgeColor',[0.5 0.5 0.5]); hold on
scatter3(M1r(:,1), M1r(:,2), M1r(:,3),'k.'); hold on
line(UnitScale*[0 ea1(1) ],UnitScale*[0 ea1(2) ],UnitScale*[0 ea1(3) ],'Color','r'); hold on
line(UnitScale*[0 eb1(1) ],UnitScale*[0 eb1(2) ],UnitScale*[0 eb1(3) ],'Color','g'); hold on
line(UnitScale*[0 ec1(1) ],UnitScale*[0 ec1(2) ],UnitScale*[0 ec1(3) ],'Color','b'); hold on
line(UnitScale*[0 ea0(1) ],UnitScale*[0 ea0(2) ],UnitScale*[0 ea0(3) ],'Color','r','LineStyle','--'); hold on
line(UnitScale*[0 eb0(1) ],UnitScale*[0 eb0(2) ],UnitScale*[0 eb0(3) ],'Color','g','LineStyle','--'); hold on
line(UnitScale*[0 ec0(1) ],UnitScale*[0 ec0(2) ],UnitScale*[0 ec0(3) ],'Color','b','LineStyle','--'); hold on
axis equal
legend('M_0^t','M_1^r^i^i','e^a_1','e^b_1','e^c_1','e^a_0','e^b_0','e^c_0');
xlabel('x');
ylabel('y');
zlabel('z');
title('Second Rotation')

%close all
%% Scale moving mesh

Pb1 = M1r(i_Pb1,:);
Pb0 = M0t(i_Pb0,:);

% scale = norm(Pb0)/norm(Pb1);
scale = 1;

M1s = M1r*scale;

figure('Position',[plot_s plot_y2 plot_w plot_h],'NumberTitle','off')
scatter3(M0t(:,1), M0t(:,2), M0t(:,3),'MarkerEdgeColor','r','MarkerFaceColor','r'); hold on
scatter3(M1s(:,1), M1s(:,2), M1s(:,3),'k.'); hold on
axis equal
%     legend('M_0^s','M_1^s')%,'e^a_1','e^b_1','e^c_1','e^a_0','e^b_0','e^c_0');
xlabel('x');
ylabel('y');
zlabel('z');
%     title('Scaled Meshes')
axis off
%     set(gcf,'renderer','painters')

%% Hausdorff Distance Calculation

M0f = M0t;
M1f = M1s;

Hausdorff = zeros(N0,1);
i_Hausdorff = zeros(N0,1);
for i = 1:N0
    d = zeros(N1,1);
    for j = 1:N1
        d(j) = norm(M0f(i,:)-M1f(j,:));
    end
    [Hausdorff(i), i_Hausdorff(i)] = min(d);
end

Hausdorff_mean = mean(Hausdorff);
Hausdorff_std = std(Hausdorff);
Hausdorff_median = median(Hausdorff);
Hausdorff_mode = mode(Hausdorff);

figure('Position',[plot_s+plot_w plot_y2 plot_w plot_h],'NumberTitle','off')
histogram(Hausdorff,ceil(sqrt(N0)));
title('Distribution of Hausdorff Distances')
xlabel('Hausdorff Distances [mm]')
ylabel('Count')

% close all

% figure('Position',[plot_s+2*plot_w plot_y2 plot_w plot_h],'NumberTitle','off')
figure('Position',[plot_s plot_y2 plot_w plot_h],'NumberTitle','off')
[f,xi] = ksdensity(Hausdorff);
plot(xi,f); hold on
A = 0;
B = max(xlim);
C = min(ylim);
D = max(ylim);
line([Hausdorff_mean Hausdorff_mean],[0 D],'Color','k','LineStyle','--'); hold on
line([Hausdorff_median Hausdorff_median],[0 D],'Color','r','LineStyle','-.'); hold on
line([Hausdorff_mean+Hausdorff_std Hausdorff_mean+Hausdorff_std],[0 D],'Color','b','LineStyle',':');
line([Hausdorff_mean-Hausdorff_std Hausdorff_mean-Hausdorff_std],[0 D],'Color','b','LineStyle',':');
xlabel('Hausdorff Distance [mm]')
axis([A B C D])
ylabel('PDF')
title(['Hausdorff Distribution Statistics (n = ',num2str(N0),')'])
legend(['PDF Estimation'],['Mean (',num2str(round(Hausdorff_mean,3)),')'],['Median (',num2str(round(Hausdorff_median,3)),')'],['StdDev (',num2str(round(Hausdorff_std,3)),')'],'Location','northeast')
clear A B C D
M1a = zeros(size(i_Hausdorff,1),3);
for i = 1:size(i_Hausdorff,1)
    M1a(i,:) = M1s(i_Hausdorff(i),:);
end; clear i



%% Plot minimal error meshes

% figure('Position',[plot_s+3*plot_w plot_y2 plot_w plot_h],'NumberTitle','off')
% figure('Position',[plot_s+plot_w plot_y2 plot_w plot_h],'NumberTitle','off')
View = [-41 38];
View = [0 90];
Color = 'viridis';
%     ColorHausdorff = 'No';

if strcmp(ColorHausdorff,'Yes')
    figure()
    scatter3(M0t(:,1), M0t(:,2), M0t(:,3),[],Hausdorff/t,'filled'); hold on
    colorbar
    oldmap = colormap(Color);
    colormap(oldmap);
    caxis([0 45.4])
    %             colormap( flipud(oldmap))
    axis equal
    set(gca,'Color',[0.8 0.8 0.8])
    if strcmp(ShowAxes,'Yes')
        xlabel('x');
        ylabel('y');
        zlabel('z');
    else
        axis off
    end
    view(View)
    %         figure()
    %             scatter3(M1a(:,1), M1a(:,2), M1a(:,3),[],Hausdorff,'.'); hold on
    %             colorbar
    %             oldmap = colormap(Color);
    %             colormap(oldmap);
    % %             colormap( flipud(oldmap))
    %             axis equal
    %             if strcmp(ShowAxes,'Yes')
    %                 xlabel('x');
    %                 ylabel('y');
    %                 zlabel('z');
    %             else
    %                 axis off
    %             end
    %             view(View)
else
    x0 = M0t(:,1);
    y0 = M0t(:,2);
    z0 = M0t(:,3);
    [xq0,yq0] = meshgrid(min(x0):1:max(x0), min(y0):1:max(y0));
    zq0 = griddata(x0,y0,z0,xq0,yq0,'linear');
    
    x1 = M1a(:,1);
    y1 = M1a(:,2);
    z1 = M1a(:,3);
    [xq1,yq1] = meshgrid(min(x1):1:max(x1), min(y1):1:max(y1));
    zq1 = griddata(x1,y1,z1,xq1,yq1,'linear');
    
    figure()
    scatter3(M0t(:,1), M0t(:,2), M0t(:,3),'r.'); hold on
    scatter3(M1a(:,1), M1a(:,2), M1a(:,3),'b.'); hold on
    mesh(xq0,yq0,zq0,'EdgeColor','none','FaceColor','m','FaceAlpha',0.75); hold on
    mesh(xq1,yq1,zq1,'EdgeColor','none','FaceColor','c','FaceAlpha',0.75);
    legend('Bar and Hinge','Laser Scan')
    view(View)
    axis equal
    if strcmp(ShowAxes,'Yes')
        xlabel('x');
        ylabel('y');
        zlabel('z');
    else
        axis off
    end
    view(View)
end

%% Data Point Function

function output_txt = ShowIndex(obj,event_obj)

% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');

% Import x and y
x = get(get(event_obj,'Target'),'XData');
y = get(get(event_obj,'Target'),'YData');
z = get(get(event_obj,'Target'),'ZData');

% Find index
index_x = find(x == pos(1));
index_y = find(y == pos(2));
index_z = find(z == pos(3));
index = intersect(intersect(index_x,index_y),index_z);

% Set output text
if length(pos) > 2
    output_txt = {['X: ',num2str(pos(1),4)], ...
        ['Y: ',num2str(pos(2),4)], ...
        ['Z: ',num2str(pos(3),4)], ...
        ['Index: ', num2str(index)]};
else
    output_txt = {['X: ',num2str(pos(1),4)], ...
        ['Y: ',num2str(pos(2),4)], ...
        ['Index: ', num2str(index)]};
end

end