%% Hardik Patil - Mar 2021 - Full Planing Hull %%

function [InputData] = Geometry_Planing(InputData)
number_divisions = InputData.numberDivisions;
Load_Magnitude = InputData.loadMagnitude;
Plot_Nodes = InputData.plotNodes;
syms u
y1 = -8.96110745294953e-16*u^6 + 1.71802291404577e-12*u^5 - 7.15943472458927e-10*u^4 - 1.07657203662384e-07*u^3 + 9.57675494987319e-05*u^2 - 0.0115561298669055*u + 0.193869195797025;
y2 = -1.00429866111023e-15*u^6 + 3.26050178033549e-12*u^5 - 4.09228922161364e-09*u^4 + 2.15858000996365e-06*u^3 - 0.000471392453018989*u^2 + 0.0394612376507102*u + 159.754291457347;
W = InputData.lengthCurve + InputData.lengthStraight;
a = 5;

%% Calculate the mesh nodes
No_Nodes = 5*number_divisions+3;
No_Nodes_Row = number_divisions+1;
Nodes = zeros(No_Nodes,3);
x_i = 7.76;

for i = 1:No_Nodes
    
    if i <= No_Nodes_Row
        
        Nodes(i,:) = [x_i subs(y2,u,x_i)+a 0];
        x_i = x_i + W/number_divisions;
        
    elseif i >= No_Nodes_Row+1 && i <= 2*No_Nodes_Row-1
        
        if i == No_Nodes_Row+1
            x_i = 7.76;
        end
        
        Nodes(i,:) = [x_i subs(y2,u,x_i) 0];
        x_i = x_i + W/number_divisions;
        
    elseif i >= 2*No_Nodes_Row && i <= 3*No_Nodes_Row-1

        if i == 2*No_Nodes_Row
            x_i = 0;
        end
        
        Nodes(i,:) = [x_i subs(y1,u,x_i) 0];
        x_i = x_i + W/number_divisions;
    
    elseif i >= 3*No_Nodes_Row && i <= 4*No_Nodes_Row-2

        if i == 3*No_Nodes_Row
            x_i = 7.76;
        end
        
        Nodes(i,:) = [x_i subs(y2,u,x_i) 0];
        x_i = x_i + W/number_divisions;

    elseif i >= 4*No_Nodes_Row-1

        if i == 4*No_Nodes_Row-1
            x_i = 7.76;
        end
        
        Nodes(i,:) = [x_i subs(y2,u,x_i)+a 0];
        x_i = x_i + W/number_divisions;

    end
end

Nodes_Fold = [(No_Nodes_Row+1:1:2*No_Nodes_Row-1) 3*No_Nodes_Row-1 (2*No_Nodes_Row:1:3*No_Nodes_Row-1) (3*No_Nodes_Row:1:4*No_Nodes_Row-2) 3*No_Nodes_Row-1];
Nodes_Valley = (2*No_Nodes_Row:1:3*No_Nodes_Row-1)';

%% Boundary conditions
edge = [];

for i = 2*No_Nodes_Row+1:1:3*No_Nodes_Row-1
    edge = [edge;i 0 0 1];
end

Supports = [edge;
            2*No_Nodes_Row 1 1 1;
            2*No_Nodes_Row+10 1 1 1];

Loads = [No_Nodes_Row+1+(number_divisions/2) 0 0 0;
         3*No_Nodes_Row+(number_divisions/2) 0 0 0];

% Loads = [No_Nodes_Row+1 0 0 0;
%          No_Nodes_Row+1+(number_divisions/2) 0 0 0;
%          3*No_Nodes_Row 0 0 0;
%          3*No_Nodes_Row+(number_divisions/2) 0 0 0];
     
%% Calculate the nodal connections (panels)

No_Panels = (number_divisions*2)*4-2;
N1 = zeros(No_Panels,1);
N2 = N1; N3 = N1; N4 = N1;
sign = -1;
panel_count = 0;

for i = 1:(No_Panels+6)/2
    N1(i) = i;
    N2(i) = N1(i) + No_Nodes_Row;
    N3(i) = N1(i) + No_Nodes_Row + 1;
    N4(i) = N1(i) + 1;
    
    if i < No_Nodes_Row-1

        if sign == 1
            panel_count = panel_count+1;
            Panels{panel_count} = [N1(i) N2(i) N3(i)];
            panel_count = panel_count+1;
            Panels{panel_count} = [N1(i) N3(i) N4(i)];
            sign = -sign;
        elseif sign == -1
            panel_count = panel_count+1;
            Panels{panel_count} = [N1(i) N2(i) N4(i)];
            panel_count = panel_count+1;
            Panels{panel_count} = [N2(i) N3(i) N4(i)];
            sign = -sign;
        end
        
    elseif i == No_Nodes_Row-1
        
        panel_count = panel_count+1;
        Panels{panel_count} = [N1(i) N2(i) N3(i)+number_divisions];
        panel_count = panel_count+1;
        Panels{panel_count} = [N1(i) N3(i)+number_divisions N4(i)];
        
    elseif i > No_Nodes_Row && i < 2*No_Nodes_Row-1
        
        if sign == 1
            panel_count = panel_count+1;
            Panels{panel_count} = [N1(i) N2(i)-1 N3(i)-1];
            panel_count = panel_count+1;
            Panels{panel_count} = [N1(i) N3(i)-1 N4(i)];
            sign = -sign;
        elseif sign == -1
            panel_count = panel_count+1;
            Panels{panel_count} = [N1(i) N2(i)-1 N4(i)];
            panel_count = panel_count+1;
            Panels{panel_count} = [N2(i)-1 N3(i)-1 N4(i)];
            sign = -sign;
        end

    elseif i == 2*No_Nodes_Row-1
        
        panel_count = panel_count + 1;
        Panels{panel_count} = [N1(i) N2(i)-1 N2(i)];
        
    elseif i >= 2*No_Nodes_Row && i < 3*No_Nodes_Row-2
        
        if sign == 1
            panel_count = panel_count+1;
            Panels{panel_count} = [N1(i) N2(i) N3(i)];
            panel_count = panel_count+1;
            Panels{panel_count} = [N1(i) N3(i) N4(i)];
            sign = -sign;
        elseif sign == -1
            panel_count = panel_count+1;
            Panels{panel_count} = [N1(i) N2(i) N4(i)];
            panel_count = panel_count+1;
            Panels{panel_count} = [N2(i) N3(i) N4(i)];
            sign = -sign;
        end
    
    elseif i == 3*No_Nodes_Row-2
    
        panel_count = panel_count + 1;
        Panels{panel_count} = [N1(i) N2(i) N4(i)];
    
    elseif i >= 3*No_Nodes_Row && i < 4*No_Nodes_Row-2

        if sign == 1
            panel_count = panel_count+1;
            Panels{panel_count} = [N1(i) N2(i)-1 N3(i)-1];
            panel_count = panel_count+1;
            Panels{panel_count} = [N1(i) N3(i)-1 N4(i)];
            sign = -sign;
        elseif sign == -1
            panel_count = panel_count+1;
            Panels{panel_count} = [N1(i) N2(i)-1 N4(i)];
            panel_count = panel_count+1;
            Panels{panel_count} = [N2(i)-1 N3(i)-1 N4(i)];
            sign = -sign;
        end
        
    elseif i == 4*No_Nodes_Row-2
        
        panel_count = panel_count+1;
        Panels{panel_count} = [N1(i) N2(i)-1 N1(i)-number_divisions];
        panel_count = panel_count+1;
        Panels{panel_count} = [N2(i)-1 N3(i)-1 N1(i)-number_divisions];
        
    end
end; clear i Add
clear N1 N2 N3 N4


%% Plot the nodes and mesh
InputData.nodes = Nodes;
InputData.panels = Panels;
InputData.supports = Supports;
InputData.loads = Loads;
InputData.foldNodes = Nodes_Fold;
InputData.valleyNodes = Nodes_Valley;

if strcmpi(InputData.plotNodes,'yes')
    figure()
    for i = 1:No_Panels
        line([InputData.nodes(InputData.panels{i}(1),1), InputData.nodes(InputData.panels{i}(2),1)], [InputData.nodes(InputData.panels{i}(1),2), InputData.nodes(InputData.panels{i}(2),2)],'Color',[140 140 140]./255); hold on
        line([InputData.nodes(InputData.panels{i}(2),1), InputData.nodes(InputData.panels{i}(3),1)], [InputData.nodes(InputData.panels{i}(2),2), InputData.nodes(InputData.panels{i}(3),2)],'Color',[140 140 140]./255); hold on
        line([InputData.nodes(InputData.panels{i}(3),1), InputData.nodes(InputData.panels{i}(1),1)], [InputData.nodes(InputData.panels{i}(3),2), InputData.nodes(InputData.panels{i}(1),2)],'Color',[140 140 140]./255); hold on
    end; clear i
    scatter(InputData.nodes(:,1),InputData.nodes(:,2),'k.'); hold on
    for i = 1:size(InputData.supports,1)
        if InputData.supports(i,2) == 1
            scatter(InputData.nodes(InputData.supports(i,1),1),InputData.nodes(InputData.supports(i,1),2),80,'^','filled','MarkerFaceColor',[65 69 136]./255); hold on
        end
        if InputData.supports(i,3) == 1
            scatter(InputData.nodes(InputData.supports(i,1),1),InputData.nodes(InputData.supports(i,1),2),60,'>','filled','MarkerFaceColor',[43 120 142]./255); hold on
        end
        if InputData.supports(i,4) == 1
            scatter(InputData.nodes(InputData.supports(i,1),1),InputData.nodes(InputData.supports(i,1),2),50,'v','filled','MarkerFaceColor',[31 168 132]./255); hold on
        end
    end
    for i = 1:length(InputData.foldNodes)-1
        if ~(mod(i,InputData.numberDivisions+1) == 0)
            if ~isempty(find(InputData.foldNodes(i) == InputData.valleyNodes,1))
                line([InputData.nodes(InputData.foldNodes(i),1), ...
                    InputData.nodes(InputData.foldNodes(i+1),1)],...
                    [InputData.nodes(InputData.foldNodes(i),2),...
                    InputData.nodes(InputData.foldNodes(i+1),2)],...
                    'Color',[67 24 83]./255,'LineWidth',2); hold on
            else
                line([InputData.nodes(InputData.foldNodes(i),1),...
                    InputData.nodes(InputData.foldNodes(i+1),1)],...
                    [InputData.nodes(InputData.foldNodes(i),2),...
                    InputData.nodes(InputData.foldNodes(i+1),2)],...
                    'Color',[31 168 132]./255,'LineWidth',2); hold on
            end
        end
    end
    scatter(InputData.nodes(InputData.loads(:,1),1),InputData.nodes(InputData.loads(:,1),2),'o','filled','MarkerFaceColor',[129 196 85]./255); hold on
    axis equal
%     axis off
end
end