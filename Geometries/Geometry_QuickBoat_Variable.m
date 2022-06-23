%% Hardik Patil - Nov 2021 - Full Planing Hull %%

function [InputData] = Geometry_QuickBoat_Variable(InputData)
number_divisions = InputData.numberDivisions;
Load_Magnitude = InputData.loadMagnitude;
Plot_Nodes = InputData.plotNodes;
h_tip = InputData.tipHeight;              % Default is 71.887
L_top = InputData.lengthStraight;         % length of the straight part on top
L_bot = InputData.lengthStraight;         % length of the straight part on bottom
L_curve = InputData.lengthCurve;

L = L_bot + L_curve;             % total width of sheet (default = 954)
W = InputData.panelWidth;        % panel width (default = 162)

syms u
y_bot = (h_tip / (L_curve)^2) * (u - L_bot)^2;           % equation of the bot curve
y_top = ((h_tip-W) / (L_curve)^2) * (u - L_top)^2 + W;      % equation of the top curve

%% Calculate the mesh nodes
No_Nodes = 5*number_divisions+1;
No_Nodes_Row = number_divisions+1;
Nodes = zeros(No_Nodes,3);
x_i = 7.76;

for i = 1:No_Nodes
    if i < No_Nodes_Row % upto 40
        if i == 1
            x_i = 12.07;
        end
        Nodes(i,:) = [x_i h_tip 0];
        x_i = x_i + (L-12.07)/number_divisions;
    
    elseif i >= No_Nodes_Row && i < 2*No_Nodes_Row-1 % 41 to 80
        if i == No_Nodes_Row
            x_i = 7.76;
        end
        if x_i < L_top
            Nodes(i,:) = [x_i, W, 0];
            x_i = x_i + L/number_divisions;
        else
            Nodes(i,:) = [x_i, subs(y_top,u,x_i), 0];
            x_i = x_i + L/number_divisions;
        end
        
    elseif i >= 2*No_Nodes_Row-1 && i < 3*No_Nodes_Row-1 % 81 to 121
        if i == 2*No_Nodes_Row-1
            x_i = 0;
        end
        if x_i < L_bot
            Nodes(i,:) = [x_i, 0, 0];
            x_i = x_i + L/number_divisions;
        else
            Nodes(i,:) = [x_i, subs(y_bot,u,x_i), 0];
            x_i = x_i + L/number_divisions;
        end
    
    elseif i >= 3*No_Nodes_Row-1 && i < 4*No_Nodes_Row-2 % 122 to 161
        if i == 3*No_Nodes_Row-1
            x_i = 7.76;
        end
        
        if x_i < L_top
            Nodes(i,:) = [x_i, W, 0];
            x_i = x_i + L/number_divisions;
        else
            Nodes(i,:) = [x_i, subs(y_top,u,x_i), 0];
            x_i = x_i + L/number_divisions;
        end

    elseif i >= 4*No_Nodes_Row-2 % 162 onwards
        if i == 4*No_Nodes_Row-2
            x_i = 12.07;
        end
        Nodes(i,:) = [x_i h_tip 0];
        x_i = x_i + (L-12.07)/number_divisions;

    end
end

Nodes_Fold = [(No_Nodes_Row:2*No_Nodes_Row-2) 3*No_Nodes_Row-2 (2*No_Nodes_Row-1:3*No_Nodes_Row-2) (3*No_Nodes_Row-1:4*No_Nodes_Row-3) 3*No_Nodes_Row-2];
Nodes_Valley = [(No_Nodes_Row:2*No_Nodes_Row-2) 3*No_Nodes_Row-2 (2*No_Nodes_Row-1:3*No_Nodes_Row-2) (3*No_Nodes_Row-1:4*No_Nodes_Row-3) 3*No_Nodes_Row-2]';

%% Boundary conditions
edge = zeros(number_divisions,4);

for i = 1:size(edge,1)
    edge(i,:) = [2*No_Nodes_Row-1+i 0 0 1];
end

Supports = [2*No_Nodes_Row-1 1 1 1;
            edge];

Loads = [No_Nodes_Row+1+(number_divisions/2) 0 0 0;
         3*No_Nodes_Row+(number_divisions/2) 0 0 0];
     
%% Calculate the nodal connections (panels)
No_Panels = (number_divisions*2-1)*4;
N1 = zeros(No_Panels,1);
N2 = N1; N3 = N1; N4 = N1;
sign = -1;
panel_count = 0;

for i = 1:(No_Panels+6)/2
    N1(i) = i;
    N2(i) = N1(i) + number_divisions;
    N3(i) = N1(i) + number_divisions + 1;
    N4(i) = N1(i) + 1;
    
    if i < No_Nodes_Row-1 % upto 39
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
        
    elseif i == No_Nodes_Row-1 % 40
        panel_count = panel_count+1;
        Panels{panel_count} = [N1(i) N2(i) N3(i)+number_divisions];
        
    elseif i >= No_Nodes_Row && i < 2*No_Nodes_Row-2 % 41 to 79
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

    elseif i == 2*No_Nodes_Row-2 %80
        
        panel_count = panel_count + 1;
        Panels{panel_count} = [N1(i) N2(i) N3(i)];
        
    elseif i >= 2*No_Nodes_Row-1 && i < 3*No_Nodes_Row-3 %81 to 119
        if sign == 1
            panel_count = panel_count+1;
            Panels{panel_count} = [N1(i) N2(i)+1 N3(i)+1];
            panel_count = panel_count+1;
            Panels{panel_count} = [N1(i) N3(i)+1 N4(i)];
            sign = -sign;
        elseif sign == -1
            panel_count = panel_count+1;
            Panels{panel_count} = [N1(i) N2(i)+1 N4(i)];
            panel_count = panel_count+1;
            Panels{panel_count} = [N2(i)+1 N3(i)+1 N4(i)];
            sign = -sign;
        end
    
    elseif i == 3*No_Nodes_Row-3 %120
        panel_count = panel_count + 1;
        Panels{panel_count} = [N1(i) N2(i)+1 N4(i)];
    
    elseif i >= 3*No_Nodes_Row-1 && i < 4*No_Nodes_Row-3 %122 to 160

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
        
    elseif i == 4*No_Nodes_Row-3 %161
        
        panel_count = panel_count+1;
        Panels{panel_count} = [N1(i) N2(i) N1(i)-number_divisions];
        
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
    axis off
end
end