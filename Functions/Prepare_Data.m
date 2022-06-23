function [PreprocessData] = Prepare_Data(InputData)

%% Define element types
Nn = max(cellfun(@max,InputData.panels)); 
Ptri = InputData.panels';
Triglraw = cell2mat(Ptri);
Triglraw = sort(Triglraw,2);
Trigl = unique(Triglraw ,'rows');

%% Formulate connectivity matrix
Comm = sparse(Nn,size(Trigl,1));
for i=1:size(Trigl,1) 
    Comm(Trigl(i,:),i) = true;
end

%% Search for fold lines
Ge = Comm'*Comm;
[mf, me] = find(triu(Ge==2)); % triangular meshes that share two common nodes
Fold = [];
Bend = [];
for i=1:length(mf)
    [link,ia,ib] = intersect(Trigl(mf(i),:),Trigl(me(i),:));
    oftpa = setdiff(1:3,ia);
    oftpb = setdiff(1:3,ib);
    
    Nd_Fld_arr=[[InputData.foldNodes';0],[0;InputData.foldNodes']];
    if max(ismember(Nd_Fld_arr,link([1,2]),'rows'))
        Fold = [Fold; [link,Trigl(mf(i),oftpa),Trigl(me(i),oftpb)]];
    elseif max(ismember(Nd_Fld_arr,link([2,1]),'rows'))
        Fold = [Fold; [link,Trigl(mf(i),oftpa),Trigl(me(i),oftpb)]];
    else
        Bend = [Bend; [link,Trigl(mf(i),oftpa),Trigl(me(i),oftpb)]];
    end
end
fdandbd = sort(Fold(:,1:2),2);
onlybd = sort(Bend(:,1:2),2);
[~,ibd] = intersect(fdandbd,onlybd,'rows');
Fold(ibd,:) = [];

%% Search for boundaries
Edge = sort([Trigl(:,1) Trigl(:,2); Trigl(:,2) Trigl(:,3); Trigl(:,3) Trigl(:,1)],2);
[u,~,n] = unique(Edge ,'rows');
counts = accumarray(n(:), 1);
Bdry = u(counts==1,:);

%% Find folding hinges and boundaries, return final triangulation
Bars = [Bend(:,1:2);Fold(:,1:2);Bdry];

%% Direction calculations
Ne = size(Bars,1);
Nn = size(InputData.nodes,1);
D = [InputData.nodes(Bars(:,2),1)-InputData.nodes(Bars(:,1),1), ...
     InputData.nodes(Bars(:,2),2)-InputData.nodes(Bars(:,1),2),...
     InputData.nodes(Bars(:,2),3)-InputData.nodes(Bars(:,1),3)];
Li = sqrt(D(:,1).^2+D(:,2).^2+D(:,3).^2);
D = [D(:,1)./Li D(:,2)./Li D(:,3)./Li];
B = sparse(repmat((1:Ne)',1,6),[3*Bars(:,1)-2 3*Bars(:,1)-1 3*Bars(:,1),...
           3*Bars(:,2)-2 3*Bars(:,2)-1 3*Bars(:,2)],[D -D],Ne,3*Nn);
B = -B;

if size(InputData.supports,1) == 0
    rs = []; 
else
    rs = [reshape([InputData.supports(:,1)*3-2,InputData.supports(:,1)*3-1,InputData.supports(:,1)*3]',[],1),...
          reshape(InputData.supports(:,2:4)',[],1)];
    rs(rs(:,2)==0,:)=[]; rs = rs(:,1);
end 

%% Calculate bar areas
% Shear coefficients
    a0 = 1.25; a1 = 0; a2 = 0.5;
       
% Effective Bar Areas
    barAreas = zeros(size(Bars,1),1);
    for i = 1:length(InputData.panels) % Loop through panels
        % Find nodes on the ith panel
        for j = 1:3
            Index(j) = InputData.panels{i}(j);
            Point{j} = InputData.nodes(Index(j),:);
        end; clear j
        % Calculate the area of the ith panel
        panelArea = 1/2*(norm(cross((Point{2}-Point{1}),(Point{3}-Point{1}))));
        % Calculate the length of the bars surrounding the ith panel
        for j = 1:3
            L(j) = norm(Point{j} - Point{j+1 - floor(j/3)*3});
        end; clear j
        Lengths = sort([L(1); L(2); L(3)]);
        AspectRatio(i) = Lengths(2)/Lengths(1); 
        % Calculate bar cross-sectional areas in the ith panel
        for j = 1:3
            if L(j) == Lengths(3)
                A(j) = (a2*AspectRatio(i)^2 + a1*AspectRatio(i) + a0)*panelArea*InputData.thickness/(L(j));
            else
                A(j) = panelArea*InputData.thickness/(L(j));
            end
        end; clear j
        % Assign panel bars to final bar area indexing
        for j = 1:3
            endCondition = 0;
            k = 0;
            while endCondition == 0
                k = k + 1;
                find1 = ~isempty(find(Bars(k,:)==Index(j),1));
                find2 = ~isempty(find(Bars(k,:)==Index(j+1 - floor(j/3)*3),1));
                finder = find1 + find2;
                if finder == 2
                    Bar_No(j) = k;
                    endCondition = 1;
                end
            end
        end; clear j
        Panel_Bar_Areas = zeros(size(Bars,1),1);
        for j = 1:3
            Panel_Bar_Areas(Bar_No(j),1) = A(j);
        end
        barAreas = barAreas + Panel_Bar_Areas;
    end
    PreprocessData.barAreas = barAreas;
%% Calculate initial fold angles from mesh

pf0 = zeros(size(Fold,1),1); 
for i = 1:size(Fold,1)
    pf0(i) = FoldKe(InputData.nodes,Fold(i,:)); 
end

%% Calculate initial bending angles from mesh
pb0 = zeros(size(Bend,1),1); 
for i = 1:size(Bend,1)
    pb0(i) = FoldKe(InputData.nodes,Bend(i,:)); 
end

%% Generate bending stiffness array 
bendingLengths = zeros(size(Bend,1),1); % Initialize bending hinge length array
PreprocessData.Kb = bendingLengths; % Initialize bending stiffness array
for ang=1:size(Bend,1)
    node1=InputData.nodes(Bend(ang,1),:);   
    node2=InputData.nodes(Bend(ang,2),:);
    node3=InputData.nodes(Bend(ang,3),:);   
    node4=InputData.nodes(Bend(ang,4),:);
    bendingLengths(ang)=norm(node1-node2);
    A1=1/2*norm(cross(node2-node1,node3-node1));
    A2=1/2*norm(cross(node2-node1,node4-node1));
    PreprocessData.Kb(ang)=1/6*bendingLengths(ang)^2*InputData.elasticModulus*InputData.thickness^3/(A1+A2);
end

%% Generate folding stiffness array
PreprocessData.mountainFolds = [];
PreprocessData.valleyFolds = [];
ValNd_mod=[[InputData.valleyNodes; 0],[0; InputData.valleyNodes]];

for ang=1:size(Fold,1)
    % Separate folds into mountain and valley folds
    if max(ismember(ValNd_mod,Fold(ang,[1,2]),'rows'))
            PreprocessData.valleyFolds=[PreprocessData.valleyFolds;ang];
    elseif max(ismember(ValNd_mod,Fold(ang,[2,1]),'rows'))
            PreprocessData.valleyFolds=[PreprocessData.valleyFolds;ang];
    else
        PreprocessData.mountainFolds=[PreprocessData.mountainFolds;ang];
    end
    % Calculate fold lengths for fold stiffness expression
    node1=InputData.nodes(Fold(ang,1),:);
    node2=InputData.nodes(Fold(ang,2),:);
    foldLengths(ang)=norm(node1-node2);
end; clear ang

% Fold stiffness expression
PreprocessData.Kf = [(foldLengths(1:InputData.numberDivisions)')./InputData.lengthScale1*InputData.elasticModulus*InputData.thickness^3/12;
                     (foldLengths(InputData.numberDivisions+1:2*InputData.numberDivisions)')./InputData.lengthScale2*InputData.elasticModulus*InputData.thickness^3/12;
                     (foldLengths(2*InputData.numberDivisions+1:3*InputData.numberDivisions)')./InputData.lengthScale1*InputData.elasticModulus*InputData.thickness^3/12];

%% Spread load onto DOF matrix (3*n by 1)
if ~isempty(InputData.loads)
    m = size(InputData.nodes,1);
    FD = zeros(3*m,1);
    indp = InputData.loads(:,1);
    FD(3*indp-2) = InputData.loads(:,2); % X-load
    FD(3*indp-1) = InputData.loads(:,3); % Y-load
    FD(3*indp) = InputData.loads(:,4);   % Z-load
    PreprocessData.load = FD;  % Puts the load back into the input options?
end

PreprocessData.trussBars = Bars; % Bar nodes
PreprocessData.Trigl = Trigl; % What's the difference between Trigl and Panels?
PreprocessData.B = B; %
PreprocessData.L = Li;
PreprocessData.fixedDOFs = unique(rs);
PreprocessData.foldingHinges = Fold;
PreprocessData.bendingHinges = Bend;
PreprocessData.initialFoldAngles = pf0;
PreprocessData.initialBendAngles = pb0;
PreprocessData.aspectRatio = AspectRatio;
% PreprocessData.onlyBendingHinges = onlybd;

for i = 1:size(PreprocessData.foldingHinges,1)
        n1 = PreprocessData.foldingHinges(i,1);
        n2 = PreprocessData.foldingHinges(i,2);
        if n2 - n1 ~= 1
            PreprocessData.foldingHinges(i,1) = n2;
            PreprocessData.foldingHinges(i,2) = n1;
        end
end
%PreprocessData.foldingHinges = sortrows(PreprocessData.foldingHinges,1);
end