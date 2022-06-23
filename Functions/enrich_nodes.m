function enriched_nodes = enrich_nodes(enrich_factor, InputData, PostprocessData)
%{
Takes in bar and hinge nodes and adds additional nodes along the length
of the non-diagonal bars. The enriched nodes can be used for mesh
comparison using the Hausdorff distance calculator.
%}

num_nodes_per_row = InputData.numberDivisions + 1;
enriched_nodes = PostprocessData.deformedNodes{end};
num_spaces = size(InputData.nodes, 1)/(num_nodes_per_row) - 1;

for i = 1:num_nodes_per_row*num_spaces
    r{i} = PostprocessData.deformedNodes{end}(i + num_nodes_per_row,:) ...
        - PostprocessData.deformedNodes{end}(i,:);
    for j = 1:enrich_factor
        p{j} = PostprocessData.deformedNodes{end}(i,:) ...
            + r{i}/enrich_factor * j;
        enriched_nodes = [enriched_nodes; p{j}];
    end
end

enriched_nodes = unique(enriched_nodes, 'rows');

figure()
scatter3(enriched_nodes(:, 1), ...
    enriched_nodes(:, 2), ...
    enriched_nodes(:, 3), ...
    '.k')
axis equal
end