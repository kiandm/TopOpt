function plotDensityMesh(nodes, enodes, x)
    % simple patch coloring by x
    clf;
    hold on;
    colormap(flipud(gray)); % black=1 (solid), white=0 (void) if we invert later
    for e = 1:size(enodes,1)
        n  = enodes(e,:);
        xy = [nodes(n,2), nodes(n,3)];
        c = 1 - max(0,min(1,x(e)));
        c = [c c c];
        patch('Faces', [1 2 3 4], ...
              'Vertices', xy, ...
              'FaceColor', c, ...
              'EdgeColor', [0.7 0.7 0.7], ...
              'LineWidth', 0.25);
    end
    axis equal tight off;
    title(sprintf('Density (vol=%.3f)', mean(x)));
    hold off;
end
