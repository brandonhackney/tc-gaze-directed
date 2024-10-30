function searchlights = getSearchlight(vertexList, xyz, radius)
% Given a list of vertices, their coordinates, and a radius,
% define the searchlights around each vertex in terms of vertex indices
%
% Assume vertexlist is a vector of vertex indices
% Assume xyz is n rows long and 3 cols wide: x, y, z
% Assume radius is in coordinate units (i.e. not a number of vertices)

numVertices = length(vertexList);
searchlights = cell(1,numVertices); % empty cell
for i = 1:numVertices
    iCoords = xyz(i,:);
    % Define searchlight as x +- rad & y +- rad & z +- rad
        xmin = iCoords(1) - radius;
        xmax = iCoords(1) + radius;
        x = xyz(:, 1) >= xmin & xyz(:,1) <= xmax;
        ymin = iCoords(2) - radius;
        ymax = iCoords(2) + radius;
        y = xyz(:,2) >= ymin & xyz(:,2) <= ymax;
        zmin = iCoords(3) - radius;
        zmax = iCoords(3) + radius;
        z = xyz(:,3) >= zmin & xyz(:,3) <= zmax;
    searchlights{i} = vertexList(x & y & z);
    % searchlights{i} is a vector with the indices of matching verts
end

end