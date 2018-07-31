function [] = computeGlobalStress(area, cellInfo, edgeInfo)
%COMPUTEGLOBALSTRESS Summary of this function goes here
%   Detailed explanation goes here
%   From equations Sugimura - 2013 
%   "The mechanical anisotropy in a tissue promotes ordering in hexagonal 
%   cell packing", Kaoru Sugimura and Shuji Ishihara

%% Parameters notation
% N_uv      -> Global stress tensor
% Delta_uv  -> Kronecker's delta
% pressure (P) -> Pressure of the cell
% tension (T)         -> Tension of the edge
% area (A)  -> Area of the cell
% edgeLength (l)         -> Edge length
% lx        -> Component X of edge length l
% ly        -> Component Y of edge length l


N_uv = cell(2, 2);
% if u is 1, represents 'x' coordinate, otherwise -> 'y'
for u = 1:2
    % if v is 1, represents 'x' coordinate, otherwise -> 'y'
    for v = 1:2
        localStress = zeros(size(cellInfo, 1), 1);
        %Kronecker's delta
        %The function is 1 if the variables are equal, and 0 otherwise
        %So, if u = v, 1; u ~= v, 0.
        if u == v
            Delta_uv = 1;
        else
            Delta_uv = 0;
        end
        for numCell = 1:size(cellInfo, 1)
            actualEdgeInfo = edgeInfo(edgeInfo.ConnectingCell1 == numCell | edgeInfo.ConnectingCell2 == numCell, :);

            if isempty(actualEdgeInfo) == 0
                edgesStress = zeros(size(actualEdgeInfo, 1), 1);
                for numEdge = 1:size(actualEdgeInfo, 1)
                    %Checking if we have to inverse the vector
                    if actualEdgeInfo.ConnectingCell1(numEdge) ~= numCell
                        inverse = 1;
                    else
                        inverse = 0;
                    end
                    tension = actualEdgeInfo.TensionValue(numEdge);
                    if u == 1 && ~inverse
                        l_u = actualEdgeInfo.vertex2_X(numEdge) - actualEdgeInfo.vertex1_X(numEdge);
                    else
                        l_u = actualEdgeInfo.vertex1_X(numEdge) - actualEdgeInfo.vertex2_X(numEdge);
                    end
                    
                    if v == 1 && ~inverse
                        l_v = actualEdgeInfo.vertex2_Y(numEdge) - actualEdgeInfo.vertex1_Y(numEdge);
                    else
                        l_v = actualEdgeInfo.vertex1_Y(numEdge) - actualEdgeInfo.vertex2_Y(numEdge);
                    end
                    
                    edgeLength = actualEdgeInfo.EdgeLength(numEdge);
                    edgesStress(numEdge) = tension * ( l_u * l_v ) / edgeLength;
                end
                cellStress = (-area(numCell).Area * cellInfo(numCell, 2) * Delta_uv );
                %Global stress
                localStress(numCell) =  (cellStress + sum(edgesStress)) / area(numCell).Area;
            end
        end
        N_uv(u, v) = {localStress};
    end
end

disp('')
%Calculate the eigenvectors of the global stress tensor



end

