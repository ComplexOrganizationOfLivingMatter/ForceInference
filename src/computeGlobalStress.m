function [] = computeGlobalStress(areaOfCells, cellInfo, edgeInfo, validCells)
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


    numDimensions = 1;
    
    local_N_uv = cell(2, 2);
    global_N_uv = zeros(2, 2);
    % if u is 1, represents 'x' coordinate, otherwise -> 'y'
    for u = 1:2
        % if v is 1, represents 'x' coordinate, otherwise -> 'y'
        for v = 1:2
            localStress = zeros(size(cellInfo, 1), numDimensions);
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

                if isempty(actualEdgeInfo) == 0 && validCells(numCell)
                    edgesStress = zeros(size(actualEdgeInfo, 1), numDimensions);
                    for numEdge = 1:size(actualEdgeInfo, 1)
                        %Checking if we have to inverse the vector
                        if actualEdgeInfo.ConnectingCell1(numEdge) ~= numCell
                            inverse = 1;
                        else
                            inverse = 0;
                        end

                        edgesStress(numEdge, :) = calculateEdgeStress(actualEdgeInfo, numEdge, inverse, u, v);
                    end
                    cellStress(numCell) = -(areaOfCells(numCell) * cellInfo(numCell, 2) * Delta_uv );
                    %Global stress
                    localStress(numCell, :) =  (cellStress(numCell) + sum(edgesStress)) / areaOfCells(numCell);
                end
            end

            for numEdge = 1:size(edgeInfo, 1)
                globalEdgeStress(numEdge, :) = calculateEdgeStress(edgeInfo, numEdge, 0, u, v);
            end

    %         validCells = logical(zeros(size(areaOfCells)));
    %         validCells(cellStress~=0) = 1;
            globalStress = (sum(cellStress) + sum(globalEdgeStress))/sum([areaOfCells]);

            local_N_uv(u, v) = {localStress(validCells)};
            global_N_uv(u, v) = globalStress;
        end
    end


    %Calculate the eigenvectors of the global stress tensor
    eigenvectors = cell(size(localStress(validCells)));
    eigenvalues = cell(size(localStress(validCells)));
    for numCell = 1:size(localStress(validCells), 1)
        actualLocalStress = zeros(2, 2);
        for u = 1:2
            for v = 1:2
                actualLocal_N = local_N_uv{u,v};
                actualLocalStress(u, v) = actualLocal_N(numCell);
            end
        end
        [eigenvectors{numCell}, eigenvalues{numCell}] = eig(actualLocalStress);
    end
    
    disp('');
end

