function [forceInferenceValue, edgeInfo] = readDatFile( fileName, correspondingImage, newSize )
%READDATFILE Summary of this function goes here
%   Detailed explanation goes here

    fileID = fopen(fileName);

    rowFile = fgetl(fileID);
    rowFile = fgetl(fileID);
    
    edgesTensions = true;
    cellInfo = [];
    edgeInfo = [];
    
    imgLabelled = bwlabel(correspondingImage);
    
    [neighbours, sideCells] = calculateNeighbours(imgLabelled);
    [ verticesInfo ] = calculateVertices( imgLabelled, neighbours);
    allVertices = vertcat(verticesInfo.verticesPerCell{:});
    
    if isempty(newSize) == 0
        allVertices = allVertices .* newSize;
    end
    
    while ischar(rowFile)
        
        if isempty(rowFile) == 0
            lineSplitted = strsplit(rowFile);
            if isempty(lineSplitted{1})
                lineSplitted(1) = [];
            end

            if edgesTensions %Tension of the edges
                vertex1_X = lineSplitted{5};
                vertex1_Y = lineSplitted{6};
                vertex2_X = lineSplitted{7};
                vertex2_Y = lineSplitted{8};
                %Id of edge
                edgeInfo(end+1, 1:2) = [str2double(vertex1_X(2:end)), str2double(vertex1_Y(1:end-1))];
                edgeInfo(end, 3:4) = [str2double(vertex2_X(2:end)), str2double(vertex2_Y(1:end-1))];
                edgeInfo(end, 5) = pdist(vertcat(edgeInfo(end, [2 1]), edgeInfo(end, [4 3])));
                edgeInfo(end, 6) = str2double(lineSplitted{2});
                
                %First we should do a correspondance with the real vertices
                %of frusta. To do this, we calculate the neighbours in the 
                %inital frames and, then, extrapolate the initial vertices
                % of the actual image.
                [minValue, indexV1] = pdist2(allVertices, edgeInfo(end, [2 1]), 'Euclidean', 'Smallest',1);
                [minValue, indexV2] = pdist2(allVertices, edgeInfo(end, [4 3]), 'Euclidean', 'Smallest', 1);
                %realVertex1 = allVertices(indexV1, :);
                %realVertex2 = allVertices(indexV2, :);
                
                neighboursVertex1 = verticesInfo.verticesConnectCells(indexV1, :);
                neighboursVertex2 = verticesInfo.verticesConnectCells(indexV2, :);
                
%                 %Associate this edge with its cells
%                 emptyImage(edgeInfo(end, 2), edgeInfo(end, 1)) = 1;
%                 neighboursVertex1 = unique(imgLabelled.*imdilate(emptyImage, dilateShape));
%                 emptyImage(edgeInfo(end, 2), edgeInfo(end, 1)) = 0;
%                 
%                 emptyImage(edgeInfo(end, 4), edgeInfo(end, 3)) = 1;
%                 neighboursVertex2 = unique(imgLabelled.*imdilate(emptyImage, dilateShape));
%                 emptyImage(edgeInfo(end, 4), edgeInfo(end, 3)) = 0;
%                 
                cellsOfTheEdge = intersect(neighboursVertex1, neighboursVertex2);
                cellsOfTheEdge(cellsOfTheEdge == 0) = [];
                if length(cellsOfTheEdge) == 2
                    edgeInfo(end, 7:8) = cellsOfTheEdge;
                end
            else %Cell preassure
                cellInfo(end+1, :) = [str2double(lineSplitted{3})+1, str2double(lineSplitted{4})];
            end
        elseif edgesTensions == true
            edgesTensions = false;
            edgeInfo = array2table(edgeInfo, 'VariableNames',{'vertex1_X', 'vertex1_Y', 'vertex2_X', 'vertex2_Y', 'EdgeLength' 'TensionValue', 'ConnectingCell1', 'ConnectingCell2'});
        end
        
        rowFile = fgetl(fileID);
    end
    
    fclose(fileID);
    
    
    for numCell = 1:size(cellInfo, 1)
         actualTensionValues = edgeInfo.TensionValue(edgeInfo.ConnectingCell1 == numCell | edgeInfo.ConnectingCell2 == numCell);
         cellInfo(numCell, 3) = mean(actualTensionValues);
         cellInfo(numCell, 4) = std(actualTensionValues);
         cellInfo(numCell, 5) = length(actualTensionValues);
         cellInfo(numCell, 6) = sideCells(numCell);
    end
    
    cellsMorphology = regionprops(imgLabelled, 'all');
    areaOfCells = [cellsMorphology.Area];
    majorAxisOfCells = [cellsMorphology.MajorAxisLength];
    minorAxisOfCells = [cellsMorphology.MinorAxisLength];
    
    angleOfCells = [cellsMorphology.Orientation];
    
    validCells = false(size(areaOfCells, 1), 1);
    validCells(1:size(cellInfo, 1)) = (cellInfo(:, 5) ~= cellInfo(:, 6)) & isnan(cellInfo(:, 3)) == 0;
    
    computeGlobalStress(areaOfCells, cellInfo, edgeInfo, validCells);
    
    [correlation, pvalue] = corrcoef(cellInfo(:, 2:5), 'Rows', 'pairwise');
    forceInferenceValue = array2table(cellInfo, 'VariableNames', {'CellID', 'PressureValue', 'MeanTension', 'STDTension', 'NumEdgesOfTension', 'RealSides'});
    
    
end

