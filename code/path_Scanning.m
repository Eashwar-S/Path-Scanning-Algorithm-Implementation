%% Creating a sample undirected graph

clc;
clear all;
close all;
figure;
s = [1 1 2 2 3 3 4 4 4 5 6 7];
t = [2 3 4 6 4 5 5 7 6 8 7 8];
G = graph(s,t);
h = plot(G)
h.XData = [-2 -0.5 -1   0 1  1.5 2   2.5];
h.YData = [ 0 -2    2.5 0 3 -2   0.3 1.5];
G.Nodes.NodeColors = [6; 1; 1; 1; 6; 1; 1; 1];
G.Edges.EdgeColors = [1; 1; 5; 1; 1; 1; 1; 1; 1; 1; 5; 1];
h.NodeCData = G.Nodes.NodeColors;
h.EdgeCData = G.Edges.EdgeColors;
h.EdgeLabel = [2.3, 2, 3, 1.5, 3.2, 2.2, 3.8, 2.6, 2.2, 2.8, 1.8, 0.8];
h.MarkerSize = 5.5;
h.LineWidth = 2;
h.NodeFontSize = 12;
h.EdgeFontSize = 12;
%% Initializing Parameters


numNodes = 8;
numEdges = 12;
vehicleCapacity = 14;
graphMatrix = [  0 2.3 2 0 0   0 0 0;
               2.3   0 0 3 0 1.5 0 0;
                 2   0 0 3.2 2.2 0 0 0;
                 0 3 3.2 0 3.8 2.6 2.2 0;
                 0 0 2.2 3.8 0 0 0 2.8;
                 0 1.5 0 2.6 0 0 1.8 0;
                 0 0 0 2.2 0 1.8 0 0.8;
                 0 0 0 0 2.8 0 0.8 0];

requiredNodes = [2 4; 6 7];
uavsInDepotNodes = [0 2];
totalUavs = sum(uavsInDepotNodes);
numrequiredEdges = 2;
depotNodes = [1 5];
arcTraversed = false;
visitedNodes = [];
taskAllocatedtoBaseStations = [];
exploredNodes.nodes = [];
allExploredNodes = repmat(exploredNodes, numNodes,1);
path.Tour = [];
path.Cost = [];
bestPath = repmat(path,numel(depotNodes),1);
%% Task Allocation: Determining which acrs should be traversed by which base station.

depotNodesCost = zeros(numrequiredEdges, numel(depotNodes));
for i=1:numel(depotNodes)
    index = depotNodes(i);
    for j=1:numrequiredEdges
        arcTraversed = false;
        distance = 0;
        visitedNodes = [];
        while ~arcTraversed
            minCost = 0;
            for k=1:numNodes
                if graphMatrix(index,k) ~= 0 && ~ismember(k,visitedNodes)
                    if k == requiredNodes(j,1) || k == requiredNodes(j,2)
                        depotNodesCost(i,j) = distance + graphMatrix(index,k) + graphMatrix(requiredNodes(j,1), requiredNodes(j,2));
                        arcTraversed = true;
                        if numel(visitedNodes) == 0
                            if k == requiredNodes(j,1)
                                visitedNodes = [index k requiredNodes(j,2)];
                            else
                                visitedNodes = [index k requiredNodes(j,1)];
                            end
                        else
                           if k == requiredNodes(j,1)
                                visitedNodes = [visitedNodes k requiredNodes(j,2)];
                            else
                                visitedNodes = [visitedNodes k requiredNodes(j,1)];
                            end 
                        end
%                         visitedNodes
                        break;
                    else
                        if minCost == 0
                            minCost = graphMatrix(index,k) + min(abs(requiredNodes(j,1) - k), abs(requiredNodes(j,2) - k));                            
                            tempIndex = k;
                            if numel(visitedNodes) == 0
                                visitedNodes = [index];
                            end
%                             visitedNodes = [index]; 
                        else
                            if graphMatrix(index,k) + min(abs(requiredNodes(j,1) - k), abs(requiredNodes(j,2) - k)) < minCost
                                minCost = graphMatrix(index,k); %+  + min(abs(requiredNodes(j,1) - k), abs(requiredNodes(j,2) - k));                               
                                tempIndex = k;
                            end
                        end
                    end
                end
            end          
            if ~arcTraversed
                index = tempIndex;
                visitedNodes = [visitedNodes tempIndex];                  
            end
            distance = distance + minCost;
        end
    end
end

for i=1:numel(uavsInDepotNodes)
    if uavsInDepotNodes(i) == 0
        depotNodesCost(i,:) = inf;
    end
end

for i=1:numel(depotNodes)
    disp(['Allocating arc ' num2str(requiredNodes(i,1)) ' - ' num2str(requiredNodes(i,2)) ' to base station ' num2str(i) ' - node ' num2str(depotNodes(find(depotNodesCost(:,i) == min(depotNodesCost(:,i)))))])
    taskAllocatedtoBaseStations = [taskAllocatedtoBaseStations depotNodes(find(depotNodesCost(:,i) == min(depotNodesCost(:,i))))];
end
%%
% Task Allocation Visualization - Part 1
figure;
s = [1 1 2 2 3 3 4 4 4 5 6 7];
t = [2 3 4 6 4 5 5 7 6 8 7 8];
G = graph(s,t);
h = plot(G)
h.XData = [-2 -0.5 -1   0 1  1.5 2   2.5];
h.YData = [ 0 -2    2.5 0 3 -2   0.3 1.5];
G.Nodes.NodeColors = [6; 1; 1; 1; 6; 1; 1; 1];
G.Edges.EdgeColors = [1; 1; 5; 1; 1; 1; 1; 1; 1; 1; 5; 1];
h.NodeCData = G.Nodes.NodeColors;
h.EdgeCData = G.Edges.EdgeColors;
h.EdgeLabel = [2.3, 2, 3, 1.5, 3.2, 2.2, 3.8, 2.6, 2.2, 2.8, 1.8, 0.8]
h.MarkerSize = 5.5;
h.LineWidth = 2;
h.NodeFontSize = 12;
h.EdgeFontSize = 12;

hold on
s1 = [4,5]
t1 = [2,4]
G1 = digraph(s1, t1);
h1 = plot(G1,'Layout','force')
h1.XData = [-2 -0.5 -1  0 1];
h1.YData = [ 0 -2 2.5  0 3];
G1.Nodes.NodeColors = [6; 1; 1; 1; 6];
G1.Edges.EdgeColors = [10; 10];
h1.NodeCData = G1.Nodes.NodeColors;
h1.EdgeCData = G1.Edges.EdgeColors;
h1.EdgeLabel = [3, 3.8];
h1.MarkerSize = 5.5;
h1.LineWidth = 2;
h1.NodeFontSize = 12;
h1.EdgeFontSize = 12;
h1.ArrowSize = 16;
%%
% Task Allocation Visualization - Part 2

figure;
s = [1 1 2 2 3 3 4 4 4 5 6 7];
t = [2 3 4 6 4 5 5 7 6 8 7 8];
G = graph(s,t);
h = plot(G)
h.XData = [-2 -0.5 -1   0 1  1.5 2   2.5];
h.YData = [ 0 -2    2.5 0 3 -2   0.3 1.5];
G.Nodes.NodeColors = [6; 1; 1; 1; 6; 1; 1; 1];
G.Edges.EdgeColors = [1; 1; 5; 1; 1; 1; 1; 1; 1; 1; 5; 1];
h.NodeCData = G.Nodes.NodeColors;
h.EdgeCData = G.Edges.EdgeColors;
h.EdgeLabel = [2.3, 2, 3, 1.5, 3.2, 2.2, 3.8, 2.6, 2.2, 2.8, 1.8, 0.8]
h.MarkerSize = 5.5;
h.LineWidth = 2;
h.NodeFontSize = 12;
h.EdgeFontSize = 12;

hold on
s2 = [5,7,8]
t2 = [8,6,7]
G2 = digraph(s2, t2);
h2 = plot(G2,'Layout','force')
h2.XData = [-2 -0.5 -1   0 1  1.5 2   2.5];
h2.YData = [ 0 -2    2.5 0 3 -2   0.3 1.5];
G2.Nodes.NodeColors = [6; 1; 1; 1; 6; 1; 1; 1];
G2.Edges.EdgeColors = [3; 3; 3];
h2.NodeCData = G2.Nodes.NodeColors;
h2.EdgeCData = G2.Edges.EdgeColors;
h2.EdgeLabel = [2.8, 1.8, 0.8];
h2.MarkerSize = 5.5;
h2.LineWidth = 2;
h2.NodeFontSize = 12;
h2.EdgeFontSize = 12;
h2.ArrowSize = 16;
xlim([-2.3 3.2])
ylim([-2.80 3.80])
%% Path Scanning Algorithm

% Constructing feasible cycles one at a time considering the following
% optimization conditions:
% 1. Decreasing the cost taken to travel from the base station to desired
% arc.
% 2. Decreasing the cost taken to travel from the desired arc back to base
% station.

for i=1:numel(depotNodes)
    index = taskAllocatedtoBaseStations(i)
    while true
        arcTraversed = false;
        cycleFormed = false;
        distance = 0;
        visitedNodes = [];
        while ~(arcTraversed && cycleFormed)
            minCost = 0;            
            for k=1:numNodes                
                if graphMatrix(index,k) ~= 0 && ~ismember(k,visitedNodes)
                    if ~arcTraversed
                        if k == requiredNodes(i,1) || k == requiredNodes(i,2)                                                       
                            bestPath(i).Cost = distance + graphMatrix(index,k) + graphMatrix(requiredNodes(i,1), requiredNodes(i,2));
                            arcTraversed = true;
                            if numel(visitedNodes) == 0
                                if k == requiredNodes(i,1)
                                    visitedNodes = [k requiredNodes(i,2)];                                  
                                else
                                    visitedNodes = [k requiredNodes(i,1)];
                                end
                            else
                               if k == requiredNodes(i,1)
                                    visitedNodes = [visitedNodes k requiredNodes(i,2)];                                 
                                   
                                else
                                    visitedNodes = [visitedNodes k requiredNodes(i,1)];
                               end 
                            end
                            tempIndex = visitedNodes(end);
                            index = tempIndex;
                            bestPath(i).Tour = visitedNodes;                            
                            visitedNodes = [];
                            distance = 0;
                            minCost = 0;
                        else
                            if minCost == 0                                                                
                                sum(ismember(allExploredNodes(index).nodes,k))
                                if sum(ismember(allExploredNodes(index).nodes,k)) == 0
                                    minCost = graphMatrix(index,k)                            
                                    tempIndex = k;
                                    if numel(visitedNodes) == 0
                                        visitedNodes = [index];
                                    end
                                end                   
                            else
                                if graphMatrix(index,k) < minCost && sum(ismember(allExploredNodes(index).nodes,k)) == 0
                                    minCost = graphMatrix(index,k);                               
                                    tempIndex = k;
                                end
                            end
                        end                     
                    else                        
                        for m=1:numNodes
                            if graphMatrix(index,m) ~= 0 && numel(find(ismember(depotNodes,m)==1)) ~= 0                                
                                bestPath(i).Cost = bestPath(i).Cost + distance + graphMatrix(index,m);
                                cycleFormed = true;
                                visitedNodes = [visitedNodes m];                                
                                bestPath(i).Tour = [bestPath(i).Tour visitedNodes];
                                break;
                            end
                        end
                        if ~cycleFormed                            
                            if minCost == 0                                
                                if sum(ismember(allExploredNodes(index).nodes,k)) == 0
                                    minCost = graphMatrix(index,k)                            
                                    tempIndex = k;      
                                end                   
                            else
                                if graphMatrix(index,k) < minCost && sum(ismember(allExploredNodes(index).nodes,k)) == 0
                                    minCost = graphMatrix(index,k); %+  + min(abs(requiredNodes(j,1) - k), abs(requiredNodes(j,2) - k));                               
                                    tempIndex = k;
                                end
                            end                            
                        else
                            break;
                        end              
                                                           
                    end
                end               
            end          
            if ~cycleFormed
                index = tempIndex;
                visitedNodes = [visitedNodes tempIndex];          
            end            
        end       
        if bestPath(i).Cost <= vehicleCapacity
            break;
        else
            allExploredNodes(bestPath(i).Tour(1)).nodes = [allExploredNodes(bestPath(i).Tour(1)).nodes bestPath(i).Tour(2)]
            bestPath(i).Cost = [];
            bestPath(i).Tour = [];
            index = taskAllocatedtoBaseStations(i);
        end
    end
end

%%
%Adding depot nodes to the solution path and removing redundant nodes from the
%solution path
for i=1:numel(depotNodes)
    temp = [];
    numel(bestPath(i).Tour)
    if depotNodes(i) ~= bestPath(i).Tour(1)
       bestPath(i).Tour = [depotNodes(i) bestPath(i).Tour];
    end    
    for j=1:numel(bestPath(i).Tour)-1           
       if bestPath(i).Tour(j) ~= bestPath(i).Tour(j+1)
           if j == numel(bestPath(i).Tour)-1
               temp = [temp bestPath(i).Tour(j) bestPath(i).Tour(j+1)];           
           else
               temp = [temp bestPath(i).Tour(j)];
           end           
       end       
    end
    bestPath(i).Tour = temp;
end
%%
% Path 1 Visualization
figure;
s = [1 1 2 2 3 3 4 4 4 5 6 7];
t = [2 3 4 6 4 5 5 7 6 8 7 8];
G = graph(s,t);
h = plot(G)
h.XData = [-2 -0.5 -1   0 1  1.5 2   2.5];
h.YData = [ 0 -2    2.5 0 3 -2   0.3 1.5];
G.Nodes.NodeColors = [6; 1; 1; 1; 6; 1; 1; 1];
G.Edges.EdgeColors = [1; 1; 5; 1; 1; 1; 1; 1; 1; 1; 5; 1];
h.NodeCData = G.Nodes.NodeColors;
h.EdgeCData = G.Edges.EdgeColors;
h.EdgeLabel = [2.3, 2, 3, 1.5, 3.2, 2.2, 3.8, 2.6, 2.2, 2.8, 1.8, 0.8]
h.MarkerSize = 5.5;
h.LineWidth = 2;
h.NodeFontSize = 12;
h.EdgeFontSize = 12;

hold on
s1 = [2,4,5]
t1 = [1,2,4]
G1 = digraph(s1, t1);
h1 = plot(G1,'Layout','force')
h1.XData = [-2 -0.5 -1  0 1];
h1.YData = [ 0 -2 2.5  0 3];
G1.Nodes.NodeColors = [6; 1; 1; 1; 6];
G1.Edges.EdgeColors = [3; 3; 3];
h1.NodeCData = G1.Nodes.NodeColors;
h1.EdgeCData = G1.Edges.EdgeColors;
h1.EdgeLabel = [2.3, 3, 3.8];
h1.MarkerSize = 5.5;
h1.LineWidth = 2;
h1.NodeFontSize = 12;
h1.EdgeFontSize = 12;
h1.ArrowSize = 16;
%%
% Path 2 Visualization
figure;
s = [1 1 2 2 3 3 4 4 4 5 6 7];
t = [2 3 4 6 4 5 5 7 6 8 7 8];
G = graph(s,t);
h = plot(G)
h.XData = [-2 -0.5 -1   0 1  1.5 2   2.5];
h.YData = [ 0 -2    2.5 0 3 -2   0.3 1.5];
G.Nodes.NodeColors = [6; 1; 1; 1; 6; 1; 1; 1];
G.Edges.EdgeColors = [1; 1; 5; 1; 1; 1; 1; 1; 1; 1; 5; 1];
h.NodeCData = G.Nodes.NodeColors;
h.EdgeCData = G.Edges.EdgeColors;
h.EdgeLabel = [2.3, 2, 3, 1.5, 3.2, 2.2, 3.8, 2.6, 2.2, 2.8, 1.8, 0.8]
h.MarkerSize = 5.5;
h.LineWidth = 2;
h.NodeFontSize = 12;
h.EdgeFontSize = 12;

hold on
s2 = [2,5,6,7,8]
t2 = [1,8,2,6,7]
G2 = digraph(s2, t2);
h2 = plot(G2,'Layout','force')
h2.XData = [-2 -0.5 -1   0 1  1.5 2   2.5];
h2.YData = [ 0 -2    2.5 0 3 -2   0.3 1.5];
G2.Nodes.NodeColors = [6; 1; 1; 1; 6; 1; 1; 1];
G2.Edges.EdgeColors = [3; 3; 3; 3; 3];
h2.NodeCData = G2.Nodes.NodeColors;
h2.EdgeCData = G2.Edges.EdgeColors;
h2.EdgeLabel = [2.3, 2.8, 1.5, 1.8, 0.8];
h2.MarkerSize = 5.5;
h2.LineWidth = 2;
h2.NodeFontSize = 12;
h2.EdgeFontSize = 12;
h2.ArrowSize = 16;
xlim([-2.3 3.2])
ylim([-2.80 3.80])