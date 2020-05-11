clear, clc

% Loading mesh information


entable_coarse = load('entable_ab2.txt');
entable_coarse(:,1)=[];

node_coarse = load('node_ab2.txt');
node_coarse(:,1)=[];
node_coarse(:,3)=[];



EN_coarse = length(entable_coarse);
for i = 1 : EN_coarse
    
    grid_coarse = entable_coarse(i, :);
    X_coarse = node_coarse(grid_coarse,1);
    Y_coarse = node_coarse(grid_coarse,2);
    patch(X_coarse, Y_coarse, 'b')
    
end
axis equal