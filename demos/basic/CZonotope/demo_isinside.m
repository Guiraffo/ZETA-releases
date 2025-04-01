demoID = 'Point is inside a constrained zonotope';
util.demoinit;


Z = CZonotope( zeros(2,1), [ 1, 0,  1;
                             1, 2, -1],...
                [-2, 1, -1], 2);

z{1} = [-1; 1];            
z{2} = [-1.8; -1.5];
z{3} = [0; 2];
z{4} = [-0.4; 1.5];


% Verify if the points z are inside Z
isin = zeros(length(z),1);
color = cell(length(z),1);
result = cell(length(z),1);
for j=1:4
    isin(j) = isinside(Z,z{j});
    if(isin(j))
        color{j} = 'b';
        result{j} = 'inside';
    else
        color{j} = 'r';
        result{j} = 'outside';   
    end        
end

       
% Plots Z, the points z, and displays the result in the command window
figure;
hold on;
plot(Z,'green',0.2);
for j=1:4
    plot(z{j}(1), z{j}(2), [color{j},'x']);
    disp(['The point (',num2str(z{j}(1)),',',num2str(z{j}(2)),') is ',result{j},' Z.']);    
end
xlabel('x_1');
ylabel('x_2');

            
           
