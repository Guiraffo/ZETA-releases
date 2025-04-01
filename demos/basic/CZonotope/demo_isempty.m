demoID = 'Is this constrained zonotope an empty set?';
util.demoinit;


Z1 = CZonotope( zeros(2,1), [ 1, 0,  1;
                             1, 2, -1],...
                [-2, 1, -1], 2);
           
Z2 = CZonotope( zeros(2,1), [ 1, 0,  1;
                             1, 2, -1],...
                [1, 1, 1], 4);            
            

figure;
plot(Z1,'blue',0.3);
title('Z_1');
xlabel('x_1');
ylabel('x_2');

figure;
plot(Z2,'red',0.3);
title('Z_2');
xlabel('x_1');
ylabel('x_2');

if(isempty(Z1))
    disp('The constrained zonotope Z1 is empty.');  
else
    disp('The constrained zonotope Z1 is not empty.');  
end

if(isempty(Z2))
    disp('The constrained zonotope Z2 is empty.');  
else
    disp('The constrained zonotope Z2 is not empty.');  
end

            
           

