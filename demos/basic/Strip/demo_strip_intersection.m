demoID = 'Intersection of two parallel strips';
util.demoinit;

S1 = Strip( [-1; 1], 1, 0.5);
S2 = Strip( [-1; 1], 0.5, 0.5);
S3 = Strip( [1; 1], 0.5, 0.5);

disp(S1);
disp(S2);
disp(S3);

S12int = intersection(S1,S2);
try
    S13int = intersection(S1,S3); % This is just to tr
catch
    disp('S1 and S3 are not parallel! The implemented intersection algorithm cannot proceed.');
end
     
canvas = [-1,1,0,2];

% Strips are unbounded... plotting them requires some tricks
figure
hold on;
plot(S1,'blue',0.8,canvas);  
plot(S2,'red',0.8,canvas); 
plot(S12int,'magenta',0.8,canvas); 
grid off;
box on;
axis(canvas);
xlabel('x_1');
ylabel('x_2');
legend('S_1','S_2','S_1 \cap S_2');