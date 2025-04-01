demoID = 'CZ inclusion: enclosing of the product of an interval matrix by a CZ';
util.demoinit;


%% Input variables

J = [Interval(2,3), Interval(1,1.5); Interval(0,2), Interval(-1,0.5)];

p = ones(2,1);
M = [ 2.5000,    1.2500,    1.5000,   0.1;
     -0.7500,    0.5000,    3.5000,   0.1];

A = [-2, 1, -1, 1;
      1, 1, -1, 1];
b = [2; 1];

X = CZonotope(p,M,A,b);

%% CZ inclusion

Z = inclusion(X,J);
         
         
%% Generating the vertices of the set J for plotting the product of these with X

Jarray = reshape(J,size(J,1)*size(J,2),1);

B = Jarray;
               
box_dimension = size(B,1);
c = mid(B);
G = 0.5*diag(diam(B));

nof_vertices = 2^box_dimension;
unitary_box_vertex = zeros(box_dimension,1);
box_vertices = zeros(box_dimension,nof_vertices);

% Box vertices by affine transformation of each vertex of the unitary box
for i=0:nof_vertices-1
    for j=0:box_dimension-1
        if bitand(i,2^j)~=2^j
            unitary_box_vertex(j+1) = -1;
        else
            unitary_box_vertex(j+1) = 1;
        end
    end
    box_vertices(:,i+1) = c + G*unitary_box_vertex;
end         


CZonotope_Jarray = cell(nof_vertices,1);
for j=1:nof_vertices
    CZonotope_Jarray{j} = CZonotope(reshape(box_vertices(:,j),2,2)*p, reshape(box_vertices(:,j),2,2)*M,A,b);
end

%% Figures

figure;
plot(X,'yellow',0.3);
xlabel('x_1')
ylabel('x_2')
title('Constrained zonotope X');


figure;
hold on;
plot(mid(J)*X,'yellow',0);
for j=1:length(CZonotope_Jarray)
    plot(CZonotope_Jarray{j},'yellow',0);
end
plot(Z,'magenta',0.2);
xlabel('x_1')
ylabel('x_2')
title('Inclusion of JX');