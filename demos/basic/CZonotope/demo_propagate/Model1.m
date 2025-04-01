function f = Model1(x)

% System equations 

f = [ 3*x(1) - (1/7)*x(1)^2 - (4*x(1)*x(2))/(4 + x(1));
     -2*x(2) + (3*x(1)*x(2))/(4 + x(1))               ];
            

end