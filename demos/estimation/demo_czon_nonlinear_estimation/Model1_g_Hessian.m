function g_Hessian = Model1_g_Hessian(x,v)


% Output equation (Hessian)

x1 = x(1);
x2 = x(2);

g_Hessian{1} = [ 0,           0, 0, 0;
                 0, sin(x2/2)/4, 0, 0;
                 0,           0, 0, 0;
                 0,           0, 0, 0];

g_Hessian{2} = Interval([  0, -1, 0, 0;
                          -1,  0, 0, 0;
                           0,  0, 0, 0;
                           0,  0, 0, 0]);
 

end

