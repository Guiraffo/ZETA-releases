function f_Hessian = Model1_f_Hessian(x,u,w,Ts)


% System dynamics (Hessian)

x1 = x(1);
x2 = x(2);


f_Hessian = cell(2,1);

f_Hessian{1} = [ (32*x2)/(x1 + 4)^3 - 2/7, -16/(x1 + 4)^2, 0, 0;
                           -16/(x1 + 4)^2,              0, 0, 0;
                                        0,              0, 0, 0;
                                        0,              0, 0, 0];

f_Hessian{2} = [ -(24*x2)/(x1 + 4)^3, 12/(x1 + 4)^2, 0, 0;
                       12/(x1 + 4)^2,             0, 0, 0;
                                   0,             0, 0, 0;
                                   0,             0, 0, 0];

end

