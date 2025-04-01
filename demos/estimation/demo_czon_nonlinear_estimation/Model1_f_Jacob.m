function f_Jacob = Model1_f_Jacob(x,w,u,Ts)


% System dynamics (Jacobian)

x1 = x(1);
x2 = x(2);

f_Jacob = [ 3 - (16*x2)/(x1 + 4)^2 - (2*x1)/7, 16/(x1 + 4) - 4,   1,  0;
                           (12*x2)/(x1 + 4)^2, 1 - 12/(x1 + 4),   0,  1];
          

end