function f = Model1_f(x,w,u,Ts)

% System equations 

f = [ x(1)*(3 - (1/7)*x(1) - (4*x(2))/(4 + x(1))) + w(1);
      x(2)*(-2 + (3*x(1))/(4 + x(1)))             + w(2)];
            

end