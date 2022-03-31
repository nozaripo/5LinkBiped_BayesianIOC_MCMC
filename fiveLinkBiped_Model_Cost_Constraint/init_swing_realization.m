function ceq = init_swing_realization(x0,p)

h = cos(x0(1))*p.l1 + cos(x0(2))*p.l2;

ceq = cos(x0(4))*p.l4 + cos(x0(5))*p.l5 - h;
