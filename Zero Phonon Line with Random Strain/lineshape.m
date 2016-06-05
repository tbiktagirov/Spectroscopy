function h = lineshape(t1,t2,t3,t4, x, x0, v1, v2, delta, gamma)
	%x1 = t1/(1-t1);
	%x2 = t2/(1-t2);
	%x3 = t3/(1-t3);
	%x4 = t4/(1-t4);
	%J = 1 / ((1-t1)*(1-t1)*(1-t2)*(1-t2)*(1-t3)*(1-t3)*(1-t4)*(1-t4));


	e1 = (1-t1)/t1;
	e2 = (1-t2)/t2;
	e3 = (1-t3)/t3;
	e4 = (1-t4)/t4;
	J = 1 / ((t1)*(t1)*(t2)*(t2)*(t3)*(t3)*(t4)*(t4));


	
	delta = sqrt( (e1+e2)*(e1+e2) + (e3+e4)*(e3+e4) );
	L1 = exp(-(x - x0 + delta)^2/2/delta/delta);
	L2 = exp(-(x - x0 - delta)^2/2/delta/delta);
	h = (L1+L2)/(e1*e1/v1/v1 + e2*e2/v2/v2 + e3*e3/v1/v1 + e4*e4/v2/v2 + gamma^2)^(2.5)*J; 

