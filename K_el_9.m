function K = K_el_9 (E,mu,x,y)
	% 3 гауссовы точки, 9 точки интегрирования
	g1 = -sqrt(3/5);
	g2 = 0;
	g3 = sqrt(3/5);
	w1 = 5/9;
	w2 = 8/9;
	w3 = 5/9;
 	D = E/(1-mu^2)*[1 mu 0;
					mu 1 0;
					0 0 0.5*(1-mu)];
	
	b11 = Bmatrix(g1,g1,x,y); w11 = w1*w1; 	Jac11 = Jacobian(x,y,g1,g1);	J11=det(Jac11);
	b12 = Bmatrix(g1,g2,x,y); w12 = w1*w2;  Jac12 = Jacobian(x,y,g1,g2);	J12=det(Jac12);
	b21 = Bmatrix(g2,g1,x,y); w21 = w12; 	Jac21 = Jacobian(x,y,g2,g1);	J21=det(Jac21);
	b22 = Bmatrix(g2,g2,x,y); w22 = w2*w2;  Jac22 = Jacobian(x,y,g2,g2);	J22=det(Jac22);
	b13 = Bmatrix(g1,g3,x,y); w13 = w1*w3;  Jac13 = Jacobian(x,y,g1,g3);	J13=det(Jac13);
	b31 = Bmatrix(g3,g1,x,y); w31 = w13;    Jac31 = Jacobian(x,y,g3,g1);	J31=det(Jac31);
	b23 = Bmatrix(g2,g3,x,y); w23 = w2*w3;  Jac23 = Jacobian(x,y,g2,g3);	J23=det(Jac23);
	b32 = Bmatrix(g3,g2,x,y); w32 = w23;    Jac32 = Jacobian(x,y,g3,g2);	J32=det(Jac32);
	b33 = Bmatrix(g3,g3,x,y); w33 = w3*w3;  Jac33 = Jacobian(x,y,g3,g3);	J33=det(Jac33);
	K = w11*b11'*D*b11*J11+...
	w12*b12'*D*b12*J12+...
	w21*b21'*D*b21*J21+...
	w22*b22'*D*b22*J22+...
	w13*b13'*D*b13*J13+...
	w31*b31'*D*b31*J31+...
	w23*b23'*D*b23*J23+...
	w32*b32'*D*b32*J32+...
	w33*b33'*D*b33*J33;
end