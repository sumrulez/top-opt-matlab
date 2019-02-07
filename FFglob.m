function dNglob = FFglob(z,n,x,y)
	%dN - столбец производных функций формы по порядку
    %1ая строка - производные фф N1 по x
    %2ая строка - производная фф N2 по y и т.д.
	Jac = Jacobian (x,y,z,n);
	dNloc = FFloc (z,n);
	invJac = inv(Jac);
	dNglob= [
		invJac*[ dNloc(1,1); dNloc(1,2)];...
		invJac*[ dNloc(2,1); dNloc(2,2)];...
		invJac*[ dNloc(3,1); dNloc(3,2)];...
		invJac*[ dNloc(4,1); dNloc(4,2)]
	];
	
end