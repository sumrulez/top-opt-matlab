function Jac = Jacobian(x,y,z,n)
	dNloc = FFloc (z,n);
	dxdz = x*dNloc(:,1);
	dydz = y*dNloc(:,1);
	dxdn = x*dNloc(:,2);
	dydn = y*dNloc(:,2);
	Jac = [ dxdz dydz;...
	dxdn dydn];
	end