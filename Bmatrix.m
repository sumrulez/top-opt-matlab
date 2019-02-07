function B = Bmatrix (z,n,x,y)
	dNglob = FFglob(z,n,x,y);
	B = [dNglob(1) 0 dNglob(3) 0 dNglob(5) 0 dNglob(7) 0;
			0 dNglob(2) 0 dNglob(4) 0 dNglob(6) 0 dNglob(8);
			dNglob(2) dNglob(1) dNglob(4) dNglob(3) dNglob(6) dNglob(5) dNglob(8) dNglob(7)];
end