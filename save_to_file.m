function save_to_file (U,count_node)
	
% write matrix in read4.txt
fid = fopen('results.hwascii','wt');
fprintf(fid,'%s \n%s \n%s \n%s \n%s \n', ...
	'ALTAIR ASCII FILE','$SUBCASE = 1 "MatrixBrowser_2"','$BINDING = NODE',...
	'$COLUMN_INFO=	ENTITY_ID','$RESULT_TYPE = Displacements (v)');
	
	z= 0;
		p=1;
for i = 1:count_node
      fprintf(fid,'%d %f %f %f\n',i,U(p),U(p+1),z);
	p=p+2;
end
fclose(fid);

end