function save__dens_to_file (dens,count_elem,it)
	
% write matrix in read4.txt

if it ==1 
    fid = fopen('densities.hwascii','wt');
fprintf(fid,'%s \n%s \n%s \n%s \n%s \n%s%d%s\n', ...
	'ALTAIR ASCII FILE','$SUBCASE = 1 "MatrixBrowser_2"','$BINDING = ELEMENT',...
	'$COLUMN_INFO=	ENTITY_ID','$RESULT_TYPE = Density (s)','$TIME = ',it,'.0 sec');
end

if it > 1
    fid = fopen('densities.hwascii','a+');
    fprintf(fid,'%s%d%s \n','$TIME = ',it,'.0 sec');
end
for i = 1:count_elem
      fprintf(fid,'%d %f\n',i,dens(i));
end
fclose(fid);

end