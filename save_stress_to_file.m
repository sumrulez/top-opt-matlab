function save_stress_to_file (stress,count_elem,it)

%file header
if it ==1 
    fid = fopen('stresses.hwascii','wt');
fprintf(fid,'%s \n%s \n%s \n%s \n%s \n%s%d%s\n', ...
	'ALTAIR ASCII FILE','$SUBCASE = 1 "MatrixBrowser_2"','$BINDING = ELEMENT',...
	'$COLUMN_INFO=	ENTITY_ID','$RESULT_TYPE = Stress_VM (s)','$TIME = ',it,'.0 sec');
end

if it > 1
    fid = fopen('stress.hwascii','a+');
    fprintf(fid,'%s%d%s \n','$TIME = ',it,'.0 sec');
end

for i = 1:count_elem
      fprintf(fid,'%d %f\n',i,stress(i));
end
fclose(fid);

end