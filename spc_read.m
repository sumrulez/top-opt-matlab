function R = spc_read (spc_num)
	%первая колонка: id spc
	%вторая колонка: номер узла
	%третья колонка: степени свободы в формате 123456
file = fopen ('spc.txt');
	%выделение памяти
	R = zeros (spc_num,3); 
	
for i=1:spc_num
	%R-matrix 1 - номер элемента
	%		2-5 - узелы
	% 5-17 координаты узлов xyz соответственно
	R(i,:) = fscanf (file, '%f,%f,%f,/n',3);
end
fclose(file);
end
