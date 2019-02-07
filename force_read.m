function R = force_read (force_num)
	%первая колонка: id force
	%вторая колонка: номер узла
	%третья колонка: магнитуда
	%4-6: вектор в формате d,d,d
file = fopen ('force.txt');
	%выделение памяти
	R = zeros (force_num,6); 
	
for i=1:force_num
	%R-matrix 1 - номер элемента
	%		2-5 - узелы
	% 5-17 координаты узлов xyz соответственно
	R(i,:) = fscanf (file, '%f,%f,%f,%f,%f,%f/n',6);
end
fclose(file);
end