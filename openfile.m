function R = openfile(count_elem)
% 	
% file = fopen ('mesh.txt');
	%выделение памяти
	R = zeros (count_elem,17);
	R = importdata('mesh.txt',',')
% for i=1:count_elem
% 	%R-matrix 1 - номер элемента
% 	%		2-5 - узелы
% 	% 5-17 координаты узлов xyz соответственно
% 	R(i,:) = fscanf (file, '%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,/n',17);
% end
% fclose(file);
end
