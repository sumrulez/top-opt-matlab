% x - вектор плотностей
% clear all
% close all
% clc
	count_elem = 16;
	count_node = 27;
	count_spc = 3;
	count_force = 3;
	E=2.1e5;  mu=0.3;
	R = openfile (count_elem);

for i=1:count_elem
	%координаты узлов 
	%строка - номер элемента
	x(i,:)= [R(i,6),R(i,9),R(i,12),R(i,15)];
	y(i,:)= [R(i,7),R(i,10),R(i,13),R(i,16)];
	z(i,:)= [R(i,8),R(i,11),R(i,14),R(i,17)];
end

IND = zeros (count_elem,4);
for i=1:count_elem
	IND (i,:) = [R(i,2),R(i,3),R(i,4),R(i,5)];
end


%составление глобальной матрицы жесткости
K_glob = zeros (2*count_node,2*count_node);
for i_el=1:count_elem
	
	K = K_el_9 (E,mu,x(i_el,:),y(i_el,:));
	
	for i_node = 1:4
		for j_node = 1:4
			i_con= IND (i_el,i_node);
			j_con= IND (i_el,j_node);
			K1 = K((2*i_node-1):(2*i_node), (2*j_node-1):(2*j_node));
			K_glob((2*i_con -1):(2*i_con) ,(2*j_con -1):(2*j_con)) = ...
				K_glob((2*i_con -1):(2*i_con) ,(2*j_con -1):(2*j_con)) + K1;
		end
	end
		
	
	
	
end

%закрепления
SPC = spc_read (count_spc); 
v = 10^15; % large number for constraint;
for i_spc = 1:count_spc
    spc_node = SPC(i_spc,2);

        if SPC(i_spc,3) == 123
        K_glob((2*spc_node-1),(2*spc_node-1)) = v;
		K_glob((2*spc_node),(2*spc_node)) = v;
		end
		
		if SPC(i_spc,3) == 23
		K_glob((2*spc_node),(2*spc_node)) = v;
					% K_glob((2*spc_node),:) = v;
					% K_glob(:,(2*spc_node)) = v;
		end
end

%силы
F = zeros ([count_node*2,1]);
FORCE = force_read (count_force);
for i=1:count_force
	force_node = FORCE (i,2);
    vector_lenght = sqrt (FORCE(i,4)^2+FORCE(i,5)^2+FORCE(i,6)^2);
    Fx = FORCE(i,3)*(FORCE(i,4)/vector_lenght);
    Fy = FORCE(i,3)*(FORCE(i,5)/vector_lenght);
	F (2*force_node-1) = -Fx;
    F (2*force_node) = -Fy;
    %силы почему то правильно задавать отрицательными
end

% A1 = chol (K_glob);
% y = transpose(A1)\F;
% % U = abs(A1\y);
% U = A1\y;
U = K_glob\F;
save_to_file (U,count_node);


%Функции формы
	% N=[ 0.25*(1-z)*(1-n ...
	% 0.25*(1+z)*(1-n) ...
	% 0.25*(1+z)*(1+n) ...
	% 0.25*(1-z)*(1+n)]









% F = 0.5* [0 0 0 0 0 1 0 1]';
% v=10^50;
% K_9 (1,1)=v; K_9 (2,2)=v;  K_9 (3,3)=v; K_9 (4,4)=v; 
% U = K_9\F;
% KE (1,1)=v; KE (2,2)=v;  KE (3,3)=v; KE (4,4)=v; 
% U1 = KE\F;
