function [U,max_disp,stress_VM,sigma_x,sigma_y,tau_xy] = FEsolve_stress (E,mu,dens,penal,IND,SPC,FORCE,Kloc,x,y,stress_flag)
% x - вектор плотностей

%данные модели
global count_elem count_spc count_force count_node    
%составление глобальной матрицы жесткости
K_glob = zeros (2*count_node,2*count_node);
for i_el=1:count_elem
	
	K = dens(i_el)^penal*Kloc(:,:,i_el);
	
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
fixeddof = zeros(count_spc,1);
n=1;
for i_spc = 1:count_spc
    spc_node = SPC(i_spc,2);
        if SPC(i_spc,3) == 123456
          fixeddof (n:1:n+1) = [2*spc_node-1 2*spc_node];
            n=n+2;
        end
        if SPC(i_spc,3) == 123
          fixeddof (n:1:n+1) = [2*spc_node-1 2*spc_node];
            n=n+2;
        end
		if SPC(i_spc,3) == 23456
          fixeddof (n) = [2*spc_node];
        end
        if SPC(i_spc,3) == 23
          fixeddof (n) = [2*spc_node];
		end
end
fixeddof = sort (fixeddof);
alldof = [1:1:count_node *2]';
freedof = setdiff (alldof,fixeddof);
clear alldof fixeddof
%силы
F = zeros ([count_node*2,1]);
for i=1:count_force
	force_node = FORCE (i,2);
    vector_lenght = sqrt (FORCE(i,4)^2+FORCE(i,5)^2+FORCE(i,6)^2);
    Fx = FORCE(i,3)*(FORCE(i,4)/vector_lenght);
    Fy = FORCE(i,3)*(FORCE(i,5)/vector_lenght);
	F (2*force_node-1) = Fx;
    F (2*force_node) = Fy;
end
U = zeros(count_node*2,1);
Kglob_sp = sparse (K_glob);
F_sp = sparse (F);
%вычищение памяти перед решением СЛАУ
clear K1 K F K_glob
U(freedof) = Kglob_sp(freedof,freedof)\F_sp(freedof);
max_disp = max(U);
save_to_file (U,count_node);

if stress_flag == 1
%Вычисление напряжений
strain_g = zeros (count_elem,3,9); %(кол-во элементов,кол-во компонентов напряжений, кол-во точек интегрирования)
stress_g = zeros (count_elem,3,9);
%гауссовы точки
g1 = -sqrt(3/5);
g2 = 0;
g3 = sqrt(3/5);
t_int = [g1   g1 ;
         g1  g2 ;
         g2  g1;
         g2  g2 ;
         g2  g3 ;
         g3  g2 ;
         g3  g3;
         g1  g3
         g3  g1];
D = E/(1-mu^2)*[1 mu 0;...
					mu 1 0;...
					0 0 0.5*(1-mu)];
for i_el = 1:count_elem
     crd_elem = zeros(4,2); % коорд. узлов элемента
             for k =1:4
                crd_elem(k, :) = [x(i_el,k) , y(i_el,k)];
             end
      u1 = IND (i_el,1); 
      u2 = IND (i_el,2);
      u3 = IND (i_el,3);
      u4 = IND (i_el,4);

      Ue = U([2*u1-1; 2*u1;...
          2*u2-1; 2*u2;...
          2*u3-1; 2*u3;...
          2*u4-1; 2*u4]);
      
       for j = 1:9 % по точкам интегрирования
%              J = Jacobian(crd_elem(i_el,1),crd_elem(i_el,2),t_int(j,1),t_int(j,2));% Якобиан  
%              J_inv = inv(J);
             B = Bmatrix(t_int(j,1),t_int(j,2),crd_elem(:,1)',crd_elem(:,2)'); 
             strain_g(i_el,:,j) = strain_g(i_el,:,j) + (B*Ue)';
             stress_g(i_el,:,j) = stress_g(i_el,:,j) + (D*strain_g(i_el,:,j)')';
       end
       %инварианты
       sigma_x(i_el) = sum(stress_g(i_el,1,:))/9;
       sigma_y(i_el) = sum(stress_g(i_el,2,:))/9;
       sigma_z(i_el) = 0;
       tau_xy(i_el)  = sum(stress_g(i_el,3,:))/9;
       tau_yz(i_el)  = 0;
       tau_xz(i_el)  = 0;
       %Вычисление главных напряжений
%        I1 = sigma_x  + sigma_y  +sigma_z ;
%        I2 = sigma_x *sigma_y + sigma_y * sigma_z + sigma_z*sigma_x - tau_xy^2 - tau_yz^2 - tau_xz^2;
%        I3 = sigma_x*sigma_y*sigma_z + 2*tau_xy*tau_yz*tau_xz + sigma_x*tau_yz^2-sigma_y*tau_xz^2 - sigma_z*tau_xy^2;
%        s = roots ([1 -I1 I2 -I3]);
%        principal = sort(s);
%        stress_VM(i_el) = (sqrt(2)/2) * sqrt ((principal(3)-principal(2))^2+(principal(3)-principal(1))^2+(principal(2)-principal(1))^2);
stress_VM(i_el) = sqrt (sigma_x(i_el)^2+sigma_y(i_el)^2+sigma_z(i_el)^2-sigma_x(i_el)*sigma_y(i_el)-sigma_y(i_el)*sigma_z(i_el)-sigma_z(i_el)*sigma_x(i_el)+3*tau_xy(i_el)^2+3*tau_yz(i_el)^2+3*tau_xz(i_el)^2);
end
it=1;
 save_stress_to_file (stress_VM,count_elem,it);
end