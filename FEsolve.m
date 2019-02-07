function [U,max_disp] = FEsolve (dens,penal,IND,SPC,FORCE,Kloc)
% Function solve plain stress finite element problem and returns displacements and maximum displacement 

%äàííûå ìîäåëè
global count_elem count_spc count_force count_node    
%ñîñòàâëåíèå ãëîáàëüíîé ìàòðèöû æåñòêîñòè
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

%çàêðåïëåíèÿ
fixeddof = zeros(count_spc,1);
n=1;
for i_spc = 1:count_spc
    spc_node = SPC(i_spc,2);
        if SPC(i_spc,3) == 123456
          fixeddof (n:1:n+1) = [2*spc_node-1 2*spc_node];
            n=n+2;
        end
		if SPC(i_spc,3) == 23456
          fixeddof (n) = [2*spc_node];
		end
end
fixeddof = sort (fixeddof);
alldof = [1:1:count_node *2]';
freedof = setdiff (alldof,fixeddof);
clear alldof fixeddof
%ñèëû
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
%âû÷èùåíèå ïàìÿòè ïåðåä ðåøåíèåì ÑËÀÓ
clear K1 K F K_glob
U(freedof) = Kglob_sp(freedof,freedof)\F_sp(freedof);
max_disp = max(U);
save_to_file (U,count_node);
