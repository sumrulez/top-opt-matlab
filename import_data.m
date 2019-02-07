function [SPC,FORCE,design_space,non_design,count_elem,count_spc,count_force,count_node,IND,x,y,z,M_el] = import_data

R = importdata('mesh.csv',','); R = R.data;
        SPC = importdata('spc.csv',','); SPC = SPC.data;
        FORCE = importdata('force.csv',','); FORCE = FORCE.data;  
%         масштабирование сил
%         FORCE (:,3) = 1e6 *FORCE (:,3);
        design_space = importdata('design_space.csv'); design_space = design_space.data; %R(:,1); 
        non_design = setdiff(R(:,1),design_space);
        count_elem = size (R,1);
        count_node  = max(max (R(:,[2,3,4,5])));
        count_spc = size (SPC,1);
        count_force = size (FORCE,1);
        %координаты узлов 
        %строка - номер элемента
        x(:,:)= [R(:,6),R(:,9),R(:,12),R(:,15)];
        y(:,:)= [R(:,7),R(:,10),R(:,13),R(:,16)];
        z(:,:)= [R(:,8),R(:,11),R(:,14),R(:,17)];
        
        M_el (1:count_elem) = 1; % массы элементов
        %массив номеров элементов и узлов 
        %номер строки = номер элемента
        % 4 столбца - узлы данного элемента
        IND = R(:,[2:5]);
end