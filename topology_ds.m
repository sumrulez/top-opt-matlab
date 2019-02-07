%%%% A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND, JANUARY 2000 %%%
%%%% CODE MODIFIED FOR INCREASED SPEED, September 2002, BY OLE SIGMUND %%%
function topology(volfrac,penal,rmin)
% INITIALIZEunt_spc count_force count_node

%%%%%%%
    	E=70000;  mu=0.3;
global count_elem count_node
        R = importdata('mesh.csv',','); R = R.data;
        SPC = importdata('spc.csv',','); SPC = SPC.data;
        FORCE = importdata('force.csv',','); FORCE = FORCE.data;  
        %масштабирование сил
%         FORCE (:,3) = 1e6 *FORCE (:,3);
        design_space = importdata('design_space.csv'); design_space = design_space.data;
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
  
        %массив номеров элементов и узлов 
        %номер строки = номер элемента
        % 4 столбца - узлы данного элемента
        IND = R(:,[2:5]);
        %%% центроиды элементов
        for i=1:count_elem
            centr (i,1) = 0.25*(x(i,1)+x(i,4)+x(i,2)+x(i,3)); 
            centr (i,2) = 0.25*(y(i,1)+y(i,2)+y(i,3)+y(i,4));
        end
        %%%определение соседних элементов
        % neighbor_el: строка - номер элемента, столбцы -соседние элементы
        for i = 1:size(design_space,1)
            j=1;
            for k = 1:size(design_space,1)
                fac = rmin-abs(sqrt((centr(design_space(i),1)-centr(design_space(k),1))^2+(centr(design_space(k),2)-centr(design_space(i),2))^2));
                if fac > 0 && design_space(k) ~= design_space(i)
                    neighbor_el (design_space(i),j) = design_space(k);
                    j=j+1;
                end
            end
        end
%         neighbor_el = sparse(neighbor_el);
        el_nodes = 4;
        h = 20; %толщина пластины в мм
        K = zeros (2*el_nodes,2*el_nodes,count_elem);
        for i_el = 1:count_elem % loop for elements  
             K(:,:,i_el) = h*K_el_9 (E,mu,x(i_el,:),y(i_el,:));
        end
    %инициализация плотностей
    dens(1:count_elem) = volfrac; 
    dens (non_design) = 1;
loop = 0; 
change = 1.;

%%
% START ITERATION
while change > 0.005  
  loop = loop + 1;
  dens_old = dens;
% FE-ANALYSIS
  [U,max_disp]=FEsolve(dens,penal,IND,SPC,FORCE,K);         
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  c = 0.;
  for i_el = 1:count_elem
      u1 = IND (i_el,1); 
      u2 = IND (i_el,2);
      u3 = IND (i_el,3);
      u4 = IND (i_el,4);
      KE = K (:,:,i_el);
      Ue = U([2*u1-1; 2*u1; 2*u2-1; 2*u2; 2*u3-1; 2*u3; 2*u4-1; 2*u4]);
      c = c + dens(i_el)^penal*Ue'*KE*Ue;
      dc(i_el) = -penal*dens(i_el)^(penal-1)*Ue'*KE*Ue;
  end
% FILTERING OF SENSITIVITIES
  [dc]   = check(count_elem,rmin,dens,dc,centr,neighbor_el,design_space);    
% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
  [dens(design_space)]    = OC(size(design_space,1),dens(design_space),volfrac,dc(design_space)); 
  save__dens_to_file (dens,count_elem,loop);
% PRINT RESULTS
  change = max(max(abs(dens-dens_old)));
  disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(dens(design_space))/(size(design_space,1)))...
       ' ch.: ' sprintf('%6.3f',change) ' maxdisp: ' sprintf('%6.3f',max_disp )])
% PLOT DENSITIES  
%   colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6);
end 
end
%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=OC(count_elem,dens,volfrac,dc)  
l1 = 0; l2 = 100000; move = 0.2;
while (l2-l1 > 1e-4)
  lmid = 0.5*(l2+l1);
  xnew = max(0.001,max(dens-move,min(1.,min(dens+move,dens.*abs(sqrt(-dc./lmid))))));
  if sum(sum(xnew)) - volfrac*count_elem > 0;
    l1 = lmid;
  else
    l2 = lmid;
  end
end
end
%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(count_elem,rmin,dens,dc,centr,neighbor_el,design_space)
dcn=zeros(1,count_elem);
for i = 1:size(design_space,1)
    sum=0.0; 
    for k = 1:size(neighbor_el,2)
        %нужны центроиды элементов
        if neighbor_el(design_space(i),k) ~= 0
            %расстояния вывести из итерационного расчета
            fac = rmin-abs(sqrt((centr(design_space(i),1)-centr(neighbor_el(design_space(i),k),1))^2+(centr(neighbor_el(design_space(i),k),2)-centr(design_space(i),2))^2));
            sum = sum+max(0,fac);
            dcn(design_space(i)) = dcn(design_space(i)) + max(0,fac)*dens(neighbor_el(design_space(i),k))*dc(neighbor_el(design_space(i),k));
        end
    end
    dcn(design_space(i)) = dcn(design_space(i))/(dens(design_space(i))*sum);
  end
end
%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by Ole Sigmund, Department of Solid         %
% Mechanics, Technical University of Denmark, DK-2800 Lyngby, Denmark.     %
% Please sent your comments to the author: sigmund@fam.dtu.dk              %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "A 99 line topology optimization code written in Matlab"                 %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The author reserves all rights but does not guaranty that the code is    %
% free from errors. Furthermore, he shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
