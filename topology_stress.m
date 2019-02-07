function topology_stress(flag)
%ѕараметры модели
E=2.1e5;  mu=0.3;
h = 10;%толщина

%ћинимизаци€ податливости
if flag == 1
    volfrac = 0.5;
    penal  = 3;
    rmin = 5;
end

%ћинимизаци€ массы
if flag == 2 
   penal = 3;
   rmin = 5;
   clasters_num = 20;
end
global count_elem count_spc count_force count_node
[SPC,FORCE,design_space,non_design,count_elem,count_spc,count_force,count_node,IND,x,y,z,M_el] = import_data;
centr = centri (x,y);%%% центроиды элементов
neighbor_el = neighbor(rmin,design_space,centr); %определение соседних элементов
K = element_matrixes(E,mu,h,x,y); %матрицы жесткостей элементов



%%%%%%%

    %инициализаци€ плотностей
    if flag == 1
    dens(1:count_elem) = volfrac; 
    dens (non_design) = 1;
    end
    if flag == 2
        dens(1:count_elem) = 0.3;
    end
        
loop = 0; 
change = 1.;

%% ћинимизаци€ податливости с ограничением по объему
if flag == 1 
% START ITERATION
clear M_el z
while change > 0.005  
  loop = loop + 1;
  if loop == 200
      disp ('Limit on iterations reached')
      break;
  end
  dens_old = dens;
% FE-ANALYSIS
  [U,max_disp]=FEsolve_stress(E,mu,dens,penal,IND,SPC,FORCE,K,x,y,0);         
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


%% ћинимизаци€ массы с ограничением по напр€жени€м
if flag == 2
    clear z
    D = E/(1-mu^2)*[1 mu 0;mu 1 0;0 0 0.5*(1-mu)]; %матрица упругих констант
while change > 0.005  
  loop = loop + 1;
  dens_old = dens;
% FE-ANALYSIS
  [U,max_disp,stress_VM,sigma_x,sigma_y,tau_xy]=FEsolve_stress(E,mu,dens,penal,IND,SPC,FORCE,K,x,y,1);
  dens = dens_filter (rmin,dens,centr,design_space,neighbor_el);  %фильтр плотностей
  stress_VM = stress_VM .* dens; %пенализаци€ напр€жений
  df0dx = obj_fun_deriv (rmin,centr,design_space,neighbor_el,M_el); %производные целевой функции(массы) к плотност€м
  [Clasters,sigma_PN] = stress_to_clasters_distr (clasters_num,stress_VM);   %распределение по кластерам
  %PN stress derivatives
  lambda = 0;
  dfdx = zeros (clasters_num,count_elem);
%   m    = The number of general constraints.
%   n    = The number of variables x_j.
% dfdx  = (m x n)-matrix with the derivatives of the constraint functions
%          f_i with respect to the variables x_j, calculated at xval.
%          dfdx(i,j) = the derivative of f_i with respect to x_j.
  for i = 1:clasters_num
      cl_el = nnz(Clasters(:,i));
      for c = 1:length(design_space)
          j= design_space(c);
          D3 = 0.5*dens(j)^(-0.5)*[sigma_x(j); sigma_y(j); tau_xy(j)];
          invmat = inv(dens(j)*K(:,:,j));
          for a = 1:cl_el
%               for b = 1:cl_el
%                   D1 = (sigma_PN(i)^penal)^(1/penal-1)*1/cl_el*stress_VM(b); %производна€ sigma_PN(i) по sigma_VM(a)
%                   D2(1) = 1/(2*stress_VM(b))*(2*sigma_x(b)-sigma_y(b)); % производна€ sigma_VM (a) по sigma_x(a)
%                   D2(2) = 1/(2*stress_VM(b))*(2*sigma_y(b)-sigma_x(b));% производна€ sigma_VM (a) по sigma_y(a)
%                   D2(3) = 3/stress_VM(b)*tau_xy(b);% производна€ sigma_VM (a) по tau_xy(a)
                  
%               end
                b = Clasters(a,i);
              D1 = (sigma_PN(i))^(1-penal)*1/cl_el*stress_VM(b)^(penal-1); %производна€ sigma_PN(i) по sigma_VM(b)
              D2(1) = 1/(2*stress_VM(b))*(2*sigma_x(b)-sigma_y(b)); % производна€ sigma_VM (b) по sigma_x(b)
              D2(2) = 1/(2*stress_VM(b))*(2*sigma_y(b)-sigma_x(b));% производна€ sigma_VM (b) по sigma_y(b)
              D2(3) = 3/stress_VM(b)*tau_xy(b);% производна€ sigma_VM (b) по tau_xy(b)
              lambda = D1*D2*D*Bmatrix(0,0,x(b,:),y(b,:))*invmat;
              
              D5 = sqrt(dens(b))*lambda;
              D6 = 0;
              for r = 1:length(design_space)
                  rr = design_space(r);
                fac = rmin-abs(sqrt((centr(j,1)-centr(rr,1)))^2+(centr(j,2)-centr(rr,2))^2);
                  fac = 1;
                  u1 = IND (j,1); 
                  u2 = IND (j,2);
                  u3 = IND (j,3);
                  u4 = IND (j,4);
                  Ue = U([2*u1-1; 2*u1; 2*u2-1; 2*u2; 2*u3-1; 2*u3; 2*u4-1; 2*u4]);
                  D6 = D6 + (penal*dens(b)^(penal-1)*K(:,:,rr))*(fac/rmin)*Ue;
              end
              dfdx (i,j) = dfdx (i,j)+  D1*D2*D3 - sqrt(dens(b))*lambda*D6; % i-constraint j - desvar
          end
      end
  end
% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    end
end
end
%% –аспределение напр€жений по кластерам
%ћетод Stress level
function [Clasters,sigma_PN] = stress_to_clasters (clasters_num,stress_VM)
global count_elem
n = clasters_num;
  max_stress = max (stress_VM);
  step = max_stress/n;
  low_lim = 0;
  up_lim = step;
  j=1; u=1; 
  p = 8; %эвристический показатель из статьи
  sigma_PN = zeros (1,n);
  while up_lim <= (max_stress+step/10)
      for i = 1:count_elem 
         if stress_VM(i) >= low_lim &&  stress_VM(i) <= up_lim
             Clasters (u,j) = i;
             sigma_PN (j) = sigma_PN (j) + (stress_VM (i))^p;
             u = u+1;
         end
      end
      sigma_PN (j) = ((1/(u-1))*sigma_PN(j))^(1/p);
      j = j+1;
      u=1;
      low_lim = up_lim;
      up_lim = up_lim + step;
  end
end
  
%% –аспределение напр€жений по кластерам
%ћетод Distributed stress
function [Clasters,sigma_PN] = stress_to_clasters_distr (clasters_num,stress_VM)
global count_elem
n = clasters_num;
%   max_stress = max (stress_VM);
%   step = max_stress/n;
%   low_lim = 0;
%   up_lim = step;
  j=1; u=1; 
  p = 8; %эвристический показатель из статьи
  sigma_PN = zeros (1,n);
  [B, I] = sort (stress_VM); %сортировка напр€жений. B -  напр€жени€, I - номера элементов
for i = 1:count_elem 
      Clasters (u,j) = I(i);
      sigma_PN (j) = sigma_PN (j) + (B(i))^p;
      j = j+1;
      if j>n 
          j=1;
          u = u+1;
      end
end
sigma_PN  = ((1/(u-1)).*sigma_PN).^(1/p);
end

%% ‘ильтр плотностей дл€ flag=2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dens_filtered] = dens_filter (rmin,dens,centr,design_space,neighbor_el)
global count_elem
 dens_filtered=dens;
for i = 1:size(design_space,1)
    sum=0.0; 
    d = design_space(i);
    for k = 1:size(neighbor_el,2)
        %нужны центроиды элементов
        if neighbor_el(d,k) ~= 0
            %рассто€ни€ вывести из итерационного расчета
            fac = rmin-abs(sqrt((centr(d,1)-centr(neighbor_el(d,k),1))^2+(centr(neighbor_el(d,k),2)-centr(d,2))^2));
            sum = sum+max(0,fac)/rmin;
            dens_filtered(d) =  dens_filtered(d) + fac/rmin*dens(neighbor_el(d,k));
        end
    end
   if dens_filtered(d) >= 1 
       dens_filtered(d) = 0.99;
   end
end
end
function [df0dx] = obj_fun_deriv (rmin,centr,design_space,neighbor_el,M_el)
global count_elem
df0dx = zeros(1,count_elem);
for i = 1:size(design_space,1)
    sum=0.0; 
    for k = 1:size(neighbor_el,2)
        %нужны центроиды элементов
        if neighbor_el(design_space(i),k) ~= 0
            %рассто€ни€ вывести из итерационного расчета
            fac = rmin-abs(sqrt((centr(design_space(i),1)-centr(neighbor_el(design_space(i),k),1))^2+(centr(neighbor_el(design_space(i),k),2)-centr(design_space(i),2))^2));
            df0dx(design_space(i)) =   df0dx(design_space(i)) + M_el(neighbor_el(i))*(rmin-max(0,fac))/rmin;
        end
    end
end
end

%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(count_elem,rmin,dens,dc,centr,neighbor_el,design_space)
dcn=zeros(1,count_elem);
for i = 1:size(design_space,1)
    sum=0.0; 
    for k = 1:size(neighbor_el,2)
        %нужны центроиды элементов
        if neighbor_el(design_space(i),k) ~= 0
            %рассто€ни€ вывести из итерационного расчета
            fac = rmin-abs(sqrt((centr(design_space(i),1)-centr(neighbor_el(design_space(i),k),1))^2+(centr(neighbor_el(design_space(i),k),2)-centr(design_space(i),2))^2));
            sum = sum+max(0,fac);
            dcn(design_space(i)) = dcn(design_space(i)) + max(0,fac)*dens(neighbor_el(design_space(i),k))*dc(neighbor_el(design_space(i),k));
        end
    end
    dcn(design_space(i)) = dcn(design_space(i))/(dens(design_space(i))*sum);
  end
end

