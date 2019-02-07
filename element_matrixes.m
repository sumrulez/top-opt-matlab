function K = element_matrixes(E,mu,h,x,y)
global count_elem
        el_nodes = 4;
        K = zeros (2*el_nodes,2*el_nodes,count_elem);
        for i_el = 1:count_elem % loop for elements  
             K(:,:,i_el) = h*K_el_9 (E,mu,x(i_el,:),y(i_el,:));
        end
        
        end