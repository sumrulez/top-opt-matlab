function neighbor_el = neighbor(rmin,design_space,centr)
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
end