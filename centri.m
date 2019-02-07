function centr = centri(x,y)
global count_elem 
for i=1:count_elem
            centr (i,1) = 0.25*(x(i,1)+x(i,4)+x(i,2)+x(i,3)); 
            centr (i,2) = 0.25*(y(i,1)+y(i,2)+y(i,3)+y(i,4));
        end
end