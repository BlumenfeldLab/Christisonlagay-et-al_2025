rights=[]; lefts=[];
for a=1:10
    
    
    load (['significant_vertices_map' num2str(a) '.mat'])
    
    rights=cat(3, rights, significantvertex_valuesR);
    lefts=cat(3, lefts, significantvertex_valuesL);

end


fraction_significant_right=sum(rights,3)/size(rights,3);
sign_right=ge(fraction_significant_right,.5);

fraction_significant_left=sum(lefts,3)/size(lefts,3);
sign_left=ge(fraction_significant_left,.5);