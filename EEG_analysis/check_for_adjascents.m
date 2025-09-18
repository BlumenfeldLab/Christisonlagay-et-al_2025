function [rows_that_pass full_row_entry]=check_for_adjascents(all_faces, current_face, sharedfaces)


checkmatches1=all_faces==current_face(1);
checkmatches2=all_faces==current_face(2);
checkmatches3=all_faces==current_face(3);

rows_that_pass=find(ge(sum((checkmatches1+checkmatches2+checkmatches3),2),sharedfaces));
full_row_entry=all_faces(rows_that_pass,:);
