function [normalf, FaceToElem_list] = ComputeSurfaceNormals3D(boundaries, vertices, elements)

if nargout == 1    
    normalf  = faceNormal(vertices', boundaries', elements');    
elseif nargout == 2   
    [normalf, FaceToElem_list]  = faceNormal(vertices', boundaries', elements');  
end

normalf  = normalf';
end

%% Compute the normal vectors of the boundary surfaces 
function [normals, id_e] = faceNormal(nodes, faces, elements)

% compute vector of first edges
v1 = nodes(faces(:,2),1:3) - nodes(faces(:,1),1:3);
v2 = nodes(faces(:,3),1:3) - nodes(faces(:,1),1:3);

% compute normals using cross product (nodes have same size)
normals = normalizeVector3d(cross(v1, v2, 2));

message = sprintf(['\n\n[\bWarning]\b: ComputeSurfaceNormals3D assumes that boundary elements have correct orientation,\n', ...
                  'so that surface normals are directed outward. No check on the sign is made.\n\n']);
fprintf(message);

%% check correct orientation (outward)
% Uncomment the following lines if your mesh does not have outward normals
%
% faces    = faces';
% elements = elements';
% nodes = nodes';
% elements_row = [elements(1,:) elements(2,:) elements(3,:) elements(4,:)];
% 
% id_e = zeros(size(faces,2),1);
% parfor i = 1 : size(faces, 2)
%     id_e(i)        = FaceToElement(faces(:,i), elements, elements_row);
%     centroid_el    = mean(nodes(:,elements(:,id_e(i))),2) ;
%     centroid_face  = mean(nodes(:,faces(:,i)),2) ;
%     sign_normalf   = normals(i,:)*(centroid_face - centroid_el);
%     if sign_normalf < 0
%         normals(i,:) = - normals(i,:);
%     end
% end

%% compute FaceToElem_list
if nargout > 1
    faces    = faces';
    elements = elements';
    
    %orginal node ids of P1 prism 
    elements_row = [elements(1,:) elements(2,:) elements(3,:) elements(4,:)];
    
    % find the elements to which a boundary face is attached
    id_e = zeros(size(faces,2),1);
    parfor i = 1 : size(faces, 2)
        id_e(i)        = FaceToElement(faces(:,i), elements, elements_row);
    end
end
    
end

function id_e = FaceToElement(face, elements,   elements_row)
% find the element to which a face belongs
i1 = face(1);
i2 = face(2);
i3 = face(3);

noe = size(elements,2);

% using the list of boundary elements search the global elements
ie = find(elements_row == i1);
for kk = 1 : length(ie)

    ie_tmp = ie(kk);
    if mod(ie_tmp,noe) == 0
        ie_tmp = noe;   %ie_tmp/noe;
    else
        ie_tmp = mod(ie_tmp,noe);
    end
    
    tmp = setdiff(elements(1:4,ie_tmp),[i1 i2 i3]);
    if length(tmp) == 1
        id_e  = ie_tmp;
        break;
    end
end

end

function vn = normalizeVector3d(v)

vn   = bsxfun(@rdivide, v, sqrt(sum(v.^2, 2)));

end
