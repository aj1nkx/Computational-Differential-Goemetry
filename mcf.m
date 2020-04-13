clear;
[vert,face] = read_ply('./data/bun_zipper.ply');
FV.vertices = vert;
FV.faces = face;
view_patch(FV);
n_iter = 20;
for i=1:n_iter
    step = 1e-5/i;
    % Find the patchcurvature and create a new mesh
    [Cmean,N,~,~,~,~,~] = patchcurvature(FV,false);
    FV.vertices = FV.vertices + step*repmat(Cmean,1,3).*N;
end
view_patch(FV);