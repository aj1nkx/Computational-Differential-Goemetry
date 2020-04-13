% Copyright 2010 Anand A. Joshi, David W. Shattuck and Richard M. Leahy 
% This file is part SVREG.
% 
% SVREG is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SVREG is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with BSE.  If not, see <http://www.gnu.org/licenses/>.

function [W] = surface_laplacian(sparseFV)
%Implementation of Loreta on a tesselation

%current density weighting
dim = size(sparseFV.vertices,1);
W = sparse(dim,dim);

for i=1:size(sparseFV.faces,1) %for all triangles
  %  progress(i,1,size(sparseFV.faces,1),5000,'triangles');
    
    %Get points from face
     trVertices = sparseFV.faces(i,:);
     points = sparseFV.vertices(trVertices,:);
    %Get new coord system
    r1 = points(2,:)-points(1,:);
    r = points(3,:)-points(1,:);
    projr_r1 = (r1*r')/rownorm(r1);
    r2 = sqrt(r*r' - abs(projr_r1)^2);
    if(r1*r')<0
        r2 = -r2;
    end
    pointsNew = [0 0; rownorm(r1) 0; projr_r1 r2];
    %area of triangle
    Area = rownorm(r1) * rownorm(r2) / 2;
    
    %find local ps functions
    A = [pointsNew [1 1 1]'];
    Ai = A^-1;
    
    %create matrix for local laplacian
    DPs = Ai(1:2,:);
    M = Area * DPs'*DPs;
    
    %update sparse matrix for global laplacian
    W(trVertices(1),trVertices(1)) = W(trVertices(1),trVertices(1)) + M(1,1);
    W(trVertices(1),trVertices(2)) = W(trVertices(1),trVertices(2)) + M(1,2);
    W(trVertices(1),trVertices(3)) = W(trVertices(1),trVertices(3)) + M(1,3);
    W(trVertices(2),trVertices(1)) = W(trVertices(2),trVertices(1)) + M(2,1);
    W(trVertices(2),trVertices(2)) = W(trVertices(2),trVertices(2)) + M(2,2);
    W(trVertices(2),trVertices(3)) = W(trVertices(2),trVertices(3)) + M(2,3);
    W(trVertices(3),trVertices(1)) = W(trVertices(3),trVertices(1)) + M(3,1);
    W(trVertices(3),trVertices(2)) = W(trVertices(3),trVertices(2)) + M(3,2);
    W(trVertices(3),trVertices(3)) = W(trVertices(3),trVertices(3)) + M(3,3);

%     
%     r1 = points(1,:);
%     r2 = points(2,:);
%     r3 = points(3,:);
%     a = r2 - r1;
%     b = r3 - r1;
%     ex = a/norm(a);
%     ey = b-dot(a,b)/dot(a,a)*a;
%     ey = ey/norm(ey);
%     xx = [0;dot(a,ex);dot(b,ex)];
%     yy = [0;dot(a,ey);dot(b,ey)];
%     npe = 3;
%     K = [ones(npe,1) xx yy];
%     E = inv(K);
%     Area = dot(a,ex)*dot(b,ey)/2;
%     gradE = E(2:end,1:end);
% 
%     S_loc = gradE'*gradE*Area;
%     
    
    
    
end
    

%disp('Condition Number:');
%disp(condest(W));
    
    
  
    
    







function nrm = rownorm(A);
%ROWNORM - Calculate the L2 norm of each row of a matrix
% function nrm = rownorm(A);
% calculate the Euclidean norm of each ROW in A, return as
% a column vector with same number of rows as A
%
% See also COLNORM

%<autobegin> ---------------------- 27-Jun-2005 10:45:37 -----------------------
% ------ Automatically Generated Comments Block Using AUTO_COMMENTS_PRE7 -------
%
% CATEGORY: Utility - Numeric
%
% At Check-in: $Author: Mosher $  $Revision: 16 $  $Date: 6/27/05 9:00a $
%
% This software is part of BrainStorm Toolbox Version 27-June-2005  
% 
% Principal Investigators and Developers:
% ** Richard M. Leahy, PhD, Signal & Image Processing Institute,
%    University of Southern California, Los Angeles, CA
% ** John C. Mosher, PhD, Biophysics Group,
%    Los Alamos National Laboratory, Los Alamos, NM
% ** Sylvain Baillet, PhD, Cognitive Neuroscience & Brain Imaging Laboratory,
%    CNRS, Hopital de la Salpetriere, Paris, France
% 
% See BrainStorm website at http://neuroimage.usc.edu for further information.
% 
% Copyright (c) 2005 BrainStorm by the University of Southern California
% This software distributed  under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html .
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%<autoend> ------------------------ 27-Jun-2005 10:45:37 -----------------------

% ----------------------------- Script History ---------------------------------
% 1994 by John C. Mosher
% May 6, 1994 JCM author
% 19-May-2004 JCM Comments Cleaning
% ----------------------------- Script History ---------------------------------

[m,n] = size(A);

if(n>1),            % multiple columns
  nrm = sqrt(sum([A.*conj(A)],2));
else                % A is col vector
  nrm = abs(A);
end

return












