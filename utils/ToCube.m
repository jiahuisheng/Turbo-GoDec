function v = ToCube(im,row,col)
% takes MxNx3 picture and returns (MN)x3 vector
sz = size(im);
for i=1:sz(1)
v(:,:,i) = reshape(im(i,:), [row col]);
end