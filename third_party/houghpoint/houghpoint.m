function [H,angles,offsets] = houghpoint(XY,angles,offsets,MATRIXMODE)
% HOUGHPOINT  Hough transform of points
%  [TRANSFORM,ANGLES,OFFSETS] = HOUGHPOINT(XY,ANGLES,OFFSETS) computes the
%  Hough transform of the points located at coordinates XY where XY(:,1) is
%  a list of X coordinates and XY(:,2) is a list of Y coordinates and XY is
%  a Nx2 array.  It computes the transform at the given ANGLES and OFFSETS.
%  For each point and angle, contributions are split between the adjacent
%  two offset bins by linear interpolation.  Any transforms that result in
%  offsets beyond the limits of OFFSETS are discarded, which means that the
%  most extreme OFFSETS may be underrepresented if an insufficent range is
%  considered.
%  
%  To convert a line with a given ANGLE and OFFSET to slope-intercept form
%  Y=A*X+B, use A=tan(ANGLE) and B=OFFSET/cos(ANGLE).
%  
%  Alternatively, one may input an integer for ANGLES and/or OFFSETS, in
%  which case ANGLES are uniformly selected from (-pi/2,pi/2] with ANGLES
%  steps and OFFSETS from [-max(hypot(X,Y)),max(hypot(X,Y))] with OFFSETS
%  steps (which should span the possible range of offsets).
% 
%  [TRANSFORM,ANGLES,OFFSETS] = HOUGHPOINT(XYW,ANGLES,OFFSETS) is the same
%  as HOUGHPOINT(XY...), but instead XYW is an Nx3 array and the function
%  uses the (X,Y,W) triplets where W=XYW(:,3) determines the weighting of
%  each point.  The XY syntax effectively uses W=ones(size(XY,1),1).
%  
%  [H,ANGLES,OFFSETS] = HOUGHPOINT(XY,ANGLES,OFFSETS,MATRIXMODE) when
%  MATRIXMODE evaluates as TRUE computes the Hough transform matrix H such
%  that TRANSFORM==reshape(H*W,numel(ANGLES),numel(OFFSETS)).  This is
%  slower than TRANSFORM=HOUGHPOINT(XYW...), but if many different weights
%  share the same coordinates then it may be faster to precompute H and
%  apply it to each set of weights using the above matrix-vector (or
%  matrix-matrix for multiple sets of weights) multiply.
% 
%  See also HOUGHIMG.

switch size(XY,2)
	case 2
		MODEW = false; % no weights (assume all weights==1)
	case 3 % split up XY and W
		W = XY(:,3);
		XY = XY(:,1:2);
		MODEW = true; % we have weights to apply
	otherwise
		error('First input must have size Nx2 or Nx3')
end

if nargin < 4 || isempty(MATRIXMODE) % do we produce the transformation matrix? or just compute the transform
	MATRIXMODE = false;
else
	MATRIXMODE = logical(MATRIXMODE);
end

npts = size(XY,1);

if isscalar(angles)
	angles = ((1:angles)/angles-0.5)*pi;
end
maxoff = max(sqrt(sum(XY.*XY,2))); % maximum possible offset for any angle
if isscalar(offsets)
	OFFSETSPACING = 2*maxoff/(offsets-1);
	UNIFORMOFFSETS = true;
	offsets = linspace(-maxoff,maxoff,offsets)';
else
	if ~issorted(offsets)
		error('OFFSETS must be sorted') % or else the 'oind' calculation will be wrong or error
	end
	OFFSETSPACING = (offsets(end)-offsets(1))/(numel(offsets)-1);
	UNIFORMOFFSETS = all(abs(diff(offsets)/OFFSETSPACING - 1) < eps(128)); % OFFSETS is a uniform vector
end

angles = angles(:);
offsets = offsets(:);
nangles = numel(angles);
noffsets = numel(offsets);

ptoffset = nan(npts,nangles);
for ang = 1:nangles % calculate offset for each point-angle pair
% 	avec = ones(npts,1)*[cos(angles(ang)),sin(angles(ang))];
% 	rot = ones(npts,1)*[-sin(angles(ang)),cos(angles(ang))];
% 	ptoffset(:,ang) = sum(rot.*(XY - ((XY.*avec)*ones(2)).*avec),2);
	avec = [cos(angles(ang)),sin(angles(ang))];
	rot = [-sin(angles(ang)),cos(angles(ang))];
	ptoffset(:,ang) = sum(bsxfun(@times,rot,bsxfun(@minus,XY,bsxfun(@times,repmat(sum(bsxfun(@times,XY,avec),2),[1,2]),avec))),2);
end
% avec = permute([cos(angles),sin(angles)],[3,2,1]); % fully vectorized but not as fast
% rot = permute([-sin(angles),cos(angles)],[3,2,1]); % fully vectorized but not as fast
% ptoffset = permute(sum(bsxfun(@times,rot,bsxfun(@minus,XY,bsxfun(@times,repmat(sum(bsxfun(@times,XY,avec),2),[1,2]),avec))),2),[1,3,2]); % fully vectorized but not as fast

if maxoff <= offsets(end) && maxoff <= -offsets(1)
	inrange = reshape(1:numel(ptoffset),[],1); % impossible to produce out-of-range offsets
	ptoffset = reshape(ptoffset,[],1);
else
	inrange = reshape(find(ptoffset>=offsets(1) & ptoffset<=offsets(end)),[],1); % ignore out-of-range offsets
	ptoffset = reshape(ptoffset(inrange),[],1);
end

if UNIFORMOFFSETS
	oind = floor((ptoffset+(OFFSETSPACING-offsets(1)))/OFFSETSPACING); % fast
	oind = min(oind,noffsets-1); % prevent exact matches to upper bound from causing errors
elseif exist('discretize','file')~=0
	oind = discretize(ptoffset,offsets); % fast but only in R2015a and later
else
	[~,oind] = histc(ptoffset,offsets); % slower
	oind = min(oind,noffsets-1); % prevent exact matches to upper bound from causing errors
% 	oind = sum(bsxfun(@le,offsets',ptoffset),2); % slowest and most memory
end
% ADDED BY HEINARO 30.8.2019
oind(oind == 0) = 1;

low = offsets(oind); % get offset below
high = offsets(oind+1); % get offset above
lowv = (high-ptoffset)./(high-low); % find weight for lower anchor
clearvars low high ptoffset % clear up some memory in case we're getting full

[pt,aind] = ind2sub([npts,nangles],inrange);
Hi = sub2ind([nangles,noffsets],aind,oind); % get output index for each point
clearvars aind oind inrange % clear up some memory in case we're getting full
if MODEW
	V = [W(pt).*lowv;W(pt).*(1-lowv)];
else
	V = [lowv;1-lowv];
end
if MATRIXMODE
	H = sparse([Hi;Hi+nangles],[pt;pt],V,nangles*noffsets,npts);
else
	H = reshape(accumarray([Hi;Hi+nangles],V,[nangles*noffsets,1]),nangles,noffsets);
end

end