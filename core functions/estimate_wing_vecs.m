%--------------------------------------------------------------------------
% Function to estimate the span and chord vectors for a wing. Based on the
% methods from hullAnalysis
%
%   INPUTS:
%       -wing_vox = voxel coordinates for the wing
%
%       -ref_vecs = reference vectors to determine direction of span and
%       chord for left vs right wing. It is assumed
%       that the reference vectors in this case are a list of [bodyCM(i),
%       bodyCM(i-1), wingTip(i-1)].
%
%       -wingLength = estimate of wing length. comes from data struct as
%           wingLength = 35 * params.pixPerCM / 232 ;
%
%       -span_ref (optional) = reference span against which to compare
%       -chord_ref (optional) = reference chord against which to compare
%
%   OUTPUTS:
%       -span = estimate of span vector
%       -chord = estimate of chord vector
%       -chord_alt = estimate of alternate chord vector
%       -N_vox = number of voxels in wing
%--------------------------------------------------------------------------
function [spanHat, chordHat, chordAltHat, Nvox] = ...
    estimate_wing_vecs(wing_vox, ref_vecs, wingLength, span_ref, chord_ref, wingCM)
%--------------------------------------------------------------------------
%% deal with function inputs
if ~exist('span_ref','var') 
    span_ref = [] ;
end
if ~exist('chord_ref','var')
    chord_ref = [] ;
end
if ~exist('wingCM','var')
    wingCM = mean(wing_vox) ; 
end
%--------------------------------------------------------------------------
%% params
LL = wingLength * 0.55 ; % used with farthestPoint
Nvox = size(wing_vox,1) ;
debugFlag = false ;

bodyCM = ref_vecs(1,:) ;
bodyCM_prev = ref_vecs(2,:) ;
wingTip_prev = ref_vecs(3,:) ;
%--------------------------------------------------------------------------
%% get span
if isempty(span_ref)
    % first pass: span = vector from body center of mass to distal wing tip
    farPoint    = farthestPoint(wing_vox, bodyCM, LL) ;
    spanHat = farPoint - wingCM ;
    spanHat = spanHat' ;
    spanHat = spanHat / norm(spanHat) ;
    
    % makes wing span point outward
    if dot(wingCM-bodyCM,spanHat) < 0
        spanHat = -spanHat;
    end
    
    wingTip = findWingTip(wing_vox, spanHat', wingCM);
    if (isnan(wingTip(1)))
        wingTip = farPoint ;
    end
    
    % second pass: recalculate span vector based on the refined wing tip
    spanHat = wingTip - wingCM ;
    spanHat = spanHat' ;
    spanHat = spanHat / norm(spanHat) ;
else
    pca_coeffs = pca(wing_vox);
    spanHat = pca_coeffs(:,1);
    
    if dot(spanHat,span_ref) <= 0 
        spanHat = -spanHat;
    end
    wingTip = findWingTip(wing_vox, spanHat', wingCM);
end
%----------------------------------------------------------------------
%% find chord vector
[chordHat,chordAltHat, ~, ~] = find_chords_quad(wing_vox, spanHat', ...
    wingTip, wingTip_prev, bodyCM, bodyCM_prev) ;

if ~isempty(chord_ref)
    if compare_vectors(chordHat, chord_ref) > ...
            compare_vectors(chordAltHat, chord_ref)
        tmp = chordHat ;
        chordHat = chordAltHat ;
        chordAltHat = tmp ;
    end
end
%----------------------------------------------------------------------
%% plot results?
if debugFlag
    figure ;
    hold on
    plot3(wing_vox(:,1), wing_vox(:,2), wing_vox(:,3),'k.')
    plot3((wingCM(1) + 24*[0, spanHat(1)]), (wingCM(2) + 24*[0 spanHat(2)]),...
        (wingCM(3) + 24*[0 spanHat(3)]),'b-','linewidth',4)
    plot3((wingCM(1) + 12*[0 chordHat(1)]), (wingCM(2) + 12*[0 chordHat(2)]),...
        (wingCM(3) + 12*[0 chordHat(3)]),'r-','linewidth',4)
    plot3((wingCM(1) + 12*[0 chordAltHat(1)]), (wingCM(2) + 12*[0 chordAltHat(2)]),...
        (wingCM(3) + 12*[0 chordAltHat(3)]),'r:','linewidth',4)
    legend({'Data','Span','Chord','Alt. Chord'})
    axis equal
    grid on
    
end


end