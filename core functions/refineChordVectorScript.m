% refine chord script

[newChordR, newAltChordR, newDiag1R, newDiag2R] = refineChordVector(data, 'R');
data.rightChordHats = newChordR ;
data.chord1AltHats  = newAltChordR ;
data.diag11Right    = newDiag1R ;
data.diag12Right    = newDiag2R ;

% keyboard ;

[newChordL, newAltChordL, newDiag1L, newDiag2L] = refineChordVector(data, 'L');
data.leftChordHats = newChordL ;
data.chord2AltHats = newAltChordL ;
data.diag21Left    = newDiag1L ;
data.diag22Left    = newDiag2L ;

clear newChordR newAltChordR newDiag1R newDiag2R
clear newChordL newAltChordL newDiag1L newDiag2L
