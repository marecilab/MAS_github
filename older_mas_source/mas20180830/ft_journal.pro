;; $id$
;;

; ------------------------------------------------------------------------------
;
; ft_journal.pro
;
; The FT_JOURNAL procedure stores the parameters of tracking and information
; about selected ROI planes of the data in a text file with the same name as
; the tracts file but with an extension ".info".  This file can be later viewed
; from the UNIX prompt with the command 'more [name].info'.
;
; SYNTAX
;	ft_journal, parameters, thresholdfile, vectorfile, ptr_planes
;
; ARGUMENTS
;  parameters
;	A structure carrying the information of the file path, name of the file
;	that the calculated tracts will be stored in, resolution factors in 3
;	dimensions, step size for tracking, threshold and length limits, and
;	density.
;
;  thresholdfile
;	Name of the file that will be used in the termination of tracking.
;
;  vectorfile
;	Name of the file that includes the unit vectors assumed to be tangent
;	to the fibers.
;
;  ptr_planes
;	The pointer pointing to the information about the normal vectors and
;	center coordinates of the selected planes.
;
;  Created on October 18, 2001 by E. Ozarslan
;
; ------------------------------------------------------------------------------


pro ft_journal, b, tfile, dfile, ptr_planes

  rfs='_rf'+strtrim(string(b.rfx),1)+strtrim(string(b.rfy),1)+$
									strtrim(string(b.rfz),1)
  sl=1/sqrt(b.sd)
  journal, b.fname+rfs+'.info'
  print, 'ROOT				:	'+b.path
  print, 'OUTPUT FILE			:	'+b.fname+rfs+'.flt'
  print, 'DIRECTIONS FILE		:	'+dfile+'.flt'
  print, 'THRESHOLD FILE			:	'+tfile+'.flt'
  print, 'THRESHOLD VALUE		:	'+strtrim(string(b.thr),1)
  print, 'STEP SIZE			:	'+strtrim(string(b.ss),1)
  print, 'LENGTH LIMIT			:	'+strtrim(string(b.ll),1)+' voxels'
  print, 'SEED DENSITY			:	'+strtrim(string(b.sd),1)+' per voxels'
  print, 'DISTANCE BETWEEN SEEDS		:	'+strtrim(string(sl),1)+' voxels'
  print, 'RESOLUTION FACTOR (X)		:	'+strtrim(string(b.rfx),1)
  print, 'RESOLUTION FACTOR (Y)		:	'+strtrim(string(b.rfy),1)
  print, 'RESOLUTION FACTOR (Z)		:	'+strtrim(string(b.rfz),1)
  print, ' '
  print, 'SEED PLANES			:'
  print, ' '
  print, '   	Normal Vector (Nx, Ny, Nz)   	'+$
				'	   Center Point (Rx, Ry, Rz)'
  print, *ptr_planes
  journal

end

