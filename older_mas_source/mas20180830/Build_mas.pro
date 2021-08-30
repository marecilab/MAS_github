;; $Id$
;;

;+
; NAME:
;       Build_mas
;
; PURPOSE:
;       This routine is a macro allowing the user to build and create .sav
;     files for mas and include iTool into current .sav also.
;
; CATEGORY:
;       Build
;
; CALLING SEQUENCE:
;       Build_mas
;
; OUTPUTS:
;   .sav file from the current mas_project
;
; KNOW BUGS:
; 	Not specific to MAS but rather how IDL uses the RESOLVE_ALL subroutine:
;		If there are any errors in compilation (e.g. cannot find a subroutine)
;		then you have to reset the IDLDE session by going to the Run menu. This
;		is important because the subroutine is NOT compiled and therefore will
;		be included in the uncompiled list and will NOT be resolved.
;
; EXAMPLE:
;       Build_mas
;
; MODIFICATION HISTORY:
;   Written by:  Ty Black, August, 2003
;   Copyright 2003 University of Florida. All Rights Reserved
;   HS 2001025 Added the ROI tool functionality for IDL version 6.1 or newer
;	by including the command ITRESOLVE and commenting out IDLitResolveiTools.
;-
;

PRO Build_mas
compile_opt idl2
                                ;to make sure nothing else is included with the build we have to
                                ;reset the session to make sure there are no latent procedures
                                ;or variables or objects.
;   .RESET_ALL
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;empty procedures at the end of large files. this forces them be included with the save
mas_callback
mas_redraw
MAS_SAVE
mas_sig_enhance
mas_extract_files
mas_volume
mas_background_filter_file
mas_roi
FSC_Color_file

RESOLVE_ALL , /CONTINUE_ON_ERROR

; HS 20061225
; Adding this part to include iTools in IDL version 6.1 or greater.
; Ensures that all iTools are included with this save also.

; For version 6.2 or more
IRESOLVE
; For version 6.0 or less
; IDLitResolveiTools


file = dialog_pickfile(/directory)
print, file
save, /routines, filename = file+'mas.sav'


END
