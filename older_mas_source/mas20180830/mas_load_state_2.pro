;; $Id
;; Copyright 2003 University of Florida. All Rights Reserved

;*****************************************************************************************************
;
; NAME:
;   mas_load_state_2
;
; PURPOSE:
;
;
; ARGUMENTS:
; any changes to this code will be to be made in mas_display to
; equally transfrom the matrix of floating values (fpvals) that
; displayed when mousing over image pixels
;
;
; MODIFICATION HISTORY:
;
;
;*****************************************************************************************************
pro mas_load_state_2, native_orientation=native_orientation
    COMPILE_OPT IDL2
    COMMON scan_data

    ci = project.ci

    ;;check to see if the current state needs to be loaded
    if project.procPramArray[ci].state_2 eq 0 then begin

        if (ptr_valid(project.dataArray[ci].state2)) then ptr_free, project.dataarray[ci].state2
        
        mas_load_state_1
        if (project.procpramarray[ci].state_1 eq 0) then return
        
        scan = mas_subset(project.dataArray[ci].state1)
        
        mas_remove_background_filter, scan

        if (not keyword_set(native_orientation)) then begin
            mas_zoom, scan;;;, /scale_voxels
        
            mas_rotate_flip, scan
        endif
        
        ;;this mask is done later on b/c the roi is drawn on this stage.
        mas_roi_mask, scan
        
        project.dataArray[ci].state2 = scan
        ;;after we have done the processing set the state 2 flag to true
        project.procPramArray[ci].state_2 = 1
    endif
    
    HEAP_GC

end
