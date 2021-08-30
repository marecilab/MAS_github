;; $Id$
;+
; Copyright 2003 University of Florida. All Rights Reserved
;-

; Subroutine name: mas_rotate_flip
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting
    ;
    ;Edited by BT 2008/04/07
    ;Fixed flip of orientation image

pro mas_rotate_flip , p_data
    COMPILE_OPT IDL2

    common scan_data
    sz_scan          = size(*p_data)
    rotate_Direction = project.procPramArray[project.ci].rotate_Direction
    flip_Direction   = project.procPramArray[project.ci].flip_Direction

;;   scan = ROTATE(scan, rotate_Direction)

    update_status_bar,'Rotating and Flipping'
    if sz_scan[0] eq 2 then begin
        ;;2d image just rotate and flip it.
        (*p_data) = ROTATE((*p_data), rotate_Direction)
    end

    if sz_scan[0] eq 3 then begin

        if rotate_Direction eq 1 or  rotate_Direction eq 3 then begin
            temp = (*p_data)
            (*p_data) = fltarr(sz_scan[2] ,sz_scan[1] ,sz_scan[3])
            
            for ii = 0 , sz_scan[3]-1 do begin
                ;;2d image just rotate and flip it.
                (*p_data)[*,*,ii] = ROTATE(temp[*,*,ii], rotate_Direction)
            endfor
        endif

        ;; why are we doing anything for 0 rotation anyway?
;;        if rotate_Direction eq 0 or  rotate_Direction eq 2 then begin
        if rotate_Direction eq 2 then begin
            for ii = 0 , sz_scan[3]-1 do begin
                ;;2d image just rotate and flip it.
                (*p_data)[*,*,ii] = temporary(ROTATE((*p_data)[*,*,ii], rotate_Direction))
            endfor
        endif

    endif

    if sz_scan[0] eq 4 then begin

        if rotate_Direction eq 1 or  rotate_Direction eq 3 then begin
            temp = (*p_data)
            (*p_data) = fltarr(sz_scan[2] ,sz_scan[1] ,sz_scan[3] ,sz_scan[4])

            for jj = 0 , sz_scan[4]-1 do $
              for ii = 0 , sz_scan[3]-1 do begin
                ;;2d image just rotate and flip it.
                (*p_data)[*,*,ii,jj] = ROTATE(temp[*,*,ii,jj], rotate_Direction)
            end
        endif

;;        if rotate_Direction eq 0 or  rotate_Direction eq 2 then begin
        if rotate_Direction eq 2 then begin
            for jj = 0 , sz_scan[4]-1 do $
              for ii = 0 , sz_scan[3]-1 do begin
                ;;2d image just rotate and flip it.
                (*p_data)[*,*,ii,jj] = temporary(ROTATE((*p_data)[*,*,ii,jj], rotate_Direction))
                
            end

        endif
        
    endif


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    
    if sz_scan[0] eq 5 then begin

        if rotate_Direction eq 1 or  rotate_Direction eq 3 then begin
            temp = (*p_data)
            (*p_data) = fltarr(sz_scan[2] ,sz_scan[1] ,sz_scan[3] ,sz_scan[4], sz_scan[5])
            for kk = 0, sz_scan[5]-1 do $
              for jj = 0 , sz_scan[4]-1 do $
              for ii = 0 , sz_scan[3]-1 do begin
                ;;2d image just rotate and flip it.
                (*p_data)[*,*,ii,jj,kk] = ROTATE(temp[*,*,ii,jj,kk], rotate_Direction)
            end
        endif

;;        if rotate_Direction eq 0 or  rotate_Direction eq 2 then begin
        if rotate_Direction eq 2 then begin
            for kk=0, sz_scan[5]-1 do $
              for jj = 0 , sz_scan[4]-1 do $
              for ii = 0 , sz_scan[3]-1 do begin
                ;;2d image just rotate and flip it.
                (*p_data)[*,*,ii,jj,kk] = temporary(ROTATE((*p_data)[*,*,ii,jj,kk], rotate_Direction))
                
            end

        endif

    endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    

    if sz_scan[0] gt 5 then begin
        print, 'WARNING: cannot rotate data set with more than 5 dimensions'
    endif

    ;; *pdata is "undefined" if the image is an Orientation image
    if (n_elements(*p_data) ne 0) then begin

        if flip_Direction eq 1 then (*p_data) = temporary(reverse((*p_data), 1))
        if flip_Direction eq 2 then (*p_data) = temporary(reverse((*p_data), 2))
        
    endif

    update_status_bar,''


end
