;; $Id$
;;

pro mas_remove_scan_all

    COMPILE_OPT IDL2

    COMMON scan_data
    COMMON common_widgets

    totalnumberofscans = project.ni

    for z=0, totalnumberofscans do begin

        ;;if there are no images to load then do nothing

        if project.ni eq 0 then begin
            Widget_control, WID_SCAN_LIST , SET_VALUE= '' ,SENSITIVE=0
            mas_redraw_GUI
            return
        endif

        if project.ni-1 eq project.ci $
          or (project.ci eq 0 and project.ni eq 0) then begin
            
            ;;print,'removing top'
            
            if (project.procpramarray[project.ci].ortho_tlb ne -1) then begin
                if (widget_info(project.procpramarray[project.ci].ortho_tlb, /valid_id)) then begin
                    widget_control, project.procpramarray[project.ci].ortho_tlb, /destroy
                endif
            endif
            
            ;;move the clean data on top of the data we want to remove
            project.imndArray[project.ci] = project.imndArray[project.ni]
            project.dataArray[project.ci] = project.dataArray[project.ni]
            project.procPramArray[project.ci] = project.procPramArray[project.ni]

            ;;decrement the indices
            project.ni = project.ci
            project.ci =  project.ci-1

            if project.ci lt 0 then project.ci = 0
            if project.ni lt 0 then project.ni = 0

            ;;       print,'after current index',project.ci
            ;;       print,'after next Index', project.ni

            if project.ci lt 0 then project.ci = 0

            project.scan_list = project.scan_list[0:project.ni]

            if project.ni eq 0 then begin
                Widget_control, WID_SCAN_LIST , SET_VALUE= '' ,SENSITIVE=0
            endif else begin
                update_scan_selection
            endelse

	    mas_redraw_GUI

	    endif else begin
                
                if project.ni-1 ne project.ci  then begin

                    ;;         print,'Removing beginning to middle'
                    
                    for ii=project.ci, project.ni do begin
                        
                        if (project.procpramarray[project.ci].ortho_tlb ne -1) then begin
                            if (widget_info(project.procpramarray[project.ci].ortho_tlb, /valid_id)) then begin
                                widget_control, project.procpramarray[project.ci].ortho_tlb, /destroy
                            endif
                        endif

                        project.imndArray[ii] = project.imndArray[ii+1]
                        project.dataArray[ii] = project.dataArray[ii+1]
                        project.procPramArray[ii] = project.procPramArray[ii+1]
                    endfor

                    ;;1 2 3 4 5 6
                    ;;shrink the scan_list by writing over the current index

                    project.ni = project.ni -1
                    if project.ci lt 0 then project.ci = 0
                    if project.ni lt 0 then project.ni = 0

	         ;print,'after current index',project.ci
	         ;print,'after next Index', project.ni

                    for ii=project.ci, project.ni-1 do begin
                        project.scan_list[ii] = project.scan_list[ii+1]
                    endfor
	         ;1 2 4 5 6 6
	         ;remove the extra 1 at the end of the array
                    temp = project.scan_list[0:project.ni-1]

                    project.scan_list = strarr ( project.ni )
                    project.scan_list = temp

                    update_scan_selection
                    
                endif

            endelse
            
	    if project.ni eq 0 then project.scan_Open_Flag = 0

        ENDFOR
        
        mas_redraw_GUI
        HEAP_GC

end




; Copyright 2003 University of Florida. All Rights Reserved
;*****************************************************************************************************
;
; NAME:
;   mas_remove_scan
;
; PURPOSE:
;
;
; ARGUMENTS:
;
;
; MODIFICATION HISTORY:
; HS 2006/09/30
; Fix spelling mistakes.
;
;*****************************************************************************************************
pro mas_remove_scan
    COMPILE_OPT IDL2

    ;print, 'remove scan'

    COMMON scan_data
    COMMON common_widgets

;    print,'current index',project.ci
;    print,'next Index', project.ni

    ;if there are no images to load then do nothing
    if project.ni eq 0 then return

    if project.ni-1 eq project.ci $
         or (project.ci eq 0 and project.ni eq 0) then begin

       ;print,'removing top'
       if widget_info(project.procpramarray[project.ci].ortho_tlb, /valid_id) then begin
          widget_control, project.procpramarray[project.ci].ortho_tlb , /destroy
       endif

       ;move the clean data on top of the data we want to remove
       project.imndArray[project.ci] = project.imndArray[project.ni]
        project.dataArray[project.ci] = project.dataArray[project.ni]
        project.procPramArray[project.ci] = project.procPramArray[project.ni]

       ;decrement the indixes
       project.ni = project.ci
       project.ci =  project.ci-1

        if project.ci lt 0 then project.ci = 0
        if project.ni lt 0 then project.ni = 0


;       print,'after current index',project.ci
;       print,'after next Index', project.ni

        if project.ci lt 0 then project.ci = 0

        project.scan_list = project.scan_list[0:project.ni]
        
        if project.ni eq 0 then begin
            Widget_control, WID_SCAN_LIST , SET_VALUE= '' ,SENSITIVE=0
        endif else begin
            update_scan_selection
        end

    mas_redraw_GUI

    endif else begin
       if project.ni-1 ne project.ci  then begin

;         print,'Removing beginning to middle'

          if widget_info(project.procpramarray[project.ci].ortho_tlb, /valid_id) then begin
             widget_control, project.procpramarray[project.ci].ortho_tlb , /destroy
          endif

         for ii=project.ci, project.ni do begin
          project.imndArray[ii] = project.imndArray[ii+1]


          project.dataArray[ii] = project.dataArray[ii+1]


          project.procPramArray[ii] = project.procPramArray[ii+1]


         endfor
         ;1 2 3 4 5 6
         ;shrink the scan_list by writing over the current index

         project.ni = project.ni -1
         if project.ci lt 0 then project.ci = 0
         if project.ni lt 0 then project.ni = 0

         ;print,'after current index',project.ci
         ;print,'after next Index', project.ni

         for ii=project.ci, project.ni-1 do begin
          project.scan_list[ii] = project.scan_list[ii+1]
         endfor
         ;1 2 4 5 6 6
         ;remove the extra 1 at the end of the array
         temp = project.scan_list[0:project.ni-1]

         project.scan_list = strarr ( project.ni )
         project.scan_list = temp

         update_scan_selection

       endif
    endelse

    if project.ni eq 0 then project.scan_Open_Flag = 0

    mas_redraw_GUI


    HEAP_GC

end
