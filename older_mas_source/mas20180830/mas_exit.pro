;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved

;*****************************************************************************************************
;+
; NAME:
;   mas_exit
;
; PURPOSE:
;
;   Distroys the main event controlling MAS.
;
; ARGUMENTS:
;
;   event [in]    pass in event structure for main window that is to be destroyed
;
; MODIFICATION HISTORY:
;
;-
;*****************************************************************************************************
pro mas_exit, event


	COMPILE_OPT IDL2


	CATCH, Error_status
    IF Error_status NE 0 THEN BEGIN

      CATCH, /CANCEL
    END


    COMMON scan_data
    COMMON common_widgets

    ;IF N_ELEMENTS(WID_BASE_MAIN) EQ 1 THEN $
    ;    IF widget_info( WID_BASE_MAIN,  /VALID_ID ) eq 1 THEN begin
    ;    widget_control,WID_BASE_MAIN,/destroy
    ;end

	widget_control,WID_BASE_MAIN,/destroy


    HEAP_FREE, project.procpramArray
    HEAP_FREE, project.imndArray
    HEAP_FREE, project.dataarray
    HEAP_FREE, project
    HEAP_GC

;    RESET_SESSION

end
