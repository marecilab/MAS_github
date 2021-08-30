;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved


; Subroutine name: mas_windowing
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro mas_windowing , scan, maxIntensity, minIntensity
    heap_gc
    COMPILE_OPT IDL2
    COMMON scan_data


    ;bring in the data that is necessary.
    intensityCen = project.procPramArray[project.ci].intensity_Cen
    intensityMax = project.procPramArray[project.ci].intensity_Max
    intensityMin = project.procPramArray[project.ci].intensity_Min

    dataValueMin = project.procPramArray[project.ci].min_display_thresh
    dataValueMax = project.procPramArray[project.ci].max_display_thresh
    
    ;if they didn't define maxIntensity, and or minIntensity
    if N_PARAMS()  eq 3 then begin
        dataValueMin = minIntensity
        dataValueMax = maxIntensity
    end


    scan=temporary(bytscl(scan, dataValueMin, dataValueMax ))

    scan=temporary(bytscl(scan, max=intensityMax, MIN=intensityMin , TOP =intensityCen  ))


end
