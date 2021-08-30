;; $Id$
;;
;; Dispatch routine for motion correction.
;; use_mc == 0 => No MC
;; use_mc == 1 => ``Edge Method'' motion cor.
;; use_mc == 2 => ``Mutual Info.'' motion cor.

PRO mas_motion_correct

    COMMON scan_data

    use_mc = project.procPramArray[project.ci].mc_enable

    if (use_mc eq 1) then begin
        
        mas_motion_correct_ed_GUI

    endif

    if (use_mc eq 2) then begin

        mas_motion_correct_mi_GUI
        
    endif
    
    return

END
