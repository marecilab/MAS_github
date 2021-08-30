; IDL code to unwrap phase images using FSL Prelude or in house Fourier based method

; Subroutine name: mas_phase_unwrap_fsl
; Created by: Magdoom Kulam
; Calling Information:
;
; ndim - Dimension of experiment (2D/3D)
; cdata - 2D/3D/4D array of complex data
; Bugs or Important Comments to Developers:


; Purpose of subroutine:

; Editing Information: Unwrap phase images using FSL Prelude or in house Fourier based method
; 

function mas_phase_unwrap, ndim, cdata

common scan_data
CI = project.ci
N = size(cdata,/dimension)
unwrap_result = make_array(N,value = 0, /float)
if n_elements(N) gt 3 then narray = N[3] else narray = 1

setenv, 'FSLOUTPUTTYPE=NIFTI_GZ'
spawn, '/usr/local/fsl/bin/prelude', errMsg
if errMsg eq '' and !VERSION.OS_FAMILY ne 'Windows' then begin
  unwrap_result = phase_unwrap_fsl(project.imndarray[CI].dimensions, cdata)
endif else begin
  if ndim eq 3 then begin
    for i = 0, narray-1 do begin
      phase_image = atan(cdata[*,*,*,i], /PHASE)
      phase_unwrap_3D, phase_image
      unwrap_result[*,*,*,i] = phase_image
    endfor
  endif else begin
    for i = 0, project.imndarray[CI].adim-1 do begin
      for j=0, N[3]-1 do begin
        phase_image = atan(cdata[*,*,j,i], /PHASE)
        phase_unwrap_2D, phase_image
        unwrap_result[*,*,j,i] = phase_image
      endfor
    endfor
  endelse
endelse

return,unwrap_result

end

