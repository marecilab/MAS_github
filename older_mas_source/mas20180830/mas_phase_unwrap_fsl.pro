; IDL code to unwrap phase images using FSL Prelude

; Subroutine name: mas_phase_unwrap_fsl
; Created by: Magdoom Kulam
; Calling Information:
;
; ndim - Dimension of experiment (2D/3D)
; cdata - Array of complex data
; Bugs or Important Comments to Developers:


; Purpose of subroutine: Laplacian function used by mas_phase_unwrap2d/3d

; Editing Information:
; 
function phase_unwrap_fsl, ndim, cdata

common scan_data

N = size(cdata,/dimension)
temp = make_array(N,value = 0, /float)
if n_elements(N) gt 3 then narray = N[3] else narray = 1

if ndim eq 3 then begin  
  for i = 0, narray-1 do begin
    phase = ptr_new(atan(cdata[*,*,*,i], /PHASE))
    magnitude = ptr_new(abs(cdata[*,*,*,i]))
    mas_export_nifti,data_ptr = phase, file_name = project.current_path + '/phase_data.nii'
    mas_export_nifti,data_ptr = magnitude, file_name = project.current_path + '/mag_data.nii'
    spawn, '/usr/local/fsl/bin/prelude -a mag_data.nii -p phase_data.nii -u phase_data_unwrap.nii'
    temp_nifti = mas_read_nifti(nifti_filename = project.current_path +'/phase_data_unwrap.nii.gz')
    temp[*,*,*,i] = *temp_nifti.voxel_data
    file_delete, project.current_path + '/phase_data.nii',project.current_path + '/mag_data.nii', project.current_path +'/phase_data_unwrap.nii.gz'
    ptr_free, phase, magnitude
  endfor
endif else begin
  for i = 0,narray-1 do begin
    for j = 0, N[2]-1 do begin
      phase = ptr_new(atan(cdata[*,*,j,i], /PHASE))
      magnitude = ptr_new(abs(cdata[*,*,j,i]))
      mas_export_nifti,data_ptr = phase, file_name = project.current_path + '/phase_data.nii'
      mas_export_nifti,data_ptr = magnitude, file_name = project.current_path + '/mag_data.nii'
      spawn, '/usr/local/fsl/bin/prelude -a mag_data.nii -p phase_data.nii -u phase_data_unwrap.nii'
      temp_nifti = mas_read_nifti(nifti_filename = project.current_path +'/phase_data_unwrap.nii.gz')
      temp[*,*,j,i] = *temp_nifti.voxel_data
      file_delete, project.current_path + '/phase_data.nii',project.current_path + '/mag_data.nii', project.current_path +'/phase_data_unwrap.nii.gz'
      ptr_free, phase, magnitude
    endfor
  endfor
endelse

return, reform(temp)  
end