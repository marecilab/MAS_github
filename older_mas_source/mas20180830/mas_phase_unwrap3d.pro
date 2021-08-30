; IDL code to unwrap 3D phase images

; Subroutine name: laplacian_3D
; Created by: Magdoom Kulam
; Calling Information:
;
; image - Input image
;
; Bugs or Important Comments to Developers:


; Purpose of subroutine: Laplacian function used by phase_unwrap_3D

; Editing Information:

function laplacian_3D, image

sz_image = size(image,/Dimension)
x_size = sz_image[0]
y_size = sz_image[1]
z_size = sz_image[2]

p = x_size*(findgen(x_size)/float(x_size-1)-0.5)
q = y_size*(findgen(y_size)/float(y_size-1)-0.5)
r = z_size*(findgen(z_size)/float(z_size-1)-0.5)
meshgrid,p,q, z = r
  
fft_part = shift(fft(image),x_size/2, y_size/2,z_size/2)
fft_part*=(p^2+q^2+r^2)
fft_part = shift(fft_part, -x_size/2, -y_size/2,-z_size/2)
ifft_part = fft(fft_part, /inverse)

laplacian_image = (-4*!pi^2)*ifft_part
return, real_part(laplacian_image)

end

; Subroutine name: inverse_laplacian_3D
; Created by: Magdoom Kulam
; Calling Information:
;
; image - Input image
;
; Bugs or Important Comments to Developers:


; Purpose of subroutine: Inverse Laplacian function used by phase_unwrap_3D
; Modifications :
;
function inverse_laplacian_3D, image

sz_image = size(image,/Dimension)
x_size = sz_image[0]
y_size = sz_image[1]
z_size = sz_image[2]

p = x_size*(findgen(x_size)/float(x_size-1)-0.5)
q = y_size*(findgen(y_size)/float(y_size-1)-0.5)
r = z_size*(findgen(z_size)/float(z_size-1)-0.5)
meshgrid,p,q, z = r

fft_part = shift(fft(image),x_size/2, y_size/2,z_size/2)
fft_part/= (p^2+q^2+r^2)
fft_part[where(finite(fft_part,/nan))] = 0.0
fft_part[where(finite(fft_part,/infinity))] = 0.0
fft_part = shift(fft_part, -x_size/2, -y_size/2, -z_size/2)
ifft_part = fft(fft_part, /inverse)

inv_laplacian_image = -(1/(4*!pi^2))*ifft_part
return, real_part(inv_laplacian_image)

end

; Subroutine name: phase_unwrap_3D
; Created by: Magdoom Kulam
; Calling Information:
;
; phase_image - Wrapped phase image
;
; Bugs or Important Comments to Developers:
; Paper : Marvin and Zhu, Fast phase unwrapping algorithm for interferometric applications, Optics Letters (28) 14, 2003
;
; Purpose of subroutine: To unwrap the given 3D phase image 
;
; Modifications :

pro phase_unwrap_3D, phase_image, tolerance = tol

sz_image = size(phase_image,/Dimensions)
  
; Resizing the image via reflection
phase_image_ref = make_array(2*sz_image,/double)
phase_image_ref[*,*,0:sz_image[2]-1] = [[phase_image, reverse(phase_image)],[reverse(phase_image,2), reverse(reverse(phase_image),2)]]
phase_image_ref[*,*,sz_image[2]:2*sz_image[2]-1] = reverse(phase_image_ref[*,*,0:sz_image[2]-1],3)

sin_phase_image = sin(phase_image_ref)
cos_phase_image = cos(phase_image_ref)

la_phase_image = laplacian_3D(phase_image_ref)
la_true_phase = cos_phase_image*laplacian_3D(sin_phase_image) - sin_phase_image*laplacian_3D(cos_phase_image)

if n_elements(tol) lt 1 then tol = 2*!DPI
n = round(inverse_laplacian_3D(la_true_phase - la_phase_image)/tol)
true_phase_image = phase_image_ref + tol*n

phase_unwrap_fine,true_phase_image, inverse_laplacian_3D(la_true_phase)
phase_image = true_phase_image[0:sz_image[0]-1,0:sz_image[1]-1,0:sz_image[2]-1]

end
