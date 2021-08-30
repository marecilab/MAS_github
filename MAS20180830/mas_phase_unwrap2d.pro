; IDL code to unwrap 2D phase images

; Subroutine name: laplacian_2D
; Created by: Magdoom Kulam
; Calling Information:
;
; image - Input image
;
; Bugs or Important Comments to Developers:


; Purpose of subroutine: Laplacian function used by mas_phase_unwrap2d/3d

; Editing Information:

function laplacian_2D, image

sz_image = size(image,/Dimension)
x_size = sz_image[0]
y_size = sz_image[1]

p = x_size*(findgen(x_size)/float(x_size-1)-0.5)
q = y_size*(findgen(y_size)/float(y_size-1)-0.5)
meshgrid,p,q

fft_part = shift(fft(image),x_size/2, y_size/2)
fft_part*=(p^2+q^2)
fft_part = shift(fft_part, -x_size/2, -y_size/2)
ifft_part = fft(fft_part, /inverse)

laplacian_image = (-4*!pi^2)*ifft_part
return, real_part(laplacian_image)

end

; Subroutine name: inverse_laplacian_2D
; Created by: Magdoom Kulam
; Calling Information:
;
; image - Input image
;
; Bugs or Important Comments to Developers:


; Purpose of subroutine: Inverse Laplacian function used by mas_phase_unwrap2d/3d
; Modifications :
; 
function inverse_laplacian_2D, image
  
sz_image = size(image,/Dimension)
x_size = sz_image[0]
y_size = sz_image[1]

p = x_size*(findgen(x_size)/float(x_size-1)-0.5)
q = y_size*(findgen(y_size)/float(y_size-1)-0.5)
meshgrid,p,q

fft_part = shift(fft(image),x_size/2, y_size/2)
fft_part/= (p^2+q^2)
fft_part[where(finite(fft_part,/nan))] = 0.0
fft_part[where(finite(fft_part,/infinity))] = 0.0
fft_part = shift(fft_part, -x_size/2, -y_size/2)
ifft_part = fft(fft_part, /inverse)

inv_laplacian_image = -(1/(4*!pi^2))*ifft_part
return, real_part(inv_laplacian_image)

end

; Subroutine name: phase_unwrap_2D
; Created by: Magdoom Kulam
; Calling Information:
;
; phase_image - Wrapped phase image
;
; Bugs or Important Comments to Developers:
; Paper : Marvin and Zhu, Fast phase unwrapping algorithm for interferometric applications, Optics Letters (28) 14, 2003
;
; Purpose of subroutine: To unwrap the given 2D phase image
;
; Modifications :
; 
pro phase_unwrap_2D, phase_image

sz_image = size(phase_image,/Dimensions)
; Resizing the image via reflection
phase_image = [[phase_image, reverse(phase_image)],[reverse(phase_image,2), reverse(reverse(phase_image),2)]]

sin_phase_image = sin(phase_image)
cos_phase_image = cos(phase_image)

la_phase_image = laplacian_2D(phase_image)
la_true_phase = cos_phase_image*laplacian_2D(sin_phase_image) - sin_phase_image*laplacian_2D(cos_phase_image)

n = round(inverse_laplacian_2D(la_true_phase - la_phase_image)/(2.*!pi))
true_phase_image = phase_image + 2.*!pi*n

;phase_unwrap_fine,true_phase_image, inverse_laplacian_2D(la_true_phase)
phase_image = true_phase_image[0:sz_image[0]-1,0:sz_image[1]-1]

end

; Subroutine name: phase_unwrap_fine
; Created by: Magdoom Kulam
; Calling Information:
;
; phi       - Unwrapped phase image from direct method to be refined
; phi_prime - True phase image based on the wrapped image
;
; Bugs or Important Comments to Developers:
; Paper : Marvin and Zhu, Fast phase unwrapping algorithm for interferometric applications, Optics Letters (28) 14, 2003
;
; Purpose of subroutine: To perform fine phase unwrapping as mentioned in the paper.
;
; Modifications :

pro phase_unwrap_fine, phi,phi_prime

res0 = max(abs(phi_prime-phi))
dres = res0
while abs(dres) ge 1e-4 do begin
   phi += 2*!pi*round((phi_prime-phi)/(2*!pi))
   res1 = max(abs(phi_prime-phi))
   dres = res1-res0
   res0 = res1
endwhile

end
