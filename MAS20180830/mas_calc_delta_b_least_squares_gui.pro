function cdbls_down_sample_array,array, factor

size_array  = size(array)
new_array = dblarr(size_array[1]/factor, size_array[2]/factor)

for i=0, size_array[1]/factor-1 do begin
  for j=0, size_array[2]/factor-1 do begin
    i_index = i*factor + (factor - 1)
    j_index = j*factor + (factor - 1)
    new_array[i,j] = array[i_index,j_index]
  endfor
endfor

return, new_array
end

FUNCTION simulation_dB_create_slice,pbuffer, fdim, pdim, roi_set_index, roi_name_index, object_chi, media_chi, interp_factor
common scan_data

if ptr_valid((project.roi.pROIs[roi_set_index])[roi_name_index]) eq 0 then begin
  mask_temp = *project.deltaB_mask
  mask = dblarr(fdim,pdim)
  mask[pbuffer-1:fdim-pbuffer-1,pbuffer-1:pdim-pbuffer-1] = mask_temp[pbuffer-1:fdim-pbuffer-1,pbuffer-1:pdim-pbuffer-1]
  mask = (1-mask)*255
  xdim = fdim
  ydim = fdim 
endif else begin  
  object = (*project.roi.pROIs[roi_set_index])[roi_name_index]
  object[0]->scale, [interp_factor, interp_factor]
  xdim = interp_factor*fdim
  ydim = interp_factor*pdim
  mask = ptr_new(object -> ComputeMask( dimensions =[ xdim, ydim ] , MASK_RULE=2))
  mask = *mask
  object[0]->scale, [1.0/interp_factor, 1.0/interp_factor]
endelse

slice = dblarr(xdim, ydim)
for i=0,xdim-1 do begin
  for j=0, ydim-1 do begin
    if mask[i,j] eq 255 then begin
      chi = object_chi
    endif else begin
      chi = media_chi
    endelse
    slice[i,j] = chi
  endfor
endfor


return, slice
END


FUNCTION simulation_dB_calc_map, slice, fdim, pdim, field_strength
  
  ;setting up the default values
  B0 = field_strength
  mu0 = 4*!DPI*1e-07;
  

  ; setting up the matrix points
  

  ;setting up the angle of main magnetic field vs the model
  zeta = 90.0*(!PI/180)
  cos_zeta = cos( zeta )
  sin_zeta = sin( zeta )

  samples_x = fdim
  samples_y = pdim

  sizearray = [samples_x, samples_y]

  xdim_center = float(round(samples_x/2))
  ydim_center = float(round(samples_y/2))
  
 delta_chi_2 = complexarr(ceil(samples_x), ceil(samples_y))
  Mk = fltarr(ceil(samples_x), ceil(samples_y))


  ; this is a flag to only calculate the bout once
  ;closed_solution = 0

  ; Assuming susceptibility << 1
  Mz = (slice) * (B0/mu0);
  Mk_temp = (FFT(Mz))

  Mk = shift(Mk_temp, xdim_center,ydim_center)

  ; Creating a progress bar outside of the looping through the core material sizes


  for j=0.0, samples_y-1 do begin
      for i=0.0, samples_x-1 do begin
          k_x = (float(i)-xdim_center)
        k_y = (float(j)-ydim_center)
        k_z = float(0)

        if k_x eq 0.0 and k_y eq 0.0 and k_z eq 0.0 then begin
          delta_chi_2[i,j] = 0.0
        endif else begin
          cosinepart = ((k_z*cos_zeta + k_x*sin_zeta)^2.0)/(k_z^2.0 + k_x^2.0 + k_y^2.0)
          delta_chi_2[i,j] = -(Mu0/3.0)*Mk[i,j]*(3.0*cosinepart - 1.0)
        endelse

      endfor
    endfor


    delta_chi_2 = shift(delta_chi_2, xdim_center,ydim_center)
    results = REAL_PART(FFT(delta_chi_2,/INVERSE))
    

  return, results

END

function simulate_dB_calc_r2, simulation, measured, r2_window_mask


no_mask = where(measured ne 0.0 and r2_window_mask ne 0.0)
SSerr = total((measured[no_mask] - simulation[no_mask])^2.0)
SStot = total((measured[no_mask] - mean(measured[no_mask]))^2.0)

rsquare = 1.0 - SSerr/SStot 
return, rsquare
end

function simulate_dB_calc_r2_image, simulation, measured, r2_mask, mask, r2_window_mask

no_mask = where(mask ne 0.0 and r2_window_mask ne 0.0 )
SSerr = total((measured[no_mask] - simulation[no_mask])^2.0)
SStot = total((measured[no_mask] - mean(measured[no_mask]))^2.0)

rsquare = 1.0 - SSerr/SStot 
return, rsquare
end

function cdbls_uncertain_dB, dB
  common scan_data
  sigma_mat = *project.deltaB_sigma
  sz = size(sigma_mat)
  dB_new = dblarr(sz[1],sz[2])
  for i=0,sz[1]-1 do begin
    for j=0, sz[2]-1 do begin
      mean_db = dB[i,j]
      dB_new[i,j] = sigma_mat[i,j]*randomu(x,/normal)+ mean_db
    endfor
  endfor
  return, dB_new
end
function mas_linear_fit, x, y

  sz_x = size(x)
  ln = 0.0
  ld = 0.0
  x_mean = mean(x)
  y_mean = mean(y)
  for i=0,sz_x[1]-1 do begin
    ln = ln + (x[i]-x_mean)*(y[i] - y_mean)
    ld = ld + (x[i]-x_mean)^2.0
  end
  
  slope = ln/ld
  intercept = y_mean - slope*x_mean

  result = [slope,intercept]
return, result
end


function line_data_first_order_pc, line_data, p_buffer, pdim, pfov, width
common scan_data
  pres = pfov/pdim
  x_all = make_array(pdim,1, /double, /index)*pres
  x_fit = [x_all[p_buffer:p_buffer+width,0] , x_all[pdim - p_buffer - width:pdim-p_buffer,0]]
  y_fit = [line_data[p_buffer:p_buffer+width] , line_data[pdim - 1 - p_buffer - width:pdim-p_buffer-1]]
    
  fit_result = mas_linear_fit(x_fit, y_fit)
  slope = fit_result[0]
  intercept = fit_result[1] 
  
  sz_x_all = size(x_all)
  line_data_corrected = line_data

  correction_term = slope*x_all + intercept
  line_data_corrected = line_data_corrected - correction_term
  line_data_corrected = line_data_corrected + mean(correction_term)
return, line_data_corrected
end

function image_first_order_pc, measured_map, p_buffer, pdim, pfov, fdim, ffov, width, row, column
common scan_data
  pres = pfov/pdim
  line_data = measured_map[column, *]
  x_all = make_array(pdim,1, /double, /index)*pres
  x_fit = [x_all[p_buffer:p_buffer+width,0] , x_all[pdim - p_buffer - width:pdim-p_buffer,0]]
  y_fit = [line_data[p_buffer:p_buffer+width] , line_data[pdim - 1 - p_buffer - width:pdim-p_buffer-1]]
  fit_result = mas_linear_fit(x_fit, y_fit)
  slope = fit_result[0]
  intercept = fit_result[1]
  
  ;applying vertical correction
  measured_map_corrected = measured_map
  for i = 0, fdim-1 do begin
    correction_term = slope*x_all + intercept
    measured_map_corrected[i,*] = measured_map_corrected[i,*] - correction_term[i]
    measured_map_corrected[i,*] = measured_map_corrected[i,*] + mean(correction_term)
  endfor  
  
  rres = ffov/fdim
  line_data = measured_map_corrected[*,row]
  x_all = make_array(fdim,1, /double, /index)*rres
  x_fit = [x_all[p_buffer:p_buffer+width,0] , x_all[pdim - p_buffer - width:pdim-p_buffer,0]]
  y_fit = [line_data[p_buffer:p_buffer+width] , line_data[pdim - 1 - p_buffer - width:pdim-p_buffer-1]]
  fit_result = mas_linear_fit(x_fit, y_fit)
  slope = fit_result[0]
  intercept = fit_result[1]
  
  ;applying horizontal correction
  for j = 0, pdim-1 do begin
    correction_term = slope*x_all + intercept
    measured_map_corrected[*,j] = measured_map_corrected[*,j] - correction_term[j]
    measured_map_corrected[*,j] = measured_map_corrected[*,j] + mean(correction_term)
  endfor 
  
return, measured_map_corrected
end



pro line_data_calc_relative_chi, info, output1 = output1
common scan_data
  ;Set up results structure
  relative_chi_struct = {measured_line:ptr_new(), simulated_line:ptr_new(), simulated_map:ptr_new(), $
                         chi_vector:ptr_new(), r2_vector:ptr_new()}
  
  ;Get the measured line data, determined by user, and do a first order background field correction
  column = info.chi_line_column
  CI = info.CI[0]
  pdim = project.IMNDArray[CI].pdim
  fdim = project.IMNDArray[CI].fdim
  pfov = project.IMNDArray[CI].p_fov
  pbuffer = info.p_buffer
  measured_map = *project.deltaB_result
  measured_line = measured_map[column, *]
  if info.field_correction eq 1 then begin
    measured_line_corrected = line_data_first_order_pc(measured_line, pbuffer, pdim, double(pfov), info.field_correction_width)
  endif else begin
    measured_line_corrected = measured_line
  endelse
  line_mask = (*project.deltaB_mask)[column,*]
  measured_line_corrected = measured_line_corrected*line_mask
  measured = measured_line_corrected[pbuffer-1:pdim - pbuffer -1]
  
  
  relative_chi_struct.measured_line = ptr_new(measured)
  
  ;Iterate through simulations of the B field until delta_r2 is less than the user specified stopping point 
  delta_r2 = double(1.0)
  r2 = 1.0
  stopping_point = double(info.stopping_point)
  count = 0
  r2_vector_all_temp = dblarr(1)
  chi_vector_all_temp = dblarr(1)
  ;progressbar = Obj_New('progressbar', Color='red', Text='Estimating Relative Chi', /NOCANCEL)
  ;progressbar -> Start
  while delta_r2 gt stopping_point do begin
    count = count + 1
    line_data_calc_relative_chi_iterate, info, measured, output1 = chi_vector, output2 = r2_vector
    sz_r2 = size(r2_vector)
    max_r2 = max(r2_vector)
    r2_index = where(r2_vector eq max_r2)
    if r2_index[0] eq 0 then begin
      info.start_chi = info.start_chi - 5.0E-6
      info.stop_chi = chi_vector[r2_index[0]] + 5.0E-6
    endif else if r2_index[0] eq sz_r2[1]-1 then begin
      info.start_chi = chi_vector[r2_index[0]] - 5.0E-6
      info.stop_chi = info.stop_chi + 5.0E-6
    endif else begin
      info.start_chi = chi_vector[r2_index[0] - 1]
      info.stop_chi = chi_vector[r2_index[0] + 1]
    endelse
    
    if count gt 10 then begin
      delta_r2 = stopping_point
      r2_vector_all_temp = [r2_vector_all_temp, r2_vector]
      chi_vector_all_temp = [chi_vector_all_temp, chi_vector]
    endif else begin
      delta_r2 = abs(r2 - max_r2)
      r2 = max_r2
      r2_vector_all_temp = [r2_vector_all_temp, r2_vector]
      chi_vector_all_temp = [chi_vector_all_temp, chi_vector]
    endelse 
   ;progressbar -> Update, (float(stopping_point)/float(delta_r2))*100
   endwhile
  ;progressbar -> Destroy
  r2_vector_all = r2_vector_all_temp[1:*]
  chi_vector_all = chi_vector_all_temp[1:*]
  r2_max = max(r2_vector_all)
  r2_index = where(r2_vector_all eq r2_max)
  chi_best = chi_vector_all[r2_index]
  
  interp_factor = info.model_interp
  slice = simulation_dB_create_slice(pbuffer,fdim,pdim,info.roi_set,info.roi_name,chi_best, info.media_chi, interp_factor)
  db_slice_temp = simulation_dB_calc_map(slice, interp_factor*fdim, interp_factor*pdim, info.field_strength)
  db_slice = cdbls_down_sample_array(db_slice_temp, interp_factor)
  db_mask = *project.deltaB_mask
  db_slice = db_slice*db_mask
  simulated_line = db_slice[column,*]
  simulated_line = simulated_line[pbuffer-1:pdim-pbuffer-1]
  
  
  relative_chi_struct.chi_vector = ptr_new(chi_vector_all)
  relative_chi_struct.r2_vector = ptr_new(r2_vector_all)
  relative_chi_struct.simulated_map = ptr_new(db_slice)
  relative_chi_struct.simulated_line = ptr_new(simulated_line)
  
  output1 = relative_chi_struct
end

pro line_data_calc_relative_chi_iterate, info, measured_line, output1= output1, output2 = output2
common scan_data

CI = info.CI[0]
fdim = project.IMNDArray[CI].fdim
pdim = project.IMNDArray[CI].pdim
field_strength = double(info.field_strength)
pbuffer = info.p_buffer
column = info.chi_line_column
roi_set = info.roi_set
roi_name = info.roi_name
interp_factor = info.model_interp

min_chi = double(info.start_chi)
max_chi = double(info.stop_chi)
media_chi = double(info.media_chi)
num_steps = double(info.num_steps)
chi_step = (max_chi - min_chi)/(double(num_steps)-1.0)
chi_vector = make_array(num_steps, /double, /index)*chi_step + min_chi

sz_chi = size(chi_vector)
r2_vector = dblarr(sz_chi[1])

for x=0,sz_chi[1]-1 do begin
  object_chi = chi_vector[x]
  slice = simulation_dB_create_slice(pbuffer,fdim,pdim,roi_set,roi_name,object_chi, media_chi, interp_factor)
  db_slice_temp = simulation_dB_calc_map(slice, interp_factor*fdim, interp_factor*pdim, field_strength)
  db_slice = cdbls_down_sample_array(db_slice_temp, interp_factor)
  db_mask = *project.deltaB_mask
  db_slice = db_slice*db_mask
  simulated_line = db_slice[column,*]
  simulated_line = simulated_line[pbuffer-1:pdim-pbuffer-1]
  sz_mask = size(simulated_line)
  r2_window_mask = dblarr(sz_mask[1])
  if info.rsquare_window_flag eq 1 then begin
    column = info.chi_line_column
    row = info.chi_line_row
    r2_window = info.chi_line_window
    r2_window = floor(r2_window/2.0)
    r2_window_mask[row - r2_window:row + r2_window] = 1.0
  endif else begin
    r2_window_mask[*] = 1.0
  endelse
  r2_vector[x] = simulate_dB_calc_r2(simulated_line, measured_line, r2_window_mask)
endfor


output1 = chi_vector
output2 = r2_vector
end

pro line_data_calc_relative_chi_image, info, output1 = output1
common scan_data
  ;Set up results structure
  relative_chi_struct = {measured_line:ptr_new(), simulated_line:ptr_new(), simulated_map:ptr_new(), $
                         chi_vector:ptr_new(), r2_vector:ptr_new()}
  
  column = info.chi_line_column
  row = info.chi_line_row
  CI = info.CI[0]
  pdim = project.IMNDArray[CI].pdim
  fdim = project.IMNDArray[CI].fdim
  pfov = project.IMNDArray[CI].p_fov
  ffov = project.IMNDArray[CI].f_fov
  pbuffer = info.p_buffer
  measured_map = *project.deltaB_result
  
  if info.field_correction eq 1 then begin
    measured_map_corrected = image_first_order_pc(measured_map, pbuffer, pdim, double(pfov), fdim, double(ffov), $ 
    info.field_correction_width, row, column)
  endif else begin
    measured_map_corrected = measured_map
  endelse
  
  measured_map = measured_map_corrected[pbuffer-1:fdim - pbuffer -1, pbuffer-1:pdim - pbuffer -1]
  
  
  
  
  relative_chi_struct.measured_line = ptr_new(measured_map[column, *])
  
  ;Iterate through simulations of the B field until delta_r2 is less than the user specified stopping point 
  delta_r2 = double(1.0)
  r2 = 1.0
  stopping_point = double(info.stopping_point)
  count = 0
  r2_vector_all_temp = dblarr(1)
  chi_vector_all_temp = dblarr(1)
  
  while delta_r2 gt stopping_point do begin
    count = count + 1
    line_data_calc_relative_chi_iterate_image, info, measured_map, output1 = chi_vector, output2 = r2_vector
    sz_r2 = size(r2_vector)
    max_r2 = max(r2_vector)
    r2_index = where(r2_vector eq max_r2)
    ;use r2_index[0] because we get multiple indices on last iterations where all chi are the same
    if r2_index[0] eq 0 then begin
      info.start_chi = info.start_chi - 5.0E-6
      info.stop_chi = chi_vector[r2_index[0]] + 5.0E-6
    endif else if r2_index[0] eq sz_r2[1]-1 then begin
      info.start_chi = chi_vector[r2_index[0]] - 5.0E-6
      info.stop_chi = info.stop_chi + 5.0E-6
    endif else begin
      info.start_chi = chi_vector[r2_index[0] - 1]
      info.stop_chi = chi_vector[r2_index[0] + 1]
    endelse
    
    if count gt 10 then begin
      delta_r2 = stopping_point
      r2_vector_all_temp = [r2_vector_all_temp, r2_vector]
      chi_vector_all_temp = [chi_vector_all_temp, chi_vector]
    endif else begin
      delta_r2 = abs(r2 - max_r2)
      r2 = max_r2
      r2_vector_all_temp = [r2_vector_all_temp, r2_vector]
      chi_vector_all_temp = [chi_vector_all_temp, chi_vector]
    endelse 
   endwhile
  r2_vector_all = r2_vector_all_temp[1:*]
  chi_vector_all = chi_vector_all_temp[1:*]
  r2_max = max(r2_vector_all)
  r2_index = where(r2_vector_all eq r2_max)
  chi_best = chi_vector_all[r2_index]
  
  interp_factor = info.model_interp
  slice = simulation_dB_create_slice(pbuffer,fdim,pdim,info.roi_set,info.roi_name,chi_best, info.media_chi, interp_factor)
  db_slice_temp = simulation_dB_calc_map(slice, interp_factor*fdim, interp_factor*pdim, info.field_strength)
  db_slice = cdbls_down_sample_array(db_slice_temp, interp_factor)
  db_mask = *project.deltaB_mask
  ;db_slice = db_slice*db_mask
  simulated_line = db_slice[column,*]
  simulated_line = simulated_line[pbuffer-1:pdim-pbuffer-1]
  
  
  relative_chi_struct.chi_vector = ptr_new(chi_vector_all)
  relative_chi_struct.r2_vector = ptr_new(r2_vector_all)
  relative_chi_struct.simulated_map = ptr_new(db_slice)
  relative_chi_struct.simulated_line = ptr_new(simulated_line)
  
  output1 = relative_chi_struct
end

pro line_data_calc_relative_chi_iterate_image, info, measured_map, output1= output1, output2 = output2
common scan_data

CI = info.CI[0]
fdim = project.IMNDArray[CI].fdim
pdim = project.IMNDArray[CI].pdim
field_strength = double(info.field_strength)
pbuffer = info.p_buffer
roi_set = info.roi_set
roi_name = info.roi_name
interp_factor = info.model_interp

min_chi = double(info.start_chi)
max_chi = double(info.stop_chi)
media_chi = double(info.media_chi)
num_steps = double(info.num_steps)
chi_step = (max_chi - min_chi)/(double(num_steps)-1.0)
chi_vector = make_array(num_steps, /double, /index)*chi_step + min_chi

sz_chi = size(chi_vector)
r2_vector = dblarr(sz_chi[1])

for x=0,sz_chi[1]-1 do begin
  object_chi = chi_vector[x]
  slice = simulation_dB_create_slice(pbuffer,fdim,pdim,roi_set,roi_name,object_chi, media_chi, interp_factor)
  db_slice_temp = simulation_dB_calc_map(slice, interp_factor*fdim, interp_factor*pdim, field_strength)
  db_slice = cdbls_down_sample_array(db_slice_temp, interp_factor)
  db_mask = *project.deltaB_mask
  db_slice = db_slice*db_mask
  simulated_map = db_slice[pbuffer-1:fdim-pbuffer-1,pbuffer-1:pdim-pbuffer-1]
  r2_map = *project.deltaB_r2
  r2_index = where(r2_map gt 0.8)
  r2_mask = dblarr(fdim, pdim)
  r2_mask[r2_index] = 1.0
  r2_mask = r2_mask[pbuffer-1:fdim-pbuffer-1,pbuffer-1:pdim-pbuffer-1]
  mask = *project.deltaB_mask
  mask = mask[pbuffer-1:fdim-pbuffer-1,pbuffer-1:pdim-pbuffer-1]
  sz_mask = size(mask)
  r2_window_mask = dblarr(sz_mask[1], sz_mask[2])
  if info.rsquare_window_flag eq 1 then begin
    column = info.chi_line_column
    row = info.chi_line_row
    r2_window = info.chi_line_window
    r2_window = floor(r2_window/2.0)
    r2_window_mask[column - r2_window:column + r2_window, row - r2_window:row + r2_window] = 1.0
  endif else begin
    r2_window_mask[*,*] = 1.0
  endelse
  
  r2_vector[x] = simulate_dB_calc_r2_image(simulated_map, measured_map, r2_mask, mask, r2_window_mask)
endfor


output1 = chi_vector
output2 = r2_vector
end

pro line_data_calc_relative_chi_print, num_steps, chi_vector, r2_vector
  
  min_chi = min(chi_vector)
  max_chi = max(chi_vector)
  r2 = max(r2_vector)
  r2_index = where(r2_vector eq r2)
  chi = chi_vector[r2_index]
  sz_chi = size(chi_vector)
  num_iter = sz_chi[1]/num_steps
  
  display_string_array = strarr(2)
  display_string_array[0] = 'Chi:' + string([09B]) + 'Rsquared' + string([09B]) + 'Max Chi' + string([09B]) + 'Min Chi' $
                             + string([09B]) + '# of iterations'
  display_format = '(E0.4)'
  display_string_array[1] = string(FORMAT=display_format,chi) + string([09B]) + string(FORMAT=display_format,r2) $
                            + string([09B]) + string(FORMAT=display_format,max_chi) $      
                            + string([09B]) + string(FORMAT=display_format,min_chi) $ 
                            + string([09B]) + string(FORMAT=display_format,num_iter) 
  display_stats, display_string_array, 'Chi Estimation from Line Data'


end

pro line_data_calc_relative_chi_plot, info, measured, simulated
   common scan_data
   
   PI = info.CI[0]
   pdim = project.IMNDArray[PI].pdim
   pfov = project.IMNDArray[PI].p_fov
   pres = pfov/pdim
   x_all = make_array(pdim,1, /double, /index)*pres
   x = x_all[info.p_buffer-1:pdim-info.p_buffer-1]
   iplot_identifier = 13
   
   iplot, x, measured, linestyle = 6, sym_index = 6, sym_size = 0.5,$ 
          overplot = iplot_identifier, color = [0,0,204]
   
   iplot, x, simulated, linestyle = 1, sym_index = 0, thick = 3.0,$ 
          overplot = iplot_identifier, color = [0,0,0],xtitle = 'Position (cm)', ytitle = 'Field Perturbation (uT)'  
  
end
function cdbls_SHARP_dB, dB, radius, xpad, ypad
  radius = round(radius)
  size_dB = size(dB)
  cropped_dB = dB[xpad:size_dB[1]-1-xpad, ypad:size_dB[2]-1-ypad]
  kernel = dblarr(2*radius+1, 2*radius+1)
  delta = dblarr(2*radius+1, 2*radius+1)
  count = 0.0
  for i = 0, 2*radius do begin
    for j = 0, 2*radius do begin
      x = i-radius
      y = j-radius
      rho = sqrt(x^2 + y^2)
      if rho le radius then begin
        kernel[i,j] = 1.0
        count = count + 1.0
        if x eq 0 and y eq 0 then begin
          delta[i,j] = 1.0
         endif
      endif
    endfor
  endfor
  kernel = kernel/count
convolved = convol(cropped_dB, kernel, /center)
new_kernel = delta - kernel
deconvolved = convol(convolved, new_kernel, /center)
return, deconvolved
end

pro cdbls_create_mask, threshold, complex_data, output = output
   size_CD = size(complex_data)
   mask = dblarr(size_CD[1], size_CD[2]) + 1.0
   for delay = 0,size_CD[3]-1 do begin
      mag_data = abs(complex_data[*,*,delay])
      max_mag = max(mag_data)
      index = where(mag_data lt threshold*max_mag)
      sz_index = size(index)
      if sz_index[0] ne 0 then mask[index] = 0.0
   endfor


   output = mask
end
pro cdbls_complex_division, info
  common scan_data
  CI = project.CI
  sdim_start = project.procpramArray[CI].sdim_start
  adim_start = project.procpramArray[CI].adim_start 
  fdim_start = project.procpramArray[CI].fdim_start
  pdim_start = project.procpramArray[CI].pdim_start
  slice_axis = project.procpramArray[CI].slice_axis
  
  CI_arr = info.CI
  index = where(CI_arr gt -1)
  size_CI = size(index)
  
  Tfp_array = dblarr(size_CI[1])
  progressbar = Obj_New('progressbar', Color='red', Text='Performing Complex Division',/NOCANCEL)
  progressbar -> Start
  if slice_axis eq 0 then begin
    if project.procpramarray[CI_arr[0]].state_1 eq 0 then begin
      project.CI = CI_arr[0]
      mas_load_state_1
      project.procpramarray[CI_arr[0]].state_1 = 1
    endif
    temp = reform((*project.dataArray[CI_arr[0]].state1)[*,*,sdim_start,adim_start])
  endif else if slice_axis eq 1 then begin
    if project.procpramarray[CI_arr[0]].state_1 eq 0 then begin
      project.CI = CI_arr[0]
      mas_load_state_1
      project.procpramarray[CI_arr[0]].state_1 = 1
    endif
    temp = reform((*project.dataArray[CI_arr[0]].state1)[*,pdim_start,*,adim_start])
  endif else if slice_axis eq 2 then begin
    if project.procpramarray[CI_arr[0]].state_1 eq 0 then begin
      project.CI = CI_arr[0]
      mas_load_state_1
      project.procpramarray[CI_arr[0]].state_1 = 1
    endif
    temp = reform((*project.dataArray[CI_arr[0]].state1)[fdim_start,*,*,adim_start])
  endif
  
  size_data = size(temp)
  complex_data_array = complexarr(size_data[1], size_data[2],size_CI[1])
  complex_data_array[*,*,0] = temp
  Tfp_array[0] = mas_extract_delays(project.imndarray[CI_arr[0]].file_path $
                                    , project.imndarray[CI_arr[0]].state1_load_procedure,te = project.imndarray[CI_arr[0]].echo_time)
  progressBar -> Update, (float(1)/float(size_CI[1]))*100.0
  for i=1,size_CI[1]-1 do begin
    if slice_axis eq 0 then begin
      if project.procpramarray[CI_arr[i]].state_1 eq 0 then begin
      project.CI = CI_arr[i]
      mas_load_state_1
      project.procpramarray[CI_arr[i]].state_1 = 1
      endif
    temp = reform((*project.dataArray[CI_arr[i]].state1)[*,*,sdim_start,adim_start])
    endif else if slice_axis eq 1 then begin
      if project.procpramarray[CI_arr[i]].state_1 eq 0 then begin
      project.CI = CI_arr[i]
      mas_load_state_1
      project.procpramarray[CI_arr[i]].state_1 = 1
      endif
      temp = reform((*project.dataArray[CI_arr[i]].state1)[*,pdim_start,*,adim_start])
    endif else if slice_axis eq 2 then begin
      if project.procpramarray[CI_arr[i]].state_1 eq 0 then begin
      project.CI = CI_arr[i]
      mas_load_state_1
      project.procpramarray[CI_arr[i]].state_1 = 1
      endif
      temp = reform((*project.dataArray[CI_arr[i]].state1)[fdim_start,*,*,adim_start])
    endif
    complex_data_array[*,*,i] = temp
    Tfp_array[i] = mas_extract_delays(project.imndarray[CI_arr[i]].file_path $
                                    , project.imndarray[CI_arr[i]].state1_load_procedure,te = project.imndarray[CI_arr[i]].echo_time)
  progressBar -> Update, (float(1+i)/float(size_CI[1]))*100.0
  endfor
  progressbar -> Destroy
  sorted_index = sort(Tfp_array)
  Tfp_array = Tfp_array[sorted_index]
  complex_data_array = complex_data_array[*,*,sorted_index]
  project.CI = CI
  if info.mask_flag eq 1 then begin
    cdbls_create_mask,info.threshold, complex_data_array, output = mask
    info.mask = ptr_new(mask)
    project.deltaB_mask = ptr_new(mask)
  endif
  
  complex_divided_array = complexarr(size_data[1], size_data[2],size_CI[1]-1)
  dTfp_array = dblarr(size_CI[1]-1)
  for d=1, size_CI[1]-1 do begin
    complex_divided_array[*,*,d-1] = complex_data_array[*,*,d]/complex_data_array[*,*,0]
    dTfp_array[d-1] = Tfp_array[d] - Tfp_array[0]
  end

info.dTfp_array = ptr_new(dTfp_array)
info.CD_array = ptr_new(complex_divided_array)
info.magnitude_image = ptr_new(abs(complex_data_array[*,*,0]))
end

pro cdbls_deltaB_display, info, output1 = output1, output2 = output2
  common scan_data
  dTfp = *info.dTfp_array
  CD_array = *info.CD_array
  size_CD = size(CD_array)
  
  if info.mask_flag eq 1 then begin
    mask = *info.mask
  endif else begin
    mask = dblarr(size_CD[1], size_CD[2]) + 1.0
  endelse
  
  dB = dblarr(size_CD[1], size_CD[2])
  rsquare = dblarr(size_CD[1], size_CD[2])
  SSerr_mat = dblarr(size_CD[1], size_CD[2])
  SStot_mat = dblarr(size_CD[1], size_CD[2])
  sigma_mat = dblarr(size_CD[1], size_CD[2])
  RSD = dblarr(size_CD[1], size_CD[2])
  x = transpose(dTfp)
  xbar = mean(x)
  temp = dblarr(size_CD[3])
  gyromagnetic_ratio = (info.gyromagnetic_ratio)*1E6 ;convert to units of Hz/Tesla
  Phase_array = atan(CD_array, /phase)
  if info.unwrap_flag eq 1 then begin
    setenv, 'FSLOUTPUTTYPE=NIFTI_GZ'
    spawn, '/usr/local/fsl/bin/prelude', errMsg
    if errMsg eq '' and !VERSION.OS_FAMILY ne 'Windows' then begin
      for i = 0,(size(Phase_array))[3]-1 do begin
        phase_ptr = ptr_new(Phase_array[*,*,i])
        mas_export_nifti,data_ptr = info.magnitude_image, file_name = project.current_path + '/mag_data.nii'
        mas_export_nifti,data_ptr = phase_ptr, file_name = project.current_path + '/phase_data.nii'
        spawn, '/usr/local/fsl/bin/prelude -a mag_data.nii -p phase_data.nii -u phase_data_unwrap.nii -s'
        temp = mas_read_nifti(nifti_filename = project.current_path +'/phase_data_unwrap.nii.gz')
        Phase_array[*,*,i] = *temp.voxel_data
        file_delete, project.current_path + '/phase_data.nii',project.current_path + '/mag_data.nii', project.current_path +'/phase_data_unwrap.nii.gz'
        ptr_free, phase_ptr
      endfor
    endif else begin ; Fourier based method
      for i = 0,(size(Phase_array))[3]-1 do begin
        phase_image = Phase_array[*,*,i]
        phase_unwrap_2D,phase_image
        Phase_array[*,*,i] = phase_image
      endfor
    endelse
  endif
  Phase_array/=(2.0*gyromagnetic_ratio*2.0*!DPI)
  for i=0, size_CD[1]-1 do begin
    for j=0, size_CD[2]-1 do begin
      if mask[i,j] ne 0.0 then begin
        y = transpose(reform(Phase_array[i,j,*]))
        
        ;Processing the least squares method of calculating dB
        if info.ls_flag eq 1 then begin
          ybar = mean(y)
          dB[i,j] = invert(transpose(x)##x)##transpose(x)##y
          SStot = 0.0 
          SSerr = 0.0
          SSxx = 0.0
          for n=0,size_CD[3]-1 do begin
            SStot = SStot + (y[n] - ybar)^2.0
            fn = dB[i,j]*x[n]
            SSerr = SSerr + (y[n] - fn)^2.0
            SSxx = SSxx + (x[n] - xbar)^2.0
          endfor
          SSerr_mat[i,j] = SSerr ;added for testing error behavior
          SStot_mat[i,j] = SStot
          sigma_mat[i,j] = sqrt(SSerr/(size_CD[3]-1.0))/sqrt(SSxx)
          rsquare[i,j] = 1.0 - SSerr/SStot
        endif else begin
        ;Processing by taking average of each calculated delta B
          B = dblarr(size_CD[3])
          for n = 0, size_CD[3]-1 do begin
            B[n] = y[n]/x[n]
          endfor
          db[i,j] = mean(B)
          RSD[i,j] = abs(stddev(B)/db[i,j])*100.0
          
        endelse
      
      endif
    endfor
  endfor

output1 = dB
if info.ls_flag eq 1 then begin
  output2 = rsquare
  
endif else begin
  output2 = RSD

endelse

project.deltaB_result = ptr_new(output1) 
project.deltaB_r2 = ptr_new(output2)
project.deltaB_sigma = ptr_new(sigma_mat)

end

pro cdbls_phase_display, info
  common scan_data
  dTfp = *info.dTfp_array
  CD_array = *info.CD_array
  size_CD = size(CD_array)
  phase_array = dblarr(size_CD[1]*size_CD[3], size_CD[2])
  if info.unwrap_flag eq 1 then begin
    setenv, 'FSLOUTPUTTYPE=NIFTI_GZ'
    spawn, '/usr/local/fsl/bin/prelude', errMsg
  endif  else errMsg = 'NoUnwrap'
  for delay=0, size_CD[3]-1 do begin
    phase = atan(reform(CD_array[*,*,delay]), /phase)
      if errMsg eq '' and !VERSION.OS_FAMILY ne 'Windows' then begin
        phase_ptr = ptr_new(phase)
        mas_export_nifti,data_ptr = info.magnitude_image, file_name = project.current_path + '/mag_data.nii'
        mas_export_nifti,data_ptr = phase_ptr, file_name = project.current_path + '/phase_data.nii'
        spawn, '/usr/local/fsl/bin/prelude -a mag_data.nii -p phase_data.nii -u phase_data_unwrap.nii -s'
        temp = mas_read_nifti(nifti_filename = project.current_path +'/phase_data_unwrap.nii.gz')
        phase = *temp.voxel_data
        file_delete, project.current_path + '/phase_data.nii',project.current_path + '/mag_data.nii', project.current_path +'/phase_data_unwrap.nii.gz'
        ptr_free, phase_ptr
      endif else begin ; Fourier based method
        phase_image = phase
        phase_unwrap_2D,phase_image
        phase = phase_image
      endelse 
    if info.mask_flag eq 1 then begin
      mask = *info.mask
      phase_array[delay*size_CD[1]:(delay+1)*size_CD[1]-1 ,*] = phase*mask
    endif else begin
      phase_array[delay*size_CD[1]:(delay+1)*size_CD[1]-1 ,*] = phase
    endelse
  endfor
  iimage, phase_array
  end

function cdbls_create_scan_list_array, CI_arr, num_scans, max_scans
  common scan_data
  index = where(CI_arr gt -1)
  size_CI = size(index)
  scan_list_array = strarr(size_CI[1])
  for i=0,num_scans-1 do begin
    ;Get scan name
    CI = CI_arr[i]
    scan_name_temp = project.scan_list[CI]
    scan_name_split = strsplit(scan_name_temp,PATH_SEP(),/EXTRACT)
    sz_split = size(scan_name_split)
    scan_list_array[i] = PATH_SEP() + scan_name_split[sz_split[1]-2] + PATH_SEP() + scan_name_split[sz_split[1]-1] + PATH_SEP()
  end
return, scan_list_array
end



pro mas_calc_delta_b_least_squares_gui_event, event

 ;bring in global variables
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets
    COMMON common_cdbls_widgets, $
           cdbls_scanlist, cdbls_addscan, cdbls_removescan, cdbls_threshold_text, cdbls_display_button, $
           cdbls_phase_button, cdbls_cdb_gr_text, cdbls_ls_button, cdbls_avg_button, cdbls_printtfp_button, $
           cdbls_export_row_text, cdbls_export_column_text, cdbls_export_button, cdbls_filter_toggle, $
           cdbls_chi_line_data_numsteps_text, cdbls_chi_line_data_stopcondition_text, cdbls_chi_line_data_crop_text, $
           cdbls_chi_line_data_media_text, cdbls_chi_line_data_start_text, cdbls_chi_line_data_stop_text, $
           cdbls_chi_line_data_action_display, cdbls_roi_set_label, cdbls_roi_set_droplist, cdbls_roi_label, $
           cdbls_roi_droplist, cdbls_roi_refresh_button, cdbls_chi_line_data_column_text, cdbls_chi_line_data_row_text, $  
           cdbls_chi_line_data_field_text, cdbls_plot_line_data, cdbls_display_dB_maps, cdbls_model_interp_droplist, cdbls_1st_order_toggle, $
           cdbls_1st_order_text, cdbls_chi_line_fittype_linedata, cdbls_chi_line_fittype_image, cdbls_chi_line_data_window_text, $
           cdbls_chi_line_data_window_toggle, cdbls_uncertainty_text, cdbls_uncertainty_toggle, cdbls_unwrap_toggle       
    
    ;grab parameter structure that is stored in the uvalue of the mas_multi_movie_gui
    widget_control, event.top, get_uvalue=info, /no_copy
    CI = project.CI
    
    CASE event.id OF
    
      cdbls_scanlist:BEGIN
        info.SI = event.index
      END
      
      cdbls_addscan:BEGIN
        check_add = where(info.CI eq CI)
        if check_add eq -1 then begin
          info.CI[info.num_scans] = CI
          info.sensitive_flag = 1
          info.num_scans = info.num_scans + 1
          max_scans = info.max_scans
          scan_list_array = cdbls_create_scan_list_array(info.CI, info.num_scans, max_scans)
          widget_control, cdbls_scanlist, set_list_select = info.num_scans-1, set_value = scan_list_array, sensitive = info.sensitive_flag
          info.SI = info.num_scans-1
        endif 
        
        if info.num_scans ge 1 then widget_control, cdbls_printtfp_button, sensitive = 1
        
        if info.num_scans ge 3 then begin
          widget_control, cdbls_display_button, sensitive = 1
          widget_control, cdbls_phase_button, sensitive = 1
          widget_control, cdbls_export_button, sensitive = 1
          widget_control, cdbls_chi_line_data_action_display, sensitive = 1
        endif else begin
          widget_control, cdbls_display_button, sensitive = 0
          widget_control, cdbls_phase_button, sensitive = 0
          widget_control, cdbls_export_button, sensitive = 0
          widget_control, cdbls_chi_line_data_action_display, sensitive = 0
        endelse
      END
      
      cdbls_removescan:BEGIN
       ;Removing selected scan from CI list
        SI = info.SI
        CIarr = info.CI
        size_CI = size(CIarr)
        ;Remove first entry from CI list, shift all others over
        if SI eq 0 then begin
          CI_back = CIarr[SI+1:size_CI[1]-1]
          new_CI = [CI_back, -1]
        ;Remove last entry from CI list, add -1
        endif else if SI eq size_CI[1]-1 then begin
          CI_front = CIarr[0:SI-1]
          new_CI = [CI_front, -1]
        ;Remove SI from CI list, shift values to the right over, add -1 to end
        endif else begin
          CIarr_front = CIarr[0:SI-1]
          CIarr_back = CIArr[SI+1:size_CI[1]-1]
          new_CI = [CIarr_front, CIarr_back,-1]
        endelse
        info.CI = new_CI  
        info.num_scans = info.num_scans - 1
        ;updating Scan List Widget
        if info.num_scans eq 0 then begin
          widget_control, cdbls_scanlist, set_value = '', sensitive = 0
        endif else begin
          scan_list_array = cdbls_create_scan_list_array(info.CI, info.num_scans, max_scans)
          widget_control, cdbls_scanlist, set_list_select = SI-1, set_value = scan_list_array, sensitive = 1
        endelse
        
        if info.num_scans ge 1 then begin
          widget_control, cdbls_printtfp_button, sensitive = 1
        endif else begin
          widget_control, cdbls_printtfp_button, sensitive = 1
        endelse
        
        if info.num_scans ge 3 then begin
          widget_control, cdbls_display_button, sensitive = 1
          widget_control, cdbls_phase_button, sensitive = 1
          widget_control, cdbls_export_button, sensitive = 1
          widget_control, cdbls_chi_line_data_action_display, sensitive = 1
        endif else begin
          widget_control, cdbls_display_button, sensitive = 0
          widget_control, cdbls_phase_button, sensitive = 0
          widget_control, cdbls_export_button, sensitive = 0
          widget_control, cdbls_chi_line_data_action_display, sensitive = 0
        endelse
      END
      
      cdbls_filter_toggle:BEGIN
        info.mask_flag = event.select
        widget_control, cdbls_threshold_text, sensitive = info.mask_flag
      END
      
      cdbls_unwrap_toggle: info.unwrap_flag = event.select
      
      cdbls_threshold_text: BEGIN
        widget_control, cdbls_threshold_text, get_value = thresh_temp
        info.threshold = thresh_temp
      END
      
       cdbls_cdb_gr_text:BEGIN
        widget_control, cdbls_cdb_gr_text, get_value=gr_temp
        info.gyromagnetic_ratio = gr_temp
      END
      
       cdbls_ls_button:BEGIN
        info.ls_flag = 1
      END
      
       cdbls_avg_button:BEGIN
        info.ls_flag = 0
      END
      
        cdbls_chi_line_fittype_linedata:BEGIN
        info.fit_type = 0
      END
      
       cdbls_chi_line_fittype_image:BEGIN
        info.fit_type = 1
      END
      
      cdbls_printtfp_button:BEGIN
        scan_list_array = cdbls_create_scan_list_array(info.CI, info.num_scans, max_scans)
        CI_arr = info.CI
        index = where(CI_arr gt -1)
        size_CI = size(index)
        Tfp_array = dblarr(size_CI[1])
        for i=0, size_CI[1]-1 do begin  
          Tfp_array[i] = mas_extract_delays(project.imndarray[CI_arr[i]].file_path $
                                    , project.imndarray[CI_arr[i]].state1_load_procedure,te = project.imndarray[CI_arr[i]].echo_time)
        endfor
        sorted_index = sort(Tfp_array)
        Tfp_array = Tfp_array[sorted_index]
        scan_list_array = scan_list_array[sorted_index]
        
        display_string_array = strarr(size_CI[1] + 1)
        display_string_array[0] = 'Scan Name:' + string([09B]) + 'Tfp:'
        display_format = '(E0.4)'
        for i=0,size_CI[1]-1 do begin
          display_string_array[i+1] = scan_list_array[i] + string([09B]) + string(FORMAT=display_format,Tfp_array[i])
        endfor
        display_stats, display_string_array, 'Tfp Delays'
      END
      
      cdbls_display_button:BEGIN
        widget_control, cdbls_cdb_gr_text, get_value=gr_temp
        info.gyromagnetic_ratio = gr_temp
        cdbls_complex_division, info
        cdbls_deltaB_display, info, output1 = dB, output2 = Error
        B0 = project.imndArray[project.ci].spect_bf1/42.58
        if info.ls_flag eq 1 then begin
          iimage, 1e6*dB/B0, title = 'dB least squares (ppm)'
          iimage, Error, title = 'r squared'
        endif else begin
          
          iimage, 1e6*dB/B0, title = 'dB average (ppm)'
          iimage, Error, title = 'relative standard deviation'
        endelse  
            
      END
      
       cdbls_phase_button:BEGIN
        widget_control, cdbls_cdb_gr_text, get_value=gr_temp
        info.gyromagnetic_ratio = gr_temp
        if info.complex_division_flag eq 0 then begin
          cdbls_complex_division, info
        endif
        cdbls_phase_display, info
        ;mas_calc_delta_b_least_squares_unwrapped, info
        ;mas_calc_delta_b_least_squares, info        
      END
      
        cdbls_export_column_text:BEGIN
        widget_control, cdbls_export_column_text, get_value=column
        info.column_export = column
        END
        
        cdbls_export_row_text:BEGIN
        widget_control, cdbls_export_row_text, get_value=row
        info.row_export = row
        END
        
        cdbls_chi_line_data_numsteps_text:BEGIN
        widget_control, cdbls_chi_line_data_numsteps_text, get_value=numsteps
        info.num_steps = numsteps
        END
        
        cdbls_chi_line_data_stopcondition_text:BEGIN
        widget_control, cdbls_chi_line_data_stopcondition_text, get_value=stopcondition
        info.stopping_point = stopcondition
        END
        
        cdbls_chi_line_data_crop_text:BEGIN
        widget_control, cdbls_chi_line_data_crop_text, get_value=crop
        info.p_buffer = crop
        END
        
        cdbls_chi_line_data_media_text:BEGIN
        widget_control, cdbls_chi_line_data_media_text, get_value=mediachi
        info.media_chi = mediachi
        END
        
        cdbls_chi_line_data_start_text:BEGIN
        widget_control, cdbls_chi_line_data_start_text, get_value=startchi
        info.start_chi = startchi
        END
        
        cdbls_chi_line_data_stop_text:BEGIN
        widget_control, cdbls_chi_line_data_stop_text, get_value=stopchi
        info.stop_chi = stopchi
        END
        
        cdbls_1st_order_toggle:BEGIN
        info.field_correction = event.select
        widget_control, cdbls_1st_order_text, sensitive = info.field_correction
        END
        
        cdbls_chi_line_data_window_toggle:BEGIN
        info.rsquare_window_flag = event.select
        widget_control, cdbls_chi_line_data_window_text, sensitive = info.rsquare_window_flag
        END
        
        cdbls_1st_order_text:BEGIN
        widget_control, cdbls_1st_order_text, get_value=width
        info.field_correction_width = width
        END
        
        cdbls_roi_set_droplist: BEGIN
        ;we check to see which item in droplist is selected
          roi_set_index = widget_info(cdbls_roi_set_droplist, /DROPLIST_SELECT)
          IF roi_set_index EQ 0 THEN BEGIN
            info.roi_set=roi_set_index
            widget_control, cdbls_roi_droplist,set_value='      ', sensitive=0
          ENDIF ELSE BEGIN

          ;collect ROI names in the selected set
          rois = (*project.roi.pRois[roi_set_index])
          IF     (size(rois))[0] EQ 0 THEN num_roi = 1 $
          ELSE num_roi =  (size(rois))[1]
            roi_names = strarr(num_roi)
            FOR temp=0 , num_roi-1 DO BEGIN
              rois[temp] -> GETPROPERTY, name= name
              roi_names[temp] = name
            END
            info.roi_set=roi_set_index
            widget_control, cdbls_roi_droplist, set_value=roi_names, sensitive= 1
            info.roi_name = 0
          ENDELSE
         END
        
        cdbls_roi_droplist: BEGIN
        roi_index = widget_info(cdbls_roi_droplist, /DROPLIST_SELECT)
        info.roi_name = roi_index
        END
        
        cdbls_roi_refresh_button:BEGIN
          widget_control, cdbls_roi_set_droplist,set_value=*project.roi.pdisplay_names
          widget_control, cdbls_roi_droplist,set_value='      ', sensitive=0
          info.roi_name = -1
        END
        
        cdbls_chi_line_data_column_text:BEGIN
        widget_control, cdbls_chi_line_data_column_text, get_value=column
        info.chi_line_column = column
        END
        
        cdbls_chi_line_data_row_text:BEGIN
        widget_control, cdbls_chi_line_data_row_text, get_value=row
        info.chi_line_row = row
        END
        
        cdbls_chi_line_data_window_text:BEGIN
        widget_control, cdbls_chi_line_data_window_text, get_value=window
        info.chi_line_window = window
        END
        
        cdbls_chi_line_data_field_text:BEGIN
        widget_control, cdbls_chi_line_data_field_text, get_value=field
        info.field_strength = field
        END
        
        cdbls_plot_line_data:BEGIN
        info.display_line_data_flag = event.select
        END
        
        cdbls_display_dB_maps:BEGIN
        info.display_maps_flag = event.select
        END
        
        cdbls_model_interp_droplist: BEGIN
        interp_index = widget_info(cdbls_model_interp_droplist, /DROPLIST_SELECT)
        info.model_interp = double(interp_index + 1.0)
        END
        
        cdbls_uncertainty_toggle:BEGIN
        info.uncertainty_estimation_flag = event.select
        widget_control, cdbls_uncertainty_text, sensitive = info.uncertainty_estimation_flag
        END
        
        cdbls_uncertainty_text:BEGIN
        widget_control, cdbls_uncertainty_text, get_value=num_sims
        info.uncertainty_sims = num_sims
        END
        
        cdbls_chi_line_data_action_display:BEGIN
        ;Grabbing necessary parameters for delta B calculation 
        widget_control, cdbls_cdb_gr_text, get_value=gr_temp
        info.gyromagnetic_ratio = gr_temp
        widget_control, cdbls_chi_line_data_stopcondition_text, get_value=stopcondition
        info.stopping_point = stopcondition
        widget_control, cdbls_chi_line_data_numsteps_text, get_value=numsteps
        info.num_steps = numsteps
        widget_control, cdbls_chi_line_data_stop_text, get_value=stopchi
        info.stop_chi = stopchi
        widget_control, cdbls_chi_line_data_crop_text, get_value=crop
        info.p_buffer = crop
        widget_control, cdbls_chi_line_data_media_text, get_value=mediachi
        info.media_chi = mediachi
        widget_control, cdbls_chi_line_data_start_text, get_value=startchi
        info.start_chi = startchi
        widget_control, cdbls_chi_line_data_stop_text, get_value=stopchi
        info.stop_chi = stopchi
        widget_control, cdbls_chi_line_data_column_text, get_value=column
        info.chi_line_column = column
        widget_control, cdbls_chi_line_data_row_text, get_value=row
        info.chi_line_row = row
        widget_control, cdbls_uncertainty_text, get_value=sims
        info.uncertainty_sims = sims
        
        ;Calculating the B field perturbation map
        cdbls_complex_division, info
        cdbls_deltaB_display, info, output1 = dB, output2 = Error
        project.deltaB_result = ptr_new(dB)
        
        if info.uncertainty_estimation_flag eq 0 then begin
          if info.fit_type eq 0.0 then begin
            line_data_calc_relative_chi, info, output1 = relative_chi_struct
          endif else begin
            line_data_calc_relative_chi_image, info, output1 = relative_chi_struct
            dB = *project.deltaB_result
            dB_sim = *relative_chi_struct.simulated_map
            diff_map = dB - dB_sim
          endelse
          chi_vector = *relative_chi_struct.chi_vector
          r2_vector = *relative_chi_struct.r2_vector
          line_data_calc_relative_chi_print, info.num_steps, chi_vector, r2_vector
          B0 = project.imndArray[project.ci].spect_bf1/42.58
          if info.display_maps_flag EQ 1 then begin
            iimage, (*project.deltaB_result)*1E6/B0, title = 'Measured Map (ppm)'
            iimage, (*relative_chi_struct.simulated_map)*1E6/B0, title = 'Simulated Delta B (ppm)'
          endif
          
          if info.display_line_data_flag EQ 1 then begin
            measured_line = (*relative_chi_struct.measured_line)*1E6
            simulated_line = (*relative_chi_struct.simulated_line)*1E6
            line_data_calc_relative_chi_plot, info, measured_line, simulated_line
          endif  
         endif else begin
            num_sims = info.uncertainty_sims
            chi_vector_sim = dblarr(num_sims)
            r2_vector_sim = dblarr(num_sims)
            progressbar = Obj_New('progressbar', Color='red', Text='Estimating Uncertainty',/NOCANCEL)
            progressbar -> Start
            for s=0,num_sims-1 do begin
              dB = *project.deltaB_result
              dB_new = cdbls_uncertain_dB(dB)
              project.deltaB_result = ptr_new(dB_new)
              if info.fit_type eq 0.0 then begin
                line_data_calc_relative_chi, info, output1 = relative_chi_struct
                chi_vector_single = *relative_chi_struct.chi_vector
                r2_vector_single = *relative_chi_struct.r2_vector
                r2_vector_sim[s] = max(r2_vector_single)
                r2_index = (where(r2_vector_single eq r2_vector_sim[s]))[0]
                chi_vector_sim[s] = chi_vector_single[r2_index]
              endif else begin
                line_data_calc_relative_chi_image, info, output1 = relative_chi_struct
                chi_vector_single = *relative_chi_struct.chi_vector
                r2_vector_single = *relative_chi_struct.r2_vector
                r2_vector_sim[s] = max(r2_vector_single)
                r2_index = (where(r2_vector_single eq r2_vector_sim[s]))[0]
                chi_vector_sim[s] = chi_vector_single[r2_index]
              endelse
              project.deltaB_result = ptr_new(dB)
              progressBar -> Update, (float(s+1)/float(num_sims))*100.0
            endfor
            progressbar -> Destroy
         print, mean(chi_vector_sim)
         print, stdev(chi_vector_sim)    
         endelse
        END
                       
        cdbls_export_button:BEGIN
        widget_control, cdbls_export_row_text, get_value=row
        widget_control, cdbls_export_column_text, get_value=column
        widget_control, cdbls_cdb_gr_text, get_value=gr_temp
        info.row_export = row
        info.column_export = column
        info.gyromagnetic_ratio = gr_temp
        cdbls_complex_division, info
        cdbls_deltaB_display, info, output1 = dB, output2 = Error
        x = dB[*,row]
        y = dB[column,*]
        title = 'Write Row'
        write_dat_file, x, title
        title = 'Write Column'
        write_dat_file, y, title
           
      END
      
    ENDCASE

widget_control, event.top, set_uvalue=info  
end


pro mas_calc_delta_b_least_squares_gui
    HEAP_GC
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets
    COMMON common_cdbls_widgets
          
    CI = project.CI
    ;max_number of scans allowed for fits, could be increased if needed
    max_scans_allowed = 20
    cdbls_window_base = widget_base(column=1, TITLE='Calculate Delta B', xoffset=420, /base_align_center)
    
    cdbls_note_base = widget_base(cdbls_window_base, /base_align_left)
    cdbls_note_label = widget_label(cdbls_note_base, value='Note: Please use complex image type')
    
    cdbls_scanlist_base = widget_base(cdbls_window_base, /base_align_center)
    cdbls_scanlist = widget_list(cdbls_scanlist_base, uname = "cdbls_scanlist", frame = 1, xsize = 50, ysize = 5, sensitive = 0)
    
    cdbls_scanlist_edit_base = widget_base(cdbls_window_base, column = 2)
    cdbls_addscan = widget_button(cdbls_scanlist_edit_base, value = 'Add Selected Scan')
    cdbls_removescan = widget_button(cdbls_scanlist_edit_base, value = 'Remove Selected Scan')
    
    cdbls_filter_base = widget_base(cdbls_window_base, row = 2)
    cdbls_filter_toggle_base = widget_base(cdbls_filter_base, /nonexclusive,row = 1)
    cdbls_filter_toggle = widget_button(cdbls_filter_toggle_base, value = 'Mask')
    cdbls_unwrap_toggle = widget_button(cdbls_filter_toggle_base, value = 'Enable Phase Unwrapping')
    cdbls_threshold_base = widget_base(cdbls_filter_base, row = 1, /base_align_center)
    cdbls_threshold_label = widget_label(cdbls_threshold_base, value = 'Threshold:')
    cdbls_threshold_text = widget_text(cdbls_threshold_base, value='0.05', xsize=5, /editable, sensitive = 0)   
    
    ;cdbls_cdb_base = widget_base(cdbls_window_base, row =1, /base_align_center )
    cdbls_cdb_gr_label = widget_label(cdbls_threshold_base, value='Gamma (MHz/T):')
    cdbls_cdb_gr_text = widget_text(cdbls_threshold_base, value='42.58', editable=0, xsize = 5)
    
    cdbls_fittype_base = widget_base(cdbls_window_base, /exclusive, row=1)
    cdbls_ls_button = widget_button(cdbls_fittype_base, value = 'Least Squares', uname = 'cdbls_ls_button')
    cdbls_avg_button = widget_button(cdbls_fittype_base, value = 'Average', uname = 'cdbls_ls_button')
    widget_control, cdbls_ls_button, set_button = 1
    
    cdbls_action_base = widget_base(cdbls_window_base, column = 3)
    cdbls_printtfp_button = widget_button(cdbls_action_base, value = 'Print Tfp', sensitive = 0)
    cdbls_phase_button = widget_button(cdbls_action_base, value = 'Display Phase', sensitive = 0)
    cdbls_display_button = widget_button(cdbls_action_base, value = 'Display Delta B', sensitive = 0)
    
    cdbls_export_base = widget_base(cdbls_window_base, row = 1)
    cdbls_export_column_label = widget_label(cdbls_export_base, value = 'Column to export:')
    cdbls_export_column_text = widget_text(cdbls_export_base, value = '0', xsize = 3, uname = 'cdbls_export_column_text', editable = 1)
    cdbls_export_row_label = widget_label(cdbls_export_base, value = 'Row to export:')
    cdbls_export_row_text = widget_text(cdbls_export_base, value = '0', xsize = 3, uname = 'cdbls_export_row_text', editable = 1)
    ;cdbls_export_action_base = widget_base(cdbls_window_base, column=1)
    cdbls_export_button = widget_button(cdbls_export_base, value = 'Export', uname = 'cdbls_export_button', sensitive=0)
 
    cdbls_chi_line_data_note_base = widget_base(cdbls_window_base, column = 1, /base_align_center)
    cdbls_chi_line_data_note = widget_label(cdbls_chi_line_data_note_base, value = 'Estimate Chiv from line data')
    
    cdbls_chi_line_data_base = widget_base(cdbls_window_base, row=12, frame = 2)
    
    cdbls_chi_line_data_parameter_base_a = widget_base(cdbls_chi_line_data_base, row=1, /base_align_center)
    cdbls_chi_line_data_numsteps_label = widget_label(cdbls_chi_line_data_parameter_base_a, value = '# of Steps:')
    cdbls_chi_line_data_numsteps_text = widget_text(cdbls_chi_line_data_parameter_base_a, uname = 'cdbls_chi_line_data_numsteps_text', value = '20', xsize = 4, /editable)
    cdbls_chi_line_data_stopcondition_label = widget_label(cdbls_chi_line_data_parameter_base_a, value = 'Stop condition:')
    cdbls_chi_line_data_stopcondition_text = widget_text(cdbls_chi_line_data_parameter_base_a, uname = 'cdbls_chi_line_data_stopcondition_text', value = '0.0005', xsize = 8, /editable)
    
    cdbls_chi_line_data_parameter_base_b = widget_base(cdbls_chi_line_data_base, row = 1, /base_align_center)
    cdbls_chi_line_data_crop_label = widget_label(cdbls_chi_line_data_parameter_base_b, value = 'Crop size (pixels):')
    cdbls_chi_line_data_crop_text = widget_text(cdbls_chi_line_data_parameter_base_b, uname = 'cdbls_chi_line_data_crop_text', value = '3', xsize=3, /editable)
    cdbls_chi_line_data_media_label = widget_label(cdbls_chi_line_data_parameter_base_b, value = 'Media Chi:')
    cdbls_chi_line_data_media_text = widget_text(cdbls_chi_line_data_parameter_base_b, uname = 'cdbls_chi_line_data_media_text', value = '-9.05E-6', xsize=9, /editable)
    
    cdbls_chi_line_data_parameter_base_c = widget_base(cdbls_chi_line_data_base, row = 1, /base_align_center)
    cdbls_chi_line_data_start_label = widget_label(cdbls_chi_line_data_parameter_base_c, value = 'Start Chi:')
    cdbls_chi_line_data_start_text = widget_text(cdbls_chi_line_data_parameter_base_c, uname = 'cdbls_chi_line_data_start_text', value = '-15.0E-6', xsize=9, /editable)
    cdbls_chi_line_data_stop_label = widget_label(cdbls_chi_line_data_parameter_base_c, value = 'Stop Chi:')
    cdbls_chi_line_data_stop_text = widget_text(cdbls_chi_line_data_parameter_base_c, uname = 'cdbls_chi_line_data_stop_text',  value = '-3.0E-6', xsize=9, /editable)
    
    cdbls_chi_line_data_parameter_base_d = widget_base(cdbls_chi_line_data_base, row = 1, /base_align_center)
    cdbls_roi_set_label = widget_label(cdbls_chi_line_data_parameter_base_d, value = 'Select Geo ROI Set')
    cdbls_roi_set_droplist = widget_droplist(cdbls_chi_line_data_parameter_base_d,value = *project.roi.pdisplay_names, uname='cdbls_roi_set_droplist')
    cdbls_roi_label = widget_label(cdbls_chi_line_data_parameter_base_d, value='Geo ROI')
    cdbls_roi_droplist = widget_droplist(cdbls_chi_line_data_parameter_base_d,value='      ', uvalue=' ', sensitive=0, /dynamic_resize)
    
    cdbls_chi_line_data_parameter_base_e = widget_base(cdbls_chi_line_data_base, row = 1, /base_align_center)
    cdbls_roi_refresh_button = widget_button(cdbls_chi_line_data_parameter_base_e, value = 'Refresh ROI List', uname = 'cdbls_roi_refresh_button')
    cdbls_chi_line_data_field_label = widget_label(cdbls_chi_line_data_parameter_base_e, value = 'Field (T):')
    cdbls_chi_line_data_field_text = widget_text(cdbls_chi_line_data_parameter_base_e, uname = 'cdbls_chi_line_data_field_text',  value = '11.1', xsize=5, /editable)
    cdbls_chi_line_data_parameter_base_e1 = widget_base(cdbls_chi_line_data_parameter_base_e, row = 1, /nonexclusive)
    cdbls_chi_line_data_window_toggle = widget_button(cdbls_chi_line_data_parameter_base_e1, value = 'Rsquare Window')
    
    cdbls_chi_line_data_parameter_base_e2 = widget_base(cdbls_chi_line_data_base, row = 1, /base_align_center)
    cdbls_chi_line_data_column_label = widget_label(cdbls_chi_line_data_parameter_base_e2, value = 'Column for fit:')
    cdbls_chi_line_data_column_text = widget_text(cdbls_chi_line_data_parameter_base_e2, uname = 'cdbls_chi_line_data_column_text',  value = '58', xsize=3, /editable)
    cdbls_chi_line_data_row_label = widget_label(cdbls_chi_line_data_parameter_base_e2, value = 'Row for fit:')
    cdbls_chi_line_data_row_text = widget_text(cdbls_chi_line_data_parameter_base_e2, uname = 'cdbls_chi_line_data_row_text',  value = '60', xsize=3, /editable)
    cdbls_chi_line_data_window_label = widget_label(cdbls_chi_line_data_parameter_base_e2, value = 'Rsquare Window:')
    cdbls_chi_line_data_window_text = widget_text(cdbls_chi_line_data_parameter_base_e2, uname = 'cdbls_chi_line_data_window_text',  value = '64', xsize=3, /editable, sensitive = 0)
    
    cdbls_chi_line_data_parameter_base_f = widget_base(cdbls_chi_line_data_base, row = 1, /base_align_center)
    cdbls_chi_line_fittype_label = widget_label(cdbls_chi_line_data_parameter_base_f, value = 'Select Fit Type:')
    cdbls_chi_line_data_parameter_base_f2 = widget_base(cdbls_chi_line_data_parameter_base_f, row = 1, /base_align_center, /exclusive)
    cdbls_chi_line_fittype_linedata = widget_button(cdbls_chi_line_data_parameter_base_f2, value = 'Line Data', uvalue = 'cdbls_chi_line_fittype_linedata') 
    cdbls_chi_line_fittype_image = widget_button(cdbls_chi_line_data_parameter_base_f2, value = 'Image', uvalue = 'cdbls_chi_line_fittype_image')
    widget_control, cdbls_chi_line_fittype_linedata, set_button = 1
    
    cdbls_chi_line_data_parameter_base_g = widget_base(cdbls_chi_line_data_base, row =1, /nonexclusive)
    cdbls_plot_line_data = widget_button(cdbls_chi_line_data_parameter_base_g, value = 'Display Line Data')
    cdbls_display_dB_maps = widget_button(cdbls_chi_line_data_parameter_base_g, value = 'Display Maps')
    
    cdbls_chi_line_data_parameter_base_h = widget_base(cdbls_chi_line_data_base, row =1)
    cdbls_model_interp_label = widget_label(cdbls_chi_line_data_parameter_base_h, value = 'Model Interpolation Factor:')
    cdbls_model_interp_droplist = widget_droplist(cdbls_chi_line_data_parameter_base_h, uname = 'cdbls_model_interp_droplist', value = ['1', '2', '3', '4','8','16'])
    
    cdbls_chi_line_data_parameter_base_i = widget_base(cdbls_chi_line_data_base, row =1)
    cdbls_chi_line_data_parameter_base_i2 = widget_base(cdbls_chi_line_data_parameter_base_i, row =1,/nonexclusive)
    cdbls_1st_order_toggle = widget_button(cdbls_chi_line_data_parameter_base_i2, value = '1st Order Field Correction')
    cdbls_1st_order_label = widget_label(cdbls_chi_line_data_parameter_base_i, value = '# of Pixels:')
    cdbls_1st_order_text = widget_text(cdbls_chi_line_data_parameter_base_i, uname = 'cdbls_1st_order_text', value = '10',xsize = 3, /editable, sensitive = 1)
    
    cdbls_chi_line_data_parameter_base_j = widget_base(cdbls_chi_line_data_base, row =1)
    cdbls_chi_line_data_parameter_base_j2 = widget_base(cdbls_chi_line_data_parameter_base_j, row =1,/nonexclusive)
    cdbls_uncertainty_toggle = widget_button(cdbls_chi_line_data_parameter_base_j2, value = 'Uncertainty Estimation')
    cdbls_uncertainty_label = widget_label(cdbls_chi_line_data_parameter_base_j, value = '# of Simulations:')
    cdbls_uncertainty_text = widget_text(cdbls_chi_line_data_parameter_base_j, uname = 'cdbls_uncertainty_text', value = '10',xsize = 3, /editable, sensitive = 0)
  
    cdbls_chi_line_data_action_base = widget_base(cdbls_chi_line_data_base, row = 1, /base_align_center)
    cdbls_chi_line_data_action_display = widget_button(cdbls_chi_line_data_action_base, value = 'Display Chi', uname = 'cdbls_chi_line_data_action_display', sensitive=0)
    widget_control, cdbls_window_base, /realize
    
    info = {SI:0, CI: intarr(max_scans_allowed)-1, num_scans:0, max_scans:max_scans_allowed $
    , sensitive_flag:0,  unwrap_flag:0, gyromagnetic_ratio:42.58 ,mask_flag: 0, threshold: 0.05 $
    ,complex_division_flag:0,slice:0, Tfp_array: ptr_new(), dTfp_array: ptr_new(), CD_array:ptr_new(), magnitude_image:ptr_new() $
    ,mask:ptr_new(), ls_flag:1, column_export:0, row_export:0, num_steps:20, stopping_point: 0.005, p_buffer: 3.0 $
    , media_chi:-9.05E-6, start_chi: -15.0E-6, stop_chi: -3.0E-6, chi_line_column:58, chi_line_row:60, chi_line_window: 64, roi_set:0, roi_name:0 $
    , field_strength:11.1,measured_line_corrected:ptr_new(), display_line_data_flag:0, display_maps_flag:0, model_interp:1.0 $
    , field_correction:1, field_correction_width:10.0, fit_type:0, rsquare_window_flag:0, uncertainty_estimation_flag:0, uncertainty_sims:10}
 
 
 
  widget_control, cdbls_window_base, set_uvalue=info, /no_copy
  widget_control, cdbls_1st_order_toggle, /set_button

    xmanager, 'mas_calc_delta_b_least_squares_gui', cdbls_window_base

end