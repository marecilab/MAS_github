;+
; :Description:
;    Constructor. Initializes the plot using the fid_data 
;    argument.
;
; :Params:
;    fid_data - complex data. 
;
; :Keywords:
;    xlim - the extreme value of the x-axis
;    zero - the initial zero point. NOTE that this is the INDEX into the
;           array of the element that is located at zero. For example,
;           if fid_data = [ -1, 10, 5, -32, 1 ], and zero=3, then the
;           x axis will designate zero at elemnt 3, having value -32:
;           -1     10     5     -32     1 (data values)
;           |------|------|------|------|
;           -3    -2     -1      0      1 (axis values)
;           |------|------|------|------|
;           0      1      2      3      4 (array indices)
;                   
;    initial_fid - for multi-dimensional data, this is the initial dimension
;                  to display
;    fq_window - this should be a reference to an idlgrwindow object. If 
;                not specified one will be created
;    reverse_xaxis - set this keyword to have the axis reverse
;
; :Author: btt
;-
function mas_complex_plot::init, fid_data, $
    xlim=xlim, zero=zero, initial_fid=initial_fid, $
    fq_window=fq_window, reverse_xaxis=reverse_xaxis

    size_fid = size(fid_data)
    self.x_lim = keyword_set(xlim) ? xlim : 25000.
    
    color_red = byte([255,0,0])      ;; abs
    color_blue = byte([0,0,255])     ;; real
    color_orange = byte([255,120,0]) ;; imag
    
    self.plot_color_abs  = color_red
    self.plot_color_imag = color_orange
    self.plot_color_real = color_blue
    
    if (keyword_set(reverse_xaxis)) then begin
        self.fq_reverse_xaxis = 1
    endif else begin
        self.fq_reverse_xaxis = 0
    endelse
    ;;self.fq_reverse_xaxis = 0
    
    case size_fid[0] of
    
        0: message, 'No data to analyize.'
        
        1: begin ;; data has one fid
            self.raw_fid = ptr_new(fid_data)
            self.num_fidpts = size_fid[1]
            self.initial_fid = 0
            self.num_fids = 1
        end
        
        2: begin ;; data has more than one fid
            self.raw_fid = ptr_new(fid_data)
            self.num_fids = size_fid[2]
            self.num_fidpts = size_fid[1]
            self.initial_fid = 0
            if (keyword_set(initial_fid)) then begin
                self.initial_fid = (long(initial_fid) > 0) < (self.num_fids-1)
            endif
        end
        
        else: message, 'Data has too many dimensions: '+string(size_fid[0], format='(G0)')
     
     endcase
     
     self.fq_zero_ref_index = keyword_set(zero) ? zero : 0
     self.fq_xaxis_pts = ptr_new(self->generateFreqXRange())
     self.fq_data_index_start = 0
     self.fq_data_index_end = n_elements(*self.fq_xaxis_pts)-1
     
     if (not keyword_set(fq_window)) then begin
         self.ofq_window = obj_new('idlgrwindow', dimensions=[850,300])
     endif else begin
         self.ofq_window = fq_window
     endelse
     
     self.ofq_xaxis_font = obj_new('idlgrfont', size=10)
     
     self->initFrequencyPlot, data=fid_data
     self->regenerateXAxis
     self->refreshFrequencyPlot
     
     return, self
     
end

;+
; :Description:
;    Destructor, called automatically by IDL
;
; :Author: btt
;-
pro mas_complex_plot::cleanup

    message, "Cleaning up...", /informational
    
    obj_destroy, self.ofq_plot_real
    obj_destroy, self.ofq_plot_imag
    obj_destroy, self.ofq_plot_abs
    obj_destroy, self.ofq_plot_sym_real
    obj_destroy, self.ofq_plot_sym_imag
    obj_destroy, self.ofq_plot_sym_abs
    obj_destroy, self.ofq_plot_model    
    obj_destroy, self.ofq_plot_view

    obj_destroy, self.ofq_xaxis
    obj_destroy, self.ofq_yaxis
    obj_destroy, self.ofq_zero_line
    obj_destroy, self.ofq_xaxis_model
    obj_destroy, self.ofq_interval_trendline
    obj_destroy, self.ofq_xaxis_view
    obj_destroy, self.ofq_xaxis_font

    obj_destroy, self.ofq_cursor_01
    obj_destroy, self.ofq_cursor_02
    obj_destroy, self.ofq_cursor_model
    obj_destroy, self.ofq_cursor_view

    obj_destroy, self.ofq_headsup_message_text
    obj_destroy, self.ofq_headsup_model
    obj_destroy, self.ofq_headsup_view
    obj_destroy, self.ofq_scene

    ptr_free, self.raw_fid
    ptr_free, self.fq_xaxis_pts
    ptr_free, self.fq_xaxis_pts_disp
    ptr_free, self.fq_data
    ptr_free, self.fq_data_displayed

end

;+
; :Description:
;   given the frequency (x-axis) value, this function
;   provides the index into the data array corresponding
;   to this frequence. 
;
; :Params:
;    freq frequency to retrieve
;
; :Author: btt
;-
function mas_complex_plot::convertFreqToDataIndex, freq
    if (0 and self.fq_reverse_xaxis eq 1) then begin
        return, ( (self.x_lim-freq)/self.x_lim * self.num_fidpts ) + self.fq_zero_ref_index        
    endif else begin
        return, (freq/self.x_lim * self.num_fidpts) + self.fq_zero_ref_index
    endelse
end

;+
; :Description:
;   given the frequency (x-axis) value, this function
;   provides the x pixel location in the display window
;
; :Params:
;    freq frequency to retrieve
;
; :Author: btt
;-
function mas_complex_plot::convertFreqToWindowX, freq
    self.ofq_window->getProperty, dimensions=dim
    return, freq/self.x_lim * dim[0]
end

;+
; :Description:
;    This function generates an xaxis array for the 
;    data given the current setup (xlim, zero_ref, etc)
;
;    Note that this axis is not reversed... rather it is
;    used internally to position elements in the ploy view.
; :Author: btt
;-
function mas_complex_plot::generateFreqXRange
    
    new_xrange = (findgen(self.num_fidpts)/self.num_fidpts - float(self.fq_zero_ref_index)/self.num_fidpts) * self.x_lim
    if (0 and self.fq_reverse_xaxis eq 1) then begin
        return, reverse(new_xrange)
    endif else begin
        return, new_xrange
    endelse
    
end

;+
; :Description:
;    Sets the x-axis limit to a new value.
;
; :Params:
;    val - the new xaxis limit value
;
; :Author: btt
;-
pro mas_complex_plot::setXAxislimit, val

    self.x_lim = float(val)
    self->regenerateXAxis
end


;+
; :Description:
;    This procedure regenerates the xaxis that is displayed
;    in the window. It takes into account the zero reference
;    the zoom position and the axis range. 
;    
;    This takes care of drawing the axis on the screen as well.
;    
;    Note that if you extend this object, you will need to
;    determine whether or not you need to regenerate the
;    axis depending on what you do.
;
; :Author: btt
;-
pro mas_complex_plot::regenerateXAxis
    
    axis_pts = findgen(self.num_fidpts+1)
    if (self.fq_reverse_xaxis) then begin
        axis_pts = reverse(axis_pts)
    endif
        
    ;; scale the axis to x_lim
    axis_pts = (axis_pts/self.num_fidpts) * self.x_lim

    ;; adjust the zero position
    if (self.fq_zero_ref_index ne 0) then begin
        if (self.fq_reverse_xaxis) then begin
            axis_pts -= self.x_lim - float(self.fq_zero_ref_index)/self.num_fidpts*self.x_lim  ;;axis_pts[self.fq_zero_ref_index]
        endif else begin
            axis_pts -= float(self.fq_zero_ref_index)/self.num_fidpts*self.x_lim  ;;axis_pts[self.fq_zero_ref_index]
        endelse
        
    endif

    ptr_free, self.fq_xaxis_pts_disp
    self.fq_xaxis_pts_disp = ptr_new(axis_pts)
    
    ;; find the min/max values for displaying in the window
    if (obj_valid(self.ofq_cursor_view)) then begin
        self.ofq_cursor_view->getProperty, viewplane_rect=vpr
        left  = vpr[0]
        right = (vpr[0] + vpr[2]) < (self.num_fidpts)       
    endif else begin
        left  = 0L
        right = self.num_fidpts
        return
    endelse
    
    ;; extract the portion of the axis to be displayed
    axis_portion = axis_pts[left:right]
    num_axis_pts = n_elements(axis_portion)
    
    ;; we put ticks at every "ti" distance
    ti = round(num_axis_pts * 0.025) ;;& print, ti
    om = self.ofq_cursor_model->getByName('THEMODEL')
    if (obj_valid(om)) then begin
        self.ofq_cursor_model->remove, om
        obj_destroy, om
    endif
    
    tick_height = vpr[3]*0.04
    labl_height = (tick_height * (1.0 + 0.2))
    om = obj_new('idlgrmodel', name='THEMODEL')
    om->translate, 0,0,0
    
    self.ofq_cursor_model->add, om
    ;;om->add, obj_new('idlgrpolyline', [  [0,0], [vpr[2], 0]  ])
    om->add, obj_new('idlgrpolyline', [  [vpr[0], tick_height], [vpr[0]+vpr[2], tick_height]  ])
    
    n = 0L
    
    for p = 0, n_elements(axis_portion)-1 do begin
    
        if (p mod ti eq 0) then begin
        
            if (n mod 5 eq 0) then begin
            
                ;; a major tick (higher and labelled)
                om->add, obj_new('idlgrpolyline', [[left+p, 0], [left+p, tick_height]])
                
                om->add, obj_new('idlgrtext', string(axis_portion[p], format='(F0.2)'), $
                        alignment=0.5, font=self.ofq_xaxis_font, recompute_dimensions=2, $
                        locations=[left+p, labl_height])
            
            endif else begin
            
                ;;; minor tick
                om->add, obj_new('idlgrpolyline', [[left+p, tick_height], [left+p, vpr[3]*0.025]])
                
                
            endelse 
            
            n++
            
        endif
    
    endfor
    
    
end

;+
; :Description:
;    For arrayed data, this chooses an array
;    element and refreshed the plot with the data
;    from that element
;
; :Params:
;    element
;
; :Author: btt
;-
pro mas_complex_plot::selectDataArrayElement, element


    if (element lt 0 or element gt (self.num_fids-1)) then begin
        ;; out of range
        return
    endif
    
    self.initial_fid = element
    
    data = reform((*self.fq_data)[*,element])
    
    self->setDisplayedData, data

end

;+
; :Description:
;    This swaps the currently displayed data with
;    the data argument
;    
;    See also: setData
;    
; :Params:
;    data
;
; :Author: btt
;-
pro mas_complex_plot::setDisplayedData, data

    if (n_elements(data) eq 0) then return
    
    ptr_free, self.fq_data_displayed
    
    self.fq_data_displayed = ptr_new(data)

    self.ofq_plot_real->setProperty, datay=real_part(*self.fq_data_displayed)
    self.ofq_plot_imag->setProperty, datay=imaginary(*self.fq_data_displayed)
    self.ofq_plot_abs->setProperty,  datay=abs(*self.fq_data_displayed)
    
    self->refreshFrequencyPlot

end

;+
; :Description:
;    This completely changes all of the data 
;    that is plotted.
;
; :Params:
;    data

; :Author: btt
;-
pro mas_complex_plot::setData, data

    ptr_free, self.fq_data
    ptr_free, self.fq_data_displayed
    ptr_free, self.raw_fid
    
    self.raw_fid = ptr_new(data)
    self.fq_data = ptr_new(data)
    self.fq_data_displayed = ptr_new(reform(data[*,self.initial_fid]))

    self.ofq_plot_real->setProperty, datay=real_part(*self.fq_data_displayed)
    self.ofq_plot_imag->setProperty, datay=imaginary(*self.fq_data_displayed)
    self.ofq_plot_abs->setProperty,  datay=abs(*self.fq_data_displayed)

    self->refreshFrequencyPlot
    
end

;+
; :Description:
;    This procedure handles all aspects of initial plot
;    creation. It creates all of the graphics objects,
;    and divines most of the parameters.
;
; :Keywords:
;    data -- the data to be plotted.
;
; :Author: btt
;-
pro mas_complex_plot::initFrequencyPlot, data=data

     ptr_free, self.raw_fid
     ptr_free, self.fq_data
     ptr_free, self.fq_data_displayed
     
     self.raw_fid = ptr_new(data)
     self.fq_data = ptr_new(data)
     self.fq_data_displayed = ptr_new(reform(data[*,self.initial_fid]))
     
     ;; the main scene
     self.ofq_scene      = obj_new('idlgrscene')
     
     fq_min = min([min(real_part(*self.fq_data_displayed)), $
                   min(imaginary(*self.fq_data_displayed)), $
                   min(abs(*self.fq_data_displayed))])

     fq_max = max([max(real_part(*self.fq_data_displayed)), $
                   max(imaginary(*self.fq_data_displayed)), $
                   max(abs(*self.fq_data_displayed))])
     
     self.fq_yoffset = 0
     self.fq_ymin = fq_min
     self.fq_ymax = fq_max
     self.fq_yscale = abs(fq_min - fq_max) ;abs(fq_min)+abs(fq_max)
     
     self.ofq_window->getProperty, dimensions=wdim
     
     ;; these are hidden to start.
     self.ofq_plot_sym_real = obj_new('idlgrsymbol', 0, size=[3/wdim[0]*self.x_lim, 3/wdim[1]*self.fq_yscale])
     self.ofq_plot_sym_imag = obj_new('idlgrsymbol', 0, size=[3/wdim[0]*self.x_lim, 3/wdim[1]*self.fq_yscale])
     self.ofq_plot_sym_abs  = obj_new('idlgrsymbol', 0, size=[3/wdim[0]*self.x_lim, 3/wdim[1]*self.fq_yscale])

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

     ;; plot model.
     self.ofq_plot_model = obj_new('idlgrmodel')
     self.ofq_plot_view  = obj_new('idlgrview', /transparent)
     ;; this is in data space: 0 ... x_lim ; fq_min ... fq_scale
     self.ofq_plot_view->setProperty, viewplane_rect=[0, self.fq_ymin-self.fq_yoffset, self.x_lim, self.fq_yscale]

     self.ofq_plot_real = obj_new('idlgrplot', *self.fq_xaxis_pts, $
                                  color=self.plot_color_real, $
                                  real_part(*self.fq_data_displayed), $
                                  symbol=self.ofq_plot_sym_real)

     self.ofq_plot_imag = obj_new('idlgrplot', *self.fq_xaxis_pts, $
                                  color=self.plot_color_imag, $
                                  imaginary(*self.fq_data_displayed), $
                                  symbol=self.ofq_plot_sym_imag)

     self.ofq_plot_abs = obj_new('idlgrplot', *self.fq_xaxis_pts, $
                                 color=self.plot_color_abs, $
                                  abs(*self.fq_data_displayed), $
                                  symbol=self.ofq_plot_sym_abs)
                                  
     ;; this is a trend line that is fit to the data points between 
     ;; the cursors
     self.ofq_interval_trendline = obj_new('idlgrpolyline', hide=0, data=[[0,0], [0,self.x_lim]], $
                                           linestyle=1, uvalue=[0.,0.], thick=2, color=[0,0,0], $
                                           alpha_channel=0.7)
                                           
     self.ofq_plot_model->add, [ self.ofq_interval_trendline, self.ofq_plot_real, self.ofq_plot_imag, self.ofq_plot_abs ]
     self.ofq_plot_view->add, self.ofq_plot_model

     ;; this line indicates the zero 
     self.ofq_zero_line = obj_new('idlgrpolyline', thick=1, color=[0,0,0], $
                                  [ [0,self.fq_yoffset], [self.x_lim,self.fq_yoffset] ])
     self.ofq_plot_model->add,  self.ofq_zero_line

     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      
     self.ofq_xaxis_model = obj_new('idlgrmodel')    
     
     ;; axis view is in the same space as the ploy view.
     self.ofq_xaxis_view = obj_new('idlgrview', /transparent)
     self.ofq_xaxis_view->setProperty, viewplane_rect=[0, 0, wdim[0], wdim[1]]
     self.ofq_xaxis_view->add, self.ofq_xaxis_model
     
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
     ;; cursor view holds the vertical cursors. 
     ;; cursor view is in data array space: 0 ... num_fidpts, 0 .. 100
     ;; (the yscale is arbitrarly chosen)
     self.ofq_cursor_view = obj_new('idlgrview', /transparent)
     self.ofq_cursor_view->setProperty, viewplane_rect=[0, 0, self.num_fidpts, 100]
     self.ofq_cursor_model = obj_new('idlgrmodel')
     self.ofq_marker_model = obj_new('idlgrmodel')
     
     self.ofq_cursor_01 = obj_new('idlgrpolyline', $
                                  [ [self.num_fidpts/4, 0], [self.num_fidpts/4, 100] ],$
                                  color=[255,0,0])
     self.ofq_cursor_02 = obj_new('idlgrpolyline', $
                                  [ [self.num_fidpts*0.75, 0], [self.num_fidpts*0.75, 100] ],$
                                  color=[0,255,0])
     self.fq_cursor_01_pos = self.num_fidpts*0.25
     self.fq_cursor_02_pos = self.num_fidpts*0.75
      
     self.ofq_cursor_model->add, [self.ofq_marker_model, self.ofq_cursor_01, self.ofq_cursor_02]
     self.ofq_cursor_view->add, self.ofq_cursor_model
     
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     
     ;; heads up view is in window space: 0 ... window x, 0 ... window y
     self.ofq_headsup_view = obj_new('idlgrview', /transparent, viewplane_rect=[0,0,wdim[0], wdim[1]])
     self.ofq_headsup_model = obj_new('idlgrmodel')
     fnt = obj_new('idlgrfont', size=10)
     
     self.ofq_headsup_message_text = obj_new('idlgrtext', '', locations=[0,wdim[1]-12], font=fnt)
     self.ofq_headsup_model->add, self.ofq_headsup_message_text
     self.ofq_headsup_view->add, self.ofq_headsup_model

     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     self.ofq_scene->add, self.ofq_headsup_view
     self.ofq_scene->add, self.ofq_cursor_view
     self.ofq_scene->add, self.ofq_plot_view
     self.ofq_scene->add, self.ofq_xaxis_view
     
end

pro complex_plotter::setReverseXAxis, value
    
    if (value ne 0) then begin
        self.fq_reverse_xaxis = 1
    endif else begin
        self.fq_reverse_xaxis = 0
    endelse
    
    self->regenerateXAxis
    self->refreshFrequencyPlot

end
 
;+
; :Description:
;    sets the zero data index reference.
;    regenerates the x-axis accordingly.
;
; :Params:
;    data_index - the index into the data array
;                 of the zeroth element
;
;
;
; :Author: btt
;-
pro mas_complex_plot::setZeroReference, data_index

    self.fq_zero_ref_index = data_index
    self->regenerateXAxis
    self->refreshFrequencyPlot

end

;+
; :Description:
;    Sets the y-offset of the plot. Basically just a
;    translation of just the plot in the + or - y
;    direction.
;
; :Params:
;    amt - amount to translate
; :Author: btt
;-
pro mas_complex_plot::applyYoffset, amt
    
    self.fq_yoffset = amt
    
    self.ofq_window->getProperty, dimensions=wdim
    
    self.ofq_plot_view->getProperty, viewplane_rect=vpr    
    vpr[1] = self.fq_ymin-self.fq_yoffset
    self.ofq_plot_view->setProperty, viewplane_rect=vpr
    
    self->refreshFrequencyPlot

end

;+
; :Description:
;    Applys a y-scaling factor.
;
; :Params:
;    zoom_amt - an integer pos/neg, scale amount is 2^zoom_amt
;
; :Author: btt
;-
pro mas_complex_plot::applyYZoom, zoom_amt
    
    scale = float(zoom_amt);
    if (scale lt 0) then begin
        factor = 2^abs(scale)
    endif else if (scale gt 0) then begin
        factor = 1.0/2^abs(scale)
    endif else begin
        factor = 1.0
    endelse
    
    self.ofq_window->getProperty, dimensions=wdim
    
    self.ofq_plot_view->getProperty, viewplane_rect=vpr    
    vpr[1] = (factor)*(self.fq_ymin-self.fq_yoffset)
    vpr[3] = (factor)*self.fq_yscale
    self.ofq_plot_view->setProperty, viewplane_rect=vpr
    
    self.ofq_plot_sym_real->getProperty, size=sym_size
    self.ofq_plot_sym_real->setProperty, size=[sym_size[0], 3/wdim[1]*vpr[3]]
    
    self->refreshFrequencyPlot

end

;+
; :Description:
;    Fills the plot window with the region that lies between
;    the two cursors.
;
; :Params:
;    zoom_left  -- not used
;    zoom_right  -- not used
;
; :Keywords:
;    reset -- set this keyword to reset the plot to display the
;             entire plot
;
; :Author: btt
;-
pro mas_complex_plot::applyXZoom, zoom_left, zoom_right, reset=reset

    ptr_free, self.fq_xaxis_pts
    self.fq_xaxis_pts = ptr_new(self->generateFreqXRange())

    if (keyword_set(reset)) then begin
        range_start = (*self.fq_xaxis_pts)[0]
        range_end   = (*self.fq_xaxis_pts)[self.num_fidpts-1]
    endif else begin
        ind1 = self.fq_cursor_01_pos
        ind2 = self.fq_cursor_02_pos
        min_ind = min([ind1,ind2], max=max_ind)
        
        min_pos = (*self.fq_xaxis_pts)[(min_ind-4) > 0]
        max_pos = (*self.fq_xaxis_pts)[(max_ind+4) < (self.num_fidpts-1)]
        
        range_start = min_pos
        range_end   = max_pos
    endelse
    range_width = abs(range_end - range_start)
    
    start_index = round(self->convertFreqToDataIndex(range_start))
    start_index = (start_index > 0) < (self.num_fidpts-1)
    
    end_index   = round(self->convertFreqToDataIndex(range_end))
    end_index   = (end_index > 0) < (self.num_fidpts-1)
    
    self.fq_data_index_start = start_index
    self.fq_data_index_end   = end_index
    
    ;;;;;;;;;;;
    
    self.ofq_window->getProperty, dimensions=wdim
    
    self.ofq_plot_real->setProperty, datax=*self.fq_xaxis_pts
    self.ofq_plot_imag->setProperty, datax=*self.fq_xaxis_pts
    self.ofq_plot_abs->setProperty, datax=*self.fq_xaxis_pts
    
    self.ofq_plot_view->getProperty, viewplane_rect=vpr_pl
    vpr_pl[0] = range_start
    vpr_pl[2] = range_width
    self.ofq_plot_view->setProperty, viewplane_rect=vpr_pl
    
    self.ofq_xaxis_view->getProperty, viewplane_rect=vpr_ax
    vpr_ax[0] = range_start
    vpr_ax[2] = range_width
    
    self.ofq_zero_line->setProperty, data=[[range_start,0], [range_end, 0]]
    
    self.ofq_plot_sym_real->getProperty, size=sym_size
    self.ofq_plot_sym_real->setProperty, size=[3/wdim[0]*range_width, sym_size[1]]
    
    self.ofq_cursor_view->getProperty, viewplane_rect=vpr_cs
    vpr_cs[0] = start_index
    vpr_cs[2] = end_index-start_index
    self.ofq_cursor_view->setProperty, viewplane_rect=vpr_cs
    
    self->regenerateXAxis
    self->refreshFrequencyPlot

end


;+
; :Description:
;    Retrieves the extreme values (real, imag, abs) of
;    the data begin plotted. For each keyword specified,
;    it assigns an anonymous structure:
;     { max_yval: the maximum value over the whole data,
;       max_xval: the location of the max_yval,
;       min_yval: the minimum value over the whole data,
;       min_xval: the x location of the min_yval }
;
; :Keywords:
;    real 
;    imag
;    abs
;
; :Author: btt
;-
pro mas_complex_plot::getExtremeValues, $
    real=real, imag=imag, abs=abs

    if (arg_present(real)) then begin
        max_yval = max(real_part(*self.fq_data_displayed), max_ind, $
                       min=min_yval, subscript_min=min_ind) 
        max_xval = (*self.fq_xaxis_pts_disp)[max_ind]
        min_xval = (*self.fq_xaxis_pts_disp)[min_ind]
        
        real = { max_yval: max_yval, $
                 min_yval: min_yval, $
                 max_xval: max_xval, $
                 min_xval: min_xval }
    endif
    
    if (arg_present(imag)) then begin
        max_yval = max(imaginary(*self.fq_data_displayed), max_ind, $
                       min=min_yval, subscript_min=min_ind) 
        max_xval = (*self.fq_xaxis_pts_disp)[max_ind]
        min_xval = (*self.fq_xaxis_pts_disp)[min_ind]
        
        imag = { max_yval: max_yval, $
                 min_yval: min_yval, $
                 max_xval: max_xval, $
                 min_xval: min_xval }
    endif
    
end

;+
; :Description:
;    Retrieves a data value located at POS, where POS is
;    the index into the data array.
;    
;    returns the retrieved value as a keyword
; :Params:
;    pos - the index into the data array
;
; :Keywords:
;    real     - set this to get the real part
;    imag     - set this to get the imag. part
;    abs      - set this to get the magn. 
;    freqval  - set this to get the frequency val
;    all_fids - set this to get the data for all array
;               elements (for multi-array plots)
;    which_fid - set this to the array element to get 
;                (for multi-array plots)
;
; :Author: btt
;-
pro mas_complex_plot::getDataValues, pos, $
    real=real, imag=imag, abs=abs, freqval=freqval, $
    all_fids=all_fids, which_fid=which_fid

    if (keyword_set(all_fids)) then begin
        which_fid = lindgen(self.num_fids)
    endif else if n_elements(which_fid) eq 0 then begin
        which_fid = self.initial_fid
    endif else begin
        which_fid = (which_fid > 0) < (self.num_fids-1)
    endelse
    
    if (arg_present(real)) then begin
        real = reform(real_part((*self.fq_data)[pos, which_fid]))
    endif
    
    if (arg_present(imag)) then begin
        imag = reform(imaginary((*self.fq_data)[pos, which_fid]))
    endif
    
    if (arg_present(imag)) then begin
        imag = reform(abs((*self.fq_data)[pos, which_fid]))
    endif
    
    if (arg_present(freqval)) then begin
        freqval = (*self.fq_xaxis_pts_disp)[pos]
    endif
    
end

;+
; :Description:
;    event handler for clicks inside the view window for 
;    the purpose of repositioning cursors
;
; :Params:
;    event - event with element "x" which gives the 
;            x location of the mouse click
;
; :Keywords:
;    updated (not used)
;
; :Author: btt
;-
pro mas_complex_plot::windowClickEvent, event, updated=updated

    self.ofq_window->getProperty, dimensions=dim
    self.ofq_cursor_view->getProperty, viewplane_rect=vpr
    ;;print, dim & print, vpr
    updated=0
    data_pos = float(event.x)/dim[0] * vpr[2] + vpr[0]

    ;; window 0 ... dim[0] maps to vpr[0] ... vpr[0]+vpr[2]
    if (event.press) then begin
                
        d1 = abs(self.fq_cursor_01_pos - data_pos) / vpr[2]
        d2 = abs(self.fq_cursor_02_pos - data_pos) / vpr[2]
        min_d = min([d1,d2], min_ind)

        if (min_d lt 0.006) then begin
            self.fq_cursor_moving = min_ind+1
        endif else begin
            self.fq_cursor_moving = 0
        endelse
                
    endif else if (event.release) then begin
    
        self.fq_cursor_moving = 0
        
    endif
    
    if (self.fq_cursor_moving ne 0) then begin
        self->setCursorPosition, long(data_pos), self.fq_cursor_moving
        updated=1
    endif
    
end

;+
; :Description:
;    Retrieves information about the data in the interval
;    between cursors.
;
; :Keywords:
;    min_x     - the minimum x value
;    max_x     - the maximum x value
;    npts      - number of data points in the interval
;    realmean  - mean of real values in interval
;    realstdev - stded of real values in interval
;    realmin   - min of real values in interval
;    realmax   - max of real values in interval
;    imagmean  - same as above for imaginary
;    imagstdev - same as above for imaginary
;    imagmin   - same as above for imaginary
;    imagmax   - same as above for imaginary
;    absmean   - same as above for magnitude
;    absstdev  - same as above for magnitude
;    trendline - trendline information: 
;                trendline[0] = slope, 
;                trendline[1] = y-intercept
;
; :Author: btt
;-
pro mas_complex_plot::getCursorIntervalStats, $
    min_x=min_x, max_x=max_x, npts=npts, $
    realmean=realmean, realstdev=realstdev, realmin=realmin, realmax=realmax, realsum=realsum, realint=realint, $
    imagmean=imagmean, imagstdev=imagstdev, imagmin=imagmin, imagmax=imagmax, imagsum=imagsum, imagint=imagint, $
    absmean=absmean, absstdev=absstdev, trendline=trendline

    min_data_x = min([self.fq_cursor_01_pos, self.fq_cursor_02_pos], $
                     max=max_data_x)
    min_data_x = (0 > min_data_x) < (self.num_fidpts-1)     
    max_data_x = (0 > max_data_x) < (self.num_fidpts-1)
    
    min_x = (*self.fq_xaxis_pts)[min_data_x]
    max_x = (*self.fq_xaxis_pts)[max_data_x]
    
    npts = max_data_x - min_data_x
    
    interval = (*self.fq_data_displayed)[min_data_x:max_data_x]
    
    realmin = min(real_part(interval), max=rmax)
    imagmin = min(imaginary(interval), max=imax)
    realmax=rmax & imagmax=imax
    
    if (arg_present(realmean) or arg_present(realstdev)) then begin
        momt = moment(real_part(interval), maxmoment=2)
        realmean = momt[0]
        realstdev = sqrt(momt[1])
    endif

    if (arg_present(imagmean) or arg_present(imagstdev)) then begin
        momt = moment(imaginary(interval), maxmoment=2)
        imagmean = momt[0]
        imagstdev = sqrt(momt[1])
    endif

    if (arg_present(absmean) or arg_present(absstdev)) then begin
        momt = moment(abs(interval), maxmoment=2)
        absmean = momt[0]
        absstdev = sqrt(momt[1])
    endif
    
    if (arg_present(realsum) or arg_present(imagsum)) then begin
        ;;dx = abs((*self.fq_xaxis_pts_disp)[min_data_x] - (*self.fq_xaxis_pts_disp)[min_data_x+1])
        realsum = total(real_part(interval));*dx
        imagsum = total(imaginary(interval));*dx
    endif
    
    if (arg_present(realint) or arg_present(imagint)) then begin
        realint = int_tabulated(reverse((*self.fq_xaxis_pts_disp)[min_data_x:max_data_x]), $
                                reverse(real_part(interval)))
        imagint = int_tabulated(reverse((*self.fq_xaxis_pts_disp)[min_data_x:max_data_x]), $
                                reverse(imaginary(interval)))        
    endif
    
    if (arg_present(trendline)) then begin
        ;self.ofq_interval_trendline->getProperty, uvalue=trendline
        trendline = fltarr(3)
        X = (*self.fq_xaxis_pts)[min_data_x:max_data_x]
        Y = real_part(interval)
        trendline[0:1] = linfit(X, Y, chisqr=xi2)
        trendline[2] = xi2
    endif
    
end
;+
; :Description:
;    retrieves information about the data at the cursor location.
;
; :Keywords:
;    cursor01_info - this keyword contains a structure:
;                    data_index: the index into the data array of current location
;                    freq_xval: the xvalue at the current location
;                    freq_yval: the yvalue at the current location
;    cursor02_info - same as above, but for cursor two
;
;    all_fids - not used
;
; :Author: btt
;-
pro mas_complex_plot::getFreqPos, Freq, pos

    xval = *self.fq_xaxis_pts_disp
    ind = where(xval gt Freq[0] and xval lt Freq[1])
    pos = [min(ind), max(ind)]
    
end


;+
; :Description:
;    retrieves information about the data at the cursor location.
;
; :Keywords:
;    cursor01_info - this keyword contains a structure:
;                    data_index: the index into the data array of current location
;                    freq_xval: the xvalue at the current location
;                    freq_yval: the yvalue at the current location
;    cursor02_info - same as above, but for cursor two
;    
;    all_fids - not used
;
; :Author: btt
;-
pro mas_complex_plot::getCursorInfo, $
        cursor01_info=cursor01_info, $
        cursor02_info=cursor02_info, all_fids=all_fids
   
    if (arg_present(cursor01_info)) then begin
        data_index = (0 > self.fq_cursor_01_pos) < (self.num_fidpts-1)
        xval = (*self.fq_xaxis_pts_disp)[data_index]
        cursor01_info = { data_index: data_index, $
                          freq_xval: xval, $
                          freq_yval: (*self.fq_data_displayed)[data_index] }
    endif

    if (arg_present(cursor02_info)) then begin
        data_index = (0 > self.fq_cursor_02_pos) < (self.num_fidpts-1)
        xval = (*self.fq_xaxis_pts_disp)[data_index]
        cursor02_info = { data_index: data_index, $
                          freq_xval: xval, $
                          freq_yval: (*self.fq_data_displayed)[data_index] }
    endif
        
end

;+
; :Description:
;    repositions the cursor indicated by "which" to the
;    data index indicated by "data_index".
;
; :Params:
;    data_index - the data index to move to
;    which      - the cursor to move (1 or 2)
;
; :Author: btt
;-
pro mas_complex_plot::setCursorPosition, data_index, which 

    if (which eq 1) then begin
        self.ofq_cursor_01->getProperty, data=data1
        data1[0,*] = data_index
        self.ofq_cursor_01->setProperty, data=data1
        self.fq_cursor_01_pos = data_index
    endif else if (which eq 2) then begin
        self.ofq_cursor_02->getProperty, data=data2
        data2[0,*] = data_index
        self.ofq_cursor_02->setProperty, data=data2
        self.fq_cursor_02_pos = data_index
    endif
    
    self->getCursorIntervalStats, trendline=tl
    if (total(finite(tl, /nan)) eq 0) then begin
        
        start_x = (*self.fq_xaxis_pts)[0]
        end_x = (*self.fq_xaxis_pts)[self.num_fidpts-1]
        self.ofq_interval_trendline->setProperty, uvalue=tl, $
                                data=[ [start_x, tl[1]*start_x+tl[0] ], $
                                       [end_x,   tl[1]*end_x+tl[0]]  ]
    endif
    
    self->refreshFrequencyPlot

end

;+
; :Description:
;    Sets the cursor style.
;
; :Params:
;    which - the cursor to modify
;
; :Keywords: these keywords correspond to the keywords described in
;            IDL's IDLGRPOLYLINE object
;    thick
;    color
;    alpha
;    hide
;    linestyle
;
; :Author: btt
;-
pro mas_complex_plot::setCursorStyle, which, $
        thick=thick, color=color, alpha=alpha, hide=hide, linestyle=linestyle

    if (n_elements(which) eq 0) then return
    
    if (which le 0 or which gt 2) then return

    curs = (which eq 1) ? self.ofq_cursor_01 : self.ofq_cursor_02
    
    if (n_elements(thick) ne 0) then begin
        curs->setProperty, thick=thick
    endif
    
    if (n_elements(color) ne 0) then begin
        curs->setProperty, color=color
    endif
        
    if (n_elements(alpha) ne 0) then begin
        curs->setProperty, alpha_channel=alpha
    endif

    if (n_elements(linestyle) ne 0) then begin
        curs->setProperty, linestyle=linestyle
    endif
    
    if (n_elements(hide) ne 0) then begin
        curs->setProperty, hide=hide
    endif    
   
    self.ofq_window->draw, self.ofq_scene
    
end

pro mas_complex_plot::setHideMarker, name, hide

    omarker = self.ofq_marker_model->getByName(name)
    
    if (obj_valid(omarker)) then begin

        omarker->getProperty, uvalue=text
        if (obj_valid(text)) then text->setProperty, hide=hide
        
        omarker->setProperty, hide=hide
        self->refreshFrequencyPlot
        return
            
    endif

end

pro mas_complex_plot::removeMarker, name

    omarker = self.ofq_marker_model->getByName(name)
    
    if (obj_valid(omarker)) then begin

        omarker->GetProperty, uvalue=text
        if (obj_valid(text)) then obj_destroy, text
        
        self.ofq_marker_model->remove, omarker
        obj_destroy, omarker
        self->refreshFrequencyPlot
        return
            
    endif

end

pro mas_complex_plot::addMarker, data_pos, name, color=color, linestyle=linestyle

    if (n_elements(color) ne 3) then begin
        color = [128,128,128]
    endif
    
    if (n_elements(linestyle) eq 0) then begin
        linestyle = 1
    endif
    
    markerpos = [ [data_pos, 0], [data_pos, 100] ]
    if (keyword_set(text)) then begin
        text = obj_new('idlgrtext', text, locations=[data_pos, 97], alignment=0.5, $
                       font=self.ofq_xaxis_font, recompute_dimensions=2)
    endif else begin
        text = 0
    endelse
    
    marker = obj_new('idlgrpolyline', markerpos, name=name, $
                     linestyle=linestyle, color=color, uvalue=text)

    self.ofq_marker_model->add, marker
    if (obj_valid(text)) then self.ofq_marker_model->add, text
        
    self->refreshFrequencyPlot
    
end
;+
; :Description:
;    Sets the headsup text. The headsup text is displayed 
;    at the top left corner of the window. if string is
;    an array of string, then each element is placed on
;    a new line.
;
; :Params:
;    string
;
; :Author: btt
;-
pro mas_complex_plot::setHeadsupMessageText, string

    self.ofq_window->getProperty, dimensions=wdim
    char_height = 13

    locations = fltarr(2,n_elements(string))
    
    for s = 0, n_elements(string)-1 do begin
        locations[*,s] = [0, wdim[1]-(s+1)*char_height]
    endfor

    self.ofq_headsup_message_text->setproperty, $
            strings=string, $
            locations=locations
     
    self->refreshFrequencyPlot

end

;+
; :Description:
;    Refreshes the display window.
;    
;    Basically just redraws the scene.
;
; :Author: btt
;-
pro mas_complex_plot::refreshFrequencyPlot

;    self.ofq_plot_view->getProperty, viewplane_rect=plot_vpr
;    self.ofq_xaxis_view->getProperty, viewplane_rect=axis_vpr
;    self.ofq_cursor_view->getProperty, viewplane_rect=curs_vpr
;    
;    print, "plot: ", plot_vpr
;    print, "axis: ", axis_vpr
;    print, "cursor: ", curs_vpr
;    print, '----------------------------------------------------------------------'
    
    ;;self->regenerateXAxis
    self.ofq_window->draw, self.ofq_scene

end

;+
; :Description:
;    Hides or shows the plot symbols.
;
; :Params:
;    hide - set to 1 to hide the symbols, 0 or show
;    
; :Author: btt
;-
pro mas_complex_plot::setHidePlotSymbols, hide

    if (n_elements(hide) ne 0) then begin
        if (hide) then begin
            self.ofq_plot_sym_real->setProperty, data=0
        endif else begin
            self.ofq_plot_sym_real->setProperty, data=4        
        endelse
    endif
    
    self->refreshFrequencyPlot

end

;+
; :Description:
;    Hides or shows the zero line.
;
; :Params:
;    hide - set to 1 to hide the line, 0 or show
; :Author: btt
;-
pro mas_complex_plot::setHideZeroLine, hide

    if (n_elements(hide) ne 0) then begin
        if (hide) then begin
            self.ofq_zero_line->setProperty, hide=1
        endif else begin
            self.ofq_zero_line->setProperty, hide=0
        endelse
    endif
    
    self->refreshFrequencyPlot

end

;+
; :Description:
;    Hides or shows the interval trendline
;
; :Params:
;    hide - set to 1 to hide the line, 0 or show
; :Author: btt
;-
pro mas_complex_plot::setHideIntervalTrendLine, hide

    if (n_elements(hide) ne 0) then begin
        if (hide) then begin
            self.ofq_interval_trendline->setProperty, hide=1
            self->refreshFrequencyPlot
        endif else begin
            self.ofq_interval_trendline->setProperty, hide=0
            self->setCursorPosition, self.fq_cursor_01_pos, 1
        endelse
    endif

end

;+
; :Description:
;    Hides or shows the plots.
;
; :Params:
;    hide - set to 1 to hide the line, 0 or show
;    
; :Keywords:
;   abs  - apply hide option to magnitude plot
;   real - apply hide option to real plot
;   imag - apply hide option to imaginary plot
;   
; :Author: btt
;-
pro mas_complex_plot::setHideFreqPlot, hide, $
    abs=abs, real=real, imag=imag
    
    if (n_elements(abs)) then begin
        self.ofq_plot_abs->setProperty, hide=abs
    endif

    if (n_elements(real)) then begin
        self.ofq_plot_real->setProperty, hide=real
    endif

    if (n_elements(imag)) then begin
        self.ofq_plot_imag->setProperty, hide=imag
    endif
        
    self->refreshFrequencyPlot
    
end

pro mas_complex_plot__define

    struct = { mas_complex_plot, $
               plot_color_real: bytarr(3), $
               plot_color_imag: bytarr(3), $
               plot_color_abs: bytarr(3), $
               raw_fid: ptr_new(), $
               num_fidpts: 0L, $
               initial_fid: 0L, $
               num_fids: 0L, $
               x_lim: float(0), $
               fq_xaxis_pts: ptr_new(), $
               fq_xaxis_pts_disp: ptr_new(), $
               fq_data: ptr_new(), $
               fq_data_displayed: ptr_new(), $
               fq_zero_ref_index: 0L, $
               fq_data_index_start:0L, $
               fq_data_index_end: 0L, $
               fq_cursor_01_pos: 0.0, $
               fq_cursor_02_pos: 0.0, $
               fq_cursor_moving: 0L, $
               fq_ymax: 0.0, $
               fq_ymin: 0.0, $
               fq_yscale: 0.0, $
               fq_yoffset: 0.0, $
               fq_reverse_xaxis: 0B, $
               ofq_plot_view: obj_new(), $
               ofq_plot_model: obj_new(), $
               ofq_plot_real: obj_new(), $
               ofq_plot_imag: obj_new(), $
               ofq_plot_abs: obj_new(), $
               ofq_plot_sym_real:obj_new(), $
               ofq_plot_sym_imag:obj_new(), $
               ofq_plot_sym_abs:obj_new(), $
               ofq_interval_trendline: obj_new(), $
               ofq_xaxis_view: obj_new(), $
               ofq_xaxis_model: obj_new(), $
               ofq_xaxis: obj_new(), $
               ofq_yaxis: obj_new(), $
               ofq_xaxis_font: obj_new(), $
               ofq_zero_line: obj_new(), $
               ofq_cursor_view: obj_new(), $
               ofq_cursor_model: obj_new(), $
               ofq_cursor_01: obj_new(), $
               ofq_cursor_02: obj_new(), $
               ofq_marker_model: obj_new(), $
               ofq_headsup_view: obj_new(), $
               ofq_headsup_model: obj_new(), $
               ofq_headsup_message_text: obj_new(), $
               ofq_scene: obj_new(), $
               ofq_window: obj_new(), $
               $
               void: 0 $
              }

end

