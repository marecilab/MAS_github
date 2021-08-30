;$Id: //depot/idl/IDL_64/idldir/lib/utilities/idlexobjview__define.pro#1 $
;
;  Copyright (c) 1997-2007, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
;+
; NAME:
;   IDLexObjview
;
; PURPOSE:
;   Provide a subclass of IDLexInscribingView that contains a hierarchy
;   of models dedicated to various types of transformations.  The
;   transformations can be updated via an Update method that can
;   be passed mouse events.  Also, via certain models in the hierarchy,
;   subjects in this view are automatically centered about the origin
;   when self's view volume is set via self->SetViewVolume.
;
;   Default stationary lights are created and added if none are
;   supplied by the user.
;
;   This class is deployed by the XOBJVIEW command.
;
; CREATION AND USE:
;       oView = obj_new('IDLexObjview')
;       oView->Add, oObj
;
;       oObj is a reference to an atomic graphics object or to an
;       IDLgrModel.  oObj can be an array of such references.
;
;   If non-default lights (or other stationary objects) are desired:
;
;       oView->Add, oObj, /stationary
;
;       As in the previous example, oObj is a reference to an atomic
;       graphics object, a reference to an IDLgrModel, or an array
;       of such references.  In this example, if oObj is, or contains,
;       one or more IDLgrLights then the view's default lights
;       are overridden.
;
; EXAMPLE:
;   See procedure IDLexObjview__example in file idlexobjview__define.pro.
;   Also see xobjview.pro.
;
;-
;
function IDLexObjview__contains_lights, oObj
;
;Purpose: If oObj is a reference to an IDLgrLight, or is a
;   reference to an IDLgrModel that contains IDLgrLights (either
;   directly or in child IDLgrModels), or is an array with at
;   least one such reference, then return 1b, else return 0b.
;
on_error, 2 ; Return to caller on error.

for i=0,n_elements(oObj)-1 do begin
    case 1 of
        obj_valid(oObj[i]) eq 0:
        obj_isa(oObj[i], 'IDLgrModel'): begin
            oChildren = oObj[i]->Get(/all)
            if IDLexObjview__contains_lights(oChildren) then begin
                return, 1b
                end
            end
        obj_isa(oObj[i], 'IDLgrLight'): begin
            return, 1
            end
        else:
        endcase
    end

return, 0b
end
;--------------------------------------------------------------------
function IDLexObjview::CalcRanges, $
    oDestDevice, $  ; IN: supplies GetTextDimensions() if needed.
    range, $ ; OUT: [3,2] coordinates of bounding box. Used in recursion.
    oObj_, $ ; IN (opt): object to find boundaries of. Used in recursion.
    oParent_ ; IN (opt): is PATH for GetTextDimensions, if needed. Used
             ; in recursion.
compile_opt hidden

if n_elements(oObj_) eq 0 then begin
    oObj = self.oModel0
    self.oModel1Manip->SetTarget, obj_new()
    self.oModel3Manip->SetTarget, obj_new()
    end $
else begin
    if not obj_valid(oObj_) then $
        return, 0
    oObj = oObj_
    end

oParent = n_elements(oParent_) eq 0 ? self : oParent_

case 1 of
    obj_isa(oObj, 'IDLgrLight'): return, 0
    obj_isa(oObj, 'IDLgrModel'): begin
;
;       Find the boundary of oObj's children.
;
        oChildArr = oObj->Get(/all, count=count)
        got_box = 0L
        path = [oParent, oObj]
        for i=0L, count-1 do begin ; Get each box and accumulate bounding box
            if self->CalcRanges(oDestDevice, child_range, oChildArr[i], path) $
            then begin
                if got_box then begin
                    range = [ $
                        [range[0:2] < child_range[0:2]], $
                        [range[3:5] > child_range[3:5]] $
                        ]
                    end $
                else begin
                    range = child_range
                    end
                got_box = 1
                endif
            end

        if got_box then begin        ;Got a box?  Transform it..
            oObj->GetProperty, transform=transform ;Transform it and get range
            p1 = fltarr(4,8)
            for i=0, 7 do $
                p1[0,i] = [ $
                    range[0, i and 1],     $
                    range[1, (i/2) and 1], $
                    range[2, (i/4) and 1], $
                    1.0 $
                    ]
            p1 = matrix_multiply(p1, transform, /atranspose)
            range = [ $
                [min(p1[*,0], max=xmax), min(p1[*,1], max=ymax), $
                    min(p1[*,2], max=zmax)], $
                [xmax, ymax, zmax] $
                ]
            case oObj of
                self.oModel0: self.model0_range = range
                self.oModel1: self.model1_range = range
                self.oModel2: self.model2_range = range
                self.oModel3: self.model3_range = range
                self.oStationaryModel: self.stationary_range = range
                else:
                endcase
            end

        result = got_box
        case oObj of
            self.oModel0: self.got_model0_range = got_box
            self.oModel1: self.got_model1_range = got_box
            self.oModel2: self.got_model2_range = got_box
            self.oModel3: self.got_model3_range = got_box
            self.oStationaryModel: self.got_stationary_range = got_box
            else:
            endcase

        end
    obj_isa(oObj, 'IDLgrGraphic'): begin
        if obj_isa(oObj, 'IDLgrText') then begin
            if obj_valid(oDestDevice) then begin
                if obj_valid(oParent[0]) then $
;
;                   GetTextDimensions is called for this side effect: oObj's
;                   xrange, yrange & zrange are updated.
;
                    void = oDestDevice->GetTextDimensions( $
                        oObj, $
                        path=[oParent] $
                        ) $
                else $
                    void = oDestDevice->GetTextDimensions(oObj)
                endif
            endif

        oObj->GetProperty, xrange=xrange, yrange=yrange, zrange=zrange, $
            xcoord_conv=xcc, ycoord_conv=ycc, zcoord_conv=zcc
;
;       Scale and save extrema.
;
        range = [[min(xrange * xcc[1] + xcc[0], max=xmax), $
            min(yrange * ycc[1] + ycc[0], max=ymax), $
            min(zrange * zcc[1] + zcc[0], max=zmax)], $
            [xmax, ymax, zmax]]
;
;       IDLgrContour objects can have NaN ranges if their properties
;       yield no contours (e.g. DATA_VALUES are all 0, or C_VALUE
;       is [-1] when DATA_VALUES are positive).  Test for NaN's
;       using FINITE().
;
;
        if min(finite(range)) lt 1 then $
            result = 0 $
        else $
            result = 1
        end
    else: begin ; Do not know what oObj is.
        result = 0
        end
    endcase

if n_elements(oObj_) eq 0 then begin
    self.oModel1Manip->SetTarget, self.oModel1, oDestDevice
    self.oModel3Manip->SetTarget, self.oModel3, oDestDevice
    end

return, result

end
;--------------------------------------------------------------------
function IDLexObjView::GetBounds, $
    oDestDevice, $  ; IN: Supplies GetTextDimensions() if needed.
    xrange, $       ; OUT
    yrange, $       ; OUT
    zrange, $       ; OUT
    ignore=ignore, $; IN: (opt) classes.
    skip=skip       ; IN: (opt) classes. (non-recursive).
compile_opt hidden
;
;Determine the overall bounding box for objects contained in self.
;
self.oModel1Manip->SetTarget, obj_new()
self.oModel3Manip->SetTarget, obj_new()

oModels = self->IDLexInscribingView::Get(/all, count=n_models)
for i=0L, n_models-1 do begin
    if get_obj_range( $
        oModels[i], $
        self, $
        oDestDevice, $
        model_range, $
        ignore=ignore, $
        skip=skip $
        ) $
    then begin
        if n_elements(range) ne 0 then $
            range = [ $
                [range[0:2] < model_range[0:2]], $
                [range[3:5] > model_range[3:5]] $
                ] $
        else $
            range = model_range
        end
    end

if n_elements(range) ne 0 then begin ;Extract individual ranges
    xrange = range[[0,3]]
    yrange = range[[1,4]]
    zrange = range[[2,5]]
    result = 1
    endif $
else $
    result = 0

self.oModel1Manip->SetTarget, self.oModel1, oDestDevice
self.oModel3Manip->SetTarget, self.oModel3, oDestDevice

return, result
end
;--------------------------------------------------------------------
pro IDLexObjview::Normalize, $
    oDestDevice, $
    xrange=xrange, $
    yrange=yrange, $
    zrange=zrange, $
    adjust_axes=adjust_axes ; IN
compile_opt hidden
;
;Purpose: Scale and Translate member models so that self's contents
;fit in a default view volume (-1 to 1 in x, y and z). Also, if
;keyword ADJUST_AXES is set, set the range and position of any axes
;that are immediate children of self.oModel3.
;
self.oModel2->Reset
self.oStationaryModel->Reset

oContents = self.oModel3->Get(/all, isa='IDLgrComponent', count=n_drawables)
if n_drawables eq 0 then begin
    self.x_len = 0
    self.y_len = 0
    self.z_len = 0
    self.xrange = [0, 0]
    self.zrange = [0, 0]
    self.yrange = [0, 0]
    return
    end

if keyword_set(adjust_axes) then begin
    axis_indx = where(obj_isa(oContents, 'IDLgrAxis'), n_axes)
    if n_axes gt 0 then begin
        oAxes = oContents[axis_indx]
;
;       Store some properties about the axes.
;
        tick_char_dimensions = fltarr(2, n_axes)
        title_char_dimensions = fltarr(2, n_axes)
        tick_lengths = fltarr(n_axes)
        for i=0,n_axes-1 do begin
            oAxes[i]->GetProperty, $
                ticktext=oTickText, $
                title=oTitleText, $
                ticklen=ticklen
            tick_lengths[i] = ticklen
            if obj_valid(oTickText) then begin
                oTickText->GetProperty, char_dimensions=char_dimensions
                tick_char_dimensions[0, i] = char_dimensions
                end
            if obj_valid(oTitleText) then begin
                oTitleText->GetProperty, char_dimensions=char_dimensions
                title_char_dimensions[0, i] = char_dimensions
                end
            end
;
;       Obtain the range to which we will be setting axes.
;
        self.oModel1Manip->SetTarget, obj_new()
        self.oModel3Manip->SetTarget, obj_new()
        if get_obj_range( $
            self.oModel3, $
            [self, self.oModel0, self.oModel1, self.oModel2], $
            oDestDevice, $
            range, $
            skip='IDLgrAxis' $
            ) $
        then begin
            xrng = range[0, *]
            yrng = range[1, *]
            zrng = range[2, *]
            end $
        else begin
            xrng = [-1, 1]
            yrng = [-1, 1]
            zrng = [-1, 1]
            end
        self.oModel1Manip->SetTarget, self.oModel1, oDestDevice
        self.oModel3Manip->SetTarget, self.oModel3, oDestDevice

        if n_elements(xrange) gt 0 then $
            xrng = xrange
        if n_elements(yrange) gt 0 then $
            yrng = yrange
        if n_elements(zrange) gt 0 then $
            zrng = zrng
;
;       Set axes to span our newly obtained ranges.
;
        x_len = xrng[1] - xrng[0]
        y_len = yrng[1] - yrng[0]
        z_len = zrng[1] - zrng[0]

        for i=0,n_axes-1 do begin
            oAxes[i]->GetProperty, direction=direction, range=range

            case direction of
                0: begin
                    tick_char_size_factor = [ $
                        self.x_len eq 0 ? 1. : self.x_len, $
                        self.y_len eq 0 ? 1. : self.y_len $
                        ]
                    tick_length_factor = self.y_len eq 0 ? 1 : self.y_len
                    title_char_size_factor = tick_char_size_factor
                    end
                1: begin
                    tick_char_size_factor = [ $
                        self.x_len eq 0 ? 1. : self.x_len, $
                        self.y_len eq 0 ? 1. : self.y_len $
                        ]
                    tick_length_factor = self.x_len eq 0 ? 1 : self.x_len
                    title_char_size_factor = reverse(tick_char_size_factor)
                    ;title_char_size_factor = tick_char_size_factor
                    end
                2: begin
                    tick_char_size_factor = [ $
                        self.x_len eq 0 ? 1. : self.x_len, $
                        self.z_len eq 0 ? 1. : self.z_len $
                        ]
                    tick_length_factor = self.x_len eq 0 ? 1 : self.x_len
                    title_char_size_factor = reverse(tick_char_size_factor)
                    end
                endcase

            tick_char_dimensions[0, i] = $
                tick_char_dimensions[*, i] / tick_char_size_factor
            title_char_dimensions[0, i] = $
                title_char_dimensions[*, i] / title_char_size_factor
            tick_lengths[i] = tick_lengths[i] / tick_length_factor

            case direction of
                0: if x_len gt 0 then begin
                    oAxes[i]->SetProperty, range=xrng
                    oAxes[i]->GetProperty, crange=crange
                    xrng[0] = xrng[0] < crange[0]
                    xrng[1] = xrng[1] > crange[1]
                    end
                1: if y_len gt 0 then begin
                    oAxes[i]->SetProperty, range=yrng
                    oAxes[i]->GetProperty, crange=crange
                    yrng[0] = yrng[0] < crange[0]
                    yrng[1] = yrng[1] > crange[1]
                    end
                2: if z_len gt 0 then begin
                    oAxes[i]->SetProperty, range=zrng
                    oAxes[i]->GetProperty, crange=crange
                    zrng[0] = zrng[0] < crange[0]
                    zrng[1] = zrng[1] > crange[1]
                    end
                endcase
            endfor
        x_len = xrng[1] - xrng[0]
        y_len = yrng[1] - yrng[0]
        z_len = zrng[1] - zrng[0]
;
;       Re-locate axes.
;
        for i=0,n_elements(oAxes)-1 do begin
            oAxes[i]->GetProperty, direction=direction, location=location
            location = location - [ $
                self.xrange[0], $
                self.yrange[0], $
                self.zrange[0] $
                ]
            case direction of
                0: begin
                    normalized_location = location / [ $
                        1, $
                        self.y_len eq 0 ? 1 : self.y_len, $
                        self.z_len eq 0 ? 1 : self.z_len $
                        ]
                    location = [0, yrng[0], zrng[0]]
                    location = location + $
                        [0, y_len, z_len] * normalized_location
                    end
                1: begin
                    normalized_location = location / [ $
                        self.x_len eq 0 ? 1 : self.x_len, $
                        1, $
                        self.z_len eq 0 ? 1 : self.z_len $
                        ]
                    location = [xrng[0], 0, zrng[0]]
                    location = location + $
                        [x_len, 0, z_len] * normalized_location
                    end
                2: begin
                    normalized_location = location / [ $
                        self.x_len eq 0 ? 1 : self.x_len, $
                        self.y_len eq 0 ? 1 : self.y_len, $
                        1 $
                        ]
                    location = [xrng[0], yrng[0], 0]
                    location = location + $
                        [x_len, y_len, 0] * normalized_location
                    end
                endcase
            oAxes[i]->SetProperty, location=location
            endfor

        self.xrange = xrng
        self.yrange = yrng
        self.zrange = zrng

        self.x_len = x_len
        self.y_len = y_len
        self.z_len = z_len
;
;       Adjust size of axis text and tickmarks.
;
        for i=0,n_elements(oAxes)-1 do begin
            oAxes[i]->GetProperty, $
                direction=direction, $
                ticktext=oTickText, $
                title=oTitleText

            if self.isotropic and (self.axes_are_isotropic eq 0) then begin
                x_len = x_len > y_len
                y_len = x_len
                end

            case direction of
                0: begin
                    oAxes[i]->SetProperty, ticklen=y_len*tick_lengths[i]
                    tick_char_size_factor = [ $
                        x_len gt 0 ? x_len : 1., $
                        y_len gt 0 ? y_len : 1. $
                        ]
                    title_char_size_factor = tick_char_size_factor
                    end
                1: begin
                    oAxes[i]->SetProperty, ticklen=x_len*tick_lengths[i]
                    tick_char_size_factor = [ $
                        x_len gt 0 ? x_len : 1., $
                        y_len gt 0 ? y_len : 1. $
                        ]
                    title_char_size_factor = reverse(tick_char_size_factor)
                    end
                2: begin
                    oAxes[i]->SetProperty, ticklen=x_len*tick_lengths[i]
                    tick_char_size_factor = [ $
                        x_len gt 0 ? x_len : 1., $
                        z_len gt 0 ? z_len : 1. $
                        ]
                    title_char_size_factor = reverse(tick_char_size_factor)
                    end
                endcase

            if obj_valid(oTickText) then begin
                oTickText->SetProperty, $
                    char_dimensions=$
                        tick_char_dimensions[*, i] * tick_char_size_factor
                end

            if obj_valid(oTitleText) then begin
                oTitleText->SetProperty, $
                    char_dimensions=$
                        title_char_dimensions[*, i] * title_char_size_factor
                end
            endfor
        endif
    self.axes_are_isotropic = self.isotropic
    end

if not self->CalcRanges(oDestDevice) then begin
    self.x_len = 0
    self.y_len = 0
    self.z_len = 0
    self.xrange = [0, 0]
    self.zrange = [0, 0]
    self.yrange = [0, 0]
    return
    end
;
;Scale self's contents equally in x, y and z to fit inside default
;view volume.
;
x_len = self.model0_range[0, 1] - self.model0_range[0, 0]
y_len = self.model0_range[1, 1] - self.model0_range[1, 0]
z_len = self.model0_range[2, 1] - self.model0_range[2, 0]

len = x_len > y_len > z_len

scale = len eq 0 ? 1 : (2. / len)
self.oModel0->Scale, scale, scale, scale
;
;Center self's contents about zero.
;
case 1 of
    self.got_model2_range and self.got_stationary_range: begin
        x_min = self.model2_range[0] < self.stationary_range[0]
        y_min = self.model2_range[1] < self.stationary_range[1]
        z_min = self.model2_range[2] < self.stationary_range[2]
        x_max = self.model2_range[3] > self.stationary_range[3]
        y_max = self.model2_range[4] > self.stationary_range[4]
        z_max = self.model2_range[5] > self.stationary_range[5]
        end
    self.got_model2_range: begin
        x_min = self.model2_range[0]
        y_min = self.model2_range[1]
        z_min = self.model2_range[2]
        x_max = self.model2_range[3]
        y_max = self.model2_range[4]
        z_max = self.model2_range[5]
        end
    self.got_stationary_range: begin
        x_min = self.stationary_range[0]
        y_min = self.stationary_range[1]
        z_min = self.stationary_range[2]
        x_max = self.stationary_range[3]
        y_max = self.stationary_range[4]
        z_max = self.stationary_range[5]
        end
    endcase
x_len = x_max - x_min
y_len = y_max - y_min
z_len = z_max - z_min

x_offset = -x_min - x_len / 2.
y_offset = -y_min - y_len / 2.
z_offset = -z_min - z_len / 2.

self.oModel2->Translate, x_offset, y_offset, z_offset
self.oStationaryModel->Translate, x_offset, y_offset, z_offset
;
if not self.isotropic then begin
;
;   Stretch self's contents to fit a cube.
;
    case 1 of
        x_len ge y_len and x_len ge z_len: begin
            x_scale = 1
            y_scale = y_len eq 0 ? 1 : (x_len / y_len)
            z_scale = z_len eq 0 ? 1 : (x_len / z_len)
            end
        y_len ge x_len and y_len ge z_len: begin
            y_scale = 1
            x_scale = x_len eq 0 ? 1 : (y_len / x_len)
            z_scale = z_len eq 0 ? 1 : (y_len / z_len)
            end
        z_len ge x_len and z_len ge y_len: begin
            z_scale = 1
            x_scale = x_len eq 0 ? 1 : (z_len / x_len)
            y_scale = y_len eq 0 ? 1 : (z_len / y_len)
            end
        endcase
    self.oModel2->Scale, x_scale, y_scale, z_scale
    self.oStationaryModel->Scale, x_scale, y_scale, z_scale
    end

end
;--------------------------------------------------------------------
pro IDLexObjview::Reset, $
    oDestDevice, $  ; IN:
    full=full, $    ; IN: (opt)
    adjust_axes=adjust_axes
compile_opt hidden
;
;Purpose: Undo all zoom, rotate and translates effected by self's
;Update, Zoom, Rotate, and Translate methods.  Set self's
;view volume to be suitable for viewing Normalized contents on
;the given destination device.
;
;KEYWORD:
;   FULL    If this keyword is set, Reset causes self to normalize its
;           contents via self->Normalize. If self contains many graphics
;           objects, normalization can be slow (i.e. CPU intensive).
;   ADJUST_AXES
;           If this keyword is set, perform a FULL Reset, and adjust
;           axes to fit other graphics objects in self.
;
self.oModel1->SetProperty, transform=identity(4)
self.oModel3->SetProperty, transform=identity(4)

self.oModel2->GetProperty, transform=transform
self.oStationaryModel->SetProperty, transform=transform

self->SetProperty, viewplane_rect=[-1, -1, 2, 2], zclip=[1, -1]

if self.isotropic then $
    self->PadViewplaneRect, oDestDevice

if keyword_set(full) or keyword_set(adjust_axes) then begin
    self->Normalize, oDestDevice, adjust_axes=keyword_set(adjust_axes)
    self.oModel2->Scale, self.scale[0], self.scale[1], self.scale[2]
    self.oStationaryModel->Scale, self.scale[0], self.scale[1], self.scale[2]
    end

end
;--------------------------------------------------------------------
pro IDLexObjview::RakeHeap
compile_opt hidden

obj_destroy, self.oModel0
obj_destroy, self.oModel1
obj_destroy, self.oModel2
obj_destroy, self.oModel3
obj_destroy, self.oStationaryModel
obj_destroy, self.oModel1Manip
obj_destroy, self.oModel3Manip
obj_destroy, self.oViewManip
self.oGetList->Remove, /all
self.oRemoveList->Remove, /all
obj_destroy, self.oGetList
obj_destroy, self.oRemoveList
ptr_free, self.pSelected

end
;--------------------------------------------------------------------
pro IDLexObjview::cleanup
compile_opt idl2, hidden

on_error, 2

self->IDLexObjView::RakeHeap
self->IDLexInscribingView::cleanup
end
;--------------------------------------------------------------------
function IDLexObjview::init, $
    oSubjects, $                ; Obsolete.  Use self->Add method instead.
    stationary=oStationary, $   ; Obsolete.  Use self->Add, /stationary.
    mode=mode, $                ; IN: (opt) 'rotate', 'zoom', ....
    scale=scale, $              ; IN: (opt) Scalar or 3-element array.
    debug=debug, $              ; IN: (opt)
    _extra=e
compile_opt idl2

on_error, 2

catch, error_stat
if error_stat ne 0 then begin
    catch, /cancel
    printf, $
        -2, $
        !error_state.msg_prefix, $
        !error_state.msg, $
        !error_state.sys_msg
    return, 0
    end
if keyword_set(debug) then $
   catch, /cancel

if not self->IDLexInscribingView::init(_extra=e) then begin
    message, 'Failed to init IDLexInscribingView part of self.'
    end

catch, error_stat
if error_stat ne 0 then begin
    catch, /cancel
    self->RakeHeap
    self->IDLexInscribingView::cleanup
    printf, $
        -2, $
        !error_state.msg_prefix, $
        !error_state.msg, $
        !error_state.sys_msg
    return, 0
    end
if keyword_set(debug) then $
    catch, /cancel

if n_elements(scale) eq 0 then begin
    self.scale = 1./sqrt(3)
    endif $
else begin
    self.scale = scale
    endelse
if min(self.scale) lt 0 then $
    message, 'SCALE must be a positive value.'
if self.scale[0] eq 0 then begin
    self.scale = 1./sqrt(3)
    end

self.oModel0 = obj_new('IDLgrModel')
self.oModel1 = obj_new('IDLgrModel')
self.oModel2 = obj_new('IDLgrModel')
self.oModel3 = obj_new('IDLgrModel')
self.oStationaryModel = obj_new('IDLgrModel')

self.oModel1Manip = obj_new('IDLexModelManip', $
    mode=1, $ ; Rotation
    rot_radius=1, $
    manip_show_sel=0, $
    manip_show_manip=0 $
    )
self.oModel1Manip->SetTarget, self.oModel1; (oDest arg not needed here.)

self.oModel3Manip = obj_new('IDLexModelManip', $
    mode=0, $ ; Translation
    translate=[1, 1, 1], $
    manip_show_sel=0, $
    manip_show_manip=0 $
    )
self.oModel3Manip->SetTarget, self.oModel3; (oDest arg not needed here.)

self.oViewManip = obj_new('IDLexViewManip', $
    mode=1, $ ; Zoom
    manip_show_manip=1, $
    manip_show_sel=0 $
    )
;
;Create some lights.
;
self.oDefaultLights[0] = obj_new('IDLgrLight', $
    loc=[2,2,5], $
    type=2, $ ; Directional (parallel rays).
    color=[255,255,255], $
    intensity=.5 $
    )
self.oDefaultLights[1] = obj_new('IDLgrLight', $
    type=0, $ ; Ambient.
    intensity=.5, $
    color=[255,255,255] $
    )
;
;Note: use of variable oStationary here is obsolete, and is included only
;for compatibilty with older code.
;
if not IDLexObjview__contains_lights(oStationary) then begin
    self.oStationaryModel->Add, self.oDefaultLights
    end
if n_elements(oStationary) gt 0 then begin
    self.oStationaryModel->Add, oStationary, /alias
    end
;
self->IDLgrView::Add, self.oModel0
self.oModel0->Add, [self.oModel1, self.oStationaryModel]
self.oModel1->Add, self.oModel2
self.oModel2->Add, self.oModel3
;
;Note: use of variable "oSubjects" here is obsolete, and is included only
;for compatibilty with older code.
;
if n_elements(oSubjects) gt 0 then $
    self.oModel3->Add, oSubjects, /alias
;
if n_elements(mode) gt 0 then begin
    self.mode = mode
    end $
else begin
    self.mode = 'rotate'
    end

self.pSelected = ptr_new(obj_new())
self.debug = keyword_set(debug)

self.oGetList = obj_new('IDL_Container')
self.oRemoveList = obj_new('IDL_Container')

return, 1 ; Success.
end
;--------------------------------------------------------------------
pro IDLexObjview::Zoom, factor
compile_opt hidden
;
;Purpose: Provide an interface to the same zoom behavior that is
;effected via self's Update method.
;
self->ZoomViewplaneRect, factor
self.oStationaryModel->Scale, 1./factor, 1./factor, 1
end
;--------------------------------------------------------------------
pro IDLexObjview::Rotate, axis, angle, _extra=e
compile_opt hidden
;
;Purpose: Provide an interface to the same rotation behavior that is
;effected via self's Update method.
;
self.oModel1->Rotate, axis, angle, _extra=e
end
;--------------------------------------------------------------------
pro IDLexObjview::Translate, move_x, move_y, move_z, _extra=e
compile_opt hidden
;
;Purpose: Provide an interface to the same translation behavior that is
;effected via self's Update method.
;
self.oModel3->Translate, move_x, move_y, move_z, _extra=e
end
;--------------------------------------------------------------------
pro IDLexObjview::Add, $
    oObj, $
    stationary=stationary, $
    position=position, $
    _extra=e

on_error, 2

if n_elements(oObj) eq 0 then $
    message, 'Argument is undefined.'

if size(oObj, /tname) ne 'OBJREF' then $
    message, 'Argument must be an object reference.'

if min(obj_valid(oObj)) eq 0 then $
    message, 'Argument contains an invalid object.'
;
;If you are adding stationary lights, then the default stationary lights
;are destroyed.
;
if keyword_set(stationary) then begin
    if IDLexObjview__contains_lights(oObj) then begin
        self.oStationaryModel->Remove, self.oDefaultLights
        obj_destroy, self.oDefaultLights
        end
    end
;
;The following code maintains the illusion that self is an IDL_Container
;that behaves like any other IDL_Container.  If the user adds an object
;at a specified, or default position, then the user can get that object
;from the expected position.  For example, obj -> add, oModel, position = 1
;result = obj -> get(position = 1).
;
if n_elements(position) gt 0 then begin
    internal_position = position
    if self.oRemoveList->Count() gt 0 then begin
        arr = self.oRemoveList->Get(/all)
        for i = 0, n_elements(position) -1 do begin
            if position[i] gt 0 then $
                arr = arr[0:position[i] < n_elements(arr) -1]

            if keyword_set(stationary) then begin
                void = where(arr eq self.oModel3, count)
                if obj_valid(self.oDefaultLights[0]) then $
                    count = count - 2 ; Compensate for both default lights.
                end $
            else begin
                void = where(arr eq self.oStationaryModel, count)
                end

            internal_position = position[i] - count

            if keyword_set(stationary) then begin
                self.oStationaryModel->Add, $
                    oObj[i], $
                    position=internal_position, $
                    _extra=e
                self.oRemoveList->Add, $
                    self.oStationaryModel, $
                    position=position[i], $
                    _extra=e
                end $
            else begin
                self.oModel3->Add, $
                    oObj[i], $
                    position=internal_position, $
                    _extra=e
                self.oRemoveList->Add, $
                    self.oModel3, $
                    position=position[i], $
                    _extra=e
                end
            self.oGetList->Add, oObj[i], position=position[i], _extra=e
            end
        end
    end $
else begin
    if keyword_set(stationary) then begin
        self.oStationaryModel->Add, oObj, _extra=e
        self.oRemoveList->Add, $
            replicate(self.oStationaryModel, n_elements(oObj)), $
            _extra=e
        end $
    else begin
        self.oModel3->Add, oObj, _extra=e
        self.oRemoveList->Add, $
            replicate(self.oModel3, n_elements(oObj)), $
            _extra=e
        end
    self.oGetList->Add, oObj, _extra=e
    end

end
;--------------------------------------------------------------------
function IDLexObjview::Count
compile_opt hidden

return, self.oGetList->Count()

end
;--------------------------------------------------------------------
function IDLexObjview::Get, _ref_extra=e
compile_opt hidden

return, self.oGetList->Get(_extra=e)

end
;--------------------------------------------------------------------
pro IDLexObjview::Remove, oObj, _extra=e
compile_opt hidden

if n_elements(oObj) gt 0 then begin
    oObjects = oObj
    end $
else begin
    oObjects = self.oGetList->Get(_extra=e)
    end

for i=0,n_elements(oObjects)-1 do begin
    if obj_valid(oObjects[i]) then begin
        if self.oGetList->IsContained(oObjects[i], position=position) $
        then begin
            oModel = self.oRemoveList->Get(position=position)
            oModel->Remove, oObjects[i]
            self.oGetList->Remove, position=position
            self.oRemoveList->Remove, position=position
            end
        end
    end

end
;--------------------------------------------------------------------
function IDLexObjview::IsContained, oObj, _extra=e
compile_opt hidden

return, self.oGetList->IsContained(oObj, _extra=e)

end
;--------------------------------------------------------------------
pro IDLexObjview::Move, source, destination
compile_opt hidden

self.oGetList->Move, source, destination
self.oRemoveList->Move, source, destination

end
;--------------------------------------------------------------------
function IDLexObjview::Update, event
;
;Purpose: Handle given event.  Return 1 if the event requires a
;re-draw, else return 0.  EVENT must be a structure.  Ignore
;(no-op) EVENTs that are not pertinent.
;
if not keyword_set(self.debug) then $
    on_error, 2 ; Return to caller on error.

result = 0 ;Initialize return value.
;
;Ignore non-Draw-Widget events.
;
if (tag_names(event, /structure_name) ne 'WIDGET_DRAW') then $
    return, result
;
widget_control, event.id, get_value=oWindow
case event.type of
    0: begin ; Button press.
        self.mmb_is_down = 1b
        case self.mode of
            'rotate': begin
                self.oModel1Manip->MouseDown, [event.x, event.y], oWindow
                end
            'zoom': begin
                self.oViewManip->SetTarget, self, oWindow
                self.oViewManip->MouseDown, [event.x, event.y], oWindow
                end
            'translate': begin
                self.oModel3Manip->MouseDown, [event.x, event.y], oWindow
                end
            'select': begin
                self.oModel1Manip->SetProperty, hide=1 ; Necessary?
                self.oModel3Manip->SetProperty, hide=1 ; Necessary?
                self.oViewManip->SetProperty, hide=1 ; Necessary?

                oSelected = oWindow->Select(self, [event.x, event.y])
                if obj_valid(oSelected[0]) then begin
                    *self.pSelected = oSelected
                    end $
                else begin
                    *self.pSelected = obj_new()
                    end

                self.oModel1Manip->SetProperty, hide=0
                self.oModel3Manip->SetProperty, hide=0
                self.oViewManip->SetProperty, hide=0
                end
            else:
            endcase
        end
    2: begin ; Button motion
        if self.mmb_is_down then begin
            result = 1b
            case self.mode of
                'rotate': begin
                    self.oModel1Manip->MouseTrack, $
                        [event.x, event.y], $
                        oWindow
                    end
                'zoom': begin
                    self->GetProperty, viewplane_rect=viewplane_rect
                    old_size = viewplane_rect[2]
                    self.oViewManip->MouseTrack, $
                        [event.x, event.y], $
                        oWindow
                    self->GetProperty, viewplane_rect=viewplane_rect
                    new_size = viewplane_rect[2]
                    self.oStationaryModel->Scale, $
                        new_size / old_size, $
                        new_size / old_size, $
                        1
                    end
                'translate': begin
                    self.oModel3Manip->MouseTrack, $
                        [event.x, event.y], $
                        oWindow
                    end
                else: result = 0b
                endcase
            end
        end
    1: begin ; Button release
        if self.mmb_is_down then begin
            result = 1b
            case self.mode of
                'rotate': begin
                    self.oModel1Manip->MouseUp, [event.x, event.y], oWindow
                    end
                'zoom': begin
                    self->GetProperty, viewplane_rect=viewplane_rect
                    old_size = viewplane_rect[2]
                    self.oViewManip->MouseUp, [event.x, event.y], oWindow
                    self.oViewManip->SetTarget, obj_new()
                    self->GetProperty, viewplane_rect=viewplane_rect
                    new_size = viewplane_rect[2]
                    self.oStationaryModel->Scale, $
                        new_size / old_size, $
                        new_size / old_size, $
                        1
                    end
                'translate': begin
                    self.oModel3Manip->MouseUp, [event.x, event.y], oWindow
                    end
                else: result = 0b
                endcase
            end
        self.mmb_is_down = 0b
        end
    else:
    endcase

return, result
end
;--------------------------------------------------------------------
pro IDLexObjview::GetProperty, $
    rot_target=rot_target, $        ; Obsolete.
    manip_target=manip_target, $    ; Obsolete.
    subject=subject, $              ; Obsolete.
    selected=selected, $
    _ref_extra=e
;
;Note: the ability to get properties ROT_TARGET & MANIP_TARGET violates
;an Object Oriented Programming design tenet.  That tenet is: "Avoid member
;functions that return pointers or references to members less accessible
;than themselves" (Scott Meyers, "Effective C++", 1992 Addison-Wesly
;Publishing ISBN 0-201-56364-9).  Self's manipulation and rotation targets
;are protected data members (as are all object data members in IDL), and
;thus, according to the above tenet, should not be accessible via
;GetProperty. The MANIP_TARGET and ROT_TARGET keyword functionality
;is included here, however, for compatibility with earlier versions
;of this class.
;
rot_target = self.oModel1
manip_target = self.oModel3
;
;Note: the SUBJECT property is obsolete.  Use self's Get
;method instead.  The SUBJECT property is supported here for
;compatibility with older code.
;
subject = self.oModel3->Get(/all)
;
selected = *self.pSelected
self->IDLexInscribingView::GetProperty, _extra=e
end
;--------------------------------------------------------------------
pro IDLexObjview::SetProperty, $
    subject=subject, $  ; Obsolete.
    scale=scale, $      ; Reset to this size.  Scalar or 3-element array.
    mode=mode, $        ; strings: 'zoom', 'rotate', 'translate'
    _extra=e

on_error, 2
;
;Note: the SUBJECT property is obsolete.  Instead, use self's Add
;method with /ALIAS.  The SUBJECT property is supported here for
;compatibility with older code.
;
if n_elements(subject) gt 0 then begin
    if self.oModel3->Count() gt 0 then $
        self.oModel3->Remove, /all

    self.oModel3->Add, subject, /alias
    end
;
if n_elements(mode) gt 0 then begin
    self.mode = strlowcase(mode)
    end

if n_elements(scale) gt 0 then begin
    if scale[0] eq 0 then begin
        self.scale = 1./sqrt(3)
        end $
    else begin
        self.scale = scale
        end
    end

self->IDLexInscribingView::SetProperty, _extra=e
end

;--------------------------------------------------------------------
pro IDLexObjview__define
compile_opt idl2, hidden

struct_hide, {IDLexObjview, $
    inherits IDLexInscribingView, $
    mode: '', $                 ; 'zoom', 'rotate', ...
    debug: 0b, $
    xrange: dblarr(2), $
    yrange: dblarr(2), $
    zrange: dblarr(2), $
    x_len: 0.0d, $
    y_len: 0.0d, $
    z_len: 0.0d, $
    mmb_is_down: 0b, $          ; Modal mouse button.
    axes_are_isotropic: 0b, $
    scale: dblarr(3), $
    pSelected: ptr_new(), $     ; Array of object references.
    oDefaultLights: objarr(2), $
    oStationaryModel: obj_new(), $
    oModel0: obj_new(), $
    oModel1: obj_new(), $
    oModel2: obj_new(), $
    oModel3: obj_new(), $
    oModel1Manip: obj_new(), $
    oModel3Manip: obj_new(), $
    model0_range: dblarr(3,2), $  ; [xmin, ymin, zmin, xmax, ymax, zmax]
    model1_range: dblarr(3,2), $
    model2_range: dblarr(3,2), $
    model3_range: dblarr(3,2), $
    stationary_range: dblarr(3,2), $
    got_model0_range: 0b, $
    got_model1_range: 0b, $
    got_model2_range: 0b, $
    got_model3_range: 0b, $
    got_stationary_range: 0b, $
    oViewManip: obj_new(), $
    oGetList: obj_new(), $
    oRemoveList: obj_new() $
    }
end
;--------------------------------------------------------------------
pro IDLexObjview__example
compile_opt idl2
;
;Purpose: Provide an example of how to use class IDLexObjview.  This
;example shows two views together (next to each other) in one XOBJVIEW scene.
;
oSurface1 = obj_new('IDLgrSurface', $
    beselj(shift(dist(40), 20, 20) / 2, 0) * 20, $
    color=[240, 0, 240], $
    shading=1, $
    name='Surface 1', $
    style=2 $
    )
oSurface2 = obj_new('IDLgrSurface', $
    beselj(shift(dist(40), 20, 20) / 4, 0) * 10, $
    color=[0, 240, 240], $
    shading=1, $
    name='Surface 2', $
    style=2 $
    )

oView1 = obj_new('IDLexObjview', $
    location=[0, 0], $
    dimensions=[.5, 1], $
    units=3 $ ; Normalized
    )
oView1->Add, oSurface1

oView2 = obj_new('IDLexObjview', $
    location=[.5, 0], $
    dimensions=[.5, 1], $
    units=3 $ ; Normalized
    )
oView2->Add, oSurface2
;
;Note: XOBJVIEW's ability to take view arguments is an undocumented feature.
;The feature may change or be removed in future versions of IDL.
;
xobjview, [oView1, oView2], /block

obj_destroy, [oView1, oView2]
end
