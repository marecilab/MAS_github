;$Id: //depot/idl/IDL_64/idldir/lib/utilities/idlexvolview__define.pro#1 $
;
;  Copyright (c) 1997-2007, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
;
pro IDLexVolView::Normalize, $
    oDestDevice, $
    xrange=xrange, $
    yrange=yrange, $
    zrange=zrange, $
    adjust_axes=adjust_axes
compile_opt hidden
;
;Ensure that self's z-axis ticktext is calculated automatically.
;
if keyword_set(adjust_axes) then begin
    self.oAxis[2]->SetProperty, ticktext=obj_new()
    end
;
;Normalize as usual.
;
self->IDLexObjView::Normalize, $
    oDestDevice, $
    xrange=xrange, $
    yrange=yrange, $
    zrange=zrange, $
    adjust_axes=adjust_axes
;
;If the fist z-axis tick label looks like it might clash with the last
;y-axis tick label, remove it.
;
if keyword_set(adjust_axes) then begin
    self.oAxis[2]->GetProperty, $
        ticktext=oTickText, $
        range=range, $
        tickvalues=tickvalues
    closeness_factor = .1 ; Estimate.  Made somewhat arbitrarily.
    if abs(range[0] - tickvalues[0]) $ ; A crude test.
        lt abs(tickvalues[1] - tickvalues[0]) * closeness_factor $
    then begin
        oTickText->GetProperty, strings=strings
        strings[0] = ''
        oTickText->SetProperty, strings=strings
        end
    end

end
;--------------------------------------------------------------------
; Routine to update texture mapped image display of slices
pro IDLexVolView::SampleOrthoPlane, axis
compile_opt Hidden

self.oVolume->GetProperty, $
    data0=data0, $
    rgb_table0=rgb_table0, $
    opacity_table0=opacity_table0, $
    /no_copy
sz = size(data0)
case axis OF
    0: begin
        slice = self.slice_index[0]
        do_transparent = self.transparent_x_image
        img = data0[slice, *, *]
        img = reform(img,sz[2],sz[3],/overwrite)
        fTxCoords = [[0,0],[1,0],[1,1],[0,1]]
        verts = [ $
            [slice, 0,     0    ], $
            [slice, sz[2], 0    ], $
            [slice, sz[2], sz[3]], $
            [slice, 0,     sz[3]] $
            ]
       end
    1: begin
        slice = self.slice_index[1]
        do_transparent = self.transparent_y_image
        img = data0[*, slice, *]
        img = reform(img,sz[1],sz[3],/overwrite)
        fTxCoords = [[0,0],[1,0],[1,1],[0,1]]
        verts = [[0,slice,0],[sz[1],slice,0], $
                 [sz[1],slice,sz[3]],[0,slice,sz[3]]]
       end
    2: begin
        slice = self.slice_index[2]
        do_transparent = self.transparent_z_image
        img = data0[*, *, slice]
        img = reform(img,sz[1],sz[2],/overwrite)
        fTxCoords = [[0,0],[1,0],[1,1],[0,1]]
        verts = [[0,0,slice],[sz[1],0,slice], $
                 [sz[1],sz[2],slice],[0,sz[2],slice]]
       end
endcase
self.oVolume->SetProperty, data0=data0, /no_copy
; Convert to 3xNxM or 4xNxM.
sz = size(img)

if do_transparent then begin
    rgbtab = bytarr(4, 256)
    rgbtab[0:2,*] = transpose(rgb_table0)
    rgbtab[3, *] = opacity_table0
    end $
else begin
    rgbtab = transpose(rgb_table0)
    end
rgb = rgbtab[*,img]
rgb = reform(rgb, do_transparent ? 4 : 3, sz[1], sz[2], /overwrite)

case axis of
    0: begin
        self.oXImage->SetProperty, data=rgb
        self.oSlices[0]->SetProperty, data=verts, texture_coord=fTxCoords
        end
    1: begin
        self.oYImage->SetProperty, data=rgb
        self.oSlices[1]->SetProperty, data=verts, texture_coord=fTxCoords
        end
    2: begin
        self.oZImage->SetProperty, data=rgb
        self.oSlices[2]->SetProperty, data=verts, texture_coord=fTxCoords
        end
    endcase

end

;--------------------------------------------------------------------
pro IDLexVolView::Add, oObj, _extra=e
compile_opt hidden

on_error, keyword_set(self.debug) ? 0 : 2

if min(obj_valid(oObj)) eq 0 then $
    message, 'Argument contains an invalid object.'

for i=0,n_elements(oObj)-1 do begin
    if obj_isa(oObj[i], 'IDLgrVolume') then $
        oVolume = oObj[i]
    end

if obj_valid(oVolume) then begin
    self->SetProperty, volume=oVolume
    self->SetProperty, isosurface_threshold=self.isosurface_threshold
    self->SampleOrthoPlane, 0
    self->SampleOrthoPlane, 1
    self->SampleOrthoPlane, 2
    end

self->IDLexObjview::Add, oObj, _extra=e
end
;--------------------------------------------------------------------
pro IDLexVolView::SetProperty, $
    volume=oVolume, $
;   pickable_coordinates=pickable_coordinates, $
;   pickable_values=pickable_values, $
    c_value=c_value, $
    hide_x_contour=hide_x_contour, $
    hide_y_contour=hide_y_contour, $
    hide_z_contour=hide_z_contour, $
    hide_x_image=hide_x_image, $
    hide_y_image=hide_y_image, $
    hide_z_image=hide_z_image, $
    hide_x_axis=hide_x_axis, $
    hide_y_axis=hide_y_axis, $
    hide_z_axis=hide_z_axis, $
    hide_isosurface=hide_isosurface, $
    hide_volume=hide_volume, $
    x_slice=x_slice, $
    y_slice=y_slice, $
    z_slice=z_slice, $
    stippled_isosurface=stippled_isosurface, $
    isosurface_threshold=isosurface_threshold, $
    isosurface_color=isosurface_color, $
    transparent_x_image=transparent_x_image, $
    transparent_y_image=transparent_y_image, $
    transparent_z_image=transparent_z_image, $
    color=color, $
    interpolate_slices=interpolate_slices, $
    _extra=e
compile_opt hidden

self->IDLexObjview::SetProperty, color=color, _extra=e

if n_elements(color) gt 0 then begin
    mute = 150
    threshold = 85 ; Somewhat arbitrary.
    if ct_luminance(color[0], color[1], color[2]) lt threshold then begin
        for i=0,n_elements(self.oAxis)-1 do begin
            self.oAxis[i]->SetProperty, color=[255,255,255]
            self.oXContourOutline->SetProperty, color=[255,255,255]-mute
            self.oYContourOutline->SetProperty, color=[255,255,255]-mute
            self.oZContourOutline->SetProperty, color=[255,255,255]-mute
            end
        end
    threshold = 170 ; Somewhat arbitrary.
    if ct_luminance(color[0], color[1], color[2]) gt threshold then begin
        for i=0,n_elements(self.oAxis)-1 do begin
            self.oAxis[i]->SetProperty, color=[0,0,0]
            self.oXContourOutline->SetProperty, color=[0,0,0]+mute
            self.oYContourOutline->SetProperty, color=[0,0,0]+mute
            self.oZContourOutline->SetProperty, color=[0,0,0]+mute
            end
        end
    end

if n_elements(oVolume) gt 0 then begin
    if size(oVolume, /dimensions) gt 1 then $
        message, 'VOLUME property must be a scalar'
    if not obj_valid(oVolume) then $
        message, 'specified VOLUME is not valid.'
    if not obj_isa(oVolume, 'IDLgrVolume') then $
        message, 'specified VOLUME must be an IDLgrVolume.'
    self.oVolume = oVolume
    oVolume->GetProperty, data0=data0, /no_copy
    if n_elements(data0) gt 1 then begin
        siz = size(data0, /dimensions)
        x_slice = self.slice_index[0] < siz[0]
        y_slice = self.slice_index[1] < siz[1]
        z_slice = self.slice_index[2] < siz[2]
        end
    oVolume->SetProperty, data0=data0, /no_copy
    end

if n_elements(pickable_coordinates) gt 0 then begin
    self.odataset1->SetProperty, data=pickable_coordinates
    *self.pInVectorX = reform(coordinates[0, *])
    *self.pInVectorY = reform(coordinates[1, *])
    *self.pInVectorZ = reform(coordinates[2, *])
    end

if n_elements(pickable_values) gt 0 then begin
    *self.pInVectorV = values
    end

if n_elements(c_value) ne 0 then begin
    *self.pContourValues = c_value
    self.n_contours = n_elements(c_value)
    self.oXContour->SetProperty, c_value=c_value, $
        n_levels=self.n_contours
    self.oYContour->SetProperty, c_value=c_value, $
        n_levels=self.n_contours
    self.oZContour->SetProperty, c_value=c_value, $
        n_levels=self.n_contours
    endif

if n_elements(hide_x_contour) eq 1 then begin
    self.oXContour->SetProperty, hide=keyword_set(hide_x_contour)
    self.oSlices[0]->GetProperty, hide=image_hide
    self.oXContourOutline->SetProperty, $
        hide=keyword_set(hide_x_contour) or image_hide eq 0
    endif
if n_elements(hide_y_contour) eq 1 then begin
    self.oYContour->SetProperty, hide=keyword_set(hide_y_contour)
    self.oSlices[1]->GetProperty, hide=image_hide
    self.oYContourOutline->SetProperty, $
        hide=keyword_set(hide_y_contour) or image_hide eq 0
    endif
if n_elements(hide_z_contour) eq 1 then begin
    self.oZContour->SetProperty, hide=keyword_set(hide_z_contour)
    self.oSlices[2]->GetProperty, hide=image_hide
    self.oZContourOutline->SetProperty, $
        hide=keyword_set(hide_z_contour) or image_hide eq 0
    endif

if n_elements(hide_x_image) eq 1 then begin
    self.oSlices[0]->SetProperty, hide=keyword_set(hide_x_image)
    self.oXContour->GetProperty, hide=contour_hide
    self.oXContourOutline->SetProperty, $
        hide=contour_hide or keyword_set(hide_x_image) eq 0
    endif
if n_elements(hide_y_image) eq 1 then begin
    self.oSlices[1]->SetProperty, hide=keyword_set(hide_y_image)
    self.oYContour->GetProperty, hide=contour_hide
    self.oYContourOutline->SetProperty, $
        hide=contour_hide or keyword_set(hide_y_image) eq 0
    endif
if n_elements(hide_z_image) eq 1 then begin
    self.oSlices[2]->SetProperty, hide=keyword_set(hide_z_image)
    self.oZContour->GetProperty, hide=contour_hide
    self.oZContourOutline->SetProperty, $
        hide=contour_hide or keyword_set(hide_z_image) eq 0
    endif

if n_elements(hide_volume) eq 1 then begin
    self.oVolume->SetProperty, hide=keyword_set(hide_volume)
    endif

if n_elements(x_slice) eq 1 then begin
    self.slice_index[0] = x_slice
    self->ContourOrthoPlane, 0
    self->SampleOrthoPlane, 0
    endif
if n_elements(y_slice) eq 1 then begin
    self.slice_index[1] = y_slice
    self->ContourOrthoPlane, 1
    self->SampleOrthoPlane, 1
    endif
if n_elements(z_slice) eq 1 then begin
    self.slice_index[2] = z_slice
    self->ContourOrthoPlane, 2
    self->SampleOrthoPlane, 2
    endif

if n_elements(hide_isosurface) eq 1 then begin
    self.oIsoPolygon->SetProperty, hide=keyword_set(hide_isosurface)
    endif

if n_elements(stippled_isosurface) eq 1 then begin
    fill_pattern = $
        keyword_set(stippled_isosurface) ? self.oStipple : obj_new()
    self.oIsoPolygon->SetProperty, fill_pattern=fill_pattern
    endif

if n_elements(isosurface_threshold) gt 0 then begin
    self.oVolume->GetProperty, data0=data0, /no_copy

    empty_indx = where(finite(data0) eq 0)
    If empty_indx[0] ne -1 then begin
        data0[empty_indx] = min(data0) - 1.e9
        endif
    ;IsoSurface, data0, isosurface_threshold, verts, conn
    shade_volume, data0, isosurface_threshold, verts, conn
    self.oVolume->SetProperty, data0=data0, /no_copy
    if (n_elements(verts) le 0) then begin
        verts=BYTARR(3,3)
        conn=0
        endif
    self.oIsoPolygon->SetProperty, data=verts, polygons=conn
    self.isosurface_threshold = isosurface_threshold
    endif

if n_elements(isosurface_color) gt 0 then begin
    self.oIsoPolygon->SetProperty, color=isosurface_color
    end

if n_elements(transparent_x_image) gt 0 then begin
    self.transparent_x_image = keyword_set(transparent_x_image)
    end

if n_elements(transparent_y_image) gt 0 then begin
    self.transparent_y_image = keyword_set(transparent_y_image)
    end

if n_elements(transparent_z_image) gt 0 then begin
    self.transparent_z_image = keyword_set(transparent_z_image)
    end

if n_elements(hide_x_axis) gt 0 then begin
    self.oAxis[0]->SetProperty, hide=hide_x_axis
    end

if n_elements(hide_y_axis) gt 0 then begin
    self.oAxis[1]->SetProperty, hide=hide_y_axis
    end

if n_elements(hide_z_axis) gt 0 then begin
    self.oAxis[2]->SetProperty, hide=hide_z_axis
    end

if n_elements(interpolate_slices) gt 0 then begin
    self.interpolate_slices = interpolate_slices
    for i=0,2 do begin
        self.oSlices[i]->SetProperty, texture_interp=interpolate_slices
        end
    end

end
;--------------------------------------------------------------------
pro IDLexVolView::GetProperty, $
    picked=picked, $
    volume=volume, $
    isosurface_color=isosurface_color, $
    isosurface_threshold=isosurface_threshold, $
    _ref_extra=e
compile_opt hidden

if arg_present(picked) then $
    picked = *self.pPicked

volume = self.oVolume
isosurface_threshold = self.isosurface_threshold
self.oIsoPolygon->GetProperty, color=isosurface_color
self->IDLexObjview::GetProperty, _extra=e
end
;--------------------------------------------------------------------
function IDLexVolView::Update, event
compile_opt hidden
;
;Purpose: Handle given event.  Return 1 of the event changed
;   self's public state, else return 0.  EVENT must be a structure.
;   Ignore (no-op) EVENTs that are not pertinent.
;
result = self->IDLexObjview::Update(event)
;
;Ignore non-Draw-Widget events.
;
if (tag_names(event, /structure_name) ne 'WIDGET_DRAW') then $
    return, result
;
;Handle other events.
;
if event.press eq 4 then begin ; Right mouse.
    widget_control, event.id, get_value=oWindow
    oSelected = oWindow->Select(self, [event.x, event.y])
    if obj_valid(oSelected[0]) then begin
        if oSelected[0] eq self.odataSet1 then begin
            oPicked = oWindow->Pickdata( $
                self,$
                self.odataSet1, $
                [event.x, event.y], $
                dataxyz $
                )
            distance = sqrt( $
                (dataXYZ[0] - *self.pInVectorX)^2 + $
                (dataXYZ[1] - *self.pInVectorY)^2 + $
                (dataXYZ[2] - *self.pInVectorZ)^2 $
                )
            void = min(distance, index)
            *self.pPicked = { $
                x: (*self.pinVectorX)[index], $
                y: (*self.pinVectorY)[index], $
                z: (*self.pinVectorZ)[index], $
                value: (*self.pinVectorV)[index] $
                }
            result = 1b
            print,'Value:',(*self.pinVectorV)[index]
            endif
        endif
    endif
return, result
end
;--------------------------------------------------------------------
pro IDLexVolView::ContourOrthoPlane, axis
compile_opt Hidden
;
;Purpose: generate contours on given slice.
;
self.oVolume->GetProperty, data0=data, /no_copy
sz = size(data, /dimensions)
slice = self.slice_index[axis]
case axis OF
    0: begin
        img = data[slice,*,*]
        img = reform(img,sz[1],sz[2],/overwrite)
        self.oXContour->SetProperty, geomz=slice, data=img
        self.oXContourOutline->SetProperty, data=[ $
            [0, 0, slice], $
            [0, sz[2], slice], $
            [sz[1], sz[2], slice], $
            [sz[1], 0, slice], $
            [0, 0, slice] $
            ]
        end
    1: begin
        img = data[*,slice,*]
        img = reform(img,sz[0],sz[2],/overwrite)
        self.oYContour->SetProperty, geomz=-slice, data=img
        self.oYContourOutline->SetProperty, data=[ $
            [0, 0, -slice], $
            [0, sz[2], -slice], $
            [sz[0], sz[2], -slice], $
            [sz[0], 0, -slice], $
            [0, 0, -slice] $
            ]
        end
    2: begin
        img = data[*,*,slice]
        img = reform(img,sz[0],sz[1],/overwrite)
        self.oZContour->SetProperty, geomz=slice, data=img
        self.oZContourOutline->SetProperty, data=[ $
            [0, 0, slice], $
            [sz[0], 0, slice], $
            [sz[0], sz[1], slice], $
            [0, sz[1], slice], $
            [0, 0, slice] $
            ]
        end
    ELSE:
endcase
self.oVolume->SetProperty, data0=data, /no_copy
end
;--------------------------------------------------------------------
pro IDLexVolView::cleanup
compile_opt hidden

self->IDLexVolView::CleanupNoninheritedMembers
self->IDLexObjview::cleanup
end
;--------------------------------------------------------------------
pro IDLexVolView::CleanupNoninheritedMembers
compile_opt hidden

obj_destroy, self.oIsoPolygon
obj_destroy, self.oXImage
obj_destroy, self.oYImage
obj_destroy, self.oZImage
obj_destroy, self.oSlices
obj_destroy, self.odataset1
obj_destroy, self.oXtitle
obj_destroy, self.oYtitle
obj_destroy, self.oZtitle
obj_destroy, self.oBox
obj_destroy, self.oXModel
obj_destroy, self.oYModel
obj_destroy, self.oZModel
obj_destroy, self.oXContour
obj_destroy, self.oYContour
obj_destroy, self.oZContour
obj_destroy, self.oStipple

ptr_free, self.pInVectorX
ptr_free, self.pInVectorY
ptr_free, self.pInVectorZ
ptr_free, self.pInVectorV
ptr_free, self.pPicked
ptr_free, self.pContourValues

end
;--------------------------------------------------------------------
pro IDLexVolView::ConstructGraphics
compile_opt Hidden
;
; Add in the objects which may be made transparent
; or semi-transparent.  The transparent visualization
; property depends on the order in which atoms are
; added to a model, and not just on their Z order.
;
;
;Create the isosurface object.
;
self.oIsoPolygon = obj_new('IDLgrPolygon', color=[127, 127, 127], hide=1,$
    shading=1, name='isosurface')
self.oModel3->Add, self.oIsoPolygon
;
;Create the slice image-form objects (texture mapped polygons).
;
self.oXImage = obj_new('IDLgrImage',dist(5))
self.oYImage = obj_new('IDLgrImage',dist(5))
self.oZImage = obj_new('IDLgrImage',dist(5))

depth_offset = 10 ; Somewhat arbitrary
oXPolygon = obj_new('IDLgrPolygon', color=[255,255,255], hide=1,$
    texture_map=self.oXImage, name='x slice', $
    texture_interp=self.interpolate_slices, depth_offset=depth_offset)
self.oModel3->Add, oXPolygon

oYPolygon = obj_new('IDLgrPolygon', color=[255,255,255], hide=1,$
    texture_map=self.oYImage, name='y slice', $
    texture_interp=self.interpolate_slices, depth_offset=depth_offset)
self.oModel3->Add, oYPolygon

oZPolygon = obj_new('IDLgrPolygon', color=[255,255,255], hide=1,$
    texture_map=self.oZImage, name='z slice', $
    texture_interp=self.interpolate_slices, depth_offset=depth_offset)
self.oModel3->Add, oZPolygon

self.oSlices = [oXPolygon, oYPolygon, oZPolygon]
;
;Add the polygon which contains the individual sampling points.
;
self.odataset1 = obj_new('IDLgrPolygon', Style = 0, color = [0, 255, 255], $
    Thick = 2.5)
self.oModel3->Add, self.odataset1

self.oXtitle = obj_new('IDLgrText', '<- X ->')
self.oAxis[0] = obj_new('IDLgrAxis', 0, /exact, $
    title=self.oXtitle, ticklen=.1, name='x axis')

self.oYtitle = obj_new('IDLgrText','<- Y ->')
self.oAxis[1] = obj_new('IDLgrAxis', 1, /exact, $
    title=self.oYtitle, ticklen=.1, name='y axis')

self.oZtitle = obj_new('IDLgrText','')
self.oAxis[2] = obj_new('IDLgrAxis', 2, /exact, $
    title=self.oZtitle, location=[0,1,0], ticklen=.1, name='z axis')

self.oModel3->Add, self.oAxis
;
;Create the wireframe box.
;
self.oBox = obj_new('IDLgrpolyline', $
    color=[200,200,200], $
    polyline=[5,0,1,3,2,0,5,4,5,7,6,4,2,0,4,2,1,5,2,2,6,2,3,7])
self.oModel3->Add, self.oBox
;
;Create the contour objects HC-11/24/98
;
self.oXModel = obj_new('IDLgrModel')
self.oModel3->Add, self.oXModel
self.oYModel = obj_new('IDLgrModel')
self.oModel3->Add, self.oYModel
self.oZModel = obj_new('IDLgrModel')
self.oModel3->Add, self.oZModel

self.oXContour = obj_new('IDLgrContour', $
    color = [250,40,40], /planar, $
    n_levels = self.n_contours, hide=1)
self.oXcontourOutline = obj_new('IDLgrPolyline', hide=1, color=[150,150,150])
self.oXModel->Add, [self.oXContour, self.oXContourOutline]

self.oYContour = obj_new('IDLgrContour', $
    color = [40,250,40], /planar, $
    n_levels = self.n_contours, hide=1)
self.oYcontourOutline = obj_new('IDLgrPolyline', hide=1, color=[150,150,150])
self.oYModel->Add, [self.oYContour, self.oYContourOutline]

self.oZContour = obj_new('IDLgrContour', $
    color = [40,40,250], /planar, $
    n_levels = self.n_contours, hide=1)
self.oZcontourOutline = obj_new('IDLgrPolyline', hide=1, color=[150,150,150])
self.oZModel->Add, [self.oZContour, self.oZContourOutline]

self.oXModel->Rotate, [0,0,1], 90
self.oXModel->Rotate, [0,1,0], 90
self.oYModel->Rotate, [1,0,0], 90

self.oStipple = obj_new('IDLgrPattern', 2, pattern=stipple_pattern(8))
end
;--------------------------------------------------------------------
function IDLexVolView::init, _extra=e

compile_opt hidden

if not self->IDLexObjview::init(_extra=e) then begin
    message, 'Failed to init IDLexObjview part of self.'
    end

self->ConstructGraphics
self.pPicked = ptr_new(/allocate_heap)
self.pInVectorX = ptr_new(/allocate_heap)
self.pInVectorY = ptr_new(/allocate_heap)
self.pInVectorZ = ptr_new(/allocate_heap)
self.pInVectorV = ptr_new(/allocate_heap)
self.pContourValues = ptr_new(/allocate_heap)

return, 1 ; Success.
end
;--------------------------------------------------------------------
pro IDLexVolView__define
compile_opt idl2, hidden

struct_hide, {IDLexVolView, $
    inherits IDLexObjview, $
    oVolume: obj_new(), $
    oIsoPolygon: obj_new(), $
    oXImage: obj_new(), $
    oYImage: obj_new(), $
    oZImage: obj_new(), $
    oSlices: objarr(3), $
    odataset1: obj_new(), $
    oXtitle: obj_new(), $
    oYtitle: obj_new(), $
    oZtitle: obj_new(), $
    oBox: obj_new(), $
    oXModel: obj_new(), $
    oYModel: obj_new(), $
    oZModel: obj_new(), $
    oXContour: obj_new(), $
    oYContour: obj_new(), $
    oZContour: obj_new(), $
    oXContourOutline: obj_new(), $
    oYContourOutline: obj_new(), $
    oZContourOutline: obj_new(), $
    oStipple: obj_new(), $
    oAxis: objarr(3), $
    pPicked: ptr_new(), $
    pInVectorX: ptr_new(), $
    pInVectorY: ptr_new(), $
    pInVectorZ: ptr_new(), $
    pInVectorV: ptr_new(), $
    pContourValues: ptr_new(), $
    n_contours: 0, $
    interpolate_slices: 0b, $
    transparent_x_image: 0b, $
    transparent_y_image: 0b, $
    transparent_z_image: 0b, $
    isosurface_threshold: 0.0d, $
    slice_index: lonarr(3) $
    }
end

