;$Id: //depot/idl/IDL_64/idldir/lib/utilities/idlexvolviewwid__define.pro#1 $
;
;  Copyright (c) 1997-2007, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
function IDLexVolViewWid::DefaultViewType
compile_opt hidden

return, 'IDLexVolView'
end
;--------------------------------------------------------------------
pro IDLexVolViewWid::SetSliderRanges
compile_opt hidden

oView = self.oViewGroup->Get(/Current)
oView->GetProperty, Volume=oVolume
oVolume->GetProperty, data0=data0, /no_copy
data_dims = size(data0, /dimensions)
widget_control, self.wXSlider, set_slider_max=data_dims[0] - 1
widget_control, self.wYSlider, set_slider_max=data_dims[1] - 1
widget_control, self.wZSlider, set_slider_max=data_dims[2] - 1
oVolume->SetProperty, data0=data0, /no_copy
end
;--------------------------------------------------------------------
pro IDLexVolViewWid__XRGBSlider_event, event
compile_opt hidden

widget_control, event.top, get_uvalue=pState

if event.id eq (*pState).wApply $
or event.id eq (*pState).wOK then begin
    widget_control, (*pState).wRGBSlider, get_value=value
    (*pState).oViewgroup->CallAll, 'SetProperty', $
        isosurface_color=value, $
        isa='IDLexVolview'
    if (*pState).swatch_changed then begin
        (*pState).oVolViewWid->Draw, /hourglass
        (*pState).applied = 1b
        (*pState).swatch_changed = 0b
        end
    end

if event.id eq (*pState).wCancel then begin
    (*pState).oViewgroup->CallAll, 'SetProperty', $
        isosurface_color=(*pState).isosurface_color, $
        isa='IDLexVolview'
    if (*pState).applied then begin
        (*pState).oVolViewWid->Draw, /hourglass
        end
    end

case 1 of
    (*pState).init_event: begin
        (*pState).init_event = 0b
        end
    tag_names(event, /structure_name) eq 'RGB_EVENT': begin
        (*pState).swatch_changed = 1b
        end
    else:
    endcase

if event.id eq (*pState).wOK $
or event.id eq (*pState).wCancel then begin
    WIDGET_CONTROL, event.top, /DESTROY
    end

end
;--------------------------------------------------------------------
pro IDLexVolViewWid__XRGBSlider_cleanup, wID
widget_control, wID, get_uvalue=pState

ptr_free, pState

end
;--------------------------------------------------------------------
pro IDLexVolViewWid::XRGBSlider, group_leader
compile_opt hidden

if widget_info(self.wRGBSliderTlb, /valid_id) then begin
    widget_control, self.wRGBSliderTllb, /show
    return
    end

prefix = obj_class(self) + 'IsoColor:'
tlb = WIDGET_BASE( $  ; Top-Level Base.
    /COLUMN, $
    TITLE='Isosurface Color', $
    GROUP_LEADER=group_leader, $
    UNAME=prefix + 'tlb', $
    MODAL=widget_info(group_leader, /modal) $
    )
self.wRGBSliderTlb = tlb
self.oViewgroup->CallCurrent, 'GetProperty', $
    isosurface_color=isosurface_color

wRGBSlider = cw_rgbslider( $
    tlb, $
    value=isosurface_color, $
    graphics_level=2, $
    /hls, $
    uname=prefix + 'pal_edit', $
    /frame $
    )

wRowBase = WIDGET_BASE(tlb, /ROW, /GRID)
wOK = WIDGET_BUTTON(wRowBase, VALUE='OK', UNAME=prefix + 'OK')
wCancel = WIDGET_BUTTON( $
    wRowBase, $
    value='Cancel', $
    uname=prefix + 'cancel' $
    )
wApply = widget_button(wRowBase, value='Apply')

if not widget_info(group_leader, /modal) then $
    WIDGET_CONTROL, tlb, MAP=0

WIDGET_CONTROL, tlb, /REALIZE

if not widget_info(group_leader, /modal) then begin
    tlb_geom = WIDGET_INFO(tlb, /GEOMETRY)
    leader_geom = WIDGET_INFO(group_leader, /GEOMETRY)
    DEVICE, GET_SCREEN_SIZE=screen_size

    x = leader_geom.xoffset
    x = x < (screen_size[0] - tlb_geom.scr_xsize)
    x = x > 0

    y = leader_geom.scr_ysize[0] + leader_geom.yoffset
    y = y < (screen_size[1] - tlb_geom.scr_ysize)
    y = y > 0

    WIDGET_CONTROL, tlb, TLB_SET_XOFFSET=x, TLB_SET_YOFFSET=y
    WIDGET_CONTROL, tlb, MAP=1
endif

WIDGET_CONTROL, tlb, SET_UVALUE=PTR_NEW({ $
    isosurface_color: isosurface_color, $
    oVolViewWid: self, $
    oViewgroup: self.oViewgroup, $
    applied: 0b, $  ; 1 == yes, 0 == no.
    init_event: 1b, $
    swatch_changed: 0b, $
    wRGBSlider: wRGBSlider, $
    wOK: wOK, $
    wCancel: wCancel, $
    wApply: wApply $
    })

XMANAGER, $
    'IDLexVolViewWid__XRGBSlider', $
    tlb, $
    CLEANUP='IDLexVolViewWid__XRGBSlider_Cleanup', $
    /NO_BLOCK

end
;--------------------------------------------------------------------
pro IDLexVolViewWid__XPaletteEdit_event, event
compile_opt hidden

widget_control, event.top, get_uvalue=pState

if event.id eq (*pState).wApply $
or event.id eq (*pState).wOK then begin
    widget_control, (*pState).wPaletteEditor, get_value=array
    oViews = (*pState).oViewgroup->Get(/all)
    if not obj_valid(oViews) then $
        return
    indx = where(obj_isa(oViews, 'IDLexVolView'))
    if indx[0] eq -1 then $
        return
    oViews = oViews[indx]
    for i=0,n_elements(oViews)-1 do begin
        oViews[i]->GetProperty, volume=oVolume
        oVolume->SetProperty, $
            rgb_table0=transpose(array[0:2, *]), $
            opacity_table0=array[3, *]
        end
    (*pState).oViewgroup->CallAll, 'SampleOrthoPlane', 0, isa='IDLexVolView'
    (*pState).oViewgroup->CallAll, 'SampleOrthoPlane', 1, isa='IDLexVolView'
    (*pState).oViewgroup->CallAll, 'SampleOrthoPlane', 2, isa='IDLexVolView'
    if (*pState).table_changed then begin
        (*pState).oVolViewWid->Draw, /hourglass
        (*pState).applied = 1b
        (*pState).table_changed = 0b
        end
    endif

if event.id eq (*pState).wCancel then begin
    oViews = (*pState).oViewgroup->Get(/all)
    if not obj_valid(oViews) then $
        return
    indx = where(obj_isa(oViews, 'IDLexVolView'))
    if indx[0] eq -1 then $
        return
    oViews = oViews[indx]
    for i=0,n_elements(oViews)-1 do begin
        oViews[i]->GetProperty, volume=oVolume
        oVolume->SetProperty, $
            rgb_table0=(*pState).rgb_table0, $
            opacity_table0=(*pState).opacity_table0
        endfor
    (*pState).oViewgroup->CallAll, 'SampleOrthoPlane', 0, isa='IDLexVolView'
    (*pState).oViewgroup->CallAll, 'SampleOrthoPlane', 1, isa='IDLexVolView'
    (*pState).oViewgroup->CallAll, 'SampleOrthoPlane', 2, isa='IDLexVolView'
    if (*pState).applied then begin
        (*pState).oVolViewWid->Draw, /hourglass
        end
    end

if event.id eq (*pState).wPaletteEditor then begin
    (*pState).table_changed = 1b
    end

if event.id eq (*pState).wOK $
or event.id eq (*pState).wCancel then begin
    WIDGET_CONTROL, event.top, /DESTROY
    endif

end
;--------------------------------------------------------------------
pro IDLexVolViewWid__XPaletteEdit_cleanup, wID
compile_opt hidden

widget_control, wID, get_uvalue=pState
ptr_free, pState

end
;--------------------------------------------------------------------
pro IDLexVolViewWid::XPaletteEdit, group_leader
compile_opt hidden

if widget_info(self.wPalEdTlb, /valid_id) then begin
    widget_control, self.wPalEdTlb, /show
    return
    end

prefix = obj_class(self) + 'ColorEdit:'
tlb = WIDGET_BASE( $  ; Top-Level Base.
    /COLUMN, $
    TITLE='Color Table', $
    GROUP_LEADER=group_leader, $
    UNAME=prefix + 'tlb', $
    MODAL=widget_info(group_leader, /modal) $
    )
self.wPalEdTlb = tlb
self.oViewgroup->CallCurrent, 'GetProperty', volume=oVolume
oVolume->GetProperty, $
    rgb_table0=rgb_table0, $
    opacity_table0=opacity_table0

data = [transpose(rgb_table0), transpose(opacity_table0)]
wPaletteEditor = CW_PALETTE_EDITOR( $
    tlb, $
    DATA=data, $
    UNAME=prefix + 'pal_edit', $
    /FRAME $
    )

wRowBase = WIDGET_BASE(tlb, /ROW, /GRID)
wOK = WIDGET_BUTTON(wRowBase, VALUE='OK', UNAME=prefix + 'OK')
wCancel = WIDGET_BUTTON( $
    wRowBase, $
    VALUE='Cancel', $
    UNAME=prefix + 'cancel' $
    )
wApply = widget_button(wRowBase, value='Apply')

if not widget_info(group_leader, /modal) then $
    WIDGET_CONTROL, tlb, MAP=0

WIDGET_CONTROL, tlb, /REALIZE

if not widget_info(group_leader, /modal) then begin
    tlb_geom = WIDGET_INFO(tlb, /GEOMETRY)
    leader_geom = WIDGET_INFO(group_leader, /GEOMETRY)
    DEVICE, GET_SCREEN_SIZE=screen_size

    x = leader_geom.scr_xsize[0] + leader_geom.xoffset
    x = x < (screen_size[0] - tlb_geom.scr_xsize)
    x = x > 0

    WIDGET_CONTROL, tlb, TLB_SET_XOFFSET=x
    WIDGET_CONTROL, tlb, MAP=1
endif

WIDGET_CONTROL, tlb, SET_UVALUE=PTR_NEW({ $
    rgb_table0: rgb_table0, $
    opacity_table0: opacity_table0, $
    oVolViewWid: self, $
    oViewgroup: self.oViewgroup, $
    applied: 0b, $  ; 1 == yes, 0 == no.
    table_changed: 0b, $
    wPaletteEditor: wPaletteEditor, $
    wOK: wOK, $
    wCancel: wCancel, $
    wApply: wApply $
    })

XMANAGER, $
    'IDLexVolViewWid__XPaletteEdit', $
    tlb, $
    CLEANUP='IDLexVolViewWid__XPaletteEdit_Cleanup', $
    /NO_BLOCK

end
;--------------------------------------------------------------------
pro IDLexVolViewWid::SetProperty, $
    xsize=xsize, $  ; Approximate size of self's widgets.
    ysize=ysize, $  ; Approximate size of self's widgets.
    _extra=e

compile_opt hidden

pad = 7 ; Estimate.
if n_elements(ysize) gt 0 then begin
    ys = ysize
    ys = ys - pad
    endif

if n_elements(xsize) gt 0 then begin
    xs = xsize
    panel_geometry = widget_info(self.wControlPanel, /geometry)
    xs = xs - panel_geometry.scr_xsize - pad
    endif

self->IDLexObjViewWid::SetProperty, xsize=xs, ysize=ys, _extra=e

end
;--------------------------------------------------------------------
pro IDLexVolViewWid::OnRealize
compile_opt hidden

widget_control, self.wDraw, get_value=oWindow
self.oWindow = oWindow

if self.resize_when_realized then begin
    panel_geometry = widget_info(self.wControlPanel, /geometry)
    widget_control, self.wDraw, ysize=panel_geometry.scr_ysize
    end

self.oViewgroup->CallAll, 'Reset', $
    oWindow, /adjust_axes, isa='IDLexObjview'
self.oViewgroup->CallAll, 'Rotate', [0,0,1], 30, isa='IDLexObjView'
self.oViewgroup->CallAll, 'Rotate', [1,0,0], -60, isa='IDLexObjview'

end
;--------------------------------------------------------------------
pro IDLexVolViewWid::Draw, $
    oDestination_, $
    oPicture, $
    _extra=e
compile_opt hidden

self->IDLexObjviewWid::Draw, $
    oDestination, $
    oPicture, $
    in_motion=self.hide_volume, $
    _extra=e

case 1 of
    n_elements(oDestination) eq 0: begin
        self.vol_on_screen = self.hide_volume eq 0
        end
    oDestination ne self.oWindow:
    else: begin
        self.vol_on_screen = self.hide_volume eq 0
        end
    endcase

end
;--------------------------------------------------------------------
function IDLexVolViewWid::WidgetDrawBase
compile_opt hidden
;
;This routine is intended to be called from IDLexObjview::init.
;Purpose: create a base widget to contain self's draw widget.
;Side effect: create a column base for widgets to the left
;of self's draw widget.
;
wRowBase = widget_base(self.wBase, /row, xpad=0)
self.wControlPanel = widget_base(wRowBase, /column, /frame)
return, wRowBase
end
;--------------------------------------------------------------------
function IDLexVolViewWid::HandleEvent, event
compile_opt hidden

case event.id of
    self.wXSlider: begin
        self.oViewgroup->CallAll, 'SetProperty', $
            x_slice=event.value, $
            isa='IDLexVolView'
        self.oViewgroup->CallAll, 'ContourOrthoPlane', 0, isa='IDLexVolView'
        self->Draw, /hourglass
        end
    self.wYSlider: begin
        self.oViewgroup->CallAll, 'SetProperty', $
            y_slice=event.value, $
            isa='IDLexVolView'
        self.oViewgroup->CallAll, 'ContourOrthoPlane', 1, isa='IDLexVolView'
        self->Draw, /hourglass
        end
    self.wZSlider: begin
        self.oViewgroup->CallAll, 'SetProperty', $
            z_slice=event.value, $
            isa='IDLexVolView'
        self.oViewgroup->CallAll, 'ContourOrthoPlane', 2, isa='IDLexVolView'
        self->Draw, /hourglass
        end
    self.wXContourCheckbox: begin
        self.oViewgroup->CallAll, 'SetProperty', $
            hide_x_contour=event.select eq 0, $
            isa='IDLexVolView'
        self->Draw, /hourglass
        end
    self.wYContourCheckbox: begin
        self.oViewgroup->CallAll, 'SetProperty', $
            hide_y_contour=event.select eq 0, $
            isa='IDLexVolView'
        self->Draw, /hourglass
        end
    self.wZContourCheckbox: begin
        self.oViewgroup->CallAll, 'SetProperty', $
            hide_z_contour=event.select eq 0, $
            isa='IDLexVolView'
        self->Draw, /hourglass
        end
    self.wAutoRender: begin
        self.hide_volume = event.select eq 0
        widget_control, self.wRender, sensitive=self.hide_volume
        self.oViewgroup->CallAll, 'SetProperty', $
            hide_volume=self.hide_volume, $
            isa='IDLexVolView'
        if self.hide_volume eq 0 and self.vol_on_screen eq 0 then begin
            self->Draw, /hourglass
            end
        end
    self.wRender: begin
        self.oViewgroup->CallAll, 'SetProperty', $
            hide_volume=0, $
            isa='IDLexVolView'
        self.hide_volume = 0
        self->Draw, /hourglas, in_motion=0
        self.oViewgroup->CallAll, 'SetProperty', $
            hide_volume=1, $
            isa='IDLexVolView'
        self.hide_volume = 1
        end
    self.wIsoLevel: begin
        self.oViewgroup->CallAll, 'SetProperty', $
            isosurface_threshold=( $
                (self.data_range[1] - self.data_range[0]) / 255. $
                ) * event.value + self.data_range[0], $
            isa='IDLexVolView'
        self->Draw, /hourglass
        end
    self.wIsotoggle0: begin
        self.oViewgroup->CallAll, 'SetProperty', /hide_isosurface
        self->Draw, /hourglass
        end
    self.wIsotoggle1: begin
        self.oViewgroup->CallAll, 'SetProperty', $
            hide_isosurface=0, $
            stippled_isosurface=0
        self->Draw, /hourglass
        end
    self.wIsotoggle2: begin
        self.oViewgroup->CallAll, 'SetProperty', $
            hide_isosurface=0, $
            stippled_isosurface=1
        self->Draw, /hourglass
        end
    self.wXDroplist: begin
        self.oViewgroup->CallAll, 'SetProperty', $
            hide_x_image=event.index eq 0, $
            transparent_x_image=event.index eq 2
        if event.index ne 0 then $
            self.oViewgroup->CallAll, 'SampleOrthoPlane', 0, $
                isa='IDLexVolview'
        self->Draw, /hourglass
        end
    self.wYDroplist: begin
        self.oViewgroup->CallAll, 'SetProperty', $
            hide_y_image=event.index eq 0, $
            transparent_y_image=event.index eq 2
        if event.index ne 0 then $
            self.oViewgroup->CallAll, 'SampleOrthoPlane', 1, $
                isa='IDLexVolview'
        self->Draw, /hourglass
        end
    self.wZDroplist: begin
        self.oViewgroup->CallAll, 'SetProperty', $
            hide_z_image=event.index eq 0, $
            transparent_z_image=event.index eq 2
        if event.index ne 0 then $
            self.oViewgroup->CallAll, 'SampleOrthoPlane', 2, $
                isa='IDLexVolview'
        self->Draw, /hourglass
        end
    self.wPaletteEditorButton: begin
        self->XPaletteEdit, event.top
        end
    self.wReset: begin
        self.oViewgroup->CallAll, 'Reset', $
            self.oWindow, /adjust_axes, isa='IDLexObjview'
        self.oViewgroup->CallAll, 'Rotate', [0,0,1], 30, isa='IDLexObjview'
        self.oViewgroup->CallAll, 'Rotate', [1,0,0], -60, isa='IDLexObjview'
        self->Draw, /hourglass
        return, self->Rewrite(event)
        end
    self.wIsoColorButton: begin
        self->XRGBSlider, event.top
        end
    self.wXAxisToggle: begin
        self.hide_x_axis = keyword_set(self.hide_x_axis) eq 0
        self.oViewgroup->CallAll, 'SetProperty', hide_x_axis=self.hide_x_axis
        widget_control, $
            self.wXAxisToggle, $
            set_value=(['Hide', 'Show'])[self.hide_x_axis] + ' X Axis'
        self->Draw, /hourglass
        end
    self.wYAxisToggle: begin
        self.hide_y_axis = keyword_set(self.hide_y_axis) eq 0
        self.oViewgroup->CallAll, 'SetProperty', hide_y_axis=self.hide_y_axis
        widget_control, $
            self.wYAxisToggle, $
            set_value=(['Hide', 'Show'])[self.hide_y_axis] + 'Y Axis'
        self->Draw, /hourglass
        end
    self.wZAxisToggle: begin
        self.hide_z_axis = keyword_set(self.hide_z_axis) eq 0
        self.oViewgroup->CallAll, 'SetProperty', hide_z_axis=self.hide_z_axis
        widget_control, $
            self.wZAxisToggle, $
            set_value=(['Hide', 'Show'])[self.hide_z_axis] + ' Z Axis'
        self->Draw, /hourglass
        end

    else:
    endcase

wPrintButton = self.wPrintButton
wExportImageButton = self.wExportImageButton
wVRMLButton = self.wVRMLButton
wDraw = self.wDraw

switch event.id of
    wPrintButton:
    wExportImageButton:
    wVRMLButton: begin
        if self.vol_on_screen and self.hide_volume then begin
            self.oViewgroup->CallAll, 'SetProperty', $
                hide_volume=0, $
                isa='IDLexVolView'
            end
        break
        end
    wDraw: begin
        if event.type eq 0 then begin ; Button press
            if self.mode eq 'select' then begin
                if self.vol_on_screen and self.hide_volume then begin
                    self.oViewgroup->CallAll, 'SetProperty', $
                        hide_volume=0, $
                        isa='IDLexVolView'
                    end
                end
            end
        end
    end

result = self->IDLexObjviewWid::HandleEvent(event)

switch event.id of
    wPrintButton:
    wExportImageButton:
    wVRMLButton: begin
        if self.vol_on_screen and self.hide_volume then begin
            self.oViewgroup->CallAll, 'SetProperty', $
                hide_volume=1, $
                isa='IDLexVolView'
            end
        break
        end
    wDraw: begin
        if event.type eq 0 then begin ; Button press
            if self.mode eq 'select' then begin
                if self.vol_on_screen and self.hide_volume then begin
                    self.oViewgroup->CallAll, 'SetProperty', $
                        hide_volume=1, $
                        isa='IDLexVolView'
                    end
                end
            end
        end
    end

return, result
end
;--------------------------------------------------------------------
pro IDLexVolViewWid::cleanup
compile_opt hidden

self->IDLexObjviewWid::cleanup
end
;--------------------------------------------------------------------
function IDLexVolViewWid::init, $
    wParent, $
    oSubjects, $
    draw_xsize=draw_xsize, $
    draw_ysize=draw_ysize, $
    _extra=e
compile_opt hidden

self.resize_when_realized = n_elements(draw_ysize) eq 0

if not self->IDLexObjviewWid::init( $
    wParent, $
    oSubjects, $
    draw_xsize=n_elements(draw_xsize) eq 0 ? 550 : draw_xsize, $
    draw_ysize=n_elements(draw_ysize) eq 0 ? 1   : draw_ysize, $
    _extra=e $
    ) $
then $
    message, 'failed to initialize IDLexObjviewWid part of self.'

wAxesButton = widget_button(self.wViewButton, value='Axes', /menu)
    self.wXAxisToggle = widget_button(wAxesButton, value='Hide X Axis')
    self.wYAxisToggle = widget_button(wAxesButton, value='Hide Y Axis')
    self.wZAxisToggle = widget_button(wAxesButton, value='Hide Z Axis')

wSelectBase = WIDGET_BASE(self.wControlPanel, /ROW)
    wLeftBase = WIDGET_BASE(wSelectBase,/COLUMN, UName = 'LeftBase')
        wLabel = WIDGET_LABEL(wLeftBase,VALUE='Image Planes:', $
            UName = 'PlaneLabel')
        wGuiBase2 = WIDGET_BASE(wLeftBase,/COLUMN, UName = 'GUIBase2')
            Choices = ['<Off>','Opaque','Transparent']
            self.wXDroplist = WIDGET_DROPLIST(wGuiBase2, VALUE=Choices,$
                TITLE='X:', UName = 'XDroplist')
            self.wYDroplist = WIDGET_DROPLIST(wGuiBase2, VALUE=Choices,$
                TITLE='Y:', UName = 'YDroplist')
            self.wZDroplist = WIDGET_DROPLIST(wGuiBase2, VALUE=Choices,$
                TITLE='Z:', UName = 'ZDroplist')

    ;Add Contour Plane Controls - HC-11/24/98
    wRightBase = WIDGET_BASE(wSelectBase,/COLUMN, UName = 'RightBase')
        wLabel = WIDGET_LABEL(wRightBase,VALUE='Contours:', $
            UName = 'ContourLabel')
        wRightButtons = WIDGET_BASE(wRightBase,/COLUMN,/NONEXCLUSIVE, $
            UName = 'RightButtonBase')
            self.wXContourCheckbox = WIDGET_BUTTON(wRightButtons, value='X')
            self.wYContourCheckbox = WIDGET_BUTTON(wRightButtons, value='Y')
            self.wZContourCheckbox = WIDGET_BUTTON(wRightButtons, value='Z')

self.wXSlider = WIDGET_SLIDER(self.wControlPanel, TITLE='X Plane')
self.wYSlider = WIDGET_SLIDER(self.wControlPanel, TITLE='Y Plane')
self.wZSlider = WIDGET_SLIDER(self.wControlPanel, TITLE='Z Plane')
;
;Bold assumption: if there are multiple volume views, they all
;contain volumes that are the same size.
;
indx = where(obj_isa(oSubjects, 'IDLexVolView'))
if indx[0] ne -1 then begin
    oSubjects[indx[0]]->GetProperty, volume=oVolume
    oVolume->GetProperty, data0=data0, /no_copy
    data_dims = size(data0, /dimensions)
    widget_control, self.wXSlider, set_slider_max=data_dims[0] - 1
    widget_control, self.wYSlider, set_slider_max=data_dims[1] - 1
    widget_control, self.wZSlider, set_slider_max=data_dims[2] - 1
    oVolume->SetProperty, data0=data0, /no_copy
    end
;
;Volume Rendering Controls
;
wVolBase = WIDGET_BASE(self.wControlPanel, /COLUMN, /FRAME, UName = 'VolBase')
    wLabel = WIDGET_LABEL(wVolBase, VALUE='Volume:', UName = 'VolLavel')
    self.wPaletteEditorButton = widget_button( $
        wVolBase, $
        value='Color and Opacity...' $
        )
    wRowBase = widget_base(wVolBase, /row, /fram)
        wNonExclusiveBase = widget_base(wRowBase, /nonexclusive)
            self.wAutoRender = widget_button( $
                wNonExclusiveBase, $
                value='Auto-Render' $
                )
            widget_control, self.wAutoRender, /set_button
        wColBase = widget_base(wRowBase, /col)
            self.wRender = WIDGET_BUTTON(wColBase, VALUE='Render', $
                UName = 'RenderButton', sensitive=0)
;
;Isosurface Controls
;
PolyBase = WIDGET_BASE(self.wControlPanel, /COLUMN, /FRAME, $
    UName = 'PolyBase')
    wLabel = WIDGET_LABEL(PolyBase,VALUE='IsoSurface:', UName = 'IsoLabel')
    self.wIsoColorButton = widget_button(polybase, value='Color...')
;   self.wIsoText = WIDGET_LABEL(PolyBase, Value = ' ', Scr_XSize = 100, $
;       UName = 'IsoText')
    wGuiBase2 = WIDGET_BASE(PolyBase,/COLUMN,/EXCLUSIVE,/FRAME, $
        UName = 'GUIBase2')
        self.wIsotoggle0 = WIDGET_BUTTON(wGuiBase2,VALUE='Isosurface Off',$
            UName = 'IsoOffButton', /No_Release)
        self.wIsotoggle1 = WIDGET_BUTTON(wGuiBase2,VALUE='Opaque Isosurface',$
            UName = 'IsoOnButton', /No_Release)
        self.wIsotoggle2 = WIDGET_BUTTON(wGuiBase2,$
            VALUE='Semi-transparent Isosurface', $
            UName = 'StippleButton', /No_Release)
    self.wIsoLevel = WIDGET_SLIDER(PolyBase,/SUPPRESS_VALUE, $
        TITLE='Level', MAXIMUM=255,MINIMUM=0,VALUE=128, $
        UName = 'IsoSlider')
;
; Status line.
;
;wGuiBase = WIDGET_BASE(self.TLB, /COLUMN, UName = 'AnotherGUIBase')
;wLabel = WIDGET_LABEL(wGuiBase, /FRAME, $
;   VALUE="Left Mouse: Trackball   Right Mouse: Return Value at Data Pt.", $
;   UName = 'MouseControlLabel')
;wLabel = WIDGET_LABEL(wGuiBase, VALUE=" ", /DYNAMIC_RESIZE, $
;   UName = 'StatusLabel')
;
; Set default widget items.
;
WIDGET_CONTROL, self.wXDroplist, SET_DROPLIST_SELECT=0
WIDGET_CONTROL, self.wYDroplist, SET_DROPLIST_SELECT=0
WIDGET_CONTROL, self.wZDroplist, SET_DROPLIST_SELECT=0
WIDGET_CONTROL, self.wIsotoggle0,SET_BUTTON=1

oVolViews = self.oViewgroup->Get(isa='IDLexVolview', /all, count=count)
for i=0,count-1 do begin
;
;   Bold assumption: if there are more than one volume view, the
;   ranges of their data are not too varied.
;
    oVolViews[i]->GetProperty, volume=oVolume
    oVolume->GetProperty, data0=data0, /no_copy
    vol_max = max(data0, min=vol_min, /nan)
    if i eq 0 then begin
        v_max = vol_max
        v_min = vol_min
        end $
    else begin
        v_max = v_max > vol_max
        v_min = v_min < vol_min
        end
    oVolume->SetProperty, data0=data0, /no_copy
    endfor

self.data_range = [v_min, v_max]
self.oViewgroup->CallAll, 'SetProperty', $
    isosurface_threshold=((v_max - v_min) / 255.) * 128 + v_min, $
    isa='IDLexVolView'

return, 1 ; success.
End

;--------------------------------------------------------------------
pro IDLexVolViewWid__define

compile_opt hidden
on_error, 2

struct_hide, {IDLexVolViewWid, $
    inherits IDLexObjViewWID, $
    wXSlider: 0L, $
    wYSlider: 0L, $
    wZSlider: 0L, $
    wXContourCheckbox: 0L, $
    wYContourCheckbox: 0L, $
    wZContourCheckbox: 0L, $
    wIsoLevel: 0L, $
    wIsoToggle0: 0L, $
    wIsoToggle1: 0L, $
    wIsoToggle2: 0L, $
    wControlPanel: 0L, $
    wRender: 0L, $
    wXDroplist: 0L, $
    wYDroplist: 0L, $
    wZDroplist: 0L, $
    wAutoRender: 0L, $
    wPaletteEditorButton: 0L, $
    wIsoColorButton: 0L, $
    wPalEdTlb: 0L, $
    wRGBSliderTlb: 0L, $
    wXAxisToggle: 0L, $
    wYAxisToggle: 0L, $
    wZAxisToggle: 0L, $
    hide_x_axis: 0b, $
    hide_y_axis: 0b, $
    hide_z_axis: 0b, $
    vol_on_screen: 0b, $
    hide_volume: 0b, $
    resize_when_realized: 0b, $
    data_range: dblarr(2) $
    }

end
