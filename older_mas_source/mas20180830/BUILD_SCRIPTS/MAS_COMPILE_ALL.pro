;; $Id$
;;
;; Note: commented out .compile lines indicate that the file exists
;;       but is not being used in MAS at the moment.

;;.RESET

;; make_dll, 'nnls', 'nnls_c', input_dir=dir, output_dir=dir, extra_cflags='-fast', /verbose, /show_all
;; where dir is the path to the nnls.c file.
;
;; Integrate "Changes.txt" into MAS
;; No longer necessary... 
;;mas_make_changelog

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; DEPENDENCIES FROM IDL stock Library                              ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

.compile cw_fslider.pro
.compile cw_animate.pro
.compile xmanager.pro

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Third-Party Dependencies                                         ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Fanning Consulting
;; http://www.dfanning.com/documents/programs.html#FSC_COLOR
.compile FSC_Color_file.pro
;; http://www.dfanning.com/documents/programs.html#PROGRESSBAR__DEFINE
.compile progressbar__define.pro
;; http://www.idlcoyote.com/programs/transform_volume.pro
.compile Transform_Volume.pro

;; Craig Markwardt's MPFIT
;; http://cow.physics.wisc.edu/~craigm/idl/fitting.html
.compile mpfit.pro               
.compile mpfitfun.pro

;; http://www.exelisvis.com/Support/HelpArticles/TabId/185/ArtMID/800/ArticleID/13847/Example-of-how-to-generate-an-ellipsoid-with-IDL-8-Graphics.aspx
.compile dj_ellipsoid_ng.pro

;; http://www.capca.ucalgary.ca/~wdobler/doc/idl/linspace.pro
;;.compile mas_linspace.pro

;; gdlffdicom, a replacement for IDL's dlm idlffdicom.
.compile gdlffdicom__assoc__define.pro
.compile gdlffdicom__assoc__test.pro
.compile gdlffdicom__assoc_generateuid.pro
.compile gdlffdicom__define.pro
; -> do not compile ;;.compile gdlffdicom__dictionary.pro
.compile gdlffdicom__shared__define.pro
.compile gdlffdicom__shared__test.pro
.compile gdlffdicom__test.pro
.compile gdlffdicom_copy_lun.pro
.compile gdlffdicom_date.pro
.compile gdlffdicom_time.pro
.compile gdlffdicom_trim.pro

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; These have been "appropriated" from the IDL library and altered  ;;
;; by someone at one point to suit our needs.                       ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
.compile idlexinscribingview__define.pro
.compile idlexmodelmanip__define.pro
.compile idlexobjview__define.pro
.compile idlexobjviewwid__define.pro
.compile idlexviewgroup__define.pro
.compile idlexviewmanip__define.pro
.compile idlexvolview__define.pro
.compile idlexvolviewwid__define.pro
.compile idlexwidget__define.pro
.compile xobjview.pro
.compile xobjview_rotate.pro
.compile xobjview_write_image.pro
.compile xroi.pro
.compile show_stream.pro

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; MAS-specific procedure files                                     ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
.compile Build_mas.pro

.compile mas_3d_slicer.pro
.compile mas_3d_slicer_GUI.pro
.compile mas_surface_render.pro
.compile mas_3d_surface.pro
.compile mas_3d_surface_GUI.pro
.compile mas.pro
.compile mas_callback.pro
.compile mas_adt_regress.pro
.compile mas_roi.pro
.compile mas_orient_data.pro
.compile mas_display.pro
.compile mas_line_data_window.pro
.compile mas_display_multi.pro
.compile mas_display_ortho.pro
.compile mas_display_single.pro
.compile REGION_SHRINK.pro
.compile FibersGUI_eventcb.pro
.compile ft_add_tracts.pro
.compile ft_display_adt.pro
;;.compile ft_journal.pro
.compile ft_roi.pro
.compile ft_roi_calculate.pro
.compile ft_track.pro
.compile ft_visualize.pro
.compile ft_xroi.pro
.compile func_rd_slicer.pro
.compile mas_fibers.pro
.compile view_fibers.pro
.compile view_fibers_GUI.pro
.compile view_fibers_hardi.pro
.compile view_fibers_hardi_GUI.pro
.compile mas_complex_plot__define.pro
.compile mas_display_1d_spect.pro
.compile mas_algebra.pro
.compile mas_background_filter_file.pro
.compile mas_charmed.pro
.compile mas_curvefit.pro
.compile mas_select_data.pro
.compile mas_curvefit2_core.pro
.compile mas_curvefit2_gui.pro
.compile mas_curvefit2_models.pro
.compile mas_read_philips_parrec.pro
.compile mas_DICOM.pro
.compile mas_dicom_writer.pro
.compile mas_dicom_reader.pro
.compile mas_zeropad.pro
.compile mas_phase_unwrap_fsl.pro
.compile mas_phase_unwrap2d.pro
.compile mas_phase_unwrap3d.pro
.compile mas_phase_unwrap.pro
.compile mas_flow.pro
.compile mas_Bez.pro
.compile mas_mip.pro
.compile mas_glyphscene_renderer__define.pro
.compile mas_diffusion_tools.pro
.compile mas_adt_superquadric.pro
.compile mas_odf__define.pro
.compile mas_odf_tools.pro
.compile mas_dot__define.pro
.compile mas_dot_tools.pro
.compile mas_nnls.pro
.compile mas_mow__define.pro
.compile mas_mow_tools.pro
.compile mas_exit.pro
.compile mas_filtering.pro
.compile mas_dce_processing.pro
.compile mas_hnb_ellipsoid.pro
.compile mas_hardi_tractography__define.pro
.compile mas_hardi_tractography_gui.pro
.compile mas_hardi_tractography.pro
.compile mas_image_filters.pro
.compile mas_image_statistics.pro
.compile mas_load_state_1.pro
.compile mas_load_state_2.pro
.compile mas_make_changelog.pro
.compile mas_make_sheet_display.pro
.compile mas_movie.pro
.compile mas_ucb_examples.pro
.compile mas_varian_procpar__define.pro
.compile mas_varian.pro
.compile mas_varian_glue.pro
.compile mas_extract_files.pro
.compile mas_extract_imnd.pro
.compile mas_nifti.pro
.compile mas_open.pro
.compile mas_probtrack.pro
.compile mas_read_philips_dicom_spect.pro
.compile mas_read_philips_sparsdat_spect.pro
.compile mas_redraw.pro
.compile mas_remove_scan.pro
.compile mas_rotate_flip.pro
.compile mas_save.pro
.compile mas_sig_enhance.pro
.compile mas_smooth.pro
.compile mas_subset.pro
.compile mas_tessellator.pro
.compile mas_timeseries_plot_window.pro
.compile mas_transform_interpolate.pro
.compile mas_volume.pro
.compile mas_windowing.pro
.compile mas_calc_delta_b_algebra_gui.pro
.compile mas_calc_delta_b_least_squares_gui.pro
;;.compile mas_write_files.pro
.compile mas_zoom.pro
.compile motion_correct_mi_gui.pro
.compile motion_correct_mi_poses.pro
.compile motion_correct.pro
.compile motion_correct_ed.pro
.compile motion_correct_mi.pro
.compile motion_correct_mi_vol.pro
.compile phillips_raw8echo.pro
.compile rd_slicer.pro
.compile roi3d_slicer.pro
.compile roi_dc.pro
.compile slicer3m.pro
.compile write_slicer.pro
RESOLVE_ALL
IRESOLVE

