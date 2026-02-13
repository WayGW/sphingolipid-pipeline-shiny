# -*- coding: utf-8 -*-
"""
Sphingolipid Analysis Pipeline - Streamlit Application
=======================================================

Run with: streamlit run app.py
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO
import base64
from pathlib import Path
import sys
import tempfile
import zipfile
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent))

from config.sphingolipid_species import (
    SPHINGOLIPID_PANEL, get_ceramides, get_dihydroceramides,
    get_sphingomyelins, get_sphingoid_bases, get_sphingoid_base_phosphates,
    get_saturated, get_unsaturated, get_very_long_chain, get_long_chain
)
from modules.data_processing import SphingolipidDataProcessor, ProcessedData, validate_data_quality
from modules.statistical_tests import StatisticalAnalyzer, format_analysis_report
from modules.visualization import SphingolipidVisualizer, create_summary_figure
from modules.report_generation import (
    ExcelReportGenerator, SignificancePlotter, 
    ComprehensiveAnalysisResults, format_apa_statistics,
    get_significant_differences_summary
)

st.set_page_config(page_title="Sphingolipid Analysis Pipeline", page_icon="🧬", layout="wide")


def create_excel_with_lod_highlighting(
    df: pd.DataFrame, 
    lod_rows: dict, 
    sheet_name: str = "Data"
) -> BytesIO:
    """
    Create an Excel file with cells highlighted where LOD values were replaced.
    
    Args:
        df: DataFrame to export
        lod_rows: Dict mapping column names to list of row indices with LOD replacements
        sheet_name: Name for the Excel sheet
        
    Returns:
        BytesIO buffer containing the Excel file
    """
    buf = BytesIO()
    
    with pd.ExcelWriter(buf, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name=sheet_name, index=False)
        
        # Get the worksheet
        ws = writer.sheets[sheet_name]
        
        # Import openpyxl styles for highlighting
        from openpyxl.styles import PatternFill, Font, Alignment
        
        # Yellow fill for LOD-replaced cells
        lod_fill = PatternFill(start_color="FFFF99", end_color="FFFF99", fill_type="solid")
        
        # Get column positions (1-indexed for openpyxl, +1 because row 1 is header)
        col_positions = {col: idx + 1 for idx, col in enumerate(df.columns)}
        
        # Highlight LOD cells
        for col_name, row_indices in lod_rows.items():
            if col_name in col_positions:
                col_idx = col_positions[col_name]
                for row_idx in row_indices:
                    # Excel rows are 1-indexed, and row 1 is the header
                    excel_row = row_idx + 2  # +1 for 0-indexing, +1 for header
                    cell = ws.cell(row=excel_row, column=col_idx)
                    cell.fill = lod_fill
        
        # Add a legend in a new sheet
        legend_ws = writer.book.create_sheet("Legend")
        legend_ws['A1'] = "Cell Highlighting Legend"
        legend_ws['A1'].font = Font(bold=True, size=12)
        legend_ws['A3'] = "Yellow cells"
        legend_ws['B3'] = "Below LOD - value was replaced based on LOD handling setting"
        legend_ws['A3'].fill = lod_fill
        
        # Auto-adjust column widths for main sheet
        for column in ws.columns:
            max_length = 0
            column_letter = column[0].column_letter
            for cell in column:
                try:
                    if len(str(cell.value)) > max_length:
                        max_length = len(str(cell.value))
                except:
                    pass
            adjusted_width = min(max_length + 2, 50)
            ws.column_dimensions[column_letter].width = adjusted_width
    
    buf.seek(0)
    return buf


def create_full_data_excel_with_highlighting(processed, metadata) -> BytesIO:
    """
    Create a comprehensive Excel file with all data sheets and LOD highlighting.
    """
    buf = BytesIO()
    lod_rows = processed.structure.analyte_lod_rows
    
    with pd.ExcelWriter(buf, engine='openpyxl') as writer:
        from openpyxl.styles import PatternFill, Font
        
        lod_fill = PatternFill(start_color="FFFF99", end_color="FFFF99", fill_type="solid")
        
        # Sheet 1: Concentrations (with highlighting)
        conc_df = pd.concat([metadata, processed.concentrations], axis=1)
        conc_df.to_excel(writer, sheet_name="Concentrations", index=False)
        ws_conc = writer.sheets["Concentrations"]
        
        # Get column positions for concentrations
        conc_col_positions = {col: idx + 1 for idx, col in enumerate(conc_df.columns)}
        
        # Highlight LOD cells in concentrations
        for col_name, row_indices in lod_rows.items():
            if col_name in conc_col_positions:
                col_idx = conc_col_positions[col_name]
                for row_idx in row_indices:
                    excel_row = row_idx + 2
                    cell = ws_conc.cell(row=excel_row, column=col_idx)
                    cell.fill = lod_fill
        
        # Sheet 2: Totals
        totals_df = pd.concat([metadata, processed.totals], axis=1)
        totals_df.to_excel(writer, sheet_name="Totals", index=False)
        
        # Sheet 3: Ratios
        ratios_df = pd.concat([metadata, processed.ratios], axis=1)
        ratios_df.to_excel(writer, sheet_name="Ratios", index=False)
        
        # Sheet 4: Percentages
        pct_df = pd.concat([metadata, processed.percentages], axis=1)
        pct_df.to_excel(writer, sheet_name="Percentages", index=False)
        
        # Sheet 5: LOD Summary
        n_samples = len(processed.sample_data)
        lod_summary = []
        for analyte, lod in sorted(processed.structure.analyte_lods.items()):
            count = processed.structure.analyte_lod_counts.get(analyte, 0)
            pct = (count / n_samples * 100) if n_samples > 0 else 0
            lod_summary.append({
                "Analyte": analyte,
                "LOD (ng/mL)": lod,
                "Below LOD Count": count,
                "% Replaced": round(pct, 1)
            })
        lod_summary_df = pd.DataFrame(lod_summary)
        lod_summary_df.to_excel(writer, sheet_name="LOD Summary", index=False)
        
        # Sheet 6: Legend
        legend_ws = writer.book.create_sheet("Legend")
        legend_ws['A1'] = "Cell Highlighting Legend"
        legend_ws['A1'].font = Font(bold=True, size=12)
        legend_ws['A3'] = "Yellow cells"
        legend_ws['B3'] = "Below LOD - value was replaced based on LOD handling setting"
        legend_ws['A3'].fill = lod_fill
        legend_ws['A5'] = "LOD Source"
        legend_ws['B5'] = processed.structure.lod_source.capitalize()
        legend_ws['A6'] = "Total Samples"
        legend_ws['B6'] = n_samples
        
        # Auto-adjust column widths for all sheets
        for sheet_name in writer.sheets:
            ws = writer.sheets[sheet_name]
            for column in ws.columns:
                max_length = 0
                column_letter = column[0].column_letter
                for cell in column:
                    try:
                        if len(str(cell.value)) > max_length:
                            max_length = len(str(cell.value))
                    except:
                        pass
                adjusted_width = min(max_length + 2, 50)
                ws.column_dimensions[column_letter].width = adjusted_width
    
    buf.seek(0)
    return buf


def init_session_state():
    """Initialize session state variables."""
    defaults = {
        'processed_data': None,
        'analysis_results': ComprehensiveAnalysisResults(),
        'figures': {},
        'report_generator': None,
        'stats_computed': False,
        'last_file': None,
        'last_settings': None  # Track settings that affect data/stats
    }
    for key, val in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = val


def get_data_affecting_settings(settings):
    """Extract settings that affect data processing and statistics."""
    return {
        'lod_handling': settings['lod_handling'],
        'lod_value': settings['lod_value'],
        'alpha': settings['alpha']
    }


def check_settings_changed(settings):
    """Check if data-affecting settings have changed, reset caches if so."""
    current = get_data_affecting_settings(settings)
    
    if st.session_state.last_settings != current:
        # Settings changed - need to reprocess data and recompute stats
        st.session_state.stats_computed = False
        st.session_state.figures = {}
        st.session_state.last_settings = current
        return True
    return False


def fig_to_bytes(fig, format='png', dpi=300):
    """Convert matplotlib figure to bytes."""
    buf = BytesIO()
    fig.savefig(buf, format=format, dpi=dpi, bbox_inches='tight')
    buf.seek(0)
    return buf.getvalue()


def store_figure(fig, name):
    """Store a figure in session state."""
    st.session_state.figures[name] = fig


def compute_all_statistics(processed, settings):
    """Compute all statistics once and cache."""
    if st.session_state.stats_computed:
        return st.session_state.analysis_results
    
    group_col = processed.structure.group_col
    if not group_col:
        return None
    
    # Detect factor info for two-way ANOVA
    factors = getattr(processed.structure, 'factors', {})
    n_factors = getattr(processed.structure, 'n_factors', 0)
    
    report_gen = ExcelReportGenerator(
        data=processed.sample_data,
        group_col=group_col,
        sphingolipid_cols=processed.structure.sphingolipid_cols,
        totals=processed.totals,
        ratios=processed.ratios,
        percentages=processed.percentages,
        alpha=settings['alpha'],
        factors=factors,
        n_factors=n_factors,
    )
    
    results = report_gen.run_all_statistics()
    st.session_state.analysis_results = results
    st.session_state.report_generator = report_gen
    st.session_state.stats_computed = True
    return results


def get_sig_pairs(result, max_pairs=5):
    """Extract significant pairs from result."""
    if not result or not result.posthoc_test or result.posthoc_test.pairwise_results is None:
        return []
    pairs = []
    df = result.posthoc_test.pairwise_results
    sig = df[df['significant']].sort_values('pvalue_adj' if 'pvalue_adj' in df else 'pvalue')
    for _, row in sig.head(max_pairs).iterrows():
        p = row.get('pvalue_adj', row.get('pvalue', 1.0))
        annot = '***' if p < 0.001 else ('**' if p < 0.01 else '*' if p < 0.05 else '')
        if annot:
            pairs.append((row['group1'], row['group2'], annot))
    return pairs


def generate_all_export_figures(processed, results, settings):
    """Generate all figure variations for export package."""
    figures = {}
    viz = SphingolipidVisualizer(color_palette=settings['color_palette'], style=settings['plot_style'])
    group_col = processed.structure.group_col
    
    if not group_col or not results:
        return figures
    
    data = processed.sample_data[processed.sample_data[group_col].notna()].copy()
    data = data[data[group_col].astype(str).str.lower() != 'nan']
    available_sls = processed.structure.sphingolipid_cols
    
    is_twoway = results.is_twoway
    
    if is_twoway:
        # =====================================================================
        # TWO-WAY ANOVA EXPORT FIGURES
        # =====================================================================
        fa_col = results.factor_a_col
        fb_col = results.factor_b_col
        fa_name = results.factor_a_name
        fb_name = results.factor_b_name
        
        # --- Concentrations ---
        conc_selections = {
            'top10': processed.concentrations.mean().nlargest(10).index.tolist(),
            'significant': [s for s in available_sls if s in results.twoway_individual_sl
                           and (results.twoway_individual_sl[s].twoway_result.factor_a_pvalue < settings['alpha']
                                or results.twoway_individual_sl[s].twoway_result.factor_b_pvalue < settings['alpha']
                                or results.twoway_individual_sl[s].twoway_result.interaction_pvalue < settings['alpha'])][:10],
            'ceramides': [s for s in get_ceramides() if s in available_sls][:10],
            'sphingomyelins': [s for s in get_sphingomyelins() if s in available_sls][:10],
        }
        
        for sel_name, selected in conc_selections.items():
            if not selected:
                continue
            tw_stats = {s: results.twoway_individual_sl.get(s) for s in selected 
                       if s in results.twoway_individual_sl}
            suffix = f"_{sel_name}"
            try:
                fig = viz.plot_twoway_multi_panel(
                    data=data, value_cols=selected,
                    factor_a_col=fa_col, factor_b_col=fb_col,
                    twoway_results=tw_stats, ncols=3,
                    factor_a_name=fa_name, factor_b_name=fb_name,
                    ylabel='Concentration (ng/mL)', show_points=settings['show_points'])
                figures[f'concentrations{suffix}'] = fig
            except Exception as _e:
                print(f'Export figure error: {_e}')
            try:
                fig_int = viz.plot_twoway_interaction_multi_panel(
                    data=data, value_cols=selected,
                    factor_a_col=fa_col, factor_b_col=fb_col,
                    twoway_results=tw_stats, ncols=3,
                    factor_a_name=fa_name, factor_b_name=fb_name,
                    ylabel='Concentration (ng/mL)')
                figures[f'concentrations{suffix}_interaction'] = fig_int
            except Exception as _e:
                print(f'Export figure error: {_e}')
        
        # --- Totals ---
        key_totals = ['total_all', 'total_ceramides', 'total_dihydroceramides', 'total_sphingomyelins', 
                      'sphingoid_bases', 'sphingoid_base_phosphates', 'very_long_chain', 'long_chain']
        available_totals = [t for t in key_totals if t in processed.totals.columns]
        
        if available_totals:
            totals_data = pd.concat([data[[fa_col, fb_col]].reset_index(drop=True),
                                     processed.totals[available_totals].loc[data.index].reset_index(drop=True)], axis=1)
            tw_stats = {t: results.twoway_totals.get(t) for t in available_totals 
                       if t in results.twoway_totals}
            try:
                fig = viz.plot_twoway_multi_panel(
                    data=totals_data, value_cols=available_totals,
                    factor_a_col=fa_col, factor_b_col=fb_col,
                    twoway_results=tw_stats, ncols=3,
                    factor_a_name=fa_name, factor_b_name=fb_name,
                    ylabel='Concentration (ng/mL)', show_points=settings['show_points'])
                figures['totals'] = fig
            except Exception as _e:
                print(f'Export figure error: {_e}')
            try:
                fig_int = viz.plot_twoway_interaction_multi_panel(
                    data=totals_data, value_cols=available_totals,
                    factor_a_col=fa_col, factor_b_col=fb_col,
                    twoway_results=tw_stats, ncols=3,
                    factor_a_name=fa_name, factor_b_name=fb_name,
                    ylabel='Concentration (ng/mL)')
                figures['totals_interaction'] = fig_int
            except Exception as _e:
                print(f'Export figure error: {_e}')
        
        # --- Percentages ---
        pct_cols = [col for col in processed.percentages.columns if col.endswith('_pct')]
        if pct_cols:
            top_pct = processed.percentages[pct_cols].mean().nlargest(10).index.tolist()
            display_map = {col: col.replace('_pct', '') for col in top_pct}
            display_cols = list(display_map.values())
            
            pct_data = pd.concat([data[[fa_col, fb_col]].reset_index(drop=True),
                                  processed.percentages[top_pct].loc[data.index].reset_index(drop=True)], axis=1)
            pct_data = pct_data.rename(columns=display_map)
            
            tw_stats = {}
            for pct_col in top_pct:
                disp = pct_col.replace('_pct', '')
                if pct_col in results.twoway_percentages:
                    tw_stats[disp] = results.twoway_percentages[pct_col]
            
            try:
                fig = viz.plot_twoway_multi_panel(
                    data=pct_data, value_cols=display_cols,
                    factor_a_col=fa_col, factor_b_col=fb_col,
                    twoway_results=tw_stats, ncols=3,
                    factor_a_name=fa_name, factor_b_name=fb_name,
                    ylabel='% of Total SL', show_points=settings['show_points'])
                figures['percentages_top10'] = fig
            except Exception as _e:
                print(f'Export figure error: {_e}')
            try:
                fig_int = viz.plot_twoway_interaction_multi_panel(
                    data=pct_data, value_cols=display_cols,
                    factor_a_col=fa_col, factor_b_col=fb_col,
                    twoway_results=tw_stats, ncols=3,
                    factor_a_name=fa_name, factor_b_name=fb_name,
                    ylabel='% of Total SL')
                figures['percentages_interaction'] = fig_int
            except Exception as _e:
                print(f'Export figure error: {_e}')
        
        # --- Ratios ---
        ratio_cols = [col for col in processed.ratios.columns if not processed.ratios[col].isna().all()]
        if ratio_cols:
            selected_ratios = ratio_cols[:9]
            ratio_data = pd.concat([data[[fa_col, fb_col]].reset_index(drop=True),
                                    processed.ratios[selected_ratios].loc[data.index].reset_index(drop=True)], axis=1)
            tw_stats = {r: results.twoway_ratios.get(r) for r in selected_ratios 
                       if r in results.twoway_ratios}
            
            try:
                fig = viz.plot_twoway_multi_panel(
                    data=ratio_data, value_cols=selected_ratios,
                    factor_a_col=fa_col, factor_b_col=fb_col,
                    twoway_results=tw_stats, ncols=3,
                    factor_a_name=fa_name, factor_b_name=fb_name,
                    ylabel='Ratio', show_points=settings['show_points'])
                figures['ratios'] = fig
            except Exception as _e:
                print(f'Export figure error: {_e}')
            try:
                fig_int = viz.plot_twoway_interaction_multi_panel(
                    data=ratio_data, value_cols=selected_ratios,
                    factor_a_col=fa_col, factor_b_col=fb_col,
                    twoway_results=tw_stats, ncols=3,
                    factor_a_name=fa_name, factor_b_name=fb_name,
                    ylabel='Ratio')
                figures['ratios_interaction'] = fig_int
            except Exception as _e:
                print(f'Export figure error: {_e}')
    
    else:
        # =====================================================================
        # ONE-WAY ANOVA EXPORT FIGURES (unchanged)
        # =====================================================================
        
        # === CONCENTRATIONS TAB FIGURES ===
        conc_selections = {
            'top10': processed.concentrations.mean().nlargest(10).index.tolist(),
            'significant': [s for s in available_sls if s in results.individual_sl_results 
                           and results.individual_sl_results[s].main_test.significant][:10],
            'ceramides': [s for s in get_ceramides() if s in available_sls][:10],
            'sphingomyelins': [s for s in get_sphingomyelins() if s in available_sls][:10],
        }
        
        for sel_name, selected in conc_selections.items():
            if not selected:
                continue
            stats_dict = {s: results.individual_sl_results.get(s) for s in selected if s in results.individual_sl_results}
            
            for log_scale in [False, True]:
                suffix = f"_{sel_name}{'_log' if log_scale else ''}"
                try:
                    fig = viz.plot_multi_panel_groups_with_stats(
                        data, selected, group_col, stats_dict, ncols=3,
                        plot_type=settings['plot_type'], log_scale=log_scale,
                        show_points=settings['show_points'])
                    figures[f'concentrations{suffix}'] = fig
                except Exception as _e:
                    print(f'Export figure error: {_e}')
        
        # === TOTALS TAB FIGURES ===
        totals_combined = pd.concat([data[[group_col]], processed.totals.loc[data.index]], axis=1)
        
        key_totals = ['total_all', 'total_ceramides', 'total_dihydroceramides', 'total_sphingomyelins', 
                      'sphingoid_bases', 'sphingoid_base_phosphates', 'very_long_chain', 'long_chain']
        available_totals = [t for t in key_totals if t in totals_combined.columns]
        
        if available_totals:
            totals_stats = {t: results.totals_results.get(t) for t in available_totals if t in results.totals_results}
            for log_scale in [False, True]:
                suffix = '_log' if log_scale else ''
                try:
                    fig = viz.plot_multi_panel_groups_with_stats(
                        totals_combined, available_totals, group_col, totals_stats, ncols=3,
                        plot_type=settings['plot_type'], log_scale=log_scale,
                        show_points=settings['show_points'])
                    figures[f'totals{suffix}'] = fig
                except Exception as _e:
                    print(f'Export figure error: {_e}')
        
        # === PERCENTAGES TAB FIGURES ===
        pct_cols = [col for col in processed.percentages.columns if col.endswith('_pct')]
        if pct_cols:
            top_pct = processed.percentages[pct_cols].mean().nlargest(10).index.tolist()
            pct_combined = pd.concat([data[[group_col]], processed.percentages.loc[data.index, top_pct]], axis=1)
            
            display_cols = [col.replace('_pct', '') for col in top_pct]
            pct_display = pct_combined.rename(columns={old: new for old, new in zip(top_pct, display_cols)})
            
            pct_stats = {}
            for pct_col, disp_col in zip(top_pct, display_cols):
                if pct_col in results.percentages_results:
                    pct_stats[disp_col] = results.percentages_results[pct_col]
            
            try:
                fig = viz.plot_multi_panel_groups_with_stats(
                    pct_display, display_cols, group_col, pct_stats, ncols=3,
                    plot_type=settings['plot_type'], log_scale=False,
                    show_points=settings['show_points'])
                figures['percentages_top10'] = fig
            except Exception as _e:
                print(f'Export figure error: {_e}')
        
        # === RATIOS TAB FIGURES ===
        ratio_cols = [col for col in processed.ratios.columns if not processed.ratios[col].isna().all()]
        if ratio_cols:
            ratios_combined = pd.concat([data[[group_col]], processed.ratios.loc[data.index, ratio_cols]], axis=1)
            ratios_stats = {r: results.ratios_results.get(r) for r in ratio_cols if r in results.ratios_results}
            
            for log_scale in [False, True]:
                suffix = '_log' if log_scale else ''
                try:
                    fig = viz.plot_multi_panel_groups_with_stats(
                        ratios_combined, ratio_cols[:9], group_col, ratios_stats, ncols=3,
                        plot_type=settings['plot_type'], log_scale=log_scale,
                        show_points=settings['show_points'])
                    figures[f'ratios{suffix}'] = fig
                except Exception as _e:
                    print(f'Export figure error: {_e}')
    
    return figures


def create_results_zip(processed, results, figures, report_gen, metadata):
    """Create ZIP with all results including LOD-highlighted Excel files."""
    buf = BytesIO()
    with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as zf:
        # Full data Excel with LOD highlighting (primary data file)
        full_excel_buf = create_full_data_excel_with_highlighting(processed, metadata)
        zf.writestr('data/sphingolipid_all_data.xlsx', full_excel_buf.getvalue())
        
        # Also include CSVs for compatibility
        for name, df in [('concentrations', processed.concentrations), 
                         ('percentages', processed.percentages),
                         ('totals', processed.totals), 
                         ('ratios', processed.ratios)]:
            csv_buf = BytesIO()
            export_df = pd.concat([metadata, df], axis=1)
            export_df.to_csv(csv_buf, index=False)
            zf.writestr(f'data/{name}.csv', csv_buf.getvalue())
        
        # Full data CSV
        full = pd.concat([processed.sample_data, processed.totals, processed.ratios, processed.percentages], axis=1)
        csv_buf = BytesIO()
        full.to_csv(csv_buf, index=False)
        zf.writestr('data/full_analysis.csv', csv_buf.getvalue())
        
        # Excel report
        if report_gen:
            excel_buf = BytesIO()
            report_gen.save_excel_report(excel_buf)
            excel_buf.seek(0)
            zf.writestr('reports/statistical_report.xlsx', excel_buf.getvalue())
        
        # Figures
        for name, fig in figures.items():
            if fig:
                zf.writestr(f'figures/{name}.png', fig_to_bytes(fig, 'png'))
                zf.writestr(f'figures/{name}.pdf', fig_to_bytes(fig, 'pdf'))
        
        # README
        readme_content = f"""Sphingolipid Analysis Results
==============================

Generated: {datetime.now():%Y-%m-%d %H:%M}

Contents:
---------
data/
  - sphingolipid_all_data.xlsx  (⭐ Primary data file with LOD highlighting)
  - concentrations.csv
  - totals.csv
  - ratios.csv
  - percentages.csv
  - full_analysis.csv

reports/
  - statistical_report.xlsx

figures/
  - Various PNG and PDF figures

LOD Highlighting:
-----------------
In sphingolipid_all_data.xlsx, yellow cells indicate values that were
below the limit of detection (LOD) and were replaced according to the
selected LOD handling method.

LOD Source: {processed.structure.lod_source}
Total Samples: {len(processed.sample_data)}
Total Analytes: {len(processed.structure.sphingolipid_cols)}
"""
        zf.writestr('README.txt', readme_content)
    
    buf.seek(0)
    return buf.getvalue()


def render_sidebar():
    """Render sidebar settings."""
    st.sidebar.markdown("## ⚙️ Settings")
    
    # LOD Settings
    st.sidebar.markdown("### 📊 Detection Limits")
    st.sidebar.caption("LODs are auto-detected per-analyte from standard curves")
    
    default_lod = 1.0  # Fallback LOD in ng/mL
    lod_value = st.sidebar.number_input("Fallback LOD (ng/mL)", 0.0, 100.0, default_lod, step= 0.5,
                                        help="Used only when auto-detection fails for an analyte")
    
    lod_handling = st.sidebar.selectbox("Below LOD handling", [ "lod","half_lod", "zero", "half_min", "drop"],
                                        format_func=lambda x: {"lod": "LOD value", "half_lod": "LOD/2", 
                                                               "zero": "Zero", "half_min": "Min/2", "drop": "NaN"}[x],
                                        help="How to replace below-LOD values (uses per-analyte LOD)")
    
    st.sidebar.markdown("### 📈 Analysis")
    alpha = st.sidebar.slider("Significance (α)", 0.01, 0.10, 0.05, 0.01)
    plot_type = st.sidebar.selectbox("Plot type", ["box", "violin", "bar"])
    show_points = st.sidebar.checkbox("Show data points", True)
    
    # Color palette selection
    st.sidebar.markdown("### 🎨 Appearance")
    palette_options = {
        "Set2": "Set2 (Default - Soft pastels)",
        "Set1": "Set1 (Bold primary)",
        "Paired": "Paired (Light/dark pairs)",
        "Dark2": "Dark2 (Darker pastels)",
        "colorblind": "Colorblind-friendly",
        "tab10": "Tab10 (Matplotlib default)",
        "Pastel1": "Pastel1 (Very soft)",
        "Accent": "Accent (High contrast)"
    }
    color_palette = st.sidebar.selectbox("Color palette", list(palette_options.keys()),
                                         format_func=lambda x: palette_options[x])
    
    # Background style selection
    style_options = {
        "whitegrid": "White with grid",
        "white": "Clean white",
        "darkgrid": "Gray with grid",
        "ticks": "White with ticks"
    }
    plot_style = st.sidebar.selectbox("Plot background", list(style_options.keys()),
                                      format_func=lambda x: style_options[x])
    
    return {'lod_handling': lod_handling, 'alpha': alpha, 'plot_type': plot_type,
            'show_points': show_points, 'lod_value': lod_value,
            'color_palette': color_palette, 'plot_style': plot_style}


def render_concentrations_tab(processed, settings):
    """Render concentrations tab."""
    st.markdown("### Individual Sphingolipid Concentrations")
    
    viz = SphingolipidVisualizer(color_palette=settings['color_palette'], style=settings['plot_style'])
    group_col = processed.structure.group_col
    results = st.session_state.analysis_results
    
    if not group_col:
        st.warning("No group column detected.")
        return
    
    data = processed.sample_data[processed.sample_data[group_col].notna()].copy()
    data = data[data[group_col].astype(str).str.lower() != 'nan']
    available_sls = processed.structure.sphingolipid_cols
    
    # Check if two-way ANOVA
    is_twoway = results is not None and results.is_twoway
    
    col1, col2 = st.columns([1, 2])
    with col1:
        quick = st.selectbox("Quick select", ["Top 10", "Significant", "Ceramides", "Sphingomyelins", "Custom"])
    
    if quick == "Top 10":
        selected = processed.concentrations.mean().nlargest(10).index.tolist()
    elif quick == "Significant":
        if is_twoway:
            selected = [s for s in available_sls if s in results.twoway_individual_sl
                       and (results.twoway_individual_sl[s].twoway_result.factor_a_pvalue < settings['alpha']
                            or results.twoway_individual_sl[s].twoway_result.factor_b_pvalue < settings['alpha']
                            or results.twoway_individual_sl[s].twoway_result.interaction_pvalue < settings['alpha'])][:10]
        else:
            selected = [s for s in available_sls if s in results.individual_sl_results 
                       and results.individual_sl_results[s].main_test.significant][:10]
        if not selected:
            st.info("No significant individual sphingolipids found.")
            selected = processed.concentrations.mean().nlargest(5).index.tolist()
    elif quick == "Ceramides":
        selected = [s for s in get_ceramides() if s in available_sls][:10]
    elif quick == "Sphingomyelins":
        selected = [s for s in get_sphingomyelins() if s in available_sls][:10]
    else:
        with col2:
            selected = st.multiselect("Select sphingolipids", available_sls, available_sls[:5])
    
    log_scale = st.checkbox("Log10 scale", key="conc_log")
    
    if selected:
        if is_twoway:
            # === TWO-WAY ANOVA PLOTS ===
            fa_col = results.factor_a_col
            fb_col = results.factor_b_col
            fa_name = results.factor_a_name
            fb_name = results.factor_b_name
            
            st.markdown(f"**Two-way ANOVA**: {fa_name} x {fb_name}")
            
            tw_stats = {s: results.twoway_individual_sl.get(s) for s in selected 
                       if s in results.twoway_individual_sl}
            
            fig = viz.plot_twoway_multi_panel(
                data=data, value_cols=selected,
                factor_a_col=fa_col, factor_b_col=fb_col,
                twoway_results=tw_stats, ncols=3,
                factor_a_name=fa_name, factor_b_name=fb_name,
                ylabel='Concentration (ng/mL)', show_points=settings['show_points']
            )
            st.pyplot(fig)
            store_figure(fig, 'concentrations_twoway')
            plt.close(fig)
            
            with st.expander("Interaction Plots"):
                fig_int = viz.plot_twoway_interaction_multi_panel(
                    data=data, value_cols=selected,
                    factor_a_col=fa_col, factor_b_col=fb_col,
                    twoway_results=tw_stats, ncols=3,
                    factor_a_name=fa_name, factor_b_name=fb_name,
                    ylabel='Concentration (ng/mL)'
                )
                st.pyplot(fig_int)
                store_figure(fig_int, 'concentrations_interaction')
                plt.close(fig_int)
            
            with st.expander("Two-Way ANOVA Summary"):
                from modules.report_generation import get_twoway_differences_summary
                st.dataframe(get_twoway_differences_summary(tw_stats), hide_index=True)
        else:
            # === ONE-WAY ANOVA PLOTS (unchanged) ===
            stats_dict = {s: results.individual_sl_results.get(s) for s in selected if s in results.individual_sl_results}
            
            fig = viz.plot_multi_panel_groups_with_stats(data, selected, group_col, stats_dict, 
                                                         ncols=3, plot_type=settings['plot_type'], log_scale=log_scale,
                                                         show_points=settings['show_points'], ylabel='Concentration (ng/mL)')
            st.pyplot(fig)
            store_figure(fig, f'concentrations{"_log" if log_scale else ""}')
            plt.close(fig)
            
            fig_other = viz.plot_multi_panel_groups_with_stats(data, selected, group_col, stats_dict,
                                                               ncols=3, plot_type=settings['plot_type'], log_scale=not log_scale,
                                                               show_points=settings['show_points'], ylabel='Concentration (ng/mL)')
            store_figure(fig_other, f'concentrations{"_log" if not log_scale else ""}')
            plt.close(fig_other)
            
            with st.expander("Statistical Summary"):
                st.dataframe(get_significant_differences_summary(stats_dict), hide_index=True)


def render_totals_tab(processed, settings):
    """Render totals tab."""
    st.markdown("### Aggregate Totals")
    
    viz = SphingolipidVisualizer(color_palette=settings['color_palette'], style=settings['plot_style'])
    group_col = processed.structure.group_col
    results = st.session_state.analysis_results
    
    if not group_col:
        return
    
    data = processed.sample_data[processed.sample_data[group_col].notna()].copy()
    is_twoway = results is not None and results.is_twoway
    
    key_totals = ['total_all', 'total_ceramides', 'total_dihydroceramides', 'total_sphingomyelins', 
                  'sphingoid_bases', 'sphingoid_base_phosphates', 'very_long_chain', 'long_chain']
    
    log_scale = st.checkbox("Log10 scale", key="totals_log")
    
    if is_twoway:
        fa_col = results.factor_a_col
        fb_col = results.factor_b_col
        combined = pd.concat([data[[fa_col, fb_col]], processed.totals.loc[data.index]], axis=1)
        available = [t for t in key_totals if t in combined.columns]
        
        tw_stats = {t: results.twoway_totals.get(t) for t in available if t in results.twoway_totals}
        
        fig = viz.plot_twoway_multi_panel(
            data=combined, value_cols=available,
            factor_a_col=fa_col, factor_b_col=fb_col,
            twoway_results=tw_stats, ncols=3,
            factor_a_name=results.factor_a_name, factor_b_name=results.factor_b_name,
            ylabel='Concentration (ng/mL)', show_points=settings['show_points']
        )
        st.pyplot(fig)
        store_figure(fig, 'totals_twoway')
        plt.close(fig)
        
        with st.expander("Two-Way ANOVA Summary"):
            from modules.report_generation import get_twoway_differences_summary
            st.dataframe(get_twoway_differences_summary(tw_stats), hide_index=True)
    else:
        combined = pd.concat([data[[group_col]], processed.totals.loc[data.index]], axis=1)
        available = [t for t in key_totals if t in combined.columns]
        stats_dict = {t: results.totals_results.get(t) for t in available if t in results.totals_results}
        
        fig = viz.plot_multi_panel_groups_with_stats(combined, available, group_col, stats_dict,
                                                     ncols=3, plot_type=settings['plot_type'], log_scale=log_scale,
                                                     show_points=settings['show_points'], ylabel='Concentration (ng/mL)')
        st.pyplot(fig)
        store_figure(fig, f'totals{"_log" if log_scale else ""}')
        plt.close(fig)
        
        fig_other = viz.plot_multi_panel_groups_with_stats(combined, available, group_col, stats_dict,
                                                           ncols=3, plot_type=settings['plot_type'], log_scale=not log_scale,
                                                           show_points=settings['show_points'], ylabel='Concentration (ng/mL)')
        store_figure(fig_other, f'totals{"_log" if not log_scale else ""}')
        plt.close(fig_other)
        
        with st.expander("Statistical Summary"):
            rows = []
            for t in available:
                res = results.totals_results.get(t)
                if res:
                    rows.append({
                        'Total': t.replace('_', ' ').title(),
                        'P-value': f"{res.main_test.pvalue:.4f}",
                        'Significant': 'Yes' if res.main_test.significant else '',
                        'APA': format_apa_statistics(res)
                    })
            st.dataframe(pd.DataFrame(rows), hide_index=True)


def render_percentages_tab(processed, settings):
    """Render percentages tab - comparing % composition across groups."""
    st.markdown("### Sphingolipid Pool Composition (%)")
    st.caption("Compare the percentage each sphingolipid contributes to the total pool across groups")
    
    viz = SphingolipidVisualizer(color_palette=settings['color_palette'], style=settings['plot_style'])
    group_col = processed.structure.group_col
    results = st.session_state.analysis_results
    
    if not group_col:
        st.warning("No group column detected.")
        return
    
    data = processed.sample_data[processed.sample_data[group_col].notna()].copy()
    data = data[data[group_col].astype(str).str.lower() != 'nan']
    percentages = processed.percentages.loc[data.index].copy()
    
    # Get available percentage columns (remove _pct suffix for display)
    pct_cols = [c for c in percentages.columns if c.endswith('_pct')]
    sl_names = [c.replace('_pct', '') for c in pct_cols]
    
    # Check if two-way ANOVA
    is_twoway = results is not None and results.is_twoway
    
    col1, col2 = st.columns([1, 2])
    with col1:
        quick = st.selectbox("Quick select", 
                            ["Top 10 by mean %", "Significant", "Ceramides", "Dihydroceramides", 
                             "Sphingomyelins", "Very Long Chain", "Custom"],
                            key="pct_quick")
    
    if quick == "Top 10 by mean %":
        top_pcts = percentages[pct_cols].mean().nlargest(10).index.tolist()
        selected_pct_cols = top_pcts
    elif quick == "Significant":
        if is_twoway:
            selected_pct_cols = [c for c in pct_cols if c in results.twoway_percentages
                               and (results.twoway_percentages[c].twoway_result.factor_a_pvalue < settings['alpha']
                                    or results.twoway_percentages[c].twoway_result.factor_b_pvalue < settings['alpha']
                                    or results.twoway_percentages[c].twoway_result.interaction_pvalue < settings['alpha'])][:10]
        else:
            selected_pct_cols = [c for c in pct_cols if c in results.percentages_results 
                               and results.percentages_results[c].main_test.significant][:10]
        if not selected_pct_cols:
            st.info("No significant percentage differences found. Showing top 10 by mean %.")
            selected_pct_cols = percentages[pct_cols].mean().nlargest(10).index.tolist()
    elif quick == "Ceramides":
        ceramide_sls = get_ceramides()
        selected_pct_cols = [f'{sl}_pct' for sl in ceramide_sls if f'{sl}_pct' in pct_cols][:10]
    elif quick == "Dihydroceramides":
        dhc_sls = get_dihydroceramides()
        selected_pct_cols = [f'{sl}_pct' for sl in dhc_sls if f'{sl}_pct' in pct_cols][:10]
    elif quick == "Sphingomyelins":
        sm_sls = get_sphingomyelins()
        selected_pct_cols = [f'{sl}_pct' for sl in sm_sls if f'{sl}_pct' in pct_cols][:10]
    elif quick == "Very Long Chain":
        vlc_sls = get_very_long_chain()
        selected_pct_cols = [f'{sl}_pct' for sl in vlc_sls if f'{sl}_pct' in pct_cols][:10]
    else:
        with col2:
            selected_sls = st.multiselect("Select sphingolipids", sl_names, sl_names[:5], key="pct_custom")
            selected_pct_cols = [f'{sl}_pct' for sl in selected_sls]
    
    if not selected_pct_cols:
        st.warning("No sphingolipids selected.")
        return
    
    # Show significant findings summary
    if is_twoway:
        sig_pcts = [c for c in selected_pct_cols if c in results.twoway_percentages
                   and (results.twoway_percentages[c].twoway_result.factor_a_pvalue < settings['alpha']
                        or results.twoway_percentages[c].twoway_result.factor_b_pvalue < settings['alpha']
                        or results.twoway_percentages[c].twoway_result.interaction_pvalue < settings['alpha'])]
    else:
        sig_pcts = [c for c in selected_pct_cols if c in results.percentages_results 
                   and results.percentages_results[c].main_test.significant]
    if sig_pcts:
        sig_names = [c.replace('_pct', '') for c in sig_pcts]
        st.success(f"**Significant differences:** {', '.join(sig_names)}")
    
    # Combine data for plotting
    plot_df = pd.concat([data[[group_col]].reset_index(drop=True), 
                         percentages[selected_pct_cols].reset_index(drop=True)], axis=1)
    
    # Rename columns for cleaner display (remove _pct suffix)
    display_cols = {col: col.replace('_pct', '') for col in selected_pct_cols}
    plot_df_display = plot_df.rename(columns=display_cols)
    display_names = list(display_cols.values())
    
    # =========================================================================
    # SECTION 1: Statistical comparison
    # =========================================================================
    st.markdown("#### Statistical Comparison by Group")
    
    if is_twoway:
        # === TWO-WAY ANOVA PLOTS ===
        fa_col = results.factor_a_col
        fb_col = results.factor_b_col
        fa_name = results.factor_a_name
        fb_name = results.factor_b_name
        
        st.markdown(f"**Two-way ANOVA**: {fa_name} x {fb_name}")
        
        tw_stats = {c.replace('_pct', ''): results.twoway_percentages.get(c) 
                   for c in selected_pct_cols if c in results.twoway_percentages}
        
        # Need display-named columns for plotting
        pct_plot_data = pd.concat([data[[fa_col, fb_col]].reset_index(drop=True),
                                   percentages[selected_pct_cols].reset_index(drop=True)], axis=1)
        pct_plot_data = pct_plot_data.rename(columns=display_cols)
        
        try:
            fig = viz.plot_twoway_multi_panel(
                data=pct_plot_data, value_cols=display_names,
                factor_a_col=fa_col, factor_b_col=fb_col,
                twoway_results=tw_stats, ncols=3,
                factor_a_name=fa_name, factor_b_name=fb_name,
                ylabel='% of Total SL', show_points=settings['show_points']
            )
            st.pyplot(fig)
            store_figure(fig, 'percentages_twoway')
            plt.close(fig)
        except Exception as e:
            st.warning(f"Could not generate two-way percentage plot: {e}")
        
        try:
            with st.expander("Interaction Plots"):
                fig_int = viz.plot_twoway_interaction_multi_panel(
                    data=pct_plot_data, value_cols=display_names,
                    factor_a_col=fa_col, factor_b_col=fb_col,
                    twoway_results=tw_stats, ncols=3,
                    factor_a_name=fa_name, factor_b_name=fb_name,
                    ylabel='% of Total SL'
                )
                st.pyplot(fig_int)
                store_figure(fig_int, 'percentages_interaction')
                plt.close(fig_int)
        except Exception as e:
            st.warning(f"Could not generate interaction plot: {e}")
        
        with st.expander("Two-Way ANOVA Summary"):
            from modules.report_generation import get_twoway_differences_summary
            tw_stats_orig = {c: results.twoway_percentages.get(c) 
                           for c in selected_pct_cols if c in results.twoway_percentages}
            if tw_stats_orig:
                st.dataframe(get_twoway_differences_summary(tw_stats_orig), hide_index=True)
            else:
                st.info("No two-way ANOVA results available for selected percentages.")
    else:
        # === ONE-WAY ANOVA PLOTS ===
        # Get stats results with display names
        stats_dict = {}
        for pct_col in selected_pct_cols:
            display_name = pct_col.replace('_pct', '')
            if pct_col in results.percentages_results:
                stats_dict[display_name] = results.percentages_results[pct_col]
        
        fig = viz.plot_multi_panel_groups_with_stats(
            plot_df_display, display_names, group_col, stats_dict,
            ncols=3, plot_type=settings['plot_type'], log_scale=False,
            show_points=settings['show_points'], ylabel='% of Total SL'
        )
        st.pyplot(fig)
        store_figure(fig, 'percentages')
        plt.close(fig)
    
    # =========================================================================
    # SECTION 2-4: Composition Charts (Pie, Bars, Stacked)
    # For two-way designs, use combined factorial groups (e.g., "Aged-Female")
    # =========================================================================
    all_pct_cols = percentages[pct_cols].mean().nlargest(15).index.tolist()
    
    if is_twoway:
        fa_col_comp = results.factor_a_col
        fb_col_comp = results.factor_b_col
        fa_name_comp = results.factor_a_name
        fb_name_comp = results.factor_b_name
        # Create combined group column for composition charts
        comp_group_col = '_factorial_group_'
        comp_data = data.copy()
        comp_data[comp_group_col] = comp_data[fa_col_comp].astype(str) + ' - ' + comp_data[fb_col_comp].astype(str)
        pie_df = pd.concat([comp_data[[comp_group_col]].reset_index(drop=True),
                            percentages[all_pct_cols].reset_index(drop=True)], axis=1)
        comp_title_suffix = f" ({fa_name_comp} x {fb_name_comp})"
    else:
        comp_group_col = group_col
        pie_df = pd.concat([data[[group_col]].reset_index(drop=True), 
                            percentages[all_pct_cols].reset_index(drop=True)], axis=1)
        comp_title_suffix = ""
    
    st.markdown("---")
    st.markdown("#### Pool Composition - Pie Charts")
    st.caption("Visual breakdown of sphingolipid pool for each group (top contributors)")
    
    fig_pie = viz.plot_composition_pie_charts(
        pie_df, comp_group_col, all_pct_cols,
        title='Sphingolipid Pool Composition by Group' + comp_title_suffix,
        top_n=10, other_threshold=2.0
    )
    st.pyplot(fig_pie)
    store_figure(fig_pie, 'percentages_pie')
    plt.close(fig_pie)
    
    st.markdown("---")
    st.markdown("#### Pool Composition - Bar Comparison")
    st.caption("Side-by-side comparison of sphingolipid percentages across groups")
    
    fig_hbar = viz.plot_composition_horizontal_bars(
        pie_df, comp_group_col, all_pct_cols,
        title='Sphingolipid Pool Composition Comparison' + comp_title_suffix,
        top_n=15, show_values=True
    )
    st.pyplot(fig_hbar)
    store_figure(fig_hbar, 'percentages_bars')
    plt.close(fig_hbar)
    
    st.markdown("---")
    st.markdown("#### Stacked Composition View")
    st.caption("Full sphingolipid pool breakdown for each group")
    
    fig_stacked = viz.plot_composition_stacked_horizontal(
        pie_df, comp_group_col, all_pct_cols,
        title='Complete Sphingolipid Pool Composition' + comp_title_suffix,
        top_n=12
    )
    st.pyplot(fig_stacked)
    store_figure(fig_stacked, 'percentages_stacked')
    plt.close(fig_stacked)
    
    # =========================================================================
    # SECTION 5: Statistical Summary Table
    # =========================================================================
    with st.expander("📊 Statistical Summary"):
        if is_twoway:
            from modules.report_generation import get_twoway_differences_summary
            tw_stats_all = {c: results.twoway_percentages.get(c) 
                          for c in selected_pct_cols if c in results.twoway_percentages}
            if tw_stats_all:
                st.dataframe(get_twoway_differences_summary(tw_stats_all), hide_index=True)
        else:
            rows = []
            for pct_col in selected_pct_cols:
                sl_name = pct_col.replace('_pct', '')
                res = results.percentages_results.get(pct_col)
                if res:
                    n_sig = 0
                    if res.posthoc_test and res.posthoc_test.pairwise_results is not None:
                        n_sig = res.posthoc_test.pairwise_results['significant'].sum()
                    rows.append({
                        'Sphingolipid': sl_name,
                        'Test': res.main_test.test_type.value,
                        'P-value': f"{res.main_test.pvalue:.4f}",
                        'Significant': '✔' if res.main_test.significant else '',
                        'Effect Size': f"{res.main_test.effect_size:.3f}" if res.main_test.effect_size else 'N/A',
                        'Sig. Pairs': n_sig,
                        'APA': format_apa_statistics(res)
                    })
            if rows:
                st.dataframe(pd.DataFrame(rows), hide_index=True)
    
    # Group mean percentages table
    with st.expander("📈 Group Mean Percentages"):
        summary_data = []
        for pct_col in selected_pct_cols:
            sl_name = pct_col.replace('_pct', '')
            for group in plot_df[group_col].unique():
                group_data = plot_df[plot_df[group_col] == group][pct_col]
                summary_data.append({
                    'Sphingolipid': sl_name,
                    'Group': group,
                    'Mean %': f"{group_data.mean():.2f}",
                    'SD': f"{group_data.std():.2f}",
                    'Median %': f"{group_data.median():.2f}"
                })
        
        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            # Pivot for easier reading
            pivot_df = summary_df.pivot_table(
                index='Sphingolipid', 
                columns='Group', 
                values='Mean %',
                aggfunc='first'
            )
            st.dataframe(pivot_df, use_container_width=True)
    
    # Raw data view
    with st.expander("📋 View percentage data"):
        st.dataframe(plot_df, use_container_width=True)


def render_ratios_tab(processed, settings):
    """Render ratios tab with multi-panel display."""
    st.markdown("### Clinical Ratios")
    
    viz = SphingolipidVisualizer(color_palette=settings['color_palette'], style=settings['plot_style'])
    group_col = processed.structure.group_col
    results = st.session_state.analysis_results
    
    if not group_col:
        st.warning("No group column detected.")
        return
    
    data = processed.sample_data[processed.sample_data[group_col].notna()].copy()
    data = data[data[group_col].astype(str).str.lower() != 'nan']
    combined = pd.concat([data[[group_col]], processed.ratios.loc[data.index]], axis=1)
    available = [c for c in processed.ratios.columns if not combined[c].isna().all()]
    
    if not available:
        st.warning("No ratio data available.")
        return
    
    # Check if two-way ANOVA
    is_twoway = results is not None and results.is_twoway
    
    # Quick select options
    col1, col2 = st.columns([1, 2])
    with col1:
        quick = st.selectbox("Quick select", 
                            ["All ratios", "Significant only", "Key ratios", "Custom"],
                            key="ratio_quick")
    
    # Define key clinical ratios for sphingolipids
    key_ratios = [
        # Individual species ratios
        'C16_to_C24_Cer', 'C24_1_to_C24_0_Cer', 
        # Enzyme activity markers
        'Cer_to_DHC_C16', 'Cer_to_DHC_C24', 'SM_to_Cer_C16', 'S1P_to_Sph',
        # Chain length comparisons
        'very_long_to_long', 'very_long_to_medium', 'long_to_medium', 
        'long_to_short', 'medium_to_short',
        # Class comparisons
        'total_Cer_to_SM', 'total_Cer_to_DHC', 'saturated_to_unsaturated'
    ]
    
    if quick == "All ratios":
        selected_ratios = available
    elif quick == "Significant only":
        if is_twoway:
            selected_ratios = [r for r in available if r in results.twoway_ratios
                             and (results.twoway_ratios[r].twoway_result.factor_a_pvalue < settings['alpha']
                                  or results.twoway_ratios[r].twoway_result.factor_b_pvalue < settings['alpha']
                                  or results.twoway_ratios[r].twoway_result.interaction_pvalue < settings['alpha'])]
        else:
            selected_ratios = [r for r in available if r in results.ratios_results 
                             and results.ratios_results[r].main_test.significant]
        if not selected_ratios:
            st.info("No significant ratio differences found. Showing all ratios.")
            selected_ratios = available
    elif quick == "Key ratios":
        selected_ratios = [r for r in key_ratios if r in available]
        if not selected_ratios:
            selected_ratios = available[:6]
    else:
        with col2:
            selected_ratios = st.multiselect("Select ratios", available, 
                                             default=available[:6] if len(available) >= 6 else available,
                                             key="ratio_custom")
    
    if not selected_ratios:
        st.warning("No ratios selected.")
        return
    
    # Show significant findings
    if is_twoway:
        sig_ratios = [r for r in selected_ratios if r in results.twoway_ratios
                     and (results.twoway_ratios[r].twoway_result.factor_a_pvalue < settings['alpha']
                          or results.twoway_ratios[r].twoway_result.factor_b_pvalue < settings['alpha']
                          or results.twoway_ratios[r].twoway_result.interaction_pvalue < settings['alpha'])]
    else:
        sig_ratios = [r for r in selected_ratios if r in results.ratios_results 
                     and results.ratios_results[r].main_test.significant]
    if sig_ratios:
        st.success(f"**Significant differences:** {', '.join(sig_ratios)}")
    
    log_scale = st.checkbox("Log₁₀ scale", key="ratio_log")
    
    if is_twoway:
        # === TWO-WAY ANOVA PLOTS ===
        fa_col = results.factor_a_col
        fb_col = results.factor_b_col
        fa_name = results.factor_a_name
        fb_name = results.factor_b_name
        
        st.markdown(f"**Two-way ANOVA**: {fa_name} x {fb_name}")
        
        tw_stats = {r: results.twoway_ratios.get(r) for r in selected_ratios 
                   if r in results.twoway_ratios}
        
        # Build data with factor columns
        ratio_plot_data = pd.concat([data[[fa_col, fb_col]].reset_index(drop=True),
                                     processed.ratios[selected_ratios].loc[data.index].reset_index(drop=True)], axis=1)
        
        try:
            fig = viz.plot_twoway_multi_panel(
                data=ratio_plot_data, value_cols=selected_ratios,
                factor_a_col=fa_col, factor_b_col=fb_col,
                twoway_results=tw_stats, ncols=3,
                factor_a_name=fa_name, factor_b_name=fb_name,
                ylabel='Ratio', show_points=settings['show_points']
            )
            st.pyplot(fig)
            store_figure(fig, 'ratios_twoway')
            plt.close(fig)
        except Exception as e:
            st.warning(f"Could not generate two-way ratio plot: {e}")
        
        try:
            with st.expander("Interaction Plots"):
                fig_int = viz.plot_twoway_interaction_multi_panel(
                    data=ratio_plot_data, value_cols=selected_ratios,
                    factor_a_col=fa_col, factor_b_col=fb_col,
                    twoway_results=tw_stats, ncols=3,
                    factor_a_name=fa_name, factor_b_name=fb_name,
                    ylabel='Ratio'
                )
                st.pyplot(fig_int)
                store_figure(fig_int, 'ratios_interaction')
                plt.close(fig_int)
        except Exception as e:
            st.warning(f"Could not generate interaction plot: {e}")
        
        with st.expander("Two-Way ANOVA Summary"):
            from modules.report_generation import get_twoway_differences_summary
            if tw_stats:
                st.dataframe(get_twoway_differences_summary(tw_stats), hide_index=True)
            else:
                st.info("No two-way ANOVA results available for selected ratios.")
    else:
        # === ONE-WAY ANOVA PLOTS ===
        # Get stats for selected ratios
        stats_dict = {r: results.ratios_results.get(r) for r in selected_ratios if r in results.ratios_results}
        
        # Multi-panel figure
        fig = viz.plot_multi_panel_groups_with_stats(
            combined, selected_ratios, group_col, stats_dict,
            ncols=3, plot_type=settings['plot_type'], log_scale=log_scale,
            show_points=settings['show_points'], ylabel='Ratio'
        )
        st.pyplot(fig)
        store_figure(fig, f'ratios{"_log" if log_scale else ""}')
        plt.close(fig)
        
        # Generate other version for export
        fig_other = viz.plot_multi_panel_groups_with_stats(
            combined, selected_ratios, group_col, stats_dict,
            ncols=3, plot_type=settings['plot_type'], log_scale=not log_scale,
            show_points=settings['show_points'], ylabel='Ratio'
        )
        store_figure(fig_other, f'ratios{"_log" if not log_scale else ""}')
        plt.close(fig_other)
    
    # Statistical summary
    with st.expander("📊 Statistical Summary"):
        if is_twoway:
            from modules.report_generation import get_twoway_differences_summary
            tw_all = {r: results.twoway_ratios.get(r) for r in selected_ratios 
                     if r in results.twoway_ratios}
            if tw_all:
                st.dataframe(get_twoway_differences_summary(tw_all), hide_index=True)
        else:
            rows = []
            for r in selected_ratios:
                res = results.ratios_results.get(r)
                if res:
                    n_sig = 0
                    if res.posthoc_test and res.posthoc_test.pairwise_results is not None:
                        n_sig = res.posthoc_test.pairwise_results['significant'].sum()
                    rows.append({
                        'Ratio': r,
                        'Test': res.main_test.test_type.value,
                        'P-value': f"{res.main_test.pvalue:.4f}",
                        'Significant': '✔' if res.main_test.significant else '',
                        'Effect Size': f"{res.main_test.effect_size:.3f}" if res.main_test.effect_size else 'N/A',
                        'Sig. Pairs': n_sig,
                        'APA': format_apa_statistics(res)
                    })
            if rows:
                st.dataframe(pd.DataFrame(rows), hide_index=True)
    
    # Group means table
    with st.expander("📈 Group Mean Ratios"):
        summary_data = []
        for ratio in selected_ratios:
            for group in combined[group_col].unique():
                group_data = combined[combined[group_col] == group][ratio]
                summary_data.append({
                    'Ratio': ratio,
                    'Group': group,
                    'Mean': f"{group_data.mean():.3f}",
                    'SD': f"{group_data.std():.3f}",
                    'Median': f"{group_data.median():.3f}"
                })
        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            pivot_df = summary_df.pivot_table(index='Ratio', columns='Group', values='Mean', aggfunc='first')
            st.dataframe(pivot_df, use_container_width=True)
    
    # Raw data
    with st.expander("📋 View ratio data"):
        display_cols = [group_col] + selected_ratios
        st.dataframe(combined[display_cols], use_container_width=True)


def render_statistics_tab(processed, settings):
    """Statistics summary tab."""
    st.markdown("### Statistics Summary")
    results = st.session_state.analysis_results
    is_twoway = results is not None and results.is_twoway
    
    if is_twoway:
        # === TWO-WAY ANOVA SUMMARY ===
        st.markdown(f"**Design**: Two-way ANOVA ({results.factor_a_name} x {results.factor_b_name})")
        
        from modules.report_generation import get_twoway_differences_summary
        from modules.statistical_tests import format_twoway_apa
        
        def _count_sig_twoway(results_dict):
            count = 0
            for r in results_dict.values():
                tw = r.twoway_result
                if (tw.factor_a_pvalue < settings['alpha'] or 
                    tw.factor_b_pvalue < settings['alpha'] or 
                    tw.interaction_pvalue < settings['alpha']):
                    count += 1
            return count
        
        sig_t = _count_sig_twoway(results.twoway_totals)
        sig_s = _count_sig_twoway(results.twoway_individual_sl)
        sig_p = _count_sig_twoway(results.twoway_percentages)
        sig_r = _count_sig_twoway(results.twoway_ratios)
        
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Sig. Totals", f"{sig_t}/{len(results.twoway_totals)}")
        c2.metric("Sig. SLs", f"{sig_s}/{len(results.twoway_individual_sl)}")
        c3.metric("Sig. Percentages", f"{sig_p}/{len(results.twoway_percentages)}")
        c4.metric("Sig. Ratios", f"{sig_r}/{len(results.twoway_ratios)}")
        
        st.markdown("#### Totals")
        if results.twoway_totals:
            st.dataframe(get_twoway_differences_summary(results.twoway_totals), hide_index=True)
        
        st.markdown("#### Individual Sphingolipids (any significant effect)")
        sig_sl = {k: v for k, v in results.twoway_individual_sl.items()
                  if (v.twoway_result.factor_a_pvalue < settings['alpha'] or
                      v.twoway_result.factor_b_pvalue < settings['alpha'] or
                      v.twoway_result.interaction_pvalue < settings['alpha'])}
        if sig_sl:
            st.dataframe(get_twoway_differences_summary(sig_sl), hide_index=True)
            
            with st.expander("APA-Formatted Results"):
                for name, result in list(sig_sl.items())[:10]:
                    st.markdown(f"**{name}**: {format_twoway_apa(result)}")
        else:
            st.info("No significant individual sphingolipid differences (two-way ANOVA).")
        
        st.markdown("#### Ratios")
        if results.twoway_ratios:
            st.dataframe(get_twoway_differences_summary(results.twoway_ratios), hide_index=True)
    else:
        # === ONE-WAY SUMMARY (unchanged) ===
        sig_t = sum(1 for r in results.totals_results.values() if r.main_test.significant)
        sig_s = sum(1 for r in results.individual_sl_results.values() if r.main_test.significant)
        sig_p = sum(1 for r in results.percentages_results.values() if r.main_test.significant)
        sig_r = sum(1 for r in results.ratios_results.values() if r.main_test.significant)
        
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Sig. Totals", f"{sig_t}/{len(results.totals_results)}")
        c2.metric("Sig. SLs", f"{sig_s}/{len(results.individual_sl_results)}")
        c3.metric("Sig. Percentages", f"{sig_p}/{len(results.percentages_results)}")
        c4.metric("Sig. Ratios", f"{sig_r}/{len(results.ratios_results)}")
        
        st.markdown("#### Totals")
        st.dataframe(get_significant_differences_summary(results.totals_results), hide_index=True)
        
        st.markdown("#### Individual Sphingolipids (Significant)")
        sig_sl = {k: v for k, v in results.individual_sl_results.items() if v.main_test.significant}
        if sig_sl:
            st.dataframe(get_significant_differences_summary(sig_sl), hide_index=True)
        else:
            st.info("No significant individual sphingolipid differences.")
        
        st.markdown("#### Percentages (Significant)")
        sig_pct = {k: v for k, v in results.percentages_results.items() if v.main_test.significant}
        if sig_pct:
            sig_pct_display = {k.replace('_pct', ''): v for k, v in sig_pct.items()}
            st.dataframe(get_significant_differences_summary(sig_pct_display), hide_index=True)
        else:
            st.info("No significant percentage differences.")


@st.fragment
def render_export_tab(processed, settings):
    """Export tab."""
    st.markdown("### Export Results")
    results = st.session_state.analysis_results
    report_gen = st.session_state.report_generator
    
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("#### Individual Downloads")
    
    # Get metadata columns (Sample_ID and Group)
        id_cols = []
        if processed.structure.sample_id_col and processed.structure.sample_id_col in processed.sample_data.columns:
            id_cols.append(processed.structure.sample_id_col)
        if processed.structure.group_col and processed.structure.group_col in processed.sample_data.columns:
            id_cols.append(processed.structure.group_col)
    
        metadata = processed.sample_data[id_cols] if id_cols else pd.DataFrame(index=processed.sample_data.index)
        
        # Concentrations with LOD highlighting (Excel)
        conc_df = pd.concat([metadata, processed.concentrations], axis=1)
        conc_buf = create_excel_with_lod_highlighting(
            conc_df, 
            processed.structure.analyte_lod_rows,
            sheet_name="Concentrations"
        )
        st.download_button(
            "📥 Concentrations (Excel) ⭐", 
            conc_buf.getvalue(), 
            "sphingolipid_concentrations.xlsx", 
            "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            key="download_concentrations",
            help="Yellow highlighting shows cells where LOD values were replaced"
        )
        
        # Other exports as CSV (no LOD replacement tracking needed)
        for name, df in [("Totals", processed.totals),
                         ("Ratios", processed.ratios), 
                         ("Percentages", processed.percentages)]:
            export_df = pd.concat([metadata, df], axis=1)
            st.download_button(
                f"📥 {name} (CSV)", 
                export_df.to_csv(index=False), 
                f"sphingolipid_{name.lower()}.csv", 
                "text/csv",
                key=f"download_{name.lower()}"
            )
    
    with col2:
        st.markdown("#### Complete Reports")
        
        # Full data Excel with all sheets and highlighting
        full_excel_buf = create_full_data_excel_with_highlighting(processed, metadata)
        st.download_button(
            "📥 All Data (Excel) ⭐", 
            full_excel_buf.getvalue(),
            f"sphingolipid_all_data_{datetime.now():%Y%m%d}.xlsx",
            "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            key="download_full_excel",
            help="Multi-sheet Excel with Concentrations, Totals, Ratios, Percentages, LOD Summary, and Legend"
        )
        
        if report_gen:
            buf = BytesIO()
            report_gen.save_excel_report(buf)
            buf.seek(0)
            st.download_button("📥 Statistical Report (Excel)", buf.getvalue(),
                              f"statistical_report_{datetime.now():%Y%m%d}.xlsx",
                              "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                              key="download_excel_report")
        
        # Generate all figures and ZIP
        with st.spinner("Generating all figures for export..."):
            # Generate comprehensive figure set for all dropdown options
            all_figures = generate_all_export_figures(processed, results, settings)
            # Merge with any existing figures (e.g., from viewing tabs)
            combined_figures = {**st.session_state.figures, **all_figures}
            zip_bytes = create_results_zip(processed, results, combined_figures, report_gen, metadata)
        
        st.download_button("📥 Complete Package (ZIP)", zip_bytes, 
                          f"sphingolipid_analysis_{datetime.now():%Y%m%d_%H%M}.zip", 
                          "application/zip",
                          key="download_zip")
        
        st.caption(f"📊 Package includes {len(combined_figures)} figures covering all analysis options")
        st.caption("⭐ = Contains yellow cell highlighting for LOD-replaced values")
    
    st.markdown("---")
    st.markdown("#### Summary Figure")
    if report_gen and st.button("🎨 Generate Summary Figure", key="gen_summary_fig"):
        plotter = SignificancePlotter()
        fig = plotter.plot_multi_panel_with_significance(
            processed.sample_data, processed.structure.group_col, report_gen,
            ['Total_All_Sphingolipids', 'Total_Ceramides', 'Total_Dihydroceramides', 'Total_Sphingomyelins'], 2, settings['plot_type'])
        st.pyplot(fig)
        store_figure(fig, 'summary')
        c1, c2 = st.columns(2)
        c1.download_button("📥 PNG", fig_to_bytes(fig), "summary.png", "image/png", key="download_summary_png")
        c2.download_button("📥 PDF", fig_to_bytes(fig, 'pdf'), "summary.pdf", "application/pdf", key="download_summary_pdf")
        plt.close(fig)


def main():
    """Main entry point."""
    init_session_state()
    settings = render_sidebar()
    
    st.markdown("# 🧬 Sphingolipid Analysis Pipeline")
    uploaded = st.file_uploader("Upload Excel/ODS file", ['xlsx', 'xls', 'ods'])
    
    if uploaded:
        # Check if file changed
        file_changed = st.session_state.last_file != uploaded.name
        if file_changed:
            st.session_state.stats_computed = False
            st.session_state.figures = {}
            st.session_state.last_file = uploaded.name
            st.session_state.last_settings = None  # Reset settings tracking for new file
        
        # Check if settings changed (LOD handling, alpha, etc.)
        settings_changed = check_settings_changed(settings)
        
        # Reprocess data if file or settings changed
        if file_changed or settings_changed or st.session_state.processed_data is None:
            with st.spinner("Processing data..." + (" (settings changed)" if settings_changed else "")):
                with tempfile.NamedTemporaryFile(delete=False, suffix=Path(uploaded.name).suffix) as tmp:
                    tmp.write(uploaded.getvalue())
                    tmp_path = tmp.name
                
                try:
                    processor = SphingolipidDataProcessor(lod_handling=settings['lod_handling'], lod_value=settings['lod_value'])
                    processed = processor.load_and_process(tmp_path)
                    st.session_state.processed_data = processed
                    if settings_changed:
                        st.success("✅ Data reprocessed with new settings!")
                    else:
                        st.success("✅ Data loaded!")
                except Exception as e:
                    st.error(f"Error: {e}")
                    return
    
    if st.session_state.processed_data:
        processed = st.session_state.processed_data
        
        quality = validate_data_quality(processed)
        c1, c2, c3 = st.columns(3)
        c1.metric("Samples", quality['n_samples'])
        c2.metric("Sphingolipids", quality['n_sphingolipids_detected'])
        c3.metric("Groups", quality.get('n_groups', 'N/A'))
        
        # Show LOD detection info
        lod_display = {"lod": "LOD value", "half_lod": "LOD/2", "zero": "Zero", "half_min": "Min/2", "drop": "NaN"}
        lod_source = processed.structure.lod_source
        n_lods = len(processed.structure.analyte_lods)
        
        if lod_source == "standards":
            st.caption(f"📊 LODs: **Auto-detected** from standard curves ({n_lods} analytes) | Below LOD → **{lod_display[settings['lod_handling']]}** | α = **{settings['alpha']}**")
        else:
            st.caption(f"📊 LODs: **Default** ({settings['lod_value']} ng/mL) | Below LOD → **{lod_display[settings['lod_handling']]}** | α = **{settings['alpha']}**")
        
        # Show two-way factor info if detected
        n_factors = getattr(processed.structure, 'n_factors', 0)
        if n_factors >= 2:
            factors = processed.structure.factors
            factor_names = list(factors.keys())
            factor_source = getattr(processed.structure, 'factor_source', 'unknown')
            st.info(f"**Two-Way Design Detected**: {factor_names[0]} x {factor_names[1]} "
                    f"(from {factor_source.replace('_', ' ')})")
        
        # Expandable section to view detected LODs
        with st.expander("🔬 View Per-Analyte LODs"):
            if processed.structure.analyte_lods:
                n_samples = quality['n_samples']
                lod_counts = processed.structure.analyte_lod_counts
                
                lod_data = []
                for analyte, lod in sorted(processed.structure.analyte_lods.items()):
                    count = lod_counts.get(analyte, 0)
                    pct = (count / n_samples * 100) if n_samples > 0 else 0
                    lod_data.append({
                        "Analyte": analyte,
                        "LOD (ng/mL)": lod,
                        "Below LOD": count,
                        "% Replaced": f"{pct:.1f}%"
                    })
                
                lod_df = pd.DataFrame(lod_data)
                
                # Highlight rows where many samples were replaced
                st.dataframe(lod_df, hide_index=True, use_container_width=True)
                
                # Show warning if any analytes have high replacement rates
                high_replacement = [d["Analyte"] for d in lod_data if float(d["% Replaced"].rstrip('%')) > 50]
                if high_replacement:
                    st.warning(f"⚠️ **High replacement rate (>50%):** {', '.join(high_replacement[:5])}{'...' if len(high_replacement) > 5 else ''}")
        
        with st.spinner("Computing statistics..."):
            compute_all_statistics(processed, settings)
        
        tabs = st.tabs(["📈 Concentrations", "📊 Totals", "📉 Percentages", "🔢 Ratios", "📊 Statistics", "💾 Export"])
        with tabs[0]: render_concentrations_tab(processed, settings)
        with tabs[1]: render_totals_tab(processed, settings)
        with tabs[2]: render_percentages_tab(processed, settings)
        with tabs[3]: render_ratios_tab(processed, settings)
        with tabs[4]: render_statistics_tab(processed, settings)
        with tabs[5]: render_export_tab(processed, settings)
    
    else:
        # Show expected input format when no data is loaded
        st.markdown("---")
        st.markdown("## 📋 Expected Input Format")
        
        st.markdown("### LC-MS sheet")
        st.markdown("""
        - **Rows:** First few are Std curves, followed by samples
        - **Columns:** "Data Filename" (Std curve and sample names), followed by sphingolipid species
        """)


        st.markdown("**Example:**")
        
        # Create example dataframe
        example_lcms_df = pd.DataFrame({
            'Data Filename': ['Std 1 ngmL', 'Std 3  ngmL', 'Std 10 ngmL'],
            'C16 Cer': [1.3, 3.0, 9.9],
            'C24-0 Cer': [1.1, 2.9, 10.1],
            'C16-SM': ['----', 3.1, 10.0],
            'S-d18-1': ['----', '----', 10.2],
            '...': ['...', '...', '...']
        })
        st.dataframe(example_lcms_df, hide_index=True, use_container_width=False)

        st.markdown("### Sample sheet")
        st.markdown("""
        - **Rows:** Samples
        - **Columns:** Sphingolipid species (matching panel names)
        - **First column(s):** Sample ID, Group/Type
        - **Values:** Concentrations (typically ng/mL)
        - **Below LOD:** Can be "-----", "LOD", "BLQ", "ND", etc.
        """)
        
        st.markdown("**Example:**")
        
        # Create example dataframe
        example_df = pd.DataFrame({
            'Type': ['Aged-1', 'Aged-2', 'Young-1'],
            'Sample_ID': ['S001', 'S002', 'S003'],
            'C16 Cer': [125.4, 142.8, 98.6],
            'C24-0 Cer': [312.5, 287.9, '-----'],
            'C16-SM': [1520.3, 1380.7, 1245.2],
            'S-d18-1': [15.2, 18.4, 12.8],
            '...': ['...', '...', '...']
        })
        st.dataframe(example_df, hide_index=True, use_container_width=False)
        
        st.info("""💡 **Auto-Detection Features:**
- Sphingolipid columns are auto-detected from your data
- Group assignments are inferred automatically  
- **Per-analyte LODs** are extracted from standard curves in the "LC-MS data" sheet (Std rows)
""")


if __name__ == "__main__":
    main()
