"""
Data Processing Module for Sphingolipid Analysis
=================================================

Handles:
1. Loading data from Excel/ODS files
2. Detecting and parsing data structure
3. Auto-detecting per-analyte LOD from standard curves
4. Data cleaning (handling LOD, missing values)
5. Normalization and calculations
6. Aggregation by sphingolipid groups
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any
from dataclasses import dataclass, field
import re
import warnings

# Import sphingolipid configuration
import sys
sys.path.insert(0, str(Path(__file__).parent.parent))
from config.sphingolipid_species import (
    SPHINGOLIPID_PANEL, ANALYSIS_GROUPS, CLINICAL_RATIOS,
    get_all_species, validate_columns, get_species_info,
    SphingoClass, ChainLength
)


@dataclass
class DataStructureInfo:
    """Information about detected data structure."""
    n_rows: int
    n_cols: int
    sample_id_col: Optional[str] = None
    group_col: Optional[str] = None
    other_metadata_cols: List[str] = field(default_factory=list)
    sphingolipid_cols: List[str] = field(default_factory=list)
    unrecognized_cols: List[str] = field(default_factory=list)
    has_standards: bool = False
    standard_rows: List[int] = field(default_factory=list)
    data_start_row: int = 0
    units: str = "ng/mL"
    sheet_used: Optional[str] = None  # Which sheet the data was loaded from
    analyte_lods: Dict[str, float] = field(default_factory=dict)  # Per-analyte LODs
    lod_source: str = "default"  # "standards" or "default"


@dataclass
class ProcessedData:
    """Container for processed sphingolipid data."""
    # Core data
    raw_data: pd.DataFrame
    sample_data: pd.DataFrame  # Samples only, no standards
    
    # Structure info
    structure: DataStructureInfo
    
    # Calculated values
    concentrations: pd.DataFrame  # Individual sphingolipid concentrations
    percentages: pd.DataFrame     # % of total for each sphingolipid
    totals: pd.DataFrame          # Group totals (ceramides, SM, etc.)
    ratios: pd.DataFrame          # Clinical/research ratios
    
    # Summary statistics by group
    group_summaries: Optional[Dict[str, pd.DataFrame]] = None


class SphingolipidDataProcessor:
    """
    Process sphingolipid LC-MS data from Excel files.
    
    Usage:
        processor = SphingolipidDataProcessor()
        processed = processor.load_and_process("data.xlsx")
    """
    
    # Common patterns for identifying columns
    SAMPLE_ID_PATTERNS = [
        r'^sample[\s_-]*(id|name|#)?$',
        r'^name$',
        r'^id$',
        r'^specimen$',
    ]
    
    GROUP_PATTERNS = [
        r'^(group|type|category|class|condition|treatment)$',
        r'^sample[\s_-]*type$',
    ]
    
    STANDARD_PATTERNS = [
        r'std[\s_-]*\d+',
        r'standard',
        r'cal[\s_-]*\d+',
        r'calibr',
        r'qc[\s_-]*\d*',
        r'blank',
    ]
    
    LOD_VALUES = ['-----', '----', '---', 'LOD', 'BLQ', 'ND', 'N/D', '<LOD', '<LOQ', 'BLOQ', '']
    
    def __init__(
        self,
        lod_handling: str = "half_lod",  # "zero", "lod", "half_lod", "half_min", "drop"
        lod_value: float = 0.1,  # Default/fallback LOD value
        custom_sphingolipids: Optional[Dict] = None
    ):
        """
        Initialize the processor.
        
        Args:
            lod_handling: How to handle below-LOD values
                - "zero": Replace with 0
                - "lod": Replace with the analyte's LOD value
                - "half_lod": Replace with half the analyte's LOD value
                - "half_min": Replace with half the minimum detected value
                - "drop": Keep as NaN
            lod_value: Default/fallback LOD value when auto-detection fails
            custom_sphingolipids: Additional sphingolipid definitions to merge with panel
        """
        self.lod_handling = lod_handling
        self.lod_value = lod_value  # Fallback LOD
        
        # Merge custom sphingolipids if provided
        self.sphingolipid_panel = SPHINGOLIPID_PANEL.copy()
        if custom_sphingolipids:
            self.sphingolipid_panel.update(custom_sphingolipids)
    
    # Preferred sheet names for processed data (in order of preference)
    PREFERRED_SHEET_PATTERNS = [
        r'serum[\s_-]*c',
        r'plasma[\s_-]*c', 
        r'processed',
        r'corrected',
        r'final',
        r'results',
        r'dry[\s_-]*(content|fecal)',
        r'liver',
    ]
    
    # Patterns that indicate the LC-MS raw data sheet with standards
    LCMS_DATA_PATTERNS = [
        r'lc[\s_-]*ms[\s_-]*data',
        r'^data$',
        r'raw[\s_-]*data',
        r'standards',
    ]
    
    def load_file(
        self, 
        filepath: Union[str, Path],
        sheet_name: Union[str, int, None] = None
    ) -> Tuple[pd.DataFrame, str]:
        """
        Load data from Excel/ODS file.
        
        Args:
            filepath: Path to the file
            sheet_name: Specific sheet to load. If None, auto-detects best sheet.
            
        Returns:
            Tuple of (DataFrame, sheet_name_used)
        """
        filepath = Path(filepath)
        
        if filepath.suffix.lower() == '.ods':
            engine = 'odf'
        elif filepath.suffix.lower() in ['.xlsx', '.xlsm']:
            engine = 'openpyxl'
        elif filepath.suffix.lower() == '.xls':
            engine = 'xlrd'
        else:
            engine = None
        
        # Get list of sheet names
        xlsx = pd.ExcelFile(filepath, engine=engine)
        available_sheets = xlsx.sheet_names
        
        # Determine which sheet to use
        if sheet_name is not None:
            # User specified a sheet
            selected_sheet = sheet_name
        else:
            # Auto-detect: look for preferred sheet names
            selected_sheet = self._find_best_sheet(available_sheets)
        
        df = pd.read_excel(filepath, sheet_name=selected_sheet, header=None, engine=engine)
        
        # Return both the dataframe and which sheet was used
        sheet_used = selected_sheet if isinstance(selected_sheet, str) else available_sheets[selected_sheet]
        return df, sheet_used
    
    def _find_lcms_data_sheet(self, sheet_names: List[str]) -> Optional[str]:
        """Find the LC-MS data sheet that contains standard curves."""
        for pattern in self.LCMS_DATA_PATTERNS:
            for sheet in sheet_names:
                if re.search(pattern, sheet, re.IGNORECASE):
                    return sheet
        
        # Often the first sheet is the LC-MS data
        if sheet_names:
            first_lower = sheet_names[0].lower()
            if 'lc' in first_lower or 'data' in first_lower or 'ms' in first_lower:
                return sheet_names[0]
        
        return None
    
    def _find_best_sheet(self, sheet_names: List[str]) -> Union[str, int]:
        """
        Find the best sheet to use based on naming patterns.
        
        Prefers sheets like 'Serum C' over raw 'LC-MS data'.
        Avoids sheets that look like analysis output (Overview, Report, etc.)
        """
        # Sheets to avoid (analysis output sheets or raw data)
        avoid_patterns = [
            r'overview',
            r'report',
            r'summary',
            r'lc[\s_-]*ms[\s_-]*data',  # Raw LC-MS data sheet
        ]
        
        # Filter out sheets to avoid
        valid_sheets = []
        for sheet in sheet_names:
            should_avoid = any(re.search(p, sheet, re.IGNORECASE) for p in avoid_patterns)
            if not should_avoid:
                valid_sheets.append(sheet)
        
        # If no valid sheets remain, use all sheets
        if not valid_sheets:
            valid_sheets = sheet_names
        
        # Check for preferred patterns among valid sheets
        for pattern in self.PREFERRED_SHEET_PATTERNS:
            for sheet in valid_sheets:
                if re.search(pattern, sheet, re.IGNORECASE):
                    return sheet
        
        # If only one valid sheet, use it
        if len(valid_sheets) == 1:
            return valid_sheets[0]
        
        # If multiple sheets and no match, prefer second sheet 
        # (often first is raw LC-MS, second is processed)
        if len(valid_sheets) >= 2:
            # But check if first sheet looks like raw data
            first_lower = valid_sheets[0].lower()
            if 'lc-ms' in first_lower or 'raw' in first_lower or 'data' in first_lower:
                return valid_sheets[1]
        
        # Default to first valid sheet
        return valid_sheets[0] if valid_sheets else 0
    
    def get_available_sheets(self, filepath: Union[str, Path]) -> List[str]:
        """Get list of available sheet names in a file."""
        filepath = Path(filepath)
        
        if filepath.suffix.lower() == '.ods':
            engine = 'odf'
        elif filepath.suffix.lower() in ['.xlsx', '.xlsm']:
            engine = 'openpyxl'
        elif filepath.suffix.lower() == '.xls':
            engine = 'xlrd'
        else:
            engine = None
        
        xlsx = pd.ExcelFile(filepath, engine=engine)
        return xlsx.sheet_names
    
    def _extract_lods_from_standards(
        self, 
        filepath: Union[str, Path],
        sphingolipid_cols: List[str]
    ) -> Dict[str, float]:
        """
        Extract per-analyte LOD from standard curve rows in LC-MS data sheet.
        
        Looks for rows with "Std X ng/mL" pattern and determines the lowest
        successful standard (with a valid numeric reading) for each analyte.
        
        Args:
            filepath: Path to Excel file
            sphingolipid_cols: List of sphingolipid column names to look for
            
        Returns:
            Dict mapping analyte names to their LOD values (in ng/mL)
        """
        analyte_lods = {}
        
        try:
            # Get available sheets
            sheets = self.get_available_sheets(filepath)
            
            # Find LC-MS data sheet
            lcms_sheet = self._find_lcms_data_sheet(sheets)
            if lcms_sheet is None:
                return analyte_lods
            
            # Load the LC-MS data sheet
            filepath = Path(filepath)
            if filepath.suffix.lower() == '.ods':
                engine = 'odf'
            elif filepath.suffix.lower() in ['.xlsx', '.xlsm']:
                engine = 'openpyxl'
            elif filepath.suffix.lower() == '.xls':
                engine = 'xlrd'
            else:
                engine = None
            
            df = pd.read_excel(filepath, sheet_name=lcms_sheet, header=None, engine=engine)
            
            # Find header row (row containing analyte names)
            header_row = self._find_header_row(df)
            if header_row < 0:
                return analyte_lods
            
            # Set column names
            df.columns = df.iloc[header_row].astype(str).str.strip()
            df = df.iloc[header_row + 1:].reset_index(drop=True)
            
            # Find the column containing "Data Filename" or standard identifiers
            std_col = None
            for col in df.columns[:5]:  # Check first few columns
                col_str = str(col).lower()
                if 'data' in col_str or 'filename' in col_str or 'name' in col_str:
                    # Check if this column has Std values
                    if df[col].astype(str).str.contains('Std', case=False, na=False).any():
                        std_col = col
                        break
            
            if std_col is None:
                # Try first column
                first_col = df.columns[0]
                if df[first_col].astype(str).str.contains('Std', case=False, na=False).any():
                    std_col = first_col
            
            if std_col is None:
                return analyte_lods
            
            # Parse standard rows and extract concentrations
            # Pattern matches: "Std 1 ng/mL", "Std 3  ng/mL", "Std 10 ng/mL", etc.
            std_pattern = re.compile(r'Std\s*(\d+(?:\.\d+)?)\s*(?:ng/?mL)?', re.IGNORECASE)
            
            std_rows = []  # List of (row_index, concentration)
            for idx in df.index:
                val = str(df.loc[idx, std_col])
                match = std_pattern.search(val)
                if match:
                    conc = float(match.group(1))
                    std_rows.append((idx, conc))
            
            if not std_rows:
                return analyte_lods
            
            # Sort by concentration (lowest first)
            std_rows.sort(key=lambda x: x[1])
            
            # For each analyte column, find the lowest concentration with a valid value
            for col in df.columns:
                col_str = str(col).strip()
                
                # Check if this column is an analyte we care about
                if col_str not in sphingolipid_cols and not self._is_analyte_column(col_str):
                    continue
                
                # Find lowest successful standard for this analyte
                for idx, conc in std_rows:
                    val = df.loc[idx, col]
                    if self._is_valid_measurement(val):
                        analyte_lods[col_str] = conc
                        break
            
            return analyte_lods
            
        except Exception as e:
            warnings.warn(f"Could not extract LODs from standards: {e}")
            return analyte_lods
    
    def _is_valid_measurement(self, val) -> bool:
        """Check if a value is a valid measurement (not blank, not '-----', not NaN)."""
        if pd.isna(val):
            return False
        
        if isinstance(val, (int, float)):
            return not np.isnan(val) and val > 0
        
        if isinstance(val, str):
            val_clean = val.strip()
            if val_clean in self.LOD_VALUES or val_clean.lower() in [v.lower() for v in self.LOD_VALUES]:
                return False
            try:
                num = float(val_clean)
                return num > 0
            except ValueError:
                return False
        
        return False
    
    def _is_analyte_column(self, col: str) -> bool:
        """Check if column name looks like an analyte."""
        if not col or pd.isna(col):
            return False
        
        col_str = str(col).strip()
        
        # Check if it's in the sphingolipid panel
        if col_str in self.sphingolipid_panel:
            return True
        
        # Common patterns for sphingolipid analyte columns
        analyte_patterns = [
            r'C\d+.*Cer',     # Ceramides (C16 Cer, C24-0 Cer, etc.)
            r'C\d+.*SM',      # Sphingomyelins
            r'C\d+.*DHC',     # Dihydroceramides
            r'S1P',           # Sphingosine-1-phosphate
            r'S-d\d+',        # Sphingoid bases
            r'GLU',           # Glucosylceramide
            r'3KDHS',         # 3-Ketosphinganine
            r'C\d+.*CP',      # Ceramide-1-phosphates
        ]
        
        return any(re.search(p, col_str, re.IGNORECASE) for p in analyte_patterns)
    
    def detect_structure(self, df: pd.DataFrame) -> DataStructureInfo:
        """
        Auto-detect the structure of the data file.
        
        Identifies:
        - Header row
        - Sample ID and group columns
        - Sphingolipid columns
        - Standard/calibration rows
        """
        info = DataStructureInfo(n_rows=len(df), n_cols=len(df.columns))
        
        # Find header row (row containing sphingolipid names)
        header_row = self._find_header_row(df)
        info.data_start_row = header_row + 1
        
        # Get column names
        if header_row >= 0:
            headers = df.iloc[header_row].astype(str).tolist()
            # Handle multi-row headers
            if header_row > 0:
                # Check if row above has units info
                prev_row = df.iloc[header_row - 1].astype(str).tolist()
                for i, (h, p) in enumerate(zip(headers, prev_row)):
                    if 'ng' in str(p).lower() or 'conc' in str(p).lower():
                        info.units = str(p)
                        break
        else:
            headers = [f"col_{i}" for i in range(len(df.columns))]
        
        # Classify columns
        known_species = set(self.sphingolipid_panel.keys())
        
        for i, col in enumerate(headers):
            col_str = str(col).strip()
            col_lower = col_str.lower()
            
            # Skip empty columns
            if not col_str or col_str == 'nan':
                continue
            
            # Check for sample ID patterns
            if any(re.match(p, col_lower) for p in self.SAMPLE_ID_PATTERNS):
                info.sample_id_col = col_str
                continue
            
            # Check for group patterns
            if any(re.match(p, col_lower) for p in self.GROUP_PATTERNS):
                info.group_col = col_str
                continue
            
            # Check if it's a known sphingolipid
            if col_str in known_species:
                info.sphingolipid_cols.append(col_str)
                continue
            
            # Check if it looks like a sphingolipid (but not in our panel)
            if self._looks_like_sphingolipid(col_str):
                info.sphingolipid_cols.append(col_str)
                continue
            
            # Otherwise it's unrecognized
            info.unrecognized_cols.append(col_str)
        
        # Find standard rows
        info.standard_rows = self._find_standard_rows(df, header_row)
        info.has_standards = len(info.standard_rows) > 0
        
        # If no group column found, try to infer from data
        if info.group_col is None:
            info.group_col = self._infer_group_column(df, headers, info)
        
        return info
    
    def _find_header_row(self, df: pd.DataFrame) -> int:
        """Find the row containing column headers."""
        known_sphingolipids = set(self.sphingolipid_panel.keys())
        
        # Check first 10 rows
        for i in range(min(10, len(df))):
            row_values = df.iloc[i].astype(str).tolist()
            matches = sum(1 for v in row_values if v.strip() in known_sphingolipids)
            
            # If we find several sphingolipid names, this is likely the header
            if matches >= 3:
                return i
            
            # Also check for patterns that look like sphingolipids
            pattern_matches = sum(1 for v in row_values if self._looks_like_sphingolipid(str(v)))
            if pattern_matches >= 5:
                return i
        
        # Default to first row
        return 0
    
    def _looks_like_sphingolipid(self, name: str) -> bool:
        """Check if a name looks like a sphingolipid."""
        if not name or pd.isna(name):
            return False
        
        name = str(name).strip()
        
        # Common sphingolipid patterns
        patterns = [
            r'^C\d+[:\-]?\d*\s*(Cer|SM|DHC|CP)',
            r'^(S|DHS|Sph)[:\-]?d\d+',
            r'^S1P',
            r'^GLU[\-]?C',
            r'^3KDHS',
            r'Ceramide',
            r'Sphingomyelin',
            r'Sphingosine',
        ]
        
        return any(re.search(p, name, re.IGNORECASE) for p in patterns)
    
    def _find_standard_rows(self, df: pd.DataFrame, header_row: int) -> List[int]:
        """Find rows that appear to be standards/calibrators."""
        standard_rows = []
        
        # Look at rows after header
        for i in range(header_row + 1, len(df)):
            row_str = ' '.join(df.iloc[i].astype(str).tolist()[:3]).lower()
            
            # Check for standard patterns
            if any(re.search(p, row_str) for p in self.STANDARD_PATTERNS):
                standard_rows.append(i)
        
        return standard_rows
    
    def _infer_group_column(
        self, 
        df: pd.DataFrame, 
        headers: List[str],
        info: DataStructureInfo
    ) -> Optional[str]:
        """Try to infer which column contains group information."""
        # Look for columns that might be groups but weren't matched by patterns
        potential_group_cols = []
        
        for i, col in enumerate(headers):
            col_str = str(col).strip()
            
            # Skip if already identified or looks like sphingolipid
            if col_str in info.sphingolipid_cols or col_str == info.sample_id_col:
                continue
            
            # Skip if empty/nan
            if not col_str or col_str.lower() == 'nan':
                continue
            
            # Get column data (skip header row)
            start_row = info.data_start_row
            if start_row < len(df) and i < len(df.columns):
                col_data = df.iloc[start_row:, i]
                
                # Check if it looks like group data
                if self._looks_like_group_column(col_data):
                    potential_group_cols.append(col_str)
        
        return potential_group_cols[0] if potential_group_cols else None
    
    def _looks_like_group_column(self, series: pd.Series) -> bool:
        """Check if a series looks like it contains group information."""
        # Drop NaN values
        non_null = series.dropna().astype(str)
        if len(non_null) < 2:
            return False
        
        # Get unique values
        unique_vals = non_null.unique()
        
        # Groups typically have few unique values (2-20)
        if len(unique_vals) < 2 or len(unique_vals) > 20:
            return False
        
        # Groups typically have repeated values
        if len(unique_vals) == len(non_null):
            return False  # All unique = probably IDs
        
        # Check that values look like group names (not numbers)
        for val in unique_vals:
            val_str = str(val).strip().lower()
            if val_str in ['nan', '', 'none']:
                continue
            try:
                float(val_str)
                return False  # Pure numbers aren't group names
            except ValueError:
                pass
        
        return True
    
    def clean_data(
        self, 
        df: pd.DataFrame, 
        structure: DataStructureInfo
    ) -> pd.DataFrame:
        """
        Clean the raw data:
        - Set proper column names
        - Remove standard rows
        - Handle LOD values (per-analyte)
        - Convert to numeric
        - Remove rows with NaN/empty group values
        """
        # Get header row and set columns
        header_row = structure.data_start_row - 1
        if header_row >= 0:
            df.columns = df.iloc[header_row].astype(str).str.strip()
            df = df.iloc[structure.data_start_row:].reset_index(drop=True)
        
        # Remove standard rows (adjust indices after removing header)
        adjusted_std_rows = [r - structure.data_start_row for r in structure.standard_rows]
        adjusted_std_rows = [r for r in adjusted_std_rows if r >= 0]
        if adjusted_std_rows:
            df = df.drop(index=adjusted_std_rows).reset_index(drop=True)
        
        # Handle LOD values and convert to numeric (using per-analyte LODs)
        for col in structure.sphingolipid_cols:
            if col in df.columns:
                # Get this analyte's LOD (or use fallback)
                analyte_lod = structure.analyte_lods.get(col, self.lod_value)
                df[col] = self._clean_numeric_column(df[col], analyte_lod)
        
        # Remove rows where group column is NaN or empty
        if structure.group_col and structure.group_col in df.columns:
            # Convert to string and check for empty/nan values
            group_series = df[structure.group_col].astype(str).str.strip()
            valid_mask = (
                df[structure.group_col].notna() & 
                (group_series != '') & 
                (group_series.str.lower() != 'nan') &
                (group_series.str.lower() != 'none')
            )
            df = df[valid_mask].reset_index(drop=True)
        
        return df
    
    def _clean_numeric_column(
        self, 
        series: pd.Series,
        analyte_lod: float
    ) -> pd.Series:
        """
        Clean a single numeric column.
        
        Handles below-LOD markers (like '-----', 'LOD', 'BLQ') based on lod_handling,
        using the analyte-specific LOD value.
        
        Args:
            series: The column data
            analyte_lod: The LOD value for this specific analyte
        """
        # Convert series to object type first to avoid downcasting issues
        series = series.astype(object)
        
        # Replace LOD marker values with None (will become NaN after to_numeric)
        lod_values_all = self.LOD_VALUES + [v.lower() for v in self.LOD_VALUES if v]
        mask = series.astype(str).str.strip().isin(lod_values_all)
        series = series.where(~mask, None)
        
        # Convert to numeric
        series = pd.to_numeric(series, errors='coerce')
        
        # Handle below-LOD values based on strategy
        if self.lod_handling == "zero":
            series = series.fillna(0)
        elif self.lod_handling == "lod":
            # Use the analyte-specific LOD value
            series = series.fillna(analyte_lod)
        elif self.lod_handling == "half_lod":
            # Use half the analyte-specific LOD value
            series = series.fillna(analyte_lod / 2)
        elif self.lod_handling == "half_min":
            min_val = series[series > 0].min()
            if pd.notna(min_val):
                series = series.fillna(min_val / 2)
            else:
                series = series.fillna(analyte_lod / 2)
        # "drop" leaves NaN
        
        return series
    
    def calculate_totals(
        self, 
        df: pd.DataFrame,
        sphingolipid_cols: List[str]
    ) -> pd.DataFrame:
        """Calculate aggregate totals for sphingolipid groups."""
        totals = pd.DataFrame(index=df.index)
        
        # Total of all sphingolipids
        totals['total_all'] = df[sphingolipid_cols].sum(axis=1)
        
        # Totals by sphingo class
        for group_name, group_species in ANALYSIS_GROUPS.items():
            available = [s for s in group_species if s in sphingolipid_cols]
            if available:
                totals[group_name] = df[available].sum(axis=1)
        
        return totals
    
    def calculate_percentages(
        self, 
        df: pd.DataFrame,
        sphingolipid_cols: List[str]
    ) -> pd.DataFrame:
        """Calculate each sphingolipid as percentage of total pool."""
        percentages = pd.DataFrame(index=df.index)
        
        total = df[sphingolipid_cols].sum(axis=1)
        
        for col in sphingolipid_cols:
            # Avoid division by zero
            pct = df[col] / total.replace(0, np.nan) * 100
            percentages[f'{col}_pct'] = pct.round(2)
        
        return percentages
    
    def calculate_ratios(
        self,
        df: pd.DataFrame,
        sphingolipid_cols: List[str],
        lod_value: float = 0.1
    ) -> pd.DataFrame:
        """Calculate clinical/research ratios."""
        ratios = pd.DataFrame(index=df.index)
        
        for ratio_name, ratio_info in CLINICAL_RATIOS.items():
            num_cols = ratio_info['numerator']
            den_cols = ratio_info['denominator']
            
            # Handle both list of column names and function calls
            if isinstance(num_cols, str):
                # It's a reference to an analysis group
                num_cols = ANALYSIS_GROUPS.get(num_cols, [])
            if isinstance(den_cols, str):
                den_cols = ANALYSIS_GROUPS.get(den_cols, [])
            
            # Get available columns
            num_available = [c for c in num_cols if c in sphingolipid_cols]
            den_available = [c for c in den_cols if c in sphingolipid_cols]
            
            if num_available and den_available:
                numerator = df[num_available].sum(axis=1)
                denominator = df[den_available].sum(axis=1)
                
                # Replace zeros with small value for ratio calculation
                denominator = denominator.replace(0, lod_value / 2)
                
                ratios[ratio_name] = numerator / denominator
        
        return ratios
    
    def calculate_group_summaries(
        self,
        df: pd.DataFrame,
        group_col: str,
        value_cols: List[str]
    ) -> Dict[str, pd.DataFrame]:
        """Calculate summary statistics for each group."""
        summaries = {}
        
        for col in value_cols:
            if col not in df.columns:
                continue
            
            summary = df.groupby(group_col)[col].agg([
                'count', 'mean', 'std', 'median', 'min', 'max'
            ]).round(4)
            
            # Add SEM
            summary['sem'] = df.groupby(group_col)[col].sem().round(4)
            
            summaries[col] = summary
        
        return summaries
    
    def load_and_process(
        self,
        filepath: Union[str, Path],
        sheet_name: Union[str, int, None] = None,
        group_col: Optional[str] = None
    ) -> ProcessedData:
        """
        Main entry point: load and fully process a sphingolipid data file.
        
        Args:
            filepath: Path to Excel/ODS file
            sheet_name: Sheet to load (name or index). If None, auto-detects best sheet.
            group_col: Override auto-detected group column
            
        Returns:
            ProcessedData object with all calculations
        """
        filepath = Path(filepath)
        
        # Load raw data (auto-detects best sheet if not specified)
        raw_df, sheet_used = self.load_file(filepath, sheet_name)
        
        # Detect structure
        structure = self.detect_structure(raw_df)
        
        # Store which sheet was used
        structure.sheet_used = sheet_used
        
        # Override group col if specified
        if group_col:
            structure.group_col = group_col
        
        # ================================================================
        # Extract per-analyte LODs from standard curves in LC-MS data sheet
        # ================================================================
        analyte_lods = self._extract_lods_from_standards(filepath, structure.sphingolipid_cols)
        
        if analyte_lods:
            structure.analyte_lods = analyte_lods
            structure.lod_source = "standards"
            print(f"✓ Auto-detected LODs for {len(analyte_lods)} analytes from standard curves")
        else:
            # Use default LOD for all analytes
            structure.analyte_lods = {col: self.lod_value for col in structure.sphingolipid_cols}
            structure.lod_source = "default"
            print(f"⚠ Using default LOD ({self.lod_value} ng/mL) for all analytes")
        
        # Clean data (now uses per-analyte LODs)
        clean_df = self.clean_data(raw_df.copy(), structure)
        
        # Extract sample data (with metadata)
        sample_df = clean_df.copy()
        
        # Get concentration data (only sphingolipid columns)
        sl_cols = [c for c in structure.sphingolipid_cols if c in sample_df.columns]
        concentrations = sample_df[sl_cols].copy()
        
        # Calculate derived values
        percentages = self.calculate_percentages(sample_df, sl_cols)
        totals = self.calculate_totals(sample_df, sl_cols)
        ratios = self.calculate_ratios(sample_df, sl_cols, lod_value=self.lod_value)
        
        # Group summaries if group column exists
        group_summaries = None
        if structure.group_col and structure.group_col in sample_df.columns:
            all_value_cols = sl_cols + list(totals.columns) + list(ratios.columns)
            
            # Merge for summary calculation
            full_df = pd.concat([sample_df, totals, ratios], axis=1)
            group_summaries = self.calculate_group_summaries(
                full_df, 
                structure.group_col,
                all_value_cols
            )
        
        return ProcessedData(
            raw_data=raw_df,
            sample_data=sample_df,
            structure=structure,
            concentrations=concentrations,
            percentages=percentages,
            totals=totals,
            ratios=ratios,
            group_summaries=group_summaries
        )
    
    def get_analysis_dataframe(
        self,
        processed: ProcessedData,
        include_percentages: bool = True,
        include_totals: bool = True,
        include_ratios: bool = True
    ) -> pd.DataFrame:
        """
        Get a complete dataframe for analysis.
        
        Combines sample data with calculated values.
        """
        parts = [processed.sample_data]
        
        if include_totals:
            parts.append(processed.totals)
        if include_percentages:
            parts.append(processed.percentages)
        if include_ratios:
            parts.append(processed.ratios)
        
        return pd.concat(parts, axis=1)


def validate_data_quality(processed: ProcessedData) -> Dict[str, Any]:
    """
    Assess data quality and return summary metrics.
    """
    quality = {
        'n_samples': len(processed.sample_data),
        'n_sphingolipids_detected': len(processed.structure.sphingolipid_cols),
        'n_groups': None,
        'missing_pct': {},
        'zero_pct': {},
        'lod_source': processed.structure.lod_source,
        'analyte_lods': processed.structure.analyte_lods,
    }
    
    # Group info
    if processed.structure.group_col:
        groups = processed.sample_data[processed.structure.group_col].unique()
        quality['n_groups'] = len(groups)
        quality['groups'] = list(groups)
    
    # Missing and zero percentages per sphingolipid
    for col in processed.structure.sphingolipid_cols:
        if col in processed.concentrations.columns:
            data = processed.concentrations[col]
            quality['missing_pct'][col] = (data.isna().sum() / len(data) * 100).round(1)
            quality['zero_pct'][col] = ((data == 0).sum() / len(data) * 100).round(1)
    
    return quality


if __name__ == "__main__":
    # Demo/test code
    print("Sphingolipid Data Processor - Ready")
    print("Example usage:")
    print("  processor = SphingolipidDataProcessor(lod_handling='half_lod')")
    print("  processed = processor.load_and_process('data.xlsx')")
    print("  print(validate_data_quality(processed))")
