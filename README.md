# Sphingolipid Analysis Pipeline

A comprehensive, user-friendly pipeline for analyzing LC-MS sphingolipid data with automated statistical testing and publication-quality visualizations.

## Features

- **Automatic data structure detection** - Recognizes 31 sphingolipid species from your LC-MS output
- **Smart statistical analysis** - Automatically tests assumptions (normality, homoscedasticity) and selects appropriate tests
- **Publication-quality figures** - Box plots, stacked bars, heatmaps, and more
- **Comprehensive calculations** - Group totals, percentages, clinical ratios
- **Export everything** - CSV files, figures in PNG/PDF, complete ZIP packages

## Installation

```bash
# Clone or download the pipeline
cd sphingolipid-analysis-pipeline

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Usage

### Web Application
[https://sphingolipid-pipeline.streamlit.app/](https://sphingolipid-analysis-pipeline-ej8tup4wgi3arlzba7uhmb.streamlit.app/)

### For Local Hosting
```bash
streamlit run app.py
```

Then open http://localhost:8501 in your browser.


## Data Format

Expected LC-MS sheet input format:
- **Rows**: First few are Std curves then samples
- **Columns**: "Data Filename" (Std curve and sample names), followed by sphingolipid species

Expected Sample Sheet input format:
- **Rows**: Samples
- **Columns**: Sphingolipid species (matching panel names)
- **First column(s)**: Sample ID, Group/Type
- **Values**: Concentrations (typically nmol/L)
- **Below LOD**: Can be "-----", "LOD", "BLQ", "ND", etc.

Example:
| Type | Sample_ID | C16 Cer | C24-0 Cer | C16-SM | S-d18-1 | ... |
|------|-----------|-----|-----|-------|-------|-----|
| Aged | 81-0210 | 37.62 | 194.1 | 231.9 | 3220.44 | ... |
| Aged | 81-0211 | 66.66 | 287.58 | 158.46 | 1534.02 | ... |
| Young | 60-677 | 45.2 | ----- | 189.3 | 2890.1 | ... |

## Recognized Sphingolipid Species

### Shingoid Bases (Free bases)
- S-d18-1 (Shingosine)
- S-d18-0 (Dihydrosphingosine (Sphinganine))
- 3KDHS (3-Keto-dihydrosphingosine)

### Sphingoid Base Phosphates
- S1P-d18-1 (Sphingosine-1-Phosphate (d18:1))
- S1P-d18-0 (Sphinganine-1-Phosphate (d18:0))

### Ceramides (Cer)
- C12 Cer (N-lauroyl-D-erythro-sphingosine)
- C14 Cer (N-myristoyl-D-erythro-sphingosine)
- C16 Cer (N-palmitoyl-D-erythro-sphingosine)
- C18-0 Cer (N-stearoyl-D-erythro-sphingosine)
- And many more....

### Ceramide-1-Phosphates
- C16-CP, C18-CP, C24-CP

### Dihydroceramides (DHC)
- C16-DHC, C18DHC, C18-1DHC, C24DHC, C24-1DHC

###Sphingomyelins
- C16-SM, C17-SM, C18-SM, ....

See `config/sphingolipid_species.py` for the complete panel definition.

## Adding New Bile Acids

Edit `config/sphingolipid_species.py`:

```python
SPHINGOLIPID_PANEL["NEW_S"] = SphingolipidSpecies(
        abbreviation="New S",
        full_name="New Sphingolipid",
        sphingo_class=SphingoClass.SPHINGOID_BASE,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.NA,
        saturation=Saturation.NA,
        molecular_weight=299.49,
        notes="Primary sphingoid base, unsaturated",
    )
```

## Statistical Tests

The pipeline automatically selects tests based on data characteristics:

| Condition | 2 Groups | >2 Groups | >=2 Independent Variables |
|-----------|----------|-----------|-----------|
| Normal + Equal variance | Independent t-test | One-way ANOVA | Two-way ANOVA |
| Normal + Unequal variance | Welch's t-test | Welch's ANOVA | ART ANOVA |
| Non-normal | Mann-Whitney U | Kruskal-Wallis | ART ANOVA | ART ANOVA |

Post-hoc tests (for >2 groups):
- Tukey HSD (parametric, equal variance)
- Games-Howell (parametric, unequal variance)
- Dunn's test (non-parametric)

## Calculated Metrics

### Group Totals
- `total_all`: Sum of all sphingolipids
- `total_ceramide`: Sum of ceramides
- `total_dihydroceramide`: Sum of dihydroceramides
- `total_sphingomyelin`: Sum of all sphingomyelins
- `total_sphingoid_base`: Sum of sphingoid bases
- `total_sphingoid_base_phosphate`: Sum of sphingoid base phosphates
- `total_ceramide_1_phosphate`: Sum of ceramide 1 phosphates

### Clinical Ratios
- chain length comparisons
- total ceramides vs total sphingomyelin
- total ceramides vs total dihydroceramides
- saturated vs unsaturated
- Individual sphingolpiid species ratios comparisons
- And more...

## Project Structure

```
sphingolipid-analysis-pipeline/
├── app.py                  # Streamlit web application
├── requirements.txt        # Python dependencies
├── README.md              # This file
├── config/
│   ├── __init__.py
│   └── sphingolioid_species.py  # Sphingolipid definitions (EDIT THIS FOR NEW BAs)
└── modules/
    ├── __init__.py
    ├── data_processing.py    # Data loading and calculations
    ├── report_generation.py  # Statistical report output
    ├── statistical_tests.py  # Automated statistical analysis
    └── visualization.py      # Publication-quality figures
```

## License

MIT License - Feel free to use and modify for research purposes.

