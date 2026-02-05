"""
Sphingolipid Species Configuration
==================================

This module defines all sphingolipid species and their classifications.
Update this file when new sphingolipids are added to the LC-MS panel.

Classification Schema:
- Class: ceramide, dihydroceramide, sphingomyelin, sphingoid_base, 
         hexosylceramide, ceramide_1_phosphate
- Sphingoid Base: d18:1 (sphingosine), d18:0 (sphinganine/dihydro)
- Chain Length: short (C12-C14), medium (C16-C18), long (C20-C22), 
                very_long (C24+)
- Saturation: saturated (X:0), monounsaturated (X:1)

Last updated: 2025-02-05
Panel source: Zhou Lab LC-MS/MS Ceramide Panel
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Set
from enum import Enum


class SphingoClass(Enum):
    """Major sphingolipid class."""
    CERAMIDE = "ceramide"
    DIHYDROCERAMIDE = "dihydroceramide"
    SPHINGOMYELIN = "sphingomyelin"
    SPHINGOID_BASE = "sphingoid_base"
    SPHINGOID_BASE_PHOSPHATE = "sphingoid_base_phosphate"
    HEXOSYLCERAMIDE = "hexosylceramide"
    CERAMIDE_1_PHOSPHATE = "ceramide_1_phosphate"


class SphingoBase(Enum):
    """Sphingoid base backbone type."""
    D18_1 = "d18:1"    # Sphingosine (unsaturated)
    D18_0 = "d18:0"    # Sphinganine/Dihydrosphingosine (saturated)
    KETO = "3-keto"    # 3-keto derivative
    NA = "not_applicable"


class ChainLength(Enum):
    """Fatty acid chain length category."""
    SHORT = "short"        # C12-C14
    MEDIUM = "medium"      # C16-C18
    LONG = "long"          # C20-C22
    VERY_LONG = "very_long"  # C24+
    NA = "not_applicable"  # For sphingoid bases without FA


class Saturation(Enum):
    """Fatty acid saturation status."""
    SATURATED = "saturated"          # X:0
    MONOUNSATURATED = "monounsaturated"  # X:1
    NA = "not_applicable"


@dataclass
class SphingolipidSpecies:
    """Definition of a single sphingolipid species."""
    abbreviation: str
    full_name: str
    sphingo_class: SphingoClass
    sphingo_base: SphingoBase
    chain_length: ChainLength
    saturation: Saturation
    carbon_number: Optional[int] = None  # FA carbon number (e.g., 16, 18, 24)
    molecular_weight: Optional[float] = None
    notes: str = ""


# =============================================================================
# SPHINGOLIPID PANEL DEFINITION
# =============================================================================
# Keys must match column names in LC-MS output exactly.

SPHINGOLIPID_PANEL: Dict[str, SphingolipidSpecies] = {
    
    # -------------------------------------------------------------------------
    # SPHINGOID BASES (Free bases)
    # -------------------------------------------------------------------------
    "S-d18-1": SphingolipidSpecies(
        abbreviation="S-d18-1",
        full_name="Sphingosine",
        sphingo_class=SphingoClass.SPHINGOID_BASE,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.NA,
        saturation=Saturation.NA,
        molecular_weight=299.49,
        notes="Primary sphingoid base, unsaturated",
    ),
    "S-d18-0": SphingolipidSpecies(
        abbreviation="S-d18-0",
        full_name="Dihydrosphingosine (Sphinganine)",
        sphingo_class=SphingoClass.SPHINGOID_BASE,
        sphingo_base=SphingoBase.D18_0,
        chain_length=ChainLength.NA,
        saturation=Saturation.NA,
        molecular_weight=301.50,
        notes="Saturated sphingoid base, precursor in de novo synthesis",
    ),
    "3KDHS": SphingolipidSpecies(
        abbreviation="3KDHS",
        full_name="3-Keto-dihydrosphingosine",
        sphingo_class=SphingoClass.SPHINGOID_BASE,
        sphingo_base=SphingoBase.KETO,
        chain_length=ChainLength.NA,
        saturation=Saturation.NA,
        molecular_weight=299.49,
        notes="First intermediate in de novo sphingolipid synthesis",
    ),
    
    # -------------------------------------------------------------------------
    # SPHINGOID BASE PHOSPHATES (Signaling molecules)
    # -------------------------------------------------------------------------
    "S1P-d18-1": SphingolipidSpecies(
        abbreviation="S1P-d18-1",
        full_name="Sphingosine-1-Phosphate (d18:1)",
        sphingo_class=SphingoClass.SPHINGOID_BASE_PHOSPHATE,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.NA,
        saturation=Saturation.NA,
        molecular_weight=379.50,
        notes="Bioactive signaling lipid, pro-survival",
    ),
    "S1P-d18-0": SphingolipidSpecies(
        abbreviation="S1P-d18-0",
        full_name="Sphinganine-1-Phosphate (d18:0)",
        sphingo_class=SphingoClass.SPHINGOID_BASE_PHOSPHATE,
        sphingo_base=SphingoBase.D18_0,
        chain_length=ChainLength.NA,
        saturation=Saturation.NA,
        molecular_weight=381.50,
        notes="Dihydro form of S1P",
    ),
    
    # -------------------------------------------------------------------------
    # CERAMIDES (Cer) - d18:1 sphingoid base
    # -------------------------------------------------------------------------
    "C12 Cer": SphingolipidSpecies(
        abbreviation="C12 Cer",
        full_name="N-lauroyl-D-erythro-sphingosine",
        sphingo_class=SphingoClass.CERAMIDE,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.SHORT,
        saturation=Saturation.SATURATED,
        carbon_number=12,
        molecular_weight=481.80,
    ),
    "C14 Cer": SphingolipidSpecies(
        abbreviation="C14 Cer",
        full_name="N-myristoyl-D-erythro-sphingosine",
        sphingo_class=SphingoClass.CERAMIDE,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.SHORT,
        saturation=Saturation.SATURATED,
        carbon_number=14,
        molecular_weight=509.85,
    ),
    "C16 Cer": SphingolipidSpecies(
        abbreviation="C16 Cer",
        full_name="N-palmitoyl-D-erythro-sphingosine",
        sphingo_class=SphingoClass.CERAMIDE,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.MEDIUM,
        saturation=Saturation.SATURATED,
        carbon_number=16,
        molecular_weight=535.88,
        notes="Most abundant ceramide in many tissues",
    ),
    "C18-0 Cer": SphingolipidSpecies(
        abbreviation="C18-0 Cer",
        full_name="N-stearoyl-D-erythro-sphingosine",
        sphingo_class=SphingoClass.CERAMIDE,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.MEDIUM,
        saturation=Saturation.SATURATED,
        carbon_number=18,
        molecular_weight=565.95,
    ),
    "C18-1 Cer": SphingolipidSpecies(
        abbreviation="C18-1 Cer",
        full_name="N-oleoyl-D-erythro-sphingosine",
        sphingo_class=SphingoClass.CERAMIDE,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.MEDIUM,
        saturation=Saturation.MONOUNSATURATED,
        carbon_number=18,
        molecular_weight=563.94,
    ),
    "C20 Cer": SphingolipidSpecies(
        abbreviation="C20 Cer",
        full_name="N-arachidoyl-D-erythro-sphingosine",
        sphingo_class=SphingoClass.CERAMIDE,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.LONG,
        saturation=Saturation.SATURATED,
        carbon_number=20,
        molecular_weight=594.01,
    ),
    "C22 Cer": SphingolipidSpecies(
        abbreviation="C22 Cer",
        full_name="N-behenoyl-D-erythro-sphingosine",
        sphingo_class=SphingoClass.CERAMIDE,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.LONG,
        saturation=Saturation.SATURATED,
        carbon_number=22,
        molecular_weight=622.06,
    ),
    "C24-0 Cer": SphingolipidSpecies(
        abbreviation="C24-0 Cer",
        full_name="N-lignoceroyl-D-erythro-sphingosine",
        sphingo_class=SphingoClass.CERAMIDE,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.VERY_LONG,
        saturation=Saturation.SATURATED,
        carbon_number=24,
        molecular_weight=648.10,
        notes="Very long chain ceramide",
    ),
    "C24-1 Cer": SphingolipidSpecies(
        abbreviation="C24-1 Cer",
        full_name="N-nervonoyl-D-erythro-sphingosine",
        sphingo_class=SphingoClass.CERAMIDE,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.VERY_LONG,
        saturation=Saturation.MONOUNSATURATED,
        carbon_number=24,
        molecular_weight=648.10,
        notes="Nervonic acid ceramide, abundant in nervous tissue",
    ),
    "C26 Cer": SphingolipidSpecies(
        abbreviation="C26 Cer",
        full_name="C26-Ceramide",
        sphingo_class=SphingoClass.CERAMIDE,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.VERY_LONG,
        saturation=Saturation.SATURATED,
        carbon_number=26,
        molecular_weight=678.17,
    ),
    
    # -------------------------------------------------------------------------
    # DIHYDROCERAMIDES (DHC) - d18:0 sphingoid base (saturated backbone)
    # -------------------------------------------------------------------------
    "C16-DHC": SphingolipidSpecies(
        abbreviation="C16-DHC",
        full_name="C16 Dihydroceramide",
        sphingo_class=SphingoClass.DIHYDROCERAMIDE,
        sphingo_base=SphingoBase.D18_0,
        chain_length=ChainLength.MEDIUM,
        saturation=Saturation.SATURATED,
        carbon_number=16,
        molecular_weight=539.917,
        notes="Precursor to C16 Cer in de novo synthesis",
    ),
    "C18DHC": SphingolipidSpecies(
        abbreviation="C18DHC",
        full_name="C18 Dihydroceramide",
        sphingo_class=SphingoClass.DIHYDROCERAMIDE,
        sphingo_base=SphingoBase.D18_0,
        chain_length=ChainLength.MEDIUM,
        saturation=Saturation.SATURATED,
        carbon_number=18,
        molecular_weight=567.97,
    ),
    "C18-1DHC": SphingolipidSpecies(
        abbreviation="C18-1DHC",
        full_name="C18:1 Dihydroceramide",
        sphingo_class=SphingoClass.DIHYDROCERAMIDE,
        sphingo_base=SphingoBase.D18_0,
        chain_length=ChainLength.MEDIUM,
        saturation=Saturation.MONOUNSATURATED,
        carbon_number=18,
        molecular_weight=565.954,
    ),
    "C24DHC": SphingolipidSpecies(
        abbreviation="C24DHC",
        full_name="C24 Dihydroceramide",
        sphingo_class=SphingoClass.DIHYDROCERAMIDE,
        sphingo_base=SphingoBase.D18_0,
        chain_length=ChainLength.VERY_LONG,
        saturation=Saturation.SATURATED,
        carbon_number=24,
        molecular_weight=749.14,
    ),
    "C24-1DHC": SphingolipidSpecies(
        abbreviation="C24-1DHC",
        full_name="C24:1 Dihydroceramide",
        sphingo_class=SphingoClass.DIHYDROCERAMIDE,
        sphingo_base=SphingoBase.D18_0,
        chain_length=ChainLength.VERY_LONG,
        saturation=Saturation.MONOUNSATURATED,
        carbon_number=24,
        molecular_weight=650.113,
    ),
    
    # -------------------------------------------------------------------------
    # SPHINGOMYELINS (SM) - Phosphocholine headgroup
    # -------------------------------------------------------------------------
    "C16-SM": SphingolipidSpecies(
        abbreviation="C16-SM",
        full_name="C16 Sphingomyelin",
        sphingo_class=SphingoClass.SPHINGOMYELIN,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.MEDIUM,
        saturation=Saturation.SATURATED,
        carbon_number=16,
        molecular_weight=703.00,
        notes="Major plasma sphingomyelin",
    ),
    "C17-SM": SphingolipidSpecies(
        abbreviation="C17-SM",
        full_name="C17 Sphingomyelin",
        sphingo_class=SphingoClass.SPHINGOMYELIN,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.MEDIUM,
        saturation=Saturation.SATURATED,
        carbon_number=17,
        molecular_weight=717.10,
        notes="Odd-chain sphingomyelin",
    ),
    "C18-SM": SphingolipidSpecies(
        abbreviation="C18-SM",
        full_name="C18 Sphingomyelin",
        sphingo_class=SphingoClass.SPHINGOMYELIN,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.MEDIUM,
        saturation=Saturation.SATURATED,
        carbon_number=18,
        molecular_weight=731.10,
    ),
    "C20-SM": SphingolipidSpecies(
        abbreviation="C20-SM",
        full_name="C20 Sphingomyelin",
        sphingo_class=SphingoClass.SPHINGOMYELIN,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.LONG,
        saturation=Saturation.SATURATED,
        carbon_number=20,
        molecular_weight=759.20,
    ),
    "C22-SM": SphingolipidSpecies(
        abbreviation="C22-SM",
        full_name="C22 Sphingomyelin",
        sphingo_class=SphingoClass.SPHINGOMYELIN,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.LONG,
        saturation=Saturation.SATURATED,
        carbon_number=22,
        molecular_weight=787.20,
    ),
    "C26-SM": SphingolipidSpecies(
        abbreviation="C26-SM",
        full_name="C26 Sphingomyelin",
        sphingo_class=SphingoClass.SPHINGOMYELIN,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.VERY_LONG,
        saturation=Saturation.SATURATED,
        carbon_number=26,
        molecular_weight=843.30,
    ),
    "C261-SM": SphingolipidSpecies(
        abbreviation="C261-SM",
        full_name="C26:1 Sphingomyelin",
        sphingo_class=SphingoClass.SPHINGOMYELIN,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.VERY_LONG,
        saturation=Saturation.MONOUNSATURATED,
        carbon_number=26,
        molecular_weight=841.30,
    ),
    
    # -------------------------------------------------------------------------
    # HEXOSYLCERAMIDES (HexCer) - Glucose/Galactose headgroup
    # -------------------------------------------------------------------------
    "GLU-C": SphingolipidSpecies(
        abbreviation="GLU-C",
        full_name="Glucosyl Ceramide",
        sphingo_class=SphingoClass.HEXOSYLCERAMIDE,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.MEDIUM,  # Typically measured as total
        saturation=Saturation.NA,
        molecular_weight=714.025,
        notes="Glucosylceramide, precursor to complex glycosphingolipids",
    ),
    
    # -------------------------------------------------------------------------
    # CERAMIDE-1-PHOSPHATES (C1P) - Signaling/Inflammatory
    # -------------------------------------------------------------------------
    "C16-CP": SphingolipidSpecies(
        abbreviation="C16-CP",
        full_name="C16 Ceramide-1-Phosphate",
        sphingo_class=SphingoClass.CERAMIDE_1_PHOSPHATE,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.MEDIUM,
        saturation=Saturation.SATURATED,
        carbon_number=16,
        molecular_weight=617.911,
        notes="Bioactive lipid, involved in inflammation",
    ),
    "C18-CP": SphingolipidSpecies(
        abbreviation="C18-CP",
        full_name="C18:1 Ceramide-1-Phosphate",
        sphingo_class=SphingoClass.CERAMIDE_1_PHOSPHATE,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.MEDIUM,
        saturation=Saturation.MONOUNSATURATED,
        carbon_number=18,
        molecular_weight=643.90,
    ),
    "C24-CP": SphingolipidSpecies(
        abbreviation="C24-CP",
        full_name="C24 Ceramide-1-Phosphate",
        sphingo_class=SphingoClass.CERAMIDE_1_PHOSPHATE,
        sphingo_base=SphingoBase.D18_1,
        chain_length=ChainLength.VERY_LONG,
        saturation=Saturation.SATURATED,
        carbon_number=24,
        molecular_weight=730.124,
    ),
}


# =============================================================================
# HELPER FUNCTIONS FOR ACCESSING CLASSIFICATIONS
# =============================================================================

def get_species_by_class(sphingo_class: SphingoClass) -> List[str]:
    """Get all sphingolipid abbreviations with a specific class."""
    return [abbr for abbr, spec in SPHINGOLIPID_PANEL.items() 
            if spec.sphingo_class == sphingo_class]


def get_species_by_base(sphingo_base: SphingoBase) -> List[str]:
    """Get all sphingolipid abbreviations with a specific sphingoid base."""
    return [abbr for abbr, spec in SPHINGOLIPID_PANEL.items() 
            if spec.sphingo_base == sphingo_base]


def get_species_by_chain_length(chain_length: ChainLength) -> List[str]:
    """Get all sphingolipid abbreviations with a specific chain length category."""
    return [abbr for abbr, spec in SPHINGOLIPID_PANEL.items() 
            if spec.chain_length == chain_length]


def get_species_by_saturation(saturation: Saturation) -> List[str]:
    """Get all sphingolipid abbreviations with a specific saturation status."""
    return [abbr for abbr, spec in SPHINGOLIPID_PANEL.items() 
            if spec.saturation == saturation]


def get_ceramides() -> List[str]:
    """Get all ceramide species."""
    return get_species_by_class(SphingoClass.CERAMIDE)


def get_dihydroceramides() -> List[str]:
    """Get all dihydroceramide species."""
    return get_species_by_class(SphingoClass.DIHYDROCERAMIDE)


def get_sphingomyelins() -> List[str]:
    """Get all sphingomyelin species."""
    return get_species_by_class(SphingoClass.SPHINGOMYELIN)


def get_sphingoid_bases() -> List[str]:
    """Get all free sphingoid base species."""
    return get_species_by_class(SphingoClass.SPHINGOID_BASE)


def get_sphingoid_base_phosphates() -> List[str]:
    """Get all sphingoid base phosphate species (S1P, etc.)."""
    return get_species_by_class(SphingoClass.SPHINGOID_BASE_PHOSPHATE)


def get_hexosylceramides() -> List[str]:
    """Get all hexosylceramide species."""
    return get_species_by_class(SphingoClass.HEXOSYLCERAMIDE)


def get_ceramide_1_phosphates() -> List[str]:
    """Get all ceramide-1-phosphate species."""
    return get_species_by_class(SphingoClass.CERAMIDE_1_PHOSPHATE)


def get_saturated() -> List[str]:
    """Get all species with saturated fatty acids."""
    return get_species_by_saturation(Saturation.SATURATED)


def get_unsaturated() -> List[str]:
    """Get all species with unsaturated fatty acids."""
    return get_species_by_saturation(Saturation.MONOUNSATURATED)


def get_very_long_chain() -> List[str]:
    """Get all very long chain (C24+) species."""
    return get_species_by_chain_length(ChainLength.VERY_LONG)


def get_long_chain() -> List[str]:
    """Get all long chain (C20-C22) species."""
    return get_species_by_chain_length(ChainLength.LONG)


def get_medium_chain() -> List[str]:
    """Get all medium chain (C16-C18) species."""
    return get_species_by_chain_length(ChainLength.MEDIUM)


def get_short_chain() -> List[str]:
    """Get all short chain (C12-C14) species."""
    return get_species_by_chain_length(ChainLength.SHORT)


def get_all_species() -> List[str]:
    """Get list of all sphingolipid abbreviations in the panel."""
    return list(SPHINGOLIPID_PANEL.keys())


def get_species_info(abbreviation: str) -> Optional[SphingolipidSpecies]:
    """Get full information for a specific sphingolipid."""
    return SPHINGOLIPID_PANEL.get(abbreviation)


def validate_columns(columns: List[str]) -> Dict[str, List[str]]:
    """
    Validate column names against the sphingolipid panel.
    
    Returns dict with 'matched', 'unmatched', and 'missing' lists.
    """
    panel_species = set(SPHINGOLIPID_PANEL.keys())
    column_set = set(columns)
    
    matched = list(panel_species & column_set)
    unmatched = list(column_set - panel_species)
    missing = list(panel_species - column_set)
    
    return {
        "matched": matched,
        "unmatched": unmatched,  # columns not in panel (metadata or unknown species)
        "missing": missing,      # panel species not in columns
    }


# =============================================================================
# PREDEFINED ANALYSIS GROUPS
# =============================================================================
# These define common groupings used in sphingolipid analysis

ANALYSIS_GROUPS = {
    # By class
    "total_ceramides": get_ceramides(),
    "total_dihydroceramides": get_dihydroceramides(),
    "total_sphingomyelins": get_sphingomyelins(),
    "total_sphingoid_bases": get_sphingoid_bases(),
    "total_sphingoid_base_phosphates": get_sphingoid_base_phosphates(),
    "total_hexosylceramides": get_hexosylceramides(),
    "total_ceramide_1_phosphates": get_ceramide_1_phosphates(),
    
    # By chain length
    "short_chain": get_short_chain(),
    "medium_chain": get_medium_chain(),
    "long_chain": get_long_chain(),
    "very_long_chain": get_very_long_chain(),
    
    # By saturation
    "saturated": get_saturated(),
    "unsaturated": get_unsaturated(),
    
    # By sphingoid base
    "d18_1_base": get_species_by_base(SphingoBase.D18_1),
    "d18_0_base": get_species_by_base(SphingoBase.D18_0),
    
    # Combined groups
    "all_cer_and_dhc": get_ceramides() + get_dihydroceramides(),
    "signaling_lipids": get_sphingoid_base_phosphates() + get_ceramide_1_phosphates(),
}


# =============================================================================
# CLINICAL/RESEARCH RATIO DEFINITIONS
# =============================================================================

CLINICAL_RATIOS = {
    # Ceramide ratios
    "C16_to_C24_Cer": {
        "numerator": ["C16 Cer"], 
        "denominator": ["C24-0 Cer"]
    },
    "C24_1_to_C24_0_Cer": {
        "numerator": ["C24-1 Cer"], 
        "denominator": ["C24-0 Cer"],
        "notes": "Desaturation index"
    },
    "C18_to_C16_Cer": {
        "numerator": ["C18-0 Cer"], 
        "denominator": ["C16 Cer"]
    },
    
    # Ceramide to Dihydroceramide (de novo synthesis activity)
    "Cer_to_DHC_C16": {
        "numerator": ["C16 Cer"], 
        "denominator": ["C16-DHC"],
        "notes": "Dihydroceramide desaturase activity marker"
    },
    "Cer_to_DHC_C24": {
        "numerator": ["C24-0 Cer"], 
        "denominator": ["C24DHC"]
    },
    
    # Sphingomyelin to Ceramide (sphingomyelinase activity)
    "SM_to_Cer_C16": {
        "numerator": ["C16-SM"], 
        "denominator": ["C16 Cer"],
        "notes": "Sphingomyelinase activity marker"
    },
    
    # Sphingosine-1-phosphate ratios
    "S1P_to_Sph": {
        "numerator": ["S1P-d18-1"], 
        "denominator": ["S-d18-1"],
        "notes": "Sphingosine kinase activity"
    },
    "S1P_d18_1_to_d18_0": {
        "numerator": ["S1P-d18-1"], 
        "denominator": ["S1P-d18-0"]
    },
    
    # Chain length ratios
    "very_long_to_long_Cer": {
        "numerator": get_very_long_chain(),
        "denominator": get_long_chain() + get_medium_chain()
    },
    
    # Saturation ratios
    "saturated_to_unsaturated": {
        "numerator": get_saturated(),
        "denominator": get_unsaturated()
    },
    
    # Total class ratios
    "total_Cer_to_SM": {
        "numerator": get_ceramides(),
        "denominator": get_sphingomyelins()
    },
    "total_Cer_to_DHC": {
        "numerator": get_ceramides(),
        "denominator": get_dihydroceramides()
    },
}


if __name__ == "__main__":
    # Print summary of the panel
    print("=" * 60)
    print("SPHINGOLIPID PANEL SUMMARY")
    print("=" * 60)
    print(f"Total species in panel: {len(SPHINGOLIPID_PANEL)}")
    print(f"\nBy class:")
    print(f"  Ceramides: {len(get_ceramides())}")
    print(f"  Dihydroceramides: {len(get_dihydroceramides())}")
    print(f"  Sphingomyelins: {len(get_sphingomyelins())}")
    print(f"  Sphingoid bases: {len(get_sphingoid_bases())}")
    print(f"  Sphingoid base phosphates: {len(get_sphingoid_base_phosphates())}")
    print(f"  Hexosylceramides: {len(get_hexosylceramides())}")
    print(f"  Ceramide-1-phosphates: {len(get_ceramide_1_phosphates())}")
    print(f"\nBy chain length:")
    print(f"  Short (C12-C14): {len(get_short_chain())}")
    print(f"  Medium (C16-C18): {len(get_medium_chain())}")
    print(f"  Long (C20-C22): {len(get_long_chain())}")
    print(f"  Very long (C24+): {len(get_very_long_chain())}")
    print(f"\nBy saturation:")
    print(f"  Saturated: {len(get_saturated())}")
    print(f"  Unsaturated: {len(get_unsaturated())}")
