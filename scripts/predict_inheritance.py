"""
AD/AR Inheritance Model — Inference Script
===========================================
Loads a pre-trained model and runs predictions. No training code.
Designed for production use with minimal dependencies.

Pipeline:
    1. Build feature vector directly from patient's HPO terms + onset
    2. Add gnomAD gene constraints
    3. Run through trained RF model
    4. Apply confidence threshold -> AD / AR / Uncertain

Usage:
    # As a module
    from predict_inheritance import load_predictor

    predictor = load_predictor('ad_ar_model.joblib')
    result = predictor.predict(
        patient_hpo_ids=['HP:0001250', 'HP:0001263'],
        gene_symbol='SCN1A',
        patient_demographics={'age': '7 years old'},
        onsets={'HP:0001250': '9 months old'},
    )
    print(result['prediction'])  # 'AD', 'AR', or 'Uncertain'

    # CLI
    python predict_inheritance.py \
        --model ad_ar_model.joblib \
        --gene SCN1A \
        --hpo HP:0001250,HP:0001263,HP:0001252 \
        --onsets HP:0001250=9m \
        --age "7 years old"

Requirements:
    pip install pyhpo pandas numpy scikit-learn joblib
"""

import re
import warnings
import numpy as np
import pandas as pd
import joblib
from pyhpo import Ontology

warnings.filterwarnings("ignore")

# ============================================================
# HPO Ontology Initialization (one-time cost)
# ============================================================
_ontology_initialized = False


def _ensure_ontology():
    global _ontology_initialized
    if not _ontology_initialized:
        _ = Ontology()
        _ontology_initialized = True


def _get_roots():
    _ensure_ontology()
    return (
        Ontology.get_hpo_object('HP:0000118'),  # phenotype root
        Ontology.get_hpo_object('HP:0000005'),  # inheritance root
    )


_TIER_IDS = None


def _get_tier_ids():
    global _TIER_IDS
    if _TIER_IDS is None:
        phenotype_root, _ = _get_roots()
        _TIER_IDS = {}
        current_level = {phenotype_root}
        for depth in range(1, 7):
            next_level = set()
            for parent in current_level:
                for child in parent.children:
                    next_level.add(child)
            _TIER_IDS[depth] = {t.id for t in next_level}
            current_level = next_level
            if not current_level:
                break
    return _TIER_IDS


# ============================================================
# Onset Parsing
# ============================================================

ONSET_KEYWORDS = {
    'Adult':      ['adult', 'late', 'elderly', 'middle age'],
    'Pediatric':  ['childhood', 'juvenile', 'infantile', 'toddler'],
    'Congenital': ['congenital', 'neonatal', 'birth', 'prenatal', 'antenatal', 'fetus']
}


def _parse_age_to_years(age_str):
    if not age_str or not isinstance(age_str, str):
        return None
    text = age_str.lower().strip()
    if 'unknown' in text:
        return None
    if any(k in text for k in ['birth', 'neonatal', 'newborn']):
        return 0.0
    if any(k in text for k in ['prenatal', 'antenatal', 'fetus', 'fetal', 'gestation']):
        return -0.1
    if 'congenital' in text:
        return 0.0
    patterns = [
        (r'(\d+\.?\d*)\s*(?:year|yr|y)', 'years'),
        (r'(\d+\.?\d*)\s*(?:month|mo|m)(?:s|\b)', 'months'),
        (r'(\d+\.?\d*)\s*(?:week|wk|w)', 'weeks'),
        (r'(\d+\.?\d*)\s*(?:day|d)', 'days'),
    ]
    for pattern, unit in patterns:
        match = re.search(pattern, text)
        if match:
            value = float(match.group(1))
            if unit == 'years': return value
            elif unit == 'months': return value / 12.0
            elif unit == 'weeks': return value / 52.0
            elif unit == 'days': return value / 365.0
    match = re.search(r'^(\d+\.?\d*)\s*$', text)
    if match:
        return float(match.group(1))
    return None


def _age_to_onset_category(age_years):
    if age_years is None:
        return "Unknown"
    if age_years < 0.08:
        return "Congenital"
    elif age_years < 18.0:
        return "Pediatric"
    else:
        return "Adult"


def _get_onset_category(raw_onset):
    if not raw_onset or not isinstance(raw_onset, str):
        return "Unknown"
    text = raw_onset.lower().strip()
    if "unknown" in text:
        return "Unknown"
    for category, keywords in ONSET_KEYWORDS.items():
        if any(k in text for k in keywords):
            return category
    age_years = _parse_age_to_years(raw_onset)
    if age_years is not None:
        return _age_to_onset_category(age_years)
    return "Unknown"


def _get_patient_onset_fallback(demographics):
    if not demographics or not isinstance(demographics, dict):
        return "Unknown"
    age_years = _parse_age_to_years(demographics.get('age'))
    if age_years is not None:
        return _age_to_onset_category(age_years)
    return "Unknown"


# ============================================================
# HPO Term Resolution & Bucketing
# ============================================================

def _resolve_hpo_term(h_id):
    _ensure_ontology()
    try:
        if hasattr(h_id, 'all_parents'):
            return h_id
        elif isinstance(h_id, (str, int)):
            return Ontology.get_hpo_object(h_id)
    except Exception:
        pass
    return None


def _get_tier_buckets(term, tier=2):
    _, inheritance_root = _get_roots()
    tier_ids = _get_tier_ids()
    lineage_ids = {p.id for p in term.all_parents} | {term.id}
    if inheritance_root.id in lineage_ids:
        return []
    if tier == 0:
        return [term.name]
    for t in range(tier, 0, -1):
        if t in tier_ids:
            matches = sorted(lineage_ids.intersection(tier_ids[t]))
            if matches:
                return [Ontology.get_hpo_object(mid).name for mid in matches]
    return []


def _extract_features(phenotype_records, tiers=(1, 3, 5),
                      use_onset=True, patient_demographics=None):
    """
    Extract multi-tier binary features from phenotype records.
    Returns dict of {feature_name: 1}.
    """
    if not isinstance(phenotype_records, list):
        return {}

    patient_onset = "Unknown"
    if patient_demographics:
        patient_onset = _get_patient_onset_fallback(patient_demographics)

    feature_set = set()

    for record in phenotype_records:
        if isinstance(record, dict):
            h_id = record.get('id') or record.get('HPO_ID')
            raw_onset = record.get('onset')
        else:
            h_id = record
            raw_onset = None

        term = _resolve_hpo_term(h_id)
        if term is None:
            continue

        onset_cat = "Unknown"
        if use_onset:
            onset_cat = _get_onset_category(raw_onset)
            if onset_cat == "Unknown" and patient_onset != "Unknown":
                onset_cat = patient_onset

        for t in tiers:
            bucket_names = _get_tier_buckets(term, tier=t)
            prefix = f"T{t}"
            for name in bucket_names:
                feature_set.add(f"{prefix}::{name}")
                if use_onset and onset_cat != "Unknown":
                    feature_set.add(f"{prefix}::{name}::Onset::{onset_cat}")

    return {k: 1 for k in feature_set}


# ============================================================
# Predictor Class
# ============================================================

class InheritancePredictor:
    """
    Production predictor for AD/AR inheritance classification.

    Builds features directly from patient's HPO terms (no disease
    matching step). Uses gnomAD gene constraints and a confidence
    threshold for three-tier output: AD / AR / Uncertain.
    """

    DEFAULT_CONFIDENCE_THRESHOLD = 0.20

    def __init__(self, model, gene_disease_index, feature_names,
                 training_tiers, training_use_onset,
                 feature_columns_kept, gnomad_lookup=None,
                 confidence_threshold=None):
        self.model = model
        self.index = gene_disease_index
        self.feature_names = feature_names
        self.tiers = training_tiers
        self.use_onset = training_use_onset
        self.feature_columns_kept = feature_columns_kept
        self.gnomad_lookup = gnomad_lookup or {}
        self.confidence_threshold = (
            confidence_threshold or self.DEFAULT_CONFIDENCE_THRESHOLD
        )

    def _classify(self, prob_ad):
        confidence = abs(prob_ad - 0.5)
        if confidence < self.confidence_threshold:
            return 'Uncertain', confidence
        elif prob_ad >= 0.5:
            return 'AD', confidence
        else:
            return 'AR', confidence

    def _build_feature_vector(self, phenotype_records, gene_symbol,
                              patient_demographics=None):
        """Build a feature vector matching the training encoding."""
        feature_dict = {col: 0.0 for col in self.feature_names}

        feat_values = _extract_features(
            phenotype_records,
            tiers=self.tiers,
            use_onset=self.use_onset,
            patient_demographics=patient_demographics,
        )

        for feat_name, val in feat_values.items():
            col = f"bucket_{feat_name}"
            if col in feature_dict:
                feature_dict[col] = val

        if 'n_hpo_terms' in feature_dict:
            feature_dict['n_hpo_terms'] = len(phenotype_records)

        # gene_role — default to unknown
        role_col = 'gene_role_unknown'
        if role_col in feature_dict:
            feature_dict[role_col] = 1

        # gnomAD constraints
        if gene_symbol in self.gnomad_lookup:
            for col, val in self.gnomad_lookup[gene_symbol].items():
                if col in feature_dict:
                    feature_dict[col] = val

        return np.array([[feature_dict[col] for col in self.feature_names]])

    def _get_gene_metadata(self, gene_symbol):
        """Get known inheritance info from the gene-disease index."""
        entries = self.index.get(gene_symbol, [])
        if not entries:
            return {
                'n_diseases': 0,
                'known_inheritance': 'unknown',
                'is_pleiotropic': False,
            }

        has_ad = any(e['is_AD'] for e in entries)
        has_ar = any(e['is_AR'] for e in entries)

        if has_ad and has_ar:
            known = 'AD+AR'
        elif has_ad:
            known = 'AD'
        elif has_ar:
            known = 'AR'
        else:
            known = 'unknown'

        return {
            'n_diseases': len(entries),
            'known_inheritance': known,
            'is_pleiotropic': len(entries) > 1,
        }

    def predict(self, patient_hpo_ids, gene_symbol,
                patient_demographics=None, onsets=None):
        """
        Predict inheritance for a single (patient, gene) pair.

        Args:
            patient_hpo_ids:      list of HPO ID strings
            gene_symbol:          str
            patient_demographics: dict with 'age' key, e.g. {'age': '7 years old'}
            onsets:               dict mapping HPO ID -> onset string,
                                  e.g. {'HP:0001250': '9 months old'}

        Returns:
            dict with:
                prediction:      'AD', 'AR', or 'Uncertain'
                probability_AD:  float (0-1)
                confidence:      float (0-0.5)
                threshold:       float
                n_diseases:      int (known OMIM diseases for this gene)
                known_inheritance: str ('AD', 'AR', 'AD+AR', or 'unknown')
                is_pleiotropic:  bool
                gene_in_omim:    bool
                gene_in_gnomad:  bool
                warning:         str or None
        """
        # Build phenotype records from patient input
        phenotype_records = []
        onset_map = onsets or {}
        for hid in patient_hpo_ids:
            record = {'id': hid}
            if hid in onset_map:
                record['onset'] = onset_map[hid]
            phenotype_records.append(record)

        # Data availability
        has_omim = gene_symbol in self.index
        has_gnomad = gene_symbol in self.gnomad_lookup

        # Build features directly from patient HPO terms
        X = self._build_feature_vector(
            phenotype_records, gene_symbol, patient_demographics
        )
        prob_ad = float(self.model.predict_proba(X)[0, 1])
        prediction, confidence = self._classify(prob_ad)

        # If gene not in gnomAD, force Uncertain (constraint features
        # default to 0.0 which misleadingly signals "not constrained")
        warning = None
        if not has_gnomad and not has_omim:
            prediction = 'Uncertain'
            confidence = 0.0
            warning = 'Gene not in OMIM or gnomAD — prediction unreliable.'
        elif not has_gnomad:
            warning = 'Gene not in gnomAD — prediction based on phenotype only.'

        meta = self._get_gene_metadata(gene_symbol)

        return {
            'prediction': prediction,
            'probability_AD': prob_ad,
            'confidence': confidence,
            'threshold': self.confidence_threshold,
            'n_diseases': meta['n_diseases'],
            'known_inheritance': meta['known_inheritance'],
            'is_pleiotropic': meta['is_pleiotropic'],
            'gene_in_omim': has_omim,
            'gene_in_gnomad': has_gnomad,
            'warning': warning,
        }

    def predict_batch(self, patient_hpo_ids, gene_symbols,
                      patient_demographics=None, onsets=None):
        """
        Predict for multiple candidate genes for the same patient.
        Returns pd.DataFrame.
        """
        results = []
        for gene in gene_symbols:
            pred = self.predict(
                patient_hpo_ids, gene,
                patient_demographics, onsets
            )
            pred['gene_symbol'] = gene
            results.append(pred)

        df = pd.DataFrame(results)
        cols = ['gene_symbol', 'prediction', 'probability_AD', 'confidence',
                'n_diseases', 'known_inheritance', 'is_pleiotropic',
                'gene_in_omim', 'gene_in_gnomad', 'warning']
        extra = [c for c in df.columns if c not in cols]
        return df[cols + extra]


# ============================================================
# Loader
# ============================================================

def load_predictor(model_path, confidence_threshold=None):
    """
    Load a trained model and return an InheritancePredictor.

    Args:
        model_path:           str, path to .joblib file from train_model.py
        confidence_threshold: float or None. Override saved threshold.

    Returns:
        InheritancePredictor instance ready for .predict() / .predict_batch()
    """
    print(f"Loading model from {model_path}...")
    saved = joblib.load(model_path)

    threshold = confidence_threshold or saved.get(
        'confidence_threshold', InheritancePredictor.DEFAULT_CONFIDENCE_THRESHOLD
    )

    predictor = InheritancePredictor(
        model=saved['model'],
        gene_disease_index=saved['gene_disease_index'],
        feature_names=saved['feature_names'],
        training_tiers=saved['training_tiers'],
        training_use_onset=saved['training_use_onset'],
        feature_columns_kept=saved['feature_columns_kept'],
        gnomad_lookup=saved.get('gnomad_lookup', {}),
        confidence_threshold=threshold,
    )

    metrics = saved.get('eval_metrics', {})
    print(f"  Version:    {saved.get('version', 'unknown')}")
    print(f"  Features:   {len(saved['feature_names'])}")
    print(f"  Tiers:      {saved['training_tiers']}")
    print(f"  Onset:      {saved['training_use_onset']}")
    print(f"  Threshold:  {threshold}")
    if metrics:
        print(f"  Test AUC:   {metrics.get('auc_test', 'N/A')}")
        print(f"  Dual AUC:   {metrics.get('auc_dual_holdout', 'N/A')}")
    print(f"  Ready.\n")

    return predictor


# ============================================================
# CLI
# ============================================================

def _parse_onsets(onset_str):
    """Parse 'HP:001=9m,HP:002=birth' into dict."""
    if not onset_str:
        return None
    onsets = {}
    for pair in onset_str.split(','):
        if '=' in pair:
            hpo, age = pair.split('=', 1)
            onsets[hpo.strip()] = age.strip()
    return onsets


TSV_COLUMNS = [
    'gene_symbol', 'prediction', 'probability_AD', 'confidence',
    'n_diseases', 'known_inheritance', 'is_pleiotropic',
    'gene_in_omim', 'gene_in_gnomad', 'warning',
]


if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser(
        description="Predict AD/AR inheritance. Outputs TSV to stdout.",
        epilog="""
Examples:
  # Single gene
  python predict_inheritance.py --model model.joblib --gene SCN1A \\
      --hpo HP:0001250,HP:0001263 --age "2 years old"

  # Multiple genes
  python predict_inheritance.py --model model.joblib \\
      --genes SCN1A,CFTR,LMNA,GNE \\
      --hpo HP:0001250,HP:0001263 --age "2 years old"

  # Pipe into downstream analysis
  python predict_inheritance.py --model model.joblib \\
      --genes SCN1A,CFTR,LMNA --hpo HP:0001250 \\
      --no-header | cut -f1,2,3

  # Read genes from file (one per line)
  python predict_inheritance.py --model model.joblib \\
      --genes-file candidates.txt --hpo HP:0001250,HP:0001263
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--model', required=True, help='Path to .joblib model')

    gene_group = parser.add_mutually_exclusive_group(required=True)
    gene_group.add_argument('--gene', help='Single gene symbol')
    gene_group.add_argument('--genes', help='Comma-separated gene symbols')
    gene_group.add_argument('--genes-file', help='File with one gene per line')

    parser.add_argument('--hpo', required=True, help='Comma-separated HPO IDs')
    parser.add_argument('--onsets', default=None,
                        help='Onset info: HP:001=9m,HP:002=birth')
    parser.add_argument('--age', default=None, help='Patient age')
    parser.add_argument('--threshold', type=float, default=None,
                        help='Confidence threshold override')
    parser.add_argument('--no-header', action='store_true',
                        help='Omit TSV header line')
    parser.add_argument('--output', default=None,
                        help='Output file (default: stdout)')
    args = parser.parse_args()

    # Suppress model loading messages when piping
    predictor = load_predictor(args.model, args.threshold)

    # Parse inputs
    hpo_ids = [h.strip() for h in args.hpo.split(',')]
    onsets = _parse_onsets(args.onsets)
    demographics = {'age': args.age} if args.age else None

    # Resolve gene list
    if args.gene:
        gene_symbols = [args.gene]
    elif args.genes:
        gene_symbols = [g.strip() for g in args.genes.split(',')]
    elif args.genes_file:
        with open(args.genes_file) as f:
            gene_symbols = [line.strip() for line in f if line.strip()]

    # Run predictions
    df = predictor.predict_batch(
        patient_hpo_ids=hpo_ids,
        gene_symbols=gene_symbols,
        patient_demographics=demographics,
        onsets=onsets,
    )

    # Output TSV
    out = open(args.output, 'w') if args.output else sys.stdout
    cols = [c for c in TSV_COLUMNS if c in df.columns]
    df[cols].to_csv(out, sep='\t', index=False, header=not args.no_header)
    if args.output:
        out.close()
        print(f"Results written to {args.output}", file=sys.stderr)
