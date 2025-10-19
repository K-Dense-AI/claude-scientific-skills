#!/usr/bin/env python3
"""
Exploratory Data Analysis Analyzer
Comprehensive data analysis tool that handles multiple file formats and generates
detailed statistical analysis, insights, and data quality reports.
"""

import os
import sys
import json
import argparse
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import normaltest, shapiro, kstest, anderson


class EDAAnalyzer:
    """Main EDA analysis engine"""

    def __init__(self, file_path: str, output_dir: Optional[str] = None):
        self.file_path = Path(file_path)
        self.output_dir = Path(output_dir) if output_dir else self.file_path.parent
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.df = None
        self.analysis_results = {}

    def load_data(self) -> pd.DataFrame:
        """Auto-detect file type and load data"""
        file_ext = self.file_path.suffix.lower()

        try:
            if file_ext == '.csv':
                self.df = pd.read_csv(self.file_path)
            elif file_ext in ['.xlsx', '.xls']:
                self.df = pd.read_excel(self.file_path)
            elif file_ext == '.json':
                self.df = pd.read_json(self.file_path)
            elif file_ext == '.parquet':
                self.df = pd.read_parquet(self.file_path)
            elif file_ext == '.tsv':
                self.df = pd.read_csv(self.file_path, sep='\t')
            elif file_ext == '.feather':
                self.df = pd.read_feather(self.file_path)
            elif file_ext == '.h5' or file_ext == '.hdf5':
                self.df = pd.read_hdf(self.file_path)
            elif file_ext == '.pkl' or file_ext == '.pickle':
                self.df = pd.read_pickle(self.file_path)
            else:
                raise ValueError(f"Unsupported file format: {file_ext}")

            print(f"✅ Successfully loaded {file_ext} file with shape {self.df.shape}")
            return self.df

        except Exception as e:
            print(f"❌ Error loading file: {str(e)}")
            sys.exit(1)

    def basic_info(self) -> Dict[str, Any]:
        """Generate basic dataset information"""
        info = {
            'rows': len(self.df),
            'columns': len(self.df.columns),
            'column_names': list(self.df.columns),
            'dtypes': self.df.dtypes.astype(str).to_dict(),
            'memory_usage_mb': self.df.memory_usage(deep=True).sum() / 1024**2
        }

        # Categorize columns by type
        info['numeric_columns'] = list(self.df.select_dtypes(include=[np.number]).columns)
        info['categorical_columns'] = list(self.df.select_dtypes(include=['object', 'category']).columns)
        info['datetime_columns'] = list(self.df.select_dtypes(include=['datetime64']).columns)
        info['boolean_columns'] = list(self.df.select_dtypes(include=['bool']).columns)

        self.analysis_results['basic_info'] = info
        return info

    def missing_data_analysis(self) -> Dict[str, Any]:
        """Analyze missing data patterns"""
        missing_counts = self.df.isnull().sum()
        missing_pct = (missing_counts / len(self.df) * 100).round(2)

        missing_info = {
            'total_missing_cells': int(self.df.isnull().sum().sum()),
            'missing_percentage': round(self.df.isnull().sum().sum() / (self.df.shape[0] * self.df.shape[1]) * 100, 2),
            'columns_with_missing': {}
        }

        for col in self.df.columns:
            if missing_counts[col] > 0:
                missing_info['columns_with_missing'][col] = {
                    'count': int(missing_counts[col]),
                    'percentage': float(missing_pct[col])
                }

        self.analysis_results['missing_data'] = missing_info
        return missing_info

    def summary_statistics(self) -> Dict[str, Any]:
        """Generate comprehensive summary statistics"""
        stats_dict = {}

        # Numeric columns
        if len(self.df.select_dtypes(include=[np.number]).columns) > 0:
            numeric_stats = self.df.describe().to_dict()
            stats_dict['numeric'] = numeric_stats

            # Additional statistics
            for col in self.df.select_dtypes(include=[np.number]).columns:
                if col not in stats_dict:
                    stats_dict[col] = {}

                data = self.df[col].dropna()
                if len(data) > 0:
                    stats_dict[col].update({
                        'skewness': float(data.skew()),
                        'kurtosis': float(data.kurtosis()),
                        'variance': float(data.var()),
                        'range': float(data.max() - data.min()),
                        'iqr': float(data.quantile(0.75) - data.quantile(0.25)),
                        'cv': float(data.std() / data.mean()) if data.mean() != 0 else np.nan
                    })

        # Categorical columns
        categorical_stats = {}
        for col in self.df.select_dtypes(include=['object', 'category']).columns:
            categorical_stats[col] = {
                'unique_values': int(self.df[col].nunique()),
                'most_common': self.df[col].mode().iloc[0] if len(self.df[col].mode()) > 0 else None,
                'most_common_freq': int(self.df[col].value_counts().iloc[0]) if len(self.df[col].value_counts()) > 0 else 0,
                'most_common_pct': float(self.df[col].value_counts(normalize=True).iloc[0] * 100) if len(self.df[col].value_counts()) > 0 else 0
            }

        if categorical_stats:
            stats_dict['categorical'] = categorical_stats

        self.analysis_results['summary_statistics'] = stats_dict
        return stats_dict

    def outlier_detection(self) -> Dict[str, Any]:
        """Detect outliers using multiple methods"""
        outliers = {}

        for col in self.df.select_dtypes(include=[np.number]).columns:
            data = self.df[col].dropna()
            if len(data) == 0:
                continue

            outliers[col] = {}

            # IQR method
            Q1 = data.quantile(0.25)
            Q3 = data.quantile(0.75)
            IQR = Q3 - Q1
            lower_bound = Q1 - 1.5 * IQR
            upper_bound = Q3 + 1.5 * IQR
            iqr_outliers = data[(data < lower_bound) | (data > upper_bound)]

            outliers[col]['iqr_method'] = {
                'count': len(iqr_outliers),
                'percentage': round(len(iqr_outliers) / len(data) * 100, 2),
                'lower_bound': float(lower_bound),
                'upper_bound': float(upper_bound)
            }

            # Z-score method (|z| > 3)
            if len(data) > 2:
                z_scores = np.abs(stats.zscore(data))
                z_outliers = data[z_scores > 3]
                outliers[col]['zscore_method'] = {
                    'count': len(z_outliers),
                    'percentage': round(len(z_outliers) / len(data) * 100, 2)
                }

        self.analysis_results['outliers'] = outliers
        return outliers

    def distribution_analysis(self) -> Dict[str, Any]:
        """Analyze distributions and test for normality"""
        distributions = {}

        for col in self.df.select_dtypes(include=[np.number]).columns:
            data = self.df[col].dropna()
            if len(data) < 8:  # Need at least 8 samples for tests
                continue

            distributions[col] = {}

            # Shapiro-Wilk test (best for n < 5000)
            if len(data) < 5000:
                try:
                    stat, p_value = shapiro(data)
                    distributions[col]['shapiro_wilk'] = {
                        'statistic': float(stat),
                        'p_value': float(p_value),
                        'is_normal': p_value > 0.05
                    }
                except:
                    pass

            # Anderson-Darling test
            try:
                result = anderson(data)
                distributions[col]['anderson_darling'] = {
                    'statistic': float(result.statistic),
                    'critical_values': result.critical_values.tolist(),
                    'significance_levels': result.significance_level.tolist()
                }
            except:
                pass

            # Distribution characteristics
            distributions[col]['characteristics'] = {
                'skewness': float(data.skew()),
                'skewness_interpretation': self._interpret_skewness(data.skew()),
                'kurtosis': float(data.kurtosis()),
                'kurtosis_interpretation': self._interpret_kurtosis(data.kurtosis())
            }

        self.analysis_results['distributions'] = distributions
        return distributions

    def correlation_analysis(self) -> Dict[str, Any]:
        """Analyze correlations between numeric variables"""
        numeric_df = self.df.select_dtypes(include=[np.number])

        if len(numeric_df.columns) < 2:
            return {}

        correlations = {}

        # Pearson correlation
        pearson_corr = numeric_df.corr(method='pearson')
        correlations['pearson'] = pearson_corr.to_dict()

        # Spearman correlation (rank-based, robust to outliers)
        spearman_corr = numeric_df.corr(method='spearman')
        correlations['spearman'] = spearman_corr.to_dict()

        # Find strong correlations (|r| > 0.7)
        strong_correlations = []
        for i in range(len(pearson_corr.columns)):
            for j in range(i + 1, len(pearson_corr.columns)):
                col1 = pearson_corr.columns[i]
                col2 = pearson_corr.columns[j]
                corr_value = pearson_corr.iloc[i, j]

                if abs(corr_value) > 0.7:
                    strong_correlations.append({
                        'variable1': col1,
                        'variable2': col2,
                        'correlation': float(corr_value),
                        'strength': self._interpret_correlation(corr_value)
                    })

        correlations['strong_correlations'] = strong_correlations

        self.analysis_results['correlations'] = correlations
        return correlations

    def data_quality_assessment(self) -> Dict[str, Any]:
        """Assess overall data quality"""
        quality = {
            'completeness': {
                'score': round((1 - self.df.isnull().sum().sum() / (self.df.shape[0] * self.df.shape[1])) * 100, 2),
                'interpretation': ''
            },
            'duplicates': {
                'count': int(self.df.duplicated().sum()),
                'percentage': round(self.df.duplicated().sum() / len(self.df) * 100, 2)
            },
            'issues': []
        }

        # Completeness interpretation
        if quality['completeness']['score'] > 95:
            quality['completeness']['interpretation'] = 'Excellent'
        elif quality['completeness']['score'] > 90:
            quality['completeness']['interpretation'] = 'Good'
        elif quality['completeness']['score'] > 80:
            quality['completeness']['interpretation'] = 'Fair'
        else:
            quality['completeness']['interpretation'] = 'Poor'

        # Identify potential issues
        if quality['duplicates']['count'] > 0:
            quality['issues'].append(f"Found {quality['duplicates']['count']} duplicate rows")

        if quality['completeness']['score'] < 90:
            quality['issues'].append("Missing data exceeds 10% threshold")

        # Check for constant columns
        constant_cols = [col for col in self.df.columns if self.df[col].nunique() == 1]
        if constant_cols:
            quality['issues'].append(f"Constant columns detected: {', '.join(constant_cols)}")
            quality['constant_columns'] = constant_cols

        # Check for high cardinality
        high_cardinality_cols = []
        for col in self.df.select_dtypes(include=['object']).columns:
            if self.df[col].nunique() > len(self.df) * 0.9:
                high_cardinality_cols.append(col)

        if high_cardinality_cols:
            quality['issues'].append(f"High cardinality columns (>90% unique): {', '.join(high_cardinality_cols)}")
            quality['high_cardinality_columns'] = high_cardinality_cols

        self.analysis_results['data_quality'] = quality
        return quality

    def generate_insights(self) -> List[str]:
        """Generate automated insights from the analysis"""
        insights = []

        # Dataset size insights
        info = self.analysis_results.get('basic_info', {})
        if info.get('rows', 0) > 1000000:
            insights.append(f"📊 Large dataset with {info['rows']:,} rows - consider sampling for faster iteration")

        # Missing data insights
        missing = self.analysis_results.get('missing_data', {})
        if missing.get('missing_percentage', 0) > 20:
            insights.append(f"⚠️ Significant missing data ({missing['missing_percentage']}%) - imputation or removal may be needed")

        # Correlation insights
        correlations = self.analysis_results.get('correlations', {})
        strong_corrs = correlations.get('strong_correlations', [])
        if len(strong_corrs) > 0:
            insights.append(f"🔗 Found {len(strong_corrs)} strong correlations - potential for feature engineering or multicollinearity")

        # Outlier insights
        outliers = self.analysis_results.get('outliers', {})
        high_outlier_cols = [col for col, data in outliers.items()
                           if data.get('iqr_method', {}).get('percentage', 0) > 5]
        if high_outlier_cols:
            insights.append(f"🎯 Columns with high outlier rates (>5%): {', '.join(high_outlier_cols)}")

        # Distribution insights
        distributions = self.analysis_results.get('distributions', {})
        skewed_cols = [col for col, data in distributions.items()
                      if abs(data.get('characteristics', {}).get('skewness', 0)) > 1]
        if skewed_cols:
            insights.append(f"📈 Highly skewed distributions detected in: {', '.join(skewed_cols)} - consider transformations")

        # Data quality insights
        quality = self.analysis_results.get('data_quality', {})
        if quality.get('duplicates', {}).get('count', 0) > 0:
            insights.append(f"🔄 {quality['duplicates']['count']} duplicate rows found - consider deduplication")

        # Categorical insights
        stats = self.analysis_results.get('summary_statistics', {})
        categorical = stats.get('categorical', {})
        imbalanced_cols = [col for col, data in categorical.items()
                          if data.get('most_common_pct', 0) > 90]
        if imbalanced_cols:
            insights.append(f"⚖️ Highly imbalanced categorical variables: {', '.join(imbalanced_cols)}")

        self.analysis_results['insights'] = insights
        return insights

    def _interpret_skewness(self, skew: float) -> str:
        """Interpret skewness value"""
        if abs(skew) < 0.5:
            return "Approximately symmetric"
        elif skew > 0.5:
            return "Right-skewed (positive skew)"
        else:
            return "Left-skewed (negative skew)"

    def _interpret_kurtosis(self, kurt: float) -> str:
        """Interpret kurtosis value"""
        if abs(kurt) < 0.5:
            return "Mesokurtic (normal-like tails)"
        elif kurt > 0.5:
            return "Leptokurtic (heavy tails)"
        else:
            return "Platykurtic (light tails)"

    def _interpret_correlation(self, corr: float) -> str:
        """Interpret correlation strength"""
        abs_corr = abs(corr)
        if abs_corr > 0.9:
            return "Very strong"
        elif abs_corr > 0.7:
            return "Strong"
        elif abs_corr > 0.5:
            return "Moderate"
        elif abs_corr > 0.3:
            return "Weak"
        else:
            return "Very weak"

    def run_full_analysis(self) -> Dict[str, Any]:
        """Run complete EDA analysis"""
        print("🔍 Starting comprehensive EDA analysis...")

        self.load_data()

        print("📊 Analyzing basic information...")
        self.basic_info()

        print("🔎 Analyzing missing data...")
        self.missing_data_analysis()

        print("📈 Computing summary statistics...")
        self.summary_statistics()

        print("🎯 Detecting outliers...")
        self.outlier_detection()

        print("📉 Analyzing distributions...")
        self.distribution_analysis()

        print("🔗 Computing correlations...")
        self.correlation_analysis()

        print("✅ Assessing data quality...")
        self.data_quality_assessment()

        print("💡 Generating insights...")
        self.generate_insights()

        print("✨ Analysis complete!")

        return self.analysis_results

    def save_results(self, format='json') -> str:
        """Save analysis results to file"""
        output_file = self.output_dir / f"eda_analysis.{format}"

        if format == 'json':
            with open(output_file, 'w') as f:
                json.dump(self.analysis_results, f, indent=2, default=str)

        print(f"💾 Results saved to: {output_file}")
        return str(output_file)


def main():
    parser = argparse.ArgumentParser(description='Perform comprehensive exploratory data analysis')
    parser.add_argument('file_path', help='Path to data file')
    parser.add_argument('-o', '--output', help='Output directory for results', default=None)
    parser.add_argument('-f', '--format', choices=['json'], default='json', help='Output format')

    args = parser.parse_args()

    analyzer = EDAAnalyzer(args.file_path, args.output)
    analyzer.run_full_analysis()
    analyzer.save_results(format=args.format)


if __name__ == '__main__':
    main()
