#!/usr/bin/env python3
"""
EDA Visualizer
Generate comprehensive visualizations for exploratory data analysis including
distribution plots, correlation heatmaps, time series, and categorical analyses.
"""

import os
import sys
import argparse
from pathlib import Path
from typing import List, Optional, Tuple
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec


class EDAVisualizer:
    """Generate comprehensive EDA visualizations"""

    def __init__(self, file_path: str, output_dir: Optional[str] = None):
        self.file_path = Path(file_path)
        self.output_dir = Path(output_dir) if output_dir else self.file_path.parent / "eda_visualizations"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.df = None

        # Set style
        sns.set_style("whitegrid")
        plt.rcParams['figure.figsize'] = (12, 8)
        plt.rcParams['figure.dpi'] = 100
        plt.rcParams['savefig.dpi'] = 300
        plt.rcParams['savefig.bbox'] = 'tight'

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

    def plot_missing_data(self) -> str:
        """Visualize missing data patterns"""
        fig, axes = plt.subplots(2, 1, figsize=(12, 8))

        # Missing data heatmap
        if self.df.isnull().sum().sum() > 0:
            # Only plot columns with missing data
            missing_cols = self.df.columns[self.df.isnull().any()].tolist()
            if missing_cols:
                sns.heatmap(self.df[missing_cols].isnull(), cbar=True, yticklabels=False,
                           cmap='viridis', ax=axes[0])
                axes[0].set_title('Missing Data Pattern', fontsize=14, fontweight='bold')
                axes[0].set_xlabel('Columns')

                # Missing data bar chart
                missing_pct = (self.df[missing_cols].isnull().sum() / len(self.df) * 100).sort_values(ascending=True)
                missing_pct.plot(kind='barh', ax=axes[1], color='coral')
                axes[1].set_title('Missing Data Percentage by Column', fontsize=14, fontweight='bold')
                axes[1].set_xlabel('Missing %')
                axes[1].set_ylabel('Columns')

                for i, v in enumerate(missing_pct):
                    axes[1].text(v + 0.5, i, f'{v:.1f}%', va='center')
            else:
                axes[0].text(0.5, 0.5, 'No missing data detected', ha='center', va='center',
                           transform=axes[0].transAxes, fontsize=14)
                axes[0].axis('off')
                axes[1].axis('off')
        else:
            axes[0].text(0.5, 0.5, 'No missing data detected', ha='center', va='center',
                       transform=axes[0].transAxes, fontsize=14)
            axes[0].axis('off')
            axes[1].axis('off')

        plt.tight_layout()
        output_path = self.output_dir / "missing_data.png"
        plt.savefig(output_path)
        plt.close()

        print(f"✅ Missing data visualization saved: {output_path}")
        return str(output_path)

    def plot_distributions(self) -> str:
        """Plot distributions for all numeric columns"""
        numeric_cols = self.df.select_dtypes(include=[np.number]).columns.tolist()

        if not numeric_cols:
            print("⚠️ No numeric columns found for distribution plots")
            return ""

        # Limit to first 20 columns if too many
        if len(numeric_cols) > 20:
            print(f"⚠️ Too many numeric columns ({len(numeric_cols)}), plotting first 20")
            numeric_cols = numeric_cols[:20]

        n_cols = min(3, len(numeric_cols))
        n_rows = (len(numeric_cols) + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 4 * n_rows))
        if n_rows == 1 and n_cols == 1:
            axes = np.array([[axes]])
        elif n_rows == 1 or n_cols == 1:
            axes = axes.reshape(n_rows, n_cols)

        for idx, col in enumerate(numeric_cols):
            row = idx // n_cols
            col_idx = idx % n_cols
            ax = axes[row, col_idx]

            data = self.df[col].dropna()

            # Create histogram with KDE
            ax.hist(data, bins=30, alpha=0.6, color='skyblue', edgecolor='black', density=True)

            # Add KDE line
            try:
                data.plot(kind='kde', ax=ax, color='red', linewidth=2)
            except:
                pass

            ax.set_title(f'{col}', fontsize=10, fontweight='bold')
            ax.set_xlabel('Value')
            ax.set_ylabel('Density')

            # Add statistics box
            stats_text = f'Mean: {data.mean():.2f}\nMedian: {data.median():.2f}\nStd: {data.std():.2f}'
            ax.text(0.98, 0.98, stats_text, transform=ax.transAxes,
                   verticalalignment='top', horizontalalignment='right',
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
                   fontsize=8)

        # Hide empty subplots
        for idx in range(len(numeric_cols), n_rows * n_cols):
            row = idx // n_cols
            col_idx = idx % n_cols
            axes[row, col_idx].axis('off')

        plt.suptitle('Distribution Analysis', fontsize=16, fontweight='bold', y=1.00)
        plt.tight_layout()

        output_path = self.output_dir / "distributions.png"
        plt.savefig(output_path)
        plt.close()

        print(f"✅ Distribution plots saved: {output_path}")
        return str(output_path)

    def plot_boxplots(self) -> str:
        """Create box plots for numeric columns to show outliers"""
        numeric_cols = self.df.select_dtypes(include=[np.number]).columns.tolist()

        if not numeric_cols:
            print("⚠️ No numeric columns found for box plots")
            return ""

        # Limit to first 20 columns if too many
        if len(numeric_cols) > 20:
            numeric_cols = numeric_cols[:20]

        n_cols = min(3, len(numeric_cols))
        n_rows = (len(numeric_cols) + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 4 * n_rows))
        if n_rows == 1 and n_cols == 1:
            axes = np.array([[axes]])
        elif n_rows == 1 or n_cols == 1:
            axes = axes.reshape(n_rows, n_cols)

        for idx, col in enumerate(numeric_cols):
            row = idx // n_cols
            col_idx = idx % n_cols
            ax = axes[row, col_idx]

            data = self.df[col].dropna()

            # Box plot with violin
            parts = ax.violinplot([data], positions=[0], widths=0.7, showmeans=True, showextrema=True)
            ax.boxplot([data], positions=[0], widths=0.3, patch_artist=True,
                      boxprops=dict(facecolor='lightblue', alpha=0.7))

            ax.set_title(f'{col}', fontsize=10, fontweight='bold')
            ax.set_ylabel('Value')
            ax.set_xticks([])

        # Hide empty subplots
        for idx in range(len(numeric_cols), n_rows * n_cols):
            row = idx // n_cols
            col_idx = idx % n_cols
            axes[row, col_idx].axis('off')

        plt.suptitle('Box Plots with Violin Plots (Outlier Detection)', fontsize=16, fontweight='bold', y=1.00)
        plt.tight_layout()

        output_path = self.output_dir / "boxplots.png"
        plt.savefig(output_path)
        plt.close()

        print(f"✅ Box plots saved: {output_path}")
        return str(output_path)

    def plot_correlation_heatmap(self) -> str:
        """Create correlation heatmap for numeric variables"""
        numeric_df = self.df.select_dtypes(include=[np.number])

        if len(numeric_df.columns) < 2:
            print("⚠️ Need at least 2 numeric columns for correlation heatmap")
            return ""

        fig, axes = plt.subplots(1, 2, figsize=(20, 8))

        # Pearson correlation
        corr_pearson = numeric_df.corr(method='pearson')
        mask = np.triu(np.ones_like(corr_pearson, dtype=bool))

        sns.heatmap(corr_pearson, mask=mask, annot=True, fmt='.2f', cmap='coolwarm',
                   center=0, square=True, linewidths=1, cbar_kws={"shrink": 0.8},
                   ax=axes[0])
        axes[0].set_title('Pearson Correlation Matrix', fontsize=14, fontweight='bold')

        # Spearman correlation
        corr_spearman = numeric_df.corr(method='spearman')
        mask = np.triu(np.ones_like(corr_spearman, dtype=bool))

        sns.heatmap(corr_spearman, mask=mask, annot=True, fmt='.2f', cmap='coolwarm',
                   center=0, square=True, linewidths=1, cbar_kws={"shrink": 0.8},
                   ax=axes[1])
        axes[1].set_title('Spearman Correlation Matrix', fontsize=14, fontweight='bold')

        plt.tight_layout()

        output_path = self.output_dir / "correlation_heatmap.png"
        plt.savefig(output_path)
        plt.close()

        print(f"✅ Correlation heatmap saved: {output_path}")
        return str(output_path)

    def plot_scatter_matrix(self) -> str:
        """Create scatter plot matrix for numeric variables"""
        numeric_df = self.df.select_dtypes(include=[np.number])

        if len(numeric_df.columns) < 2:
            print("⚠️ Need at least 2 numeric columns for scatter matrix")
            return ""

        # Limit to first 6 columns if too many (scatter matrix gets too large)
        if len(numeric_df.columns) > 6:
            print(f"⚠️ Too many columns for scatter matrix, using first 6")
            numeric_df = numeric_df.iloc[:, :6]

        fig = plt.figure(figsize=(15, 15))
        pd.plotting.scatter_matrix(numeric_df, alpha=0.6, figsize=(15, 15),
                                  diagonal='kde', hist_kwds={'bins': 20})
        plt.suptitle('Scatter Plot Matrix', fontsize=16, fontweight='bold', y=1.00)

        output_path = self.output_dir / "scatter_matrix.png"
        plt.savefig(output_path)
        plt.close()

        print(f"✅ Scatter matrix saved: {output_path}")
        return str(output_path)

    def plot_categorical_analysis(self) -> str:
        """Analyze and visualize categorical variables"""
        categorical_cols = self.df.select_dtypes(include=['object', 'category']).columns.tolist()

        if not categorical_cols:
            print("⚠️ No categorical columns found")
            return ""

        # Limit to first 12 columns if too many
        if len(categorical_cols) > 12:
            print(f"⚠️ Too many categorical columns ({len(categorical_cols)}), plotting first 12")
            categorical_cols = categorical_cols[:12]

        n_cols = min(3, len(categorical_cols))
        n_rows = (len(categorical_cols) + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(8 * n_cols, 5 * n_rows))
        if n_rows == 1 and n_cols == 1:
            axes = np.array([[axes]])
        elif n_rows == 1 or n_cols == 1:
            axes = axes.reshape(n_rows, n_cols)

        for idx, col in enumerate(categorical_cols):
            row = idx // n_cols
            col_idx = idx % n_cols
            ax = axes[row, col_idx]

            # Get top 10 categories
            value_counts = self.df[col].value_counts().head(10)

            # Create bar chart
            value_counts.plot(kind='barh', ax=ax, color='steelblue')
            ax.set_title(f'{col} (Top 10)', fontsize=11, fontweight='bold')
            ax.set_xlabel('Count')
            ax.set_ylabel('')

            # Add value labels
            for i, v in enumerate(value_counts):
                ax.text(v + max(value_counts) * 0.01, i, str(v), va='center')

        # Hide empty subplots
        for idx in range(len(categorical_cols), n_rows * n_cols):
            row = idx // n_cols
            col_idx = idx % n_cols
            axes[row, col_idx].axis('off')

        plt.suptitle('Categorical Variable Analysis', fontsize=16, fontweight='bold', y=1.00)
        plt.tight_layout()

        output_path = self.output_dir / "categorical_analysis.png"
        plt.savefig(output_path)
        plt.close()

        print(f"✅ Categorical analysis saved: {output_path}")
        return str(output_path)

    def plot_time_series(self) -> str:
        """Create time series visualizations if datetime columns exist"""
        datetime_cols = self.df.select_dtypes(include=['datetime64']).columns.tolist()

        # Also check for columns that might be dates but stored as strings
        for col in self.df.columns:
            if self.df[col].dtype == 'object':
                try:
                    pd.to_datetime(self.df[col].head(100))
                    datetime_cols.append(col)
                except:
                    pass

        if not datetime_cols:
            print("⚠️ No datetime columns found for time series analysis")
            return ""

        # Take first datetime column as index
        date_col = datetime_cols[0]
        df_temp = self.df.copy()

        if df_temp[date_col].dtype == 'object':
            df_temp[date_col] = pd.to_datetime(df_temp[date_col])

        df_temp = df_temp.sort_values(date_col)

        # Get numeric columns
        numeric_cols = df_temp.select_dtypes(include=[np.number]).columns.tolist()

        if not numeric_cols:
            print("⚠️ No numeric columns found for time series plots")
            return ""

        # Limit to first 6 numeric columns
        if len(numeric_cols) > 6:
            numeric_cols = numeric_cols[:6]

        n_rows = len(numeric_cols)
        fig, axes = plt.subplots(n_rows, 1, figsize=(14, 4 * n_rows))

        if n_rows == 1:
            axes = [axes]

        for idx, col in enumerate(numeric_cols):
            ax = axes[idx]

            # Plot time series
            ax.plot(df_temp[date_col], df_temp[col], linewidth=1, alpha=0.8)
            ax.set_title(f'{col} over Time', fontsize=12, fontweight='bold')
            ax.set_xlabel('Date')
            ax.set_ylabel(col)
            ax.grid(True, alpha=0.3)

            # Add trend line
            try:
                z = np.polyfit(range(len(df_temp)), df_temp[col].fillna(df_temp[col].mean()), 1)
                p = np.poly1d(z)
                ax.plot(df_temp[date_col], p(range(len(df_temp))), "r--", linewidth=2, alpha=0.8, label='Trend')
                ax.legend()
            except:
                pass

        plt.suptitle('Time Series Analysis', fontsize=16, fontweight='bold', y=1.00)
        plt.tight_layout()

        output_path = self.output_dir / "time_series.png"
        plt.savefig(output_path)
        plt.close()

        print(f"✅ Time series plots saved: {output_path}")
        return str(output_path)

    def generate_all_visualizations(self) -> List[str]:
        """Generate all visualizations"""
        print("🎨 Starting visualization generation...")

        self.load_data()
        generated_files = []

        print("📊 Creating missing data visualization...")
        missing_plot = self.plot_missing_data()
        if missing_plot:
            generated_files.append(missing_plot)

        print("📈 Creating distribution plots...")
        dist_plot = self.plot_distributions()
        if dist_plot:
            generated_files.append(dist_plot)

        print("📦 Creating box plots...")
        box_plot = self.plot_boxplots()
        if box_plot:
            generated_files.append(box_plot)

        print("🔥 Creating correlation heatmap...")
        corr_plot = self.plot_correlation_heatmap()
        if corr_plot:
            generated_files.append(corr_plot)

        print("🔢 Creating scatter matrix...")
        scatter_plot = self.plot_scatter_matrix()
        if scatter_plot:
            generated_files.append(scatter_plot)

        print("📊 Creating categorical analysis...")
        cat_plot = self.plot_categorical_analysis()
        if cat_plot:
            generated_files.append(cat_plot)

        print("⏱️ Creating time series plots...")
        ts_plot = self.plot_time_series()
        if ts_plot:
            generated_files.append(ts_plot)

        print(f"✨ Generated {len(generated_files)} visualizations!")
        print(f"📁 Saved to: {self.output_dir}")

        return generated_files


def main():
    parser = argparse.ArgumentParser(description='Generate comprehensive EDA visualizations')
    parser.add_argument('file_path', help='Path to data file')
    parser.add_argument('-o', '--output', help='Output directory for visualizations', default=None)

    args = parser.parse_args()

    visualizer = EDAVisualizer(args.file_path, args.output)
    visualizer.generate_all_visualizations()


if __name__ == '__main__':
    main()
