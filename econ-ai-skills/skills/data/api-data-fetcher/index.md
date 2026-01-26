---
name: api-data-fetcher
description: Fetch economic data from FRED, World Bank, and other APIs
workflow_stage: data
compatibility:
  - claude-code
  - cursor
  - codex
  - gemini-cli
author: Awesome Econ AI Community
version: 1.0.0
tags:
  - Python
  - API
  - FRED
  - World-Bank
  - data-collection
---

# API Data Fetcher

## Purpose

This skill helps economists fetch data from major economic data APIs including FRED (Federal Reserve Economic Data), World Bank, IMF, BLS, and OECD. It generates clean, documented Python code with proper error handling.

## When to Use

- Downloading macroeconomic indicators
- Building custom datasets from multiple sources
- Automating data updates for ongoing projects
- Fetching cross-country panel data

## Instructions

### Step 1: Identify Data Requirements

Ask the user:
1. What data do you need? (GDP, unemployment, inflation, etc.)
2. What time period and frequency?
3. What countries/regions?
4. Preferred output format? (CSV, DataFrame, etc.)

### Step 2: Select Appropriate API

| Data Type | Best Source | Package |
|-----------|------------|---------|
| US macro | FRED | `fredapi` |
| Global development | World Bank | `wbdata` |
| Labor statistics | BLS | `bls` |
| Cross-country | OECD | `pandasdmx` |
| Financial | Yahoo Finance | `yfinance` |

### Step 3: Generate Clean Code

Include:
- API key handling (environment variables)
- Error handling for API failures
- Data cleaning and formatting
- Documentation of series definitions

## Example Output

```python
"""
Economic Data Fetcher
=====================
Downloads macroeconomic data from FRED and World Bank APIs.
Requires: fredapi, wbdata, pandas

Setup: Set FRED_API_KEY environment variable
Get a free key from: https://fred.stlouisfed.org/docs/api/api_key.html
"""

import os
import pandas as pd
from datetime import datetime, timedelta
from typing import List, Optional, Dict

# ============================================
# FRED Data Fetcher
# ============================================

def fetch_fred_series(
    series_ids: List[str],
    start_date: str = "2000-01-01",
    end_date: Optional[str] = None,
    api_key: Optional[str] = None
) -> pd.DataFrame:
    """
    Fetch time series data from FRED.
    
    Parameters
    ----------
    series_ids : list of str
        FRED series IDs (e.g., ['GDP', 'UNRATE', 'CPIAUCSL'])
    start_date : str
        Start date in YYYY-MM-DD format
    end_date : str, optional
        End date (defaults to today)
    api_key : str, optional
        FRED API key (defaults to FRED_API_KEY env var)
    
    Returns
    -------
    pd.DataFrame
        DataFrame with date index and series as columns
    
    Example
    -------
    >>> df = fetch_fred_series(['GDP', 'UNRATE'], '2010-01-01')
    """
    try:
        from fredapi import Fred
    except ImportError:
        raise ImportError("Install fredapi: pip install fredapi")
    
    # Get API key
    api_key = api_key or os.environ.get('FRED_API_KEY')
    if not api_key:
        raise ValueError(
            "FRED API key required. Set FRED_API_KEY environment variable "
            "or pass api_key parameter. Get a key at: "
            "https://fred.stlouisfed.org/docs/api/api_key.html"
        )
    
    fred = Fred(api_key=api_key)
    end_date = end_date or datetime.now().strftime('%Y-%m-%d')
    
    # Fetch each series
    data = {}
    for series_id in series_ids:
        try:
            series = fred.get_series(
                series_id,
                observation_start=start_date,
                observation_end=end_date
            )
            data[series_id] = series
            print(f"✓ Downloaded {series_id}")
        except Exception as e:
            print(f"✗ Failed to download {series_id}: {e}")
    
    # Combine into DataFrame
    df = pd.DataFrame(data)
    df.index.name = 'date'
    
    return df


# Common FRED series for economists
FRED_SERIES = {
    # GDP and Output
    'GDP': 'Gross Domestic Product',
    'GDPC1': 'Real GDP',
    'GDPPOT': 'Real Potential GDP',
    
    # Labor Market
    'UNRATE': 'Unemployment Rate',
    'PAYEMS': 'Total Nonfarm Payrolls',
    'CIVPART': 'Labor Force Participation Rate',
    
    # Prices
    'CPIAUCSL': 'Consumer Price Index',
    'PCEPI': 'PCE Price Index',
    'CPILFESL': 'Core CPI',
    
    # Interest Rates
    'FEDFUNDS': 'Federal Funds Rate',
    'DGS10': '10-Year Treasury Rate',
    'T10Y2Y': '10Y-2Y Treasury Spread',
    
    # Money and Credit
    'M2SL': 'M2 Money Stock',
    'TOTRESNS': 'Total Reserves',
}


# ============================================
# World Bank Data Fetcher
# ============================================

def fetch_world_bank_data(
    indicators: Dict[str, str],
    countries: List[str] = ['USA', 'GBR', 'DEU', 'FRA', 'JPN'],
    start_year: int = 2000,
    end_year: Optional[int] = None
) -> pd.DataFrame:
    """
    Fetch indicator data from World Bank.
    
    Parameters
    ----------
    indicators : dict
        Dict mapping indicator codes to names
        e.g., {'NY.GDP.PCAP.CD': 'gdp_per_capita'}
    countries : list of str
        ISO 3-letter country codes
    start_year : int
        Start year
    end_year : int, optional
        End year (defaults to current year)
    
    Returns
    -------
    pd.DataFrame
        Panel data with country and year
    
    Example
    -------
    >>> indicators = {
    ...     'NY.GDP.PCAP.CD': 'gdp_per_capita',
    ...     'SP.POP.TOTL': 'population'
    ... }
    >>> df = fetch_world_bank_data(indicators, ['USA', 'GBR'])
    """
    try:
        import wbdata
    except ImportError:
        raise ImportError("Install wbdata: pip install wbdata")
    
    end_year = end_year or datetime.now().year
    
    all_data = []
    
    for indicator_code, indicator_name in indicators.items():
        try:
            # Fetch data
            data = wbdata.get_dataframe(
                {indicator_code: indicator_name},
                country=countries,
            )
            data = data.reset_index()
            all_data.append(data)
            print(f"✓ Downloaded {indicator_name}")
            
        except Exception as e:
            print(f"✗ Failed to download {indicator_name}: {e}")
    
    # Merge all indicators
    if all_data:
        df = all_data[0]
        for other_df in all_data[1:]:
            df = df.merge(other_df, on=['country', 'date'], how='outer')
        
        # Filter years
        df['year'] = pd.to_datetime(df['date']).dt.year
        df = df[(df['year'] >= start_year) & (df['year'] <= end_year)]
        
        return df
    
    return pd.DataFrame()


# Common World Bank indicators
WORLD_BANK_INDICATORS = {
    # Income and Growth
    'NY.GDP.PCAP.CD': 'GDP per capita (current US$)',
    'NY.GDP.PCAP.KD.ZG': 'GDP per capita growth (%)',
    'NY.GDP.MKTP.KD.ZG': 'GDP growth (%)',
    
    # Population
    'SP.POP.TOTL': 'Population, total',
    'SP.URB.TOTL.IN.ZS': 'Urban population (%)',
    
    # Trade
    'NE.TRD.GNFS.ZS': 'Trade (% of GDP)',
    'BX.KLT.DINV.WD.GD.ZS': 'FDI, net inflows (% of GDP)',
    
    # Human Capital
    'SE.XPD.TOTL.GD.ZS': 'Education expenditure (% of GDP)',
    'SH.XPD.CHEX.GD.ZS': 'Health expenditure (% of GDP)',
    
    # Inequality
    'SI.POV.GINI': 'Gini index',
    'SI.POV.DDAY': 'Poverty headcount ratio ($1.90/day)',
}


# ============================================
# Usage Example
# ============================================

if __name__ == "__main__":
    # Example 1: Fetch US macro data from FRED
    us_macro = fetch_fred_series(
        series_ids=['GDP', 'UNRATE', 'CPIAUCSL', 'FEDFUNDS'],
        start_date='2010-01-01'
    )
    
    print("\nUS Macro Data (FRED):")
    print(us_macro.tail())
    
    # Save to CSV
    us_macro.to_csv('data/us_macro_fred.csv')
    print("\nSaved to data/us_macro_fred.csv")
    
    # Example 2: Fetch cross-country data from World Bank
    indicators = {
        'NY.GDP.PCAP.CD': 'gdp_per_capita',
        'SP.POP.TOTL': 'population',
        'NY.GDP.MKTP.KD.ZG': 'gdp_growth'
    }
    
    cross_country = fetch_world_bank_data(
        indicators=indicators,
        countries=['USA', 'GBR', 'DEU', 'FRA', 'JPN', 'CHN', 'IND', 'BRA'],
        start_year=2000
    )
    
    print("\nCross-Country Data (World Bank):")
    print(cross_country.head(10))
    
    # Save to CSV
    cross_country.to_csv('data/cross_country_wb.csv', index=False)
    print("\nSaved to data/cross_country_wb.csv")
```

## Requirements

### Python Packages
```bash
pip install fredapi wbdata pandas
```

### API Keys
- **FRED**: Free key from https://fred.stlouisfed.org/docs/api/api_key.html
- **World Bank**: No key required
- **BLS**: Free key from https://www.bls.gov/developers/

Set environment variables:
```bash
export FRED_API_KEY="your_key_here"
```

## Best Practices

1. **Store API keys in environment variables** - never hardcode
2. **Add rate limiting** for bulk downloads
3. **Cache data locally** to avoid repeated API calls
4. **Document series definitions** from the source
5. **Check for revisions** in real-time data

## Common Pitfalls

- ❌ Hardcoding API keys in scripts
- ❌ Not handling API rate limits
- ❌ Ignoring data vintages/revisions
- ❌ Mixing data frequencies without proper handling

## References

- [FRED API Documentation](https://fred.stlouisfed.org/docs/api/)
- [World Bank Data API](https://datahelpdesk.worldbank.org/knowledgebase/topics/125589)
- [QuantEcon: Python Data Sources](https://python-programming.quantecon.org/)

## Changelog

### v1.0.0
- Initial release with FRED and World Bank support
