---
name: social-media-analytics
description: Track, analyze, and report on social media performance across Facebook, Instagram, LinkedIn, and TikTok. Computes engagement rates, identifies top-performing content, benchmarks against competitors, calculates ROI, and generates actionable growth recommendations. Produces analytics reports in Markdown or DOCX format.
allowed-tools: [Bash, Read, Write, Edit, WebSearch, WebFetch, Agent]
license: MIT license
metadata:
    skill-author: K-Dense Inc.
---

# Social Media Analytics

## Overview

This skill provides a full analytics engine for social media performance tracking and optimization. It collects metrics from Facebook, Instagram, LinkedIn, and TikTok, computes engagement KPIs, benchmarks performance against competitors and industry averages, calculates ROI, and generates actionable reports with growth recommendations.

## When to Use

- **Performance Tracking**: Monitor engagement, reach, and follower growth over time.
- **Content Optimization**: Identify top-performing content formats and topics to guide future strategy.
- **Competitor Benchmarking**: Compare your metrics against competitors to find gaps and opportunities.
- **Reporting**: Generate weekly, monthly, or quarterly analytics reports for stakeholders.
- **ROI Analysis**: Quantify the return on time, ad spend, and tooling invested in social media.
- **Strategy Planning**: Use data-driven recommendations to refine posting cadence, content mix, and audience targeting.

## Quick Start

### Collect Metrics
```bash
python scripts/collect_metrics.py \
  --platforms instagram,linkedin,tiktok \
  --period 30d \
  --output metrics.json
```

### Analyze Performance
```bash
python scripts/analyze_performance.py \
  --metrics metrics.json \
  --output analysis.json
```

### Benchmark Against Competitors
```bash
python scripts/competitor_benchmark.py \
  --brand-metrics metrics.json \
  --competitor-profiles @competitor1,@competitor2 \
  --output benchmark.json
```

### Generate Report
```bash
python scripts/generate_report.py \
  --analysis analysis.json \
  --format md \
  --output report.md
```

### Calculate ROI
```bash
python scripts/roi_calculator.py \
  --metrics metrics.json \
  --costs costs.json \
  --revenue-attributable 12000 \
  --output roi.json
```

## Core Metrics Tracked

### Awareness Metrics
- **Reach**: Unique users who saw the content.
- **Impressions**: Total number of times content was displayed.
- **Share of Voice**: Brand mentions relative to total category mentions.
- **Brand Mention Rate**: Frequency of brand mentions across platforms.

### Engagement Metrics
- **Engagement Rate**: Total interactions divided by reach or followers.
- **Comment Rate**: Comments per post relative to followers.
- **Comment Sentiment**: Ratio of positive, neutral, and negative comments.
- **Share/Repost Rate**: Shares per post relative to followers.
- **Save/Bookmark Rate**: Saves per post relative to reach.
- **Amplification Rate**: Shares divided by followers.
- **Virality Rate**: Shares divided by impressions.

### Content Metrics
- **Click-Through Rate (CTR)**: Link clicks divided by impressions.
- **Video View Rate**: Video views divided by impressions.
- **Average Watch Duration**: Mean seconds watched per video view.
- **Video Completion Rate**: Percentage of viewers who watched the entire video.
- **Story Completion Rate**: Percentage of viewers who watched all story frames.
- **Story Reply Rate**: Replies per story view.

### Growth Metrics
- **Follower Growth Rate**: Net new followers over a period relative to total.
- **Follower-to-Engagement Ratio**: Engagement rate normalized by follower count.
- **Profile Visit Rate**: Profile visits per impression.
- **Website Click Rate**: Website clicks per profile visit.

### Conversion Metrics
- **Conversion Rate**: Conversions divided by link clicks.
- **Cost Per Click (CPC)**: Ad spend divided by clicks.
- **Cost Per Engagement (CPE)**: Ad spend divided by engagements.
- **Return on Ad Spend (ROAS)**: Revenue generated divided by ad spend.
- **Content ROI**: (Revenue attributed - total cost) / total cost.

## Analysis Workflows

### 1. Performance Collection via Composio/Chrome

Metrics are collected through one of two methods:

- **Composio MCP Integration**: If Composio connectors are available for the target platforms, the skill uses them to pull analytics data directly via API.
- **Chrome Automation**: If API access is not available, the skill documents step-by-step Chrome automation workflows to extract metrics from each platform's native analytics dashboard.

The collected data is normalized into a structured JSON format shared across all downstream scripts.

### 2. Content Performance Ranking

Each post is scored using a weighted composite of engagement metrics. Posts are ranked and the top 10% and bottom 10% are flagged. Pattern analysis identifies which content attributes (format, length, topic, posting time, hashtags) correlate with high performance.

### 3. Competitor Benchmarking

Public competitor profiles are analyzed for posting frequency, engagement rates, follower growth, and content mix. A gap analysis highlights where the brand underperforms and an opportunity score quantifies the upside of closing each gap.

### 4. Report Generation (MD or DOCX)

Reports include:
- Executive summary with key takeaways
- Platform-by-platform breakdown
- Top-performing content showcase
- Trend analysis with period-over-period comparisons
- Competitor benchmark summary
- Actionable recommendations ranked by expected impact

DOCX output uses the `document-skills/docx` skill for formatting. Charts use the `scientific-visualization` skill.

### 5. ROI Calculation

All costs (time invested at an hourly rate, ad spend, tool subscriptions) are aggregated and compared against attributable revenue. ROI is broken down by platform and content type so teams can allocate resources to the highest-returning channels.

## Growth Recommendations Engine

The skill produces prioritized recommendations based on the analysis results. Each recommendation includes:

- **Action**: What to do.
- **Rationale**: Which data points support the recommendation.
- **Expected Impact**: Estimated improvement range.
- **Priority**: High, Medium, or Low based on effort-to-impact ratio.

See `references/growth_strategies.md` for the full playbook of strategies the engine draws from.

## Integration Points

- **scientific-visualization**: Generate matplotlib/seaborn charts for engagement trends, platform comparisons, and content performance distributions.
- **document-skills/docx**: Produce formatted DOCX reports with branded headers, tables, and embedded charts.
- **Composio MCP**: Pull platform analytics data when API connectors are configured.
- **Chrome Automation**: Fall back to browser-based data extraction when APIs are unavailable.
