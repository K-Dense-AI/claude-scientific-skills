---
name: beamer-presentation
description: Create academic presentations in Beamer with professional themes
workflow_stage: communication
compatibility:
  - claude-code
  - cursor
  - codex
  - gemini-cli
author: Awesome Econ AI Community
version: 1.0.0
tags:
  - LaTeX
  - Beamer
  - presentations
  - slides
---

# Beamer Presentation Creator

## Purpose

This skill helps economists create professional academic presentations using LaTeX Beamer. It provides templates for conference talks, job market presentations, and seminar presentations with proper structure and clean aesthetics.

## When to Use

- Preparing conference presentations
- Creating job market talk slides
- Making seminar/workshop presentations
- Converting a paper into presentation slides

## Instructions

### Step 1: Understand the Context

Ask the user:
1. What type of presentation? (20-min conference, 90-min seminar, job market)
2. What's the paper/project about?
3. What's the target audience expertise level?
4. Do they have specific style preferences?

### Step 2: Structure by Time

| Duration | Structure |
|----------|-----------|
| 15-20 min | Motivation (2) → Question (1) → Method (2) → Results (3-4) → Conclusion (1) |
| 45-60 min | Add literature review, more results detail, robustness |
| 90 min | Full seminar with theoretical framework, extensive empirics |

### Step 3: Follow Presentation Best Practices

- **One idea per slide**
- **Minimal text** - use bullets of 3-6 words
- **Big fonts** - minimum 20pt for content
- **Consistent colors** - use a limited palette
- **Reveal incrementally** using `\pause` or `<+->` for complex slides

## Example Output

```latex
\documentclass[aspectratio=169, 11pt]{beamer}

% ============================================
% THEME AND APPEARANCE
% ============================================

% Clean minimal theme
\usetheme{metropolis}
\usecolortheme{default}

% Or for a more traditional look:
% \usetheme{Madrid}
% \usecolortheme{whale}

% Custom colors
\definecolor{darkblue}{RGB}{0, 51, 102}
\definecolor{lightgray}{RGB}{245, 245, 245}

\setbeamercolor{frametitle}{bg=darkblue, fg=white}
\setbeamercolor{title}{fg=darkblue}
\setbeamercolor{structure}{fg=darkblue}

% Remove navigation symbols
\setbeamertemplate{navigation symbols}{}

% Frame numbers
\setbeamertemplate{footline}[frame number]

% ============================================
% PACKAGES
% ============================================

\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{tikz}
\usepackage{pgfplots}
\pgfplotsset{compat=1.17}

% ============================================
% TITLE PAGE
% ============================================

\title{The Effect of X on Y: \\Evidence from Z}
\subtitle{Short and Descriptive}
\author{Your Name}
\institute{Your University}
\date{Conference Name \\ Month Year}

\begin{document}

% Title slide
\begin{frame}[plain]
    \titlepage
\end{frame}

% ============================================
% MOTIVATION (2-3 slides)
% ============================================

\begin{frame}{Motivation: Why This Matters}
    \begin{itemize}
        \item<1-> \textbf{Big picture:} [One sentence on broad relevance]
        \item<2-> \textbf{Specific puzzle:} [What we don't know]
        \item<3-> \textbf{Stakes:} [Why should we care?]
    \end{itemize}
    
    \vspace{1em}
    
    \only<4>{
    \begin{block}{Key Statistic}
        \Large \textbf{X\%} of [outcome] can be explained by [factor]
    \end{block}
    }
\end{frame}

\begin{frame}{What We Know (and Don't Know)}
    \textbf{Previous literature:}
    \begin{itemize}
        \item Author et al. (2020): Finding 1
        \item Other Author (2019): Finding 2
    \end{itemize}
    
    \vspace{1em}
    
    \textbf{Gap we fill:}
    \begin{itemize}
        \item[\textcolor{red}{?}] [Open question our paper addresses]
    \end{itemize}
\end{frame}

% ============================================
% RESEARCH QUESTION (1 slide)
% ============================================

\begin{frame}{This Paper}
    \begin{center}
        \Large
        \textbf{Research Question:} \\[1em]
        Does [X] cause [Y]? \\[2em]
    \end{center}
    
    \textbf{Preview of findings:}
    \begin{itemize}
        \item Main result in plain language
        \item Key magnitude: [Quantitative summary]
    \end{itemize}
\end{frame}

% ============================================
% EMPIRICAL STRATEGY (2-3 slides)
% ============================================

\begin{frame}{Data}
    \textbf{Sources:}
    \begin{itemize}
        \item Dataset 1: [Description, years, N]
        \item Dataset 2: [Description, matching method]
    \end{itemize}
    
    \vspace{1em}
    
    \textbf{Sample:}
    \begin{itemize}
        \item Unit of observation: [What is an observation?]
        \item Final sample: [N] observations, [Time period]
    \end{itemize}
\end{frame}

\begin{frame}{Identification Strategy}
    \textbf{Challenge:} [Endogeneity concern in one sentence]
    
    \vspace{1em}
    
    \textbf{Solution:} We exploit [natural experiment / instrument / RDD]
    
    \vspace{1em}
    
    \textbf{Key assumption:} [Identification assumption in plain language]
    
    \begin{equation*}
        Y_{it} = \alpha + \beta \cdot \text{Treatment}_{it} + \gamma X_{it} + \mu_i + \delta_t + \varepsilon_{it}
    \end{equation*}
\end{frame}

% ============================================
% RESULTS (3-5 slides)
% ============================================

\begin{frame}{Main Result}
    \begin{center}
        \includegraphics[width=0.8\textwidth]{figures/main_result.pdf}
    \end{center}
    
    \vspace{0.5em}
    
    \textbf{Takeaway:} [One sentence interpretation]
\end{frame}

\begin{frame}{Main Result: Regression Table}
    \begin{table}
        \centering
        \small
        \begin{tabular}{lccc}
            \toprule
            & (1) & (2) & (3) \\
            & OLS & + Controls & + FE \\
            \midrule
            Treatment & 0.052*** & 0.048*** & 0.041** \\
                      & (0.012)  & (0.011)  & (0.015) \\
            \midrule
            Controls & No & Yes & Yes \\
            Fixed Effects & No & No & Yes \\
            N & 10,000 & 9,850 & 9,850 \\
            \bottomrule
        \end{tabular}
    \end{table}
    
    \textbf{Economic magnitude:} 1 SD increase in X $\rightarrow$ Y\% increase in outcome
\end{frame}

\begin{frame}{Robustness Checks}
    \begin{itemize}
        \item[\checkmark] Alternative specifications
        \item[\checkmark] Placebo tests
        \item[\checkmark] Different sample cuts
        \item[\checkmark] [Other relevant checks]
    \end{itemize}
    
    \vspace{1em}
    
    $\rightarrow$ Results robust across specifications
\end{frame}

% ============================================
% CONCLUSION (1 slide)
% ============================================

\begin{frame}{Takeaways}
    \begin{enumerate}
        \item \textbf{Finding 1:} [Main result]
        \item \textbf{Finding 2:} [Secondary result]
        \item \textbf{Implication:} [Policy/theory takeaway]
    \end{enumerate}
    
    \vspace{2em}
    
    \begin{center}
        \Large Thank you! \\[0.5em]
        \normalsize your.email@university.edu
    \end{center}
\end{frame}

% ============================================
% APPENDIX
% ============================================

\appendix

\begin{frame}[noframenumbering]{Appendix: Additional Results}
    [Backup slides for Q\&A]
\end{frame}

\end{document}
```

## Theme Recommendations

| Audience | Theme | Notes |
|----------|-------|-------|
| Academic | `metropolis` | Clean, modern, minimal |
| Conference | `Madrid` | Traditional, professional |
| Job market | `default` with custom colors | Safe, customizable |
| Policy | `CambridgeUS` | Authoritative look |

## Best Practices

1. **One message per slide** - if you need more, split it
2. **Use figures over tables** when possible
3. **Highlight key numbers** in results tables
4. **Build complex slides** incrementally with `\pause`
5. **Prepare backup slides** for anticipated questions
6. **Practice timing** - 1-2 minutes per slide max

## Common Pitfalls

- ❌ Too much text on slides
- ❌ Reading slides word-for-word
- ❌ Tables with too many columns
- ❌ Skipping the roadmap/preview
- ❌ Ending with "Questions?" instead of takeaways

## References

- [Shapiro (2019) How to Give Applied Micro Talk](https://www.brown.edu/Research/Shapiro/pdfs/applied_micro_slides.pdf)
- [Beamer User Guide](https://ctan.org/pkg/beamer)
- [Metropolis Theme](https://github.com/matze/mtheme)

## Changelog

### v1.0.0
- Initial release with conference talk template
