---
name: academic-paper-writer
description: Draft economics papers with proper structure and academic style
workflow_stage: writing
compatibility:
  - claude-code
  - cursor
  - codex
  - gemini-cli
author: Awesome Econ AI Community
version: 1.0.0
tags:
  - LaTeX
  - academic-writing
  - papers
  - economics
---

# Academic Paper Writer

## Purpose

This skill helps economists draft, structure, and polish academic papers with proper conventions for economics journals. It provides templates for different paper types and guidance on academic writing style.

## When to Use

- Starting a new research paper from scratch
- Restructuring an existing draft
- Writing specific sections (introduction, literature review, conclusion)
- Preparing papers for journal submission

## Instructions

### Step 1: Identify Paper Type

Ask the user:
1. Is this empirical or theoretical?
2. What is the target journal/audience?
3. What stage is the paper at? (outline, first draft, revision)
4. What sections need help?

### Step 2: Follow the IMRAD Structure

For empirical papers, use:
1. **Introduction** - Motivation, research question, contribution
2. **Literature Review** - Related work and positioning
3. **Data & Methods** - Sources, sample, empirical strategy
4. **Results** - Main findings with tables/figures
5. **Discussion** - Interpretation, mechanisms, limitations
6. **Conclusion** - Summary and implications

### Step 3: Apply Economics Writing Conventions

- **First paragraph** should state the research question and main finding
- **Use present tense** for established facts, past tense for your findings
- **Be precise** with causal language (effect vs. association)
- **Cite heavily** in the literature review
- **Lead with results** in the results section

## Example Output: Introduction Template

```latex
\section{Introduction}

% Hook - Why does this matter?
[TOPIC] is a fundamental question in economics, with implications for 
[POLICY AREA] and [BROADER RELEVANCE]. Despite extensive research, 
we still lack clear evidence on [SPECIFIC GAP].

% Research question
This paper asks: [RESEARCH QUESTION IN PLAIN LANGUAGE]? 
Specifically, we examine whether [PRECISE FORMULATION OF THE QUESTION].

% Preview of answer
We find that [MAIN RESULT IN ONE SENTENCE]. This effect is 
[economically significant / modest / heterogeneous], with 
[QUANTITATIVE SUMMARY: e.g., "a one standard deviation increase 
in X associated with a Y percent increase in Z"].

% Methodology (brief)
To identify this effect, we exploit [IDENTIFICATION STRATEGY: 
natural experiment / RCT / instrumental variable / RDD]. 
Our data come from [DATA SOURCE], covering [TIME PERIOD] 
and [SAMPLE SIZE] observations.

% Contribution / Related literature
Our paper contributes to several strands of literature. 
First, we extend the work of \citet{Author2020} by [EXTENSION]. 
Second, we provide new evidence on [MECHANISM/CHANNEL] that 
complements \citet{OtherAuthor2019}. Finally, our findings 
have implications for [POLICY/FUTURE RESEARCH].

% Roadmap
The remainder of the paper is organized as follows. 
Section~\ref{sec:background} provides background and reviews 
related literature. Section~\ref{sec:data} describes our data 
and empirical strategy. Section~\ref{sec:results} presents our 
main findings. Section~\ref{sec:robustness} discusses robustness 
checks. Section~\ref{sec:conclusion} concludes.
```

## Example Output: Results Section Template

```latex
\section{Results}
\label{sec:results}

% Lead with the main finding
Table~\ref{tab:main} presents our main results. Column (1) shows 
the baseline OLS specification without controls. The coefficient 
on [TREATMENT VARIABLE] is [POINT ESTIMATE] (s.e. = [SE]), 
statistically significant at the [1/5/10] percent level.

% Add controls incrementally
In column (2), we add [CONTROL SET 1]. The point estimate 
[increases/decreases slightly/remains stable] to [ESTIMATE]. 
Column (3) includes [CONTROL SET 2] and adds [FIXED EFFECTS]. 
Our preferred specification in column (4) includes [FULL CONTROLS] 
and yields [FINAL ESTIMATE].

% Interpret magnitude
To gauge economic significance, note that [INTERPRETATION]. 
A one standard deviation increase in [X] is associated with 
a [Y] percent [increase/decrease] in [OUTCOME], or roughly 
[COMPARISON TO MEAN/OTHER BENCHMARK].

% Brief mention of mechanisms/heterogeneity if relevant
Table~\ref{tab:hetero} explores heterogeneity by [DIMENSION]. 
We find that the effect is [larger/concentrated among] 
[SUBGROUP], suggesting that [INTERPRETATION].

\begin{table}[htbp]
\centering
\caption{Main Results: Effect of X on Y}
\label{tab:main}
\begin{tabular}{lcccc}
\hline\hline
 & (1) & (2) & (3) & (4) \\
 & OLS & + Controls & + FE & Preferred \\
\hline
Treatment & 0.052*** & 0.048*** & 0.041** & 0.039** \\
          & (0.012)  & (0.011)  & (0.015) & (0.016) \\
\\
Controls       & No  & Yes & Yes & Yes \\
Fixed Effects  & No  & No  & Yes & Yes \\
Cluster SE     & No  & No  & No  & Yes \\
\\
Observations   & 10,000 & 9,850 & 9,850 & 9,850 \\
R-squared      & 0.05   & 0.12  & 0.35  & 0.35  \\
\hline\hline
\multicolumn{5}{l}{\footnotesize Notes: * p<0.10, ** p<0.05, *** p<0.01.} \\
\multicolumn{5}{l}{\footnotesize Standard errors in parentheses.} \\
\end{tabular}
\end{table}
```

## Example Output: Conclusion Template

```latex
\section{Conclusion}
\label{sec:conclusion}

% Restate question and answer
This paper examined [RESEARCH QUESTION]. Using [METHOD/DATA], 
we found that [MAIN FINDING]. This result is robust to 
[ROBUSTNESS CHECKS].

% Implications
Our findings have several implications. For policy, they suggest 
that [POLICY IMPLICATION]. For theory, they provide support for 
[THEORETICAL MECHANISM] and challenge [ALTERNATIVE VIEW].

% Limitations (brief, honest)
Several limitations warrant mention. First, [LIMITATION 1: 
e.g., external validity]. Second, [LIMITATION 2: e.g., 
data constraints]. Future research could address these by 
[SUGGESTION].

% Future directions
This paper opens several avenues for future work. 
[DIRECTION 1]. [DIRECTION 2]. We hope our findings 
stimulate further research on [BROADER TOPIC].
```

## Writing Tips

### For Introductions
- **First sentence should grab attention** - not "This paper examines..."
- **State your contribution clearly** - what's new about this paper?
- **Be specific about magnitudes** - don't just say "large effect"
- **Acknowledge limitations** preemptively in the last paragraph

### For Results
- **Lead with numbers** - put the coefficient in the first sentence
- **Interpret economically** - what does a 0.05 coefficient mean?
- **Guide the reader** through tables column by column
- **Don't oversell** - distinguish statistical from economic significance

### For Conclusions
- **Don't introduce new results** - synthesize what you've shown
- **Be honest about limitations** - reviewers will find them anyway
- **End on the contribution** - remind readers why this matters

## Common Pitfalls

- ❌ Burying the main result in the middle of the paper
- ❌ Using "significant" without specifying statistical or economic
- ❌ Over-claiming causality without proper identification
- ❌ Literature review that's just a list of papers
- ❌ Conclusion that's just a summary

## References

- [Cochrane (2005) Writing Tips for PhD Students](https://www.johnhcochrane.com/research-all/writing-tips-for-phd-studentsnbsp)
- [Shapiro (2019) How to Give an Applied Micro Talk](https://www.brown.edu/Research/Shapiro/pdfs/applied_micro_slides.pdf)
- [Thomson (2011) A Guide for the Young Economist](https://mitpress.mit.edu/books/guide-young-economist)

## Changelog

### v1.0.0
- Initial release with introduction, results, and conclusion templates
