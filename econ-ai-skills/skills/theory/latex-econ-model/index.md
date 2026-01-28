---
name: latex-econ-model
description: Write and typeset economic models in LaTeX with proper notation
workflow_stage: theory
compatibility:
  - claude-code
  - cursor
  - codex
  - gemini-cli
author: Awesome Econ AI Community
version: 1.0.0
tags:
  - LaTeX
  - economic-models
  - theory
  - equations
---

# LaTeX Economic Model Writer

## Purpose

This skill helps economists write and typeset economic models in LaTeX with proper mathematical notation, consistent formatting, and academic conventions. It covers utility maximization, equilibrium conditions, dynamic programming, and game theory models.

## When to Use

- Writing the theory section of a paper
- Formalizing an economic model
- Creating lecture notes with equations
- Typesetting optimization problems and equilibrium conditions

## Instructions

### Step 1: Understand the Model Type

Ask the user:
1. What type of model? (Static optimization, dynamic programming, game theory, general equilibrium)
2. What are the agent types? (Consumer, firm, government, etc.)
3. What constraints apply?
4. Should it include proofs or just the setup?

### Step 2: Use Standard Notation

Follow economics conventions:
- Utility: $U$, $u$, or $v$
- Production: $F$, $f$, or $Y$
- Prices: $p$, $w$, $r$
- Quantities: $q$, $x$, $y$
- Time subscripts: $t$, $t+1$
- Greek letters for parameters: $\alpha$, $\beta$, $\gamma$, $\delta$

### Step 3: Structure the Model

Organize as:
1. **Environment** - Agents, goods, time horizon
2. **Preferences/Technology** - Utility functions, production functions
3. **Constraints** - Budget, resource, incentive compatibility
4. **Equilibrium concept** - Definition and conditions
5. **Solution** - First-order conditions, closed-form solutions

## Example Output

```latex
\documentclass{article}
\usepackage{amsmath, amssymb, amsthm}
\usepackage{mathtools}

% Theorem environments
\newtheorem{definition}{Definition}
\newtheorem{proposition}{Proposition}
\newtheorem{lemma}{Lemma}

% Custom commands for economics
\newcommand{\E}{\mathbb{E}}          % Expectation
\newcommand{\R}{\mathbb{R}}          % Real numbers
\newcommand{\pd}[2]{\frac{\partial #1}{\partial #2}}  % Partial derivative

\begin{document}

\section{A Simple Consumer Problem}

\subsection{Environment}

Consider a consumer who lives for two periods, $t \in \{1, 2\}$. The consumer has preferences over consumption $c_t$ represented by the utility function:
%
\begin{equation}
    U(c_1, c_2) = u(c_1) + \beta u(c_2)
\end{equation}
%
where $\beta \in (0,1)$ is the discount factor and $u(\cdot)$ is strictly increasing and strictly concave.

\subsection{Constraints}

The consumer earns income $y_1$ in period 1 and $y_2$ in period 2. She can save at gross interest rate $R = 1 + r$. The budget constraints are:
%
\begin{align}
    c_1 + s &= y_1 \label{eq:bc1}\\
    c_2 &= y_2 + Rs \label{eq:bc2}
\end{align}
%
where $s$ denotes savings. Combining \eqref{eq:bc1} and \eqref{eq:bc2} yields the intertemporal budget constraint:
%
\begin{equation}
    c_1 + \frac{c_2}{R} = y_1 + \frac{y_2}{R} \equiv W
\end{equation}

\subsection{Optimization Problem}

The consumer solves:
%
\begin{equation}
    \max_{c_1, c_2} \quad u(c_1) + \beta u(c_2) 
    \quad \text{s.t.} \quad c_1 + \frac{c_2}{R} = W
\end{equation}

\subsection{Solution}

The Lagrangian is:
%
\begin{equation}
    \mathcal{L} = u(c_1) + \beta u(c_2) + \lambda\left(W - c_1 - \frac{c_2}{R}\right)
\end{equation}

First-order conditions:
%
\begin{align}
    \pd{\mathcal{L}}{c_1} &= u'(c_1) - \lambda = 0 \\
    \pd{\mathcal{L}}{c_2} &= \beta u'(c_2) - \frac{\lambda}{R} = 0
\end{align}

Combining these yields the \textbf{Euler equation}:
%
\begin{equation}
    \boxed{u'(c_1) = \beta R \cdot u'(c_2)}
\end{equation}

\begin{proposition}[Consumption Smoothing]
If $\beta R = 1$, then $c_1^* = c_2^*$ (perfect consumption smoothing).
\end{proposition}

\begin{proof}
When $\beta R = 1$, the Euler equation becomes $u'(c_1) = u'(c_2)$. Since $u$ is strictly concave, $u'$ is strictly decreasing, which implies $c_1 = c_2$.
\end{proof}

%====================================
\section{A Firm's Dynamic Problem}
%====================================

Consider a firm that maximizes the present value of profits:
%
\begin{equation}
    \max_{\{k_{t+1}, n_t\}_{t=0}^{\infty}} \sum_{t=0}^{\infty} \beta^t \left[ F(k_t, n_t) - w_t n_t - I_t \right]
\end{equation}
%
subject to the capital accumulation equation:
%
\begin{equation}
    k_{t+1} = (1 - \delta) k_t + I_t
\end{equation}

The Bellman equation is:
%
\begin{equation}
    V(k) = \max_{k', n} \left\{ F(k, n) - wn - k' + (1-\delta)k + \beta V(k') \right\}
\end{equation}

\end{document}
```

## LaTeX Packages to Include

```latex
% Essential packages for economics papers
\usepackage{amsmath}    % Enhanced math environments
\usepackage{amssymb}    % Mathematical symbols
\usepackage{amsthm}     % Theorem environments
\usepackage{mathtools}  % Extensions to amsmath
\usepackage{bm}         % Bold math symbols
\usepackage{dsfont}     % \mathds for indicator functions
```

## Useful Custom Commands

```latex
% Expectation and probability
\newcommand{\E}{\mathbb{E}}
\newcommand{\Var}{\text{Var}}
\newcommand{\Cov}{\text{Cov}}
\newcommand{\Prob}{\mathbb{P}}

% Indicator function
\newcommand{\ind}{\mathds{1}}

% Partial derivatives
\newcommand{\pd}[2]{\frac{\partial #1}{\partial #2}}
\newcommand{\pdd}[2]{\frac{\partial^2 #1}{\partial #2^2}}

% Argmax/argmin
\DeclareMathOperator*{\argmax}{arg\,max}
\DeclareMathOperator*{\argmin}{arg\,min}

% Blackboard bold
\newcommand{\R}{\mathbb{R}}
\newcommand{\N}{\mathbb{N}}
\newcommand{\Z}{\mathbb{Z}}
```

## Best Practices

1. **Use `align` environment** for multiline equations
2. **Label important equations** with `\label{}` and reference with `\eqref{}`
3. **Use `\text{}` for words in equations** (not bare text)
4. **Box key results** with `\boxed{}`
5. **Define custom commands** for repeated notation
6. **Use consistent subscript conventions** ($t$ for time, $i$ for individuals)

## Common Pitfalls

- ❌ Using `*` for multiplication (use `\cdot` or implicit multiplication)
- ❌ Forgetting `\left(` and `\right)` for auto-sizing brackets
- ❌ Inconsistent notation across the paper
- ❌ Not aligning equations at `=` signs
- ❌ Using `$$ ... $$` instead of proper environments

## References

- [AMS-LaTeX User's Guide](https://www.ams.org/arc/resources/amslatex-about.html)
- [The Not So Short Introduction to LaTeX](https://tobi.oetiker.ch/lshort/lshort.pdf)
- [Mathpix](https://mathpix.com/) - Convert handwritten equations to LaTeX

## Changelog

### v1.0.0
- Initial release with consumer, firm, and game theory templates
