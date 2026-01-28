---
name: research-ideation
description: Generate research questions from economic phenomena
workflow_stage: ideation
compatibility:
  - claude-code
  - cursor
  - codex
  - gemini-cli
author: Awesome Econ AI Community
version: 1.0.0
tags:
  - research
  - ideation
  - hypothesis
  - brainstorming
---

# Research Ideation

## Purpose

This skill helps economists generate and refine research questions by applying economic thinking frameworks. It guides the process from observing phenomena to formulating testable hypotheses.

## When to Use

- Starting a new research project
- Brainstorming dissertation topics
- Finding angles for policy analysis
- Developing extensions of existing work

## Instructions

### Step 1: Understand the Starting Point

Ask the user:
1. Do you have a phenomenon, puzzle, or data source in mind?
2. What field/subfield of economics?
3. What's your methodological preference? (Theory, empirical, structural)
4. Any constraints? (Data access, timeline, etc.)

### Step 2: Apply Ideation Frameworks

Use these approaches:

**1. The Puzzle Approach**
- What's surprising about current patterns?
- What contradicts conventional economic wisdom?
- Where do standard models fail?

**2. The Policy Approach**
- What policies lack rigorous evaluation?
- What natural experiments remain unexploited?
- What interventions might solve important problems?

**3. The Data Approach**
- What new data sources have become available?
- What can existing data tell us that hasn't been explored?
- What linkages between datasets are possible?

**4. The Extension Approach**
- How can seminal papers be extended?
- What mechanisms remain unexplored?
- Can methods from one field apply to another?

### Step 3: Evaluate and Refine

For each idea, assess:
- **Feasibility**: Can this be done with available data/methods?
- **Contribution**: What's new here?
- **Interest**: Who cares about this question?
- **Identification**: Can effects be credibly estimated?

## Example Interaction

**User Query:** "I'm interested in labor economics and have access to LinkedIn data"

**Generated Response:**

## Research Ideas: LinkedIn Data × Labor Economics

### 1. Network Effects in Job Search
**Question:** How do professional networks affect job transitions and wage changes?

**Puzzle:** Standard search models assume random matching, but most jobs come through networks. How much do networks matter for outcomes?

**Approach:**
- Track profile updates signaling job changes
- Measure network overlap between current/future employers
- Instrument with exogenous network shocks (school alumni connections)

**Contribution:** Quantify the causal role of networks vs. unobserved ability

---

### 2. Skill Signaling and Returns
**Question:** Do workers who list specific skills on profiles earn wage premiums?

**Puzzle:** Are listed skills signals of ability, or just cheap talk? What's the return to skill acquisition vs. skill signaling?

**Approach:**
- Compare wages before/after workers add skills
- Use skill endorsements as validation
- Exploit skill assessment features as quasi-experiments

---

### 3. Gender Gaps in Self-Promotion
**Question:** Do men and women describe equivalent achievements differently?

**Puzzle:** Lab evidence shows women understate accomplishments. Does this appear in real profiles and affect outcomes?

**Approach:**
- NLP analysis of profile language (adjectives, quantification)
- Compare language for workers in same role/company
- Link to recruiter contact rates if available

---

### 4. Remote Work and Labor Market Geography
**Question:** How has remote work changed the geographic reach of job matching?

**Approach:**
- Compare pre/post-2020 job applications by distance
- Measure changes in location signals on profiles
- Study wage convergence across cities

---

## Evaluation Matrix

| Idea | Data Feasibility | Identification | Policy Relevance |
|------|-----------------|----------------|------------------|
| Network effects | High | Medium (need IV) | High |
| Skill signaling | High | Medium | Medium |
| Gender language | High | High (descriptive) | High |
| Remote geography | High | High (COVID natural experiment) | High |

## Recommended Next Steps

1. **For empirical work:** Start with Idea #4 (clearest identification)
2. **For policy impact:** Idea #1 has most implications for job matching platforms
3. **For methodological novelty:** Idea #3 combines NLP with quasi-experimental methods

---

## Frameworks for Generating Questions

### The "5 Whys" for Economics

Start with an observation and drill down:
1. GDP growth is slowing → Why?
2. Productivity is stagnant → Why?
3. Investment is low → Why?
4. Uncertainty is high → Why?
5. Policy is unpredictable → **Testable: Does policy uncertainty cause low investment?**

### The "What If" Generator

- What if [policy X] were implemented?
- What if [technology Y] became widespread?
- What if [assumption Z] were relaxed?

### The Cross-Field Pollinator

Take a method from one field and apply to another:
- IO techniques → Labor markets (how do firms choose wages like prices?)
- Finance models → Education (returns to schooling as asset pricing)
- Macro shocks → Micro outcomes (firm-level effects of exchange rate changes)

## Common Pitfalls

- ❌ Questions that are too broad ("What causes inequality?")
- ❌ Questions without clean identification ("Does education cause income?")
- ❌ Questions without data ("Were medieval peasants happy?")
- ❌ Questions already well-answered

## References

- [Shapiro (2022) How to Get Started on Research in Economics](https://www.brown.edu/Research/Shapiro/pdfs/research.pdf)
- [Angrist & Pischke on Mostly Harmless research design](https://www.mostlyharmlesseconometrics.com/)
- [AEA Research Pipelines](https://www.aeaweb.org/rfe/)

## Changelog

### v1.0.0
- Initial release with ideation frameworks
