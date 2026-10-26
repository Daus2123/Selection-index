# Workshop slide guide: selection pipeline for single and multiple locations

**Audience:** breeders, trial coordinators, and analysts who will use the Shiny app.  
**Length:** about 2 hours, including two demonstrations and one group exercise.  
**Outcome:** participants can prepare a trial workbook, explain the main result tables, and defend an ADVANCE, RETEST, or DISCARD recommendation.

This is a slide outline with speaker notes and visual directions. Use synthetic trial data in the demonstrations. Keep one set of example genotypes throughout: show them at one location first, then show the same entries across several environments. Replace the example trait names with traits familiar to the team.

## Suggested teaching example

Use 8–12 fictional genotypes, 2–3 replications, and three traits such as yield, disease score, and fruit quality. Give the traits different units and priorities: yield is **higher better**, disease score is **lower better**, and fruit quality has a moderate weight. Mark one or two entries as checks. For the multi-location example, use several environments and make one genotype strong at one location but variable elsewhere. This creates a useful discussion about performance versus stability. Add a missing genotype–environment combination to demonstrate the presence and uncertainty review. Avoid real breeding data in shared training files.

## Part 1 — Introduction: what the pipeline does (slides 1–5, 20 minutes)

### Slide 1 — Title and workshop objective (2 minutes)

**On slide:** “From trial data to selection decisions” and the two paths, Single-Location Trial and Multi-Environment Trial.

**Explain:** The app helps the team move from plot observations to a traceable genotype recommendation. The decision is the end of a chain: input data, statistical estimates, trait priorities, ranking, then decision rules. Participants should leave knowing where each number came from and when they should question it. State that the workshop uses a fictional trial, so everyone can discuss decisions freely.

**Visual:** A simple horizontal diagram: Excel workbook → analysis → evidence tables and charts → ADVANCE / RETEST / DISCARD. Put the two analysis paths inside the analysis step.

### Slide 2 — The breeding decision behind the analysis (4 minutes)

**On slide:** “Which entries deserve the next round of testing?”

**Explain:** In one location, the practical question is which genotypes performed well against the trial's priorities and checks. Across locations, the question expands: did they perform well consistently, or only in a favorable environment? A high rank is evidence for discussion, not an automatic release or advancement instruction. The breeder still considers trial quality and the program's objectives.

**Visual:** Two small field maps. The left shows one site with replicated plots. The right shows several sites with the same genotype labels. Highlight one genotype whose rank changes across sites.

**Ask the room:** “Would you advance an entry that tops one site but is poor at the others?” Use answers to introduce stability.

### Slide 3 — One app, two analysis paths (4 minutes)

**On slide:** A comparison table with columns for Single Location and Multi Location. Rows: input, main estimate, extra evidence, final ranking.

**Explain:** Both paths use the same general workflow: Data → Analyze → Results and Charts → Export. Single Location combines traits into a weighted selection index for one trial. Multi Location models genotype and environment effects, adds stability evidence, then integrates traits across environments. The outputs share the ADVANCE / RETEST / DISCARD vocabulary, making decisions easier to compare.

**Visual:** Two branching flows that rejoin at a decision board. Avoid showing every statistical method yet; this slide is the orientation map.

### Slide 4 — What the workbook must contain (6 minutes)

**On slide:** An annotated excerpt from the synthetic Excel workbook.

**Explain:** Each observation belongs to a genotype and a replication. Trait columns hold measurements. Multi Location also requires an `Environment` column. An optional entry-type column can identify checks. Special weight and direction rows tell the app which traits matter most and whether higher, lower, or a target value is desirable. Emphasize that direction is part of the breeding goal: a lower disease score must contribute favorably, even though its numeric value is smaller. The app reads the first worksheet, so the teaching example should put the data there.

**Visual:** Color-code ID, replication, environment, trait measurements, check label, weight row, and direction row in a short spreadsheet screenshot. Add a callout with one example of each direction.

**Teaching check:** Ask participants to spot one error in a deliberately flawed copy, such as a misspelled environment or a disease trait marked “Higher better.”

### Slide 5 — The common workflow and quality gate (4 minutes)

**On slide:** Data → Analyze → Results → Charts → Export, with “review input and notes” under the first and third steps.

**Explain:** Uploading is only the start. Confirm that the trial design, replication, checks, trait names, directions, and weights reflect the study. After running, read quality-control notes before the ranking. Results tables and charts answer different questions: tables give exact estimates and decisions; charts show patterns. Exports provide a record for the selection meeting. Explain that changing certain decision settings updates rankings without uploading the file again; a new trial or changed raw observations requires a new upload and run.

**Visual:** App navigation screenshot with numbered markers. If the live interface changes, use a current screenshot rather than an old mockup.

## Part 2 — Single-Location Trial (slides 6–11, 35 minutes)

### Slide 6 — Single Location: from plot data to trait evidence (5 minutes)

**On slide:** “First understand each trait.”

**Explain:** The app analyzes each trait with the selected trial model. Single Location supports CRD, RCBD, and LMM; the choice should match the design, and RCBD or LMM requires replication information. ANOVA or model results help assess whether observed genotype differences are distinguishable from trial variation. Mean comparisons, where requested, show which entries differ under LSD or Tukey grouping. Explain letters carefully: sharing a letter means the comparison does not establish a difference at the selected threshold. It does not prove the entries are biologically identical.

**Visual:** Dot plot of genotype means for one trait, with uncertainty intervals if available and compact letter groups. Mark the check entry. Put a small callout for the model selected in the demo.

### Slide 7 — Check-relative performance (4 minutes)

**On slide:** “How much better or worse than the check?”

**Explain:** A check anchors the analysis to an entry the team knows. The pipeline can use manually selected checks, detected check labels, or a fallback when none is supplied. During training, choose checks explicitly and show their names on screen. Discuss superiority over checks as a practical comparison, while reminding participants that one favorable comparison on one trait does not settle the full multi-trait decision.

**Visual:** Horizontal bars of trait difference from the check, with zero as a clear reference line. Use a second small table showing the selected checks.

**Demo cue:** Change the selected check and ask which comparisons change.

### Slide 8 — Why standardize, weight, and set direction (7 minutes)

**On slide:** An example with yield measured in large units and disease score on a small scale.

**Explain:** Raw traits cannot be added fairly when their units differ. The pipeline adjusts each trait to a common scale, orients it so desirable performance contributes positively, and applies the trait weights. A higher weight means that trait has more influence on the combined score. A target trait rewards closeness to the target rather than simply a high or low raw number. The weighted contributions table is the best place to explain *why* a genotype reached its index rank. Briefly point out that trait correlations help the team see whether traits tend to move together or conflict.

**Visual:** A three-stage worked example for two genotypes: raw values → direction-adjusted standardized values → weighted contributions. Show numbers as illustrative teaching values, not claimed app output. A stacked contribution bar is useful for the final stage.

**Ask the room:** “What happens if disease resistance becomes twice as important?” Then reveal the changed contribution or ranking in the demo.

### Slide 9 — Ranking, genetic gain, and decision rules (7 minutes)

**On slide:** “A rank tells us order; a decision rule tells us action.”

**Explain:** The combined index ranks the entries. Heritability and expected genetic gain provide context for how useful selection on a trait may be, given the data and assumptions. The final decision uses index thresholds plus priority-trait rules; an entry with a severe weakness in an important trait may need extra caution even if its overall index is attractive. Explain the program's chosen cutoffs before showing the labels. Show the exact reason or supporting evidence for any entry that prompts discussion.

**Visual:** A ranked index plot with decision colors, paired with a small decision table and a single highlighted entry's trait contributions. Use consistent colors for ADVANCE, RETEST, and DISCARD throughout the deck.

### Slide 10 — Single Location walkthrough (8 minutes)

**On slide:** Five demo checkpoints: upload, settings, trait evidence, index, export.

**Explain while demonstrating:** Upload the synthetic workbook in Data. Choose Single-Location Trial in Analyze, confirm the model and checks, then run. In Results, visit the trait ANOVA and mean comparison, index ranking, weighted contributions, heritability/gain, and final decision. In Charts, show a trait comparison and ranking chart. Change one check or decision cutoff and point out how the displayed decision responds. Download the analysis workbook and identify the sheets that support the recommendation.

**Visual:** Live app. Keep screenshots of the five checkpoints as a fallback if the app or internet is unavailable.

### Slide 11 — Single Location interpretation exercise (4 minutes)

**On slide:** Three synthetic genotypes with their index rank, check comparison, and one priority-trait concern.

**Explain:** Give participants one minute to choose an entry to advance and one to retest. Ask them to cite a trait result and the contribution or decision table, not only the rank. Summarize the discussion by distinguishing observed performance, statistical confidence, and breeding preference.

**Visual:** A compact evidence card for each genotype, with the recommendation initially hidden.

## Part 3 — Multi-Environment Trial (slides 12–18, 45 minutes)

### Slide 12 — Why the location effect matters (5 minutes)

**On slide:** “The same genotype can respond differently across environments.”

**Explain:** Environments can differ in climate, management, pressure from disease, or other conditions. A genotype's average is useful, but it can hide uneven performance. Genotype × environment interaction describes changing relative performance across environments. The purpose of the Multi Location path is to estimate performance while showing the strength and pattern of this interaction. Keep the biological interpretation tied to the trial: an entry can have broad adaptation or be promising for particular target environments.

**Visual:** A line chart with environment on the x-axis and trait performance on the y-axis. Show one consistently good genotype and one with a high peak and a low trough. Label the lines directly.

### Slide 13 — Data cleaning, coverage, and checks (6 minutes)

**On slide:** “Know what data reached the model.”

**Explain:** The MET pipeline first checks potential outliers within genotype–environment cells and then within environments. The before/after plots and row-count summary show what was removed. Removal deserves review: a genuine extreme observation can carry biological information. Next, inspect the genotype × environment presence table. Entries observed at only one environment have limited evidence for broad adaptation and are flagged with higher uncertainty. The app uses chosen checks or a fallback based on entries present in every environment, and its notes explain such gaps and fallback choices. Do not silently treat missing sites as poor performance.

**Visual:** Three panels: before/after distribution, outlier count table, and presence heatmap. Use a visible “missing” pattern distinct from low trait performance.

**Demo cue:** Point to one entry with a missing environment and read the corresponding note.

### Slide 14 — Model choice and the meaning of estimates (7 minutes)

**On slide:** A model-choice table: Standard LMM, Breeding LMM (`sommer`), and joint ANOVA (RCBD).

**Explain:** The app offers different engines to fit the multi-environment design. The mixed-model paths report BLUPs for genotype effects; the joint ANOVA path reports adjusted means, also called BLUEs. These labels identify different kinds of estimates and should stay visible when results are exported or discussed. The model incorporates genotype, environment, their interaction, and design terms as applicable. Encourage analysts to review model diagnostics and variance components before interpreting small rank differences. For the workshop, choose one primary model for the live run and explain that the alternative exists for the team's statistical workflow.

**Visual:** A plain-language table with “When used in this app,” “Estimate label,” and “What to inspect.” Add a small diagram showing genotype, environment, and interaction as contributors to observed values.

### Slide 15 — Mean performance and stability evidence (8 minutes)

**On slide:** “Performance and stability answer related questions.”

**Explain:** The genotype estimate summarizes performance across environments. Finlay–Wilkinson regression describes how each genotype responds as environments become more or less favorable. AMMI separates the main effects from interaction patterns, and its stability value can help compare consistency. GGE charts show genotype and genotype × environment patterns, including mean versus stability and which genotypes perform best in groups of environments. A biplot is a compressed view, so read it with its explained-variance information and the numeric tables. The app only produces the AMMI/GGE views when enough environments are available under the selected minimum; absence of a biplot is not a failure of the core model.

**Visual:** Three small panels using the same synthetic trait: response lines, AMMI biplot, and GGE mean-versus-stability plot. Highlight one genotype consistently across all panels. Add a short caption under each stating the question it answers.

### Slide 16 — From per-trait ranking to integrated ranking (7 minutes)

**On slide:** “One trait at a time, then all traits together.”

**Explain:** For each trait, the pipeline combines the mean estimate with stability components. The breeder can adjust the mean, Finlay–Wilkinson, and AMMI stability value weights. Then the cross-trait step combines the selected traits according to their trait weights, directions, and targets. This produces the integrated ranking used for the final decision board. Explain that changing component weights or selection percentage changes the decision view without refitting every trait model; changing raw data or the model choice requires another analysis run.

**Visual:** A layered diagram: each trait's performance and stability evidence → per-trait ranking → integrated ranking. Beside it, show a small “what if” example where increasing stability weight swaps two genotype ranks. Label the example as illustrative.

### Slide 17 — Decision board and target environments (6 minutes)

**On slide:** “Read the reason alongside the recommendation.”

**Explain:** The decision board applies the chosen selection intensity to the integrated ranking and reports ADVANCE, RETEST, or DISCARD. Show the rule thresholds beside the board so the meeting can see how many entries were selected and why. Review stability class, location coverage, and possible target environments. For an entry with limited coverage or conflicting evidence, RETEST can be the most useful outcome because it identifies the next experiment needed. Avoid presenting the board as proof that a genotype is universally superior.

**Visual:** A cropped decision board with one broad candidate, one location-specific candidate, and one uncertain candidate highlighted. Add a small map or environment list only if the actual trial has meaningful site labels.

### Slide 18 — Multi Location walkthrough and comparison (6 minutes)

**On slide:** Demo checklist: outlier review, coverage, model, stability, integrated ranking, decision.

**Explain while demonstrating:** Upload the multi-environment synthetic workbook. Select traits and a MET model, confirm checks and the minimum number of environments for biplots, and run. Read pipeline notes and QC first. Then show one trait's genotype estimates, a stability view, the integrated ranking, and the decision board. Compare the same example genotype's Single Location decision with its Multi Location decision and ask what new evidence caused the difference. Export the workbook as the record of the discussion.

**Visual:** Live app plus a prepared side-by-side comparison table. Keep static screenshots as backup.

## Part 4 — Practice and close (slides 19–20, 20 minutes)

### Slide 19 — Team selection exercise (15 minutes)

**On slide:** “Select two entries to advance and one to retest. Explain each choice.”

**Explain:** Divide participants into small groups. Give each group a one-page extract from the synthetic MET output: genotype estimates, location coverage, stability evidence, integrated ranking, decision rule table, and notes. Ask for a concise recommendation with evidence and one remaining uncertainty. Have groups report back briefly. If recommendations differ, use the difference to discuss breeding priorities and uncertainty instead of forcing agreement from rank alone.

**Visual:** A blank evidence card with fields for genotype, performance, stability, coverage, priority traits, recommendation, and next trial action.

### Slide 20 — Close: a repeatable decision routine (5 minutes)

**On slide:** Five questions: Is the workbook correct? Is the model suitable? What do trait results say? How stable is performance? Why did the decision rule assign this action?

**Explain:** Recap the sequence participants should use after every run. The Single Location path establishes a within-trial ranking. The Multi Location path revisits the same material across environments and adds stability. The best final output is a recommendation the team can trace to exported evidence. Point participants to the Results, Charts, and Export areas for their next real trial.

**Visual:** A compact checklist on top of the original workbook-to-decision flow, closing the loop with Slide 1.

## Presenter preparation

1. Prepare two synthetic Excel workbooks: one for a single location and one with an `Environment` column for several locations. Keep genotype names and trait meanings consistent across them.
2. Run both analyses before the workshop. Save screenshots and exports for every live-demo checkpoint. Confirm which charts appear for the chosen number of environments.
3. Decide the example model, checks, trait weights, directions, and decision settings in advance. Put those values on the relevant slides so participants can reproduce the run.
4. Make one slide or handout defining key terms in plain language: check, ANOVA, adjusted mean, BLUP, genotype × environment interaction, and stability. Use it as a backup reference rather than interrupting every slide with definitions.
5. Give groups an exported subset of the synthetic results for Slide 19, including notes and thresholds. The exercise should be solvable without access to the live app.

## Source notes

This guide follows the current project descriptions in [ARCHITECTURE.md](ARCHITECTURE.md), [PRD.md](PRD.md), and [README.md](../README.md). Check the current app labels and settings before converting this guide into final slides, especially if the app has changed since the training screenshots were captured.
