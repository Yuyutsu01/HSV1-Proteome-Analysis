# Phase 15 Methodology: Biological and Contextual Interpretation

## 1. Scientific Purpose & Core Question
Phase 15 provides literature-anchored biological contextualization for the computational findings established across Phases 11–14.
Core Question:
> *"What biological context can be used to interpret the computational representation, clustering, classification, and error-analysis findings?"*

## 2. Frozen Computational Evidence
- **Dataset:** Exactly 74 HSV-1 strain 17 proteins (RefSeq `NC_001806.2`).
- **Class Balance:** Immediate-Early (IE) = 5, Early = 15, Late = 54.
- **Authoritative Annotations:** `data/annotations/temporal_annotations_final.csv`.
- **Zero Model Retraining / Zero Label Modification:** All computational metrics remain frozen from Phases 11–14.

## 3. Evidence Hierarchy & Tripartite Distinction Rule
To ensure scientific rigor, all discussions strictly distinguish between:
1. **DOCUMENTED FACT:** Experimentally demonstrated molecular or biochemical property documented in peer-reviewed literature (e.g. UL36 is a 3,139 aa tegument protein).
2. **COMPUTATIONAL OBSERVATION:** Quantifiable model behavior or geometric metric (e.g. UL36 had out-of-fold accuracy = 0% and mean margin = 0.333).
3. **INTERPRETIVE HYPOTHESIS:** Plausible biological or methodological hypothesis linking sequence properties to model behavior (e.g. UL36's extreme size and multidomain architecture place it as a representation outlier).

## 4. Summary of Key Contextual Interpretations
### A. Functional Heterogeneity Within Temporal Classes
- Temporal classes are not functionally monolithic. Late proteins encompass enzymatic terminases (UL15), kinases (UL13), transcription factors (VP16/UL48), envelope glycoproteins (gB, gD, gH/gL), and structural capsid components (VP5/UL19).
### B. Difficult Proteins Case Studies (UL36 & US12)
- **UL36 (Large Tegument Protein):** 3,139 aa giant structural hub. Classified as Moderate-Margin Incorrect (Delta_p = 0.333). Its multi-functional domain architecture bridges capsid, motors, and tegument.
- **US12 (ICP47):** 88 aa small peptide blocking host TAP. Classified as High-Confidence Incorrect (Delta_p = 0.669). Lacks classical transactivation domains of canonical IE factors.
### C. Representation Disagreement Context
- 30/74 proteins exhibited disagreement between Physicochemical and ProtBERT models. ProtBERT captured complex domain topologies for multi-protein complexes, whereas Physicochemical descriptors captured global physical properties (length, pI, hydrophobicity).

## 5. Prohibited Claims & Boundaries
- Computational misclassification does NOT imply misannotation.
- Sequence associations do NOT establish molecular causal mechanisms.
- No therapeutic efficacy or clinical translation claims are made.
