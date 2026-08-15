# CLAUDE.md — Mycelium Living Repository

This repository is a **mycelium-enabled living repository**. It carries its own memory and grows smarter over time through structured traces of every action.

## Quick Orientation

1. **Read `.living/INDEX.md` first** — A one-screen knowledge map: tag clusters, most-recent entries, and a tag → entry-ID inverted index. Drill into the underlying file (`learnings.md`, `decisions.md`, `conventions.md`) only when a row matches your task. The SessionStart hook regenerates this index every fresh session — trust it.
2. **For targeted lookup:** `python3 skills/core/scripts/recall_lessons.py --living-dir .living/ --tag <tag>` — fetches only the matching entries instead of pulling the whole file. Also accepts `--id L-42`, `--since YYYY-MM-DD`.
3. **Read `ENVIRONMENTS_INSTALLATIONS.md`** — Environment setup, dependencies, and installation gotchas.
4. **Read the relevant manifest** — Each top-level directory has a descriptive manifest (`ANALYSIS_MANIFEST.md`, `DATA_MANIFEST.md`, `ALGORITHM_MANIFEST.md`, `REFERENCE_MANIFEST.md`).

## Installed Convention Packs

<!-- Updated by install_convention.py -->
<!-- Check .living/conventions/ACTIVE_CONVENTIONS.yaml for the full list -->

### Core (auto-installed)

- **robust-analysis** — Defensive execution practices: strict error handling, data validation checks, parameter sensitivity sweeps, null hypothesis testing, adversarial self-challenge. See `.living/conventions/robust-analysis/analysis-conventions.md` for the entry point.
- **report-generator** — Structured LaTeX PDF report generation. See `.living/conventions/report-generator/analysis-conventions.md` for the workflow.
- **idea-generator** — Persona-based creative ideation for new analysis directions. See `.living/conventions/idea-generator/analysis-conventions.md` for the entry point.

### Domain (opt-in)

- **bioinformatics** — RNA-seq, single-cell, and genomics workflow conventions. See `.living/conventions/bioinformatics/analysis-conventions.md`.
- **image-analysis** — Segmentation, microscopy, and quantification conventions. See `.living/conventions/image-analysis/analysis-conventions.md`.
- **skill-bridge** — Routes to external skill repos cloned in `skillpacks/` (scientific-agent-skills, bioSkills, Autonomous-Science). See `.living/conventions/skill-bridge/analysis-conventions.md`.

## Repository Structure

```
├── algorithms/         — Reusable computational methods (see ALGORITHM_MANIFEST.md)
├── analysis/           — Analytical work (see ANALYSIS_MANIFEST.md)
├── data/               — Data assets: raw (immutable), processed, metadata (see DATA_MANIFEST.md)
├── reference_material/ — External references (see REFERENCE_MANIFEST.md)
├── skillpacks/         — Inert clones of external skill repos (scientific-agent-skills, bioSkills)
├── todo/               — Future work items and ideas (see todo/TODO_REGISTRY.md)
└── .living/            — Repository memory layer
    ├── INDEX.md                    — Auto-regenerated knowledge map (read first)
    ├── decisions.md
    ├── learnings.md
    ├── conventions.md              — Repo-specific overrides
    ├── conventions/                — Installed convention packs
    │   ├── ACTIVE_CONVENTIONS.yaml
    │   ├── robust-analysis/
    │   ├── report-generator/
    │   ├── idea-generator/
    │   ├── bioinformatics/
    │   ├── image-analysis/
    │   └── skill-bridge/
    ├── generated-conventions/      — Conventions crystallized from learnings
    ├── log/                        — Append-only event log
    │   └── LOG_REGISTRY.md
    ├── findings/                   — Scientific findings by topic
    │   ├── FINDINGS_REGISTRY.md
    │   └── {topic-slug}.md
    └── outputs/                        — Derived reports and transfer logs
        └── knowledge-transfers/        — Cross-project transfer audit trail
```

## Workflow

### Before Starting Work

1. **Open `.living/INDEX.md`** — its tag clusters and recent-entries lists tell you which past learnings/decisions are likely relevant. Don't skip this — it's already loaded into context by the SessionStart hook, but the agent should re-read for full detail when starting non-trivial work.
2. **Drill in selectively** — for entries that look relevant, either:
   - Fetch just the entries: `python3 skills/core/scripts/recall_lessons.py --living-dir .living/ --tag <tag>` or `--id L-42`
   - Read the whole file (`learnings.md`/`decisions.md`) only if you need broader context
3. Read the manifest for the area you'll be working in (e.g., `ANALYSIS_MANIFEST.md`)
4. Check `.living/conventions.md` for any repo-specific overrides
5. If a domain convention is active, read its conventions in `.living/conventions/[domain]/`

### While Working

- Follow analysis conventions: every analysis gets its own folder with UPPER_SNAKE_CASE.md documentation, scripts, outputs, reports
- Follow statistical conventions: report effect sizes, confidence intervals, document assumptions
- **Follow robust-analysis conventions** (`.living/conventions/robust-analysis/`): fail loudly on unexpected data, assert shapes/types/ranges, log row counts at every step, run sensitivity analyses for every decision, test null hypotheses via permutation/bootstrap
- **Do not subset data** without explicit user confirmation and justification
- Use Python scripts for reproducible pipelines; every analysis must have a `run.sh` or `run.py` that reproduces final outputs
- This repo uses **pixi** for environment/dependency management and a **DuckDB** database for measured phenotypes — see `ENVIRONMENTS_INSTALLATIONS.md` and `scripts/db/`

### After Every Significant Action (Post-Action Hook Protocol)

**This is critical.** After completing any significant step:

1. **Update manifests**: Update the relevant manifest (`ANALYSIS_MANIFEST.md`, `DATA_MANIFEST.md`, etc.) with new/changed entries
2. **Update documentation**: Update or create the UPPER_SNAKE_CASE.md file in the affected subfolder
3. **Log decisions**: If a non-obvious choice was made, append to `.living/decisions.md`
4. **Log learnings**: If something unexpected happened, append to `.living/learnings.md` (consider promoting to conventions if the pattern recurs 3+ times)
5. **Log findings**: If the work produced a scientific finding (empirical observation, validated/invalidated hypothesis, quantitative result, or domain discovery), crystallize it to `.living/findings/{topic}.md`. Check existing topics first for consistency. Prefer broad topic names.
6. **Log todos**: If future work is identified, add items to `todo/TODO_REGISTRY.md` (create a detailed `todo/[item].md` writeup for complex items)
7. **Validate structure**: Run `validate_structure.py` to confirm repo structure is correct
8. **Crystallize conventions**: Review recent learnings — if 3+ entries share a pattern, promote to `.living/conventions.md` or a named convention pack
9. **Convention feedback**: Note whether convention pack practices were helpful or had gaps
10. **Session summary**: Write `.claude/last-session.md` with a brief summary of what was done, decisions made, and next steps

### Automated Enforcement

Mycelium hooks are installed in `.claude/settings.local.json`:

| Hook | Event | What it does |
|------|-------|--------------|
| `mycelium-health.sh` | SessionStart | Checks `.living/` health, records session timestamp |
| `mycelium-post-action.sh` | PostToolUse (Bash) | Detects code execution and directs Claude to run the full post-action protocol (debounced) |
| `mycelium-activity-tracker.sh` | PostToolUse (Edit\|Write) | Tracks file edits |
| `mycelium-read-tracker.sh` | PostToolUse (Read) | Logs `.living/` access |
| `mycelium-stop-check.sh` | Stop | Safety net — blocks session end if the post-action hook fired but `.living/` was never updated |

## Data Conventions

- `data/raw/` is **IMMUTABLE** — never modify original files
- Every dataset has metadata in `data/metadata/[dataset-name]/`
- Every dataset has a manifest entry in `data/DATA_MANIFEST.md`
- Large files are gitignored with download instructions documented

## Key Files

| File | Purpose |
|------|---------|
| `ENVIRONMENTS_INSTALLATIONS.md` | How to set up the environment |
| `.living/INDEX.md` | Auto-regenerated knowledge map — entry point for `.living/` |
| `.living/decisions.md` | Why choices were made |
| `.living/learnings.md` | Accumulated insights and gotchas |
| `.living/conventions.md` | Repo-specific convention overrides |
| `.living/conventions/ACTIVE_CONVENTIONS.yaml` | Registry of installed convention packs |
| `todo/TODO_REGISTRY.md` | Master list of future work items |
| `*/*_MANIFEST.md` | Registry of contents in each directory |