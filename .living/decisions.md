# Decision Log

Append-only log of non-obvious decisions and their rationale.

**Entry template:** copy from `skills/core/templates/decision-log-entry.md` (includes Context, Decision, Alternatives considered, Rationale, Consequences, Tags fields).

### [2026-08-15] D-1 — Completed mycelium init: install all convention packs; skip broken skillpack
- **Context**: Finishing a partially-completed mycelium `init` for this repo. Core packs are mandatory; domain packs were user-selected.
- **Decision**: Installed core (`robust-analysis`, `report-generator`, `idea-generator`) plus domain (`bioinformatics`, `image-analysis`, `skill-bridge`) convention packs. Cloned `scientific-agent-skills` and `bioSkills` into `skillpacks/`; did **not** clone `Autonomous-Science` because `github.com/arjunrajlaboratory/Autonomous-Science` no longer exists publicly.
- **Alternatives considered**: Install only core packs; skip domain packs until needed; wait for `Autonomous-Science` to be restored.
- **Rationale**: User explicitly requested all three domain packs. The dead repo cannot be cloned, and blocking on it stalls init; the `skill-bridge` pack still functions for the two successfully-cloned skill repos.
- **Consequences**: `skill-bridge` persona-routing (personas/library) is currently unavailable (see L-2). Official CLAUDE.md generated; `.living/INDEX.md` built; `validate_structure.py` passes.
- **Tags**: setup, mycelium-init, conventions, skillpacks
