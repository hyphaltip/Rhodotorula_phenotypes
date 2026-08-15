# Learnings

Append-only log of gotchas, surprises, and insights.

**Entry template:** copy from `skills/core/templates/learning-entry.md` (includes Category, What happened, Why it matters, Resolution, Tags fields). The `**Tags**:` line is consumed by `generate_index.py --summary-heuristic` to build the cluster summary in INDEX.md — use them.

### [2026-08-15] L-1 — mycelium scripts require Python >= 3.10
- **Category**: tooling
- **What happened**: `install_convention.py` (and other mycelium core scripts) use PEP 604 union syntax (`Path | None`), which throws `TypeError: unsupported operand type(s) for |` under the system default `python3` (3.9.18 on this host).
- **Why it matters**: Home-page python3 is too old; running the scripts fails loudly but misleadingly during mycelium init.
- **Resolution**: Run with `/usr/bin/python3.12` (available) or under the repo pixi env; also pass `--network-dir` explicitly because the script's auto-detect may miss the plugin marketplace path.
- **mitigation_type**: structural — add a version guard/`sys.version_info` check or a shebang in the scripts.
- **structural_mitigation_candidate**: assert `sys.version_info >= (3,10)` with a clear message at script start.
- **Tags**: tooling, mycelium, python, setup

### [2026-08-15] L-2 — skillpacks/Autonomous-Science repo is no longer public
- **Category**: tooling
- **What happened**: `git clone https://github.com/arjunrajlaboratory/Autonomous-Science.git` fails (`Repository not found`) — the repo is absent from the org's public repos. Only `scientific-agent-skills` and `bioSkills` cloned successfully.
- **Why it matters**: The `skill-bridge` convention's `skill-sources.yaml` lists this as a verified persona source (`personas/library/arjun_raj.json`); the persona routing feature is therefore currently unavailable.
- **Resolution**: Left `skillpacks/` with the two working repos. TODO: re-point `skill-bridge` `skill-sources.yaml` to a live personas source or drop the source entry.
- **mitigation_type**: ambient-awareness — monitor whether the org restores/renames the repo.
- **Tags**: tooling, skill-bridge, skillpacks, third-party
