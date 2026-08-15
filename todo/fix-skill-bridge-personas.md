# Re-point skill-bridge personas source

| Field | Value |
|-------|-------|
| **Date** | 2026-08-15 |
| **Author** | jstajich |
| **Priority** | medium |
| **Status** | open |
| **Category** | infrastructure |
| **Related analyses** | — |
| **Related data** | — |

## Description

`skillpacks/Autonomous-Science` (github.com/arjunrajlaboratory/Autonomous-Science) no longer exists publicly, so skill-bridge's persona-routing source is dead. The `skill-bridge` `skill-sources.yaml` still lists it as a verified source with `verify: personas/library/arjun_raj.json`.

## Motivation

The skill-bridge convention pack's persona-guided review feature (personas/library of 50 researcher personas) is currently unavailable. Either restore a working source or drop the entry so the pack doesn't reference a dead repo.

## Proposed Approach

- Check whether the org restored/renamed the repo; if found, update the `url` and `verify` path in `.living/conventions/skill-bridge/skill-sources.yaml`.
- If permanently gone, remove the `autonomous-science-personas` source entry and note the removal in `.living/conventions.md`.

## Acceptance Criteria

- [ ] `skill-sources.yaml` references only repos that exist (both cloned skillpacks verify OK).
- [ ] Decision recorded in `.living/decisions.md`.

## Notes

See `.living/learnings.md` entry L-2.