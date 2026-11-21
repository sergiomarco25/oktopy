# Region glossary

One row per region code of `config/schema.toml` and `config/harmonise.toml`. Papers name regions in many ways; claims
and datasets must agree on what a code means, or a claim is tested in the wrong tissue. Fill it before extracting
claims, and check every claim's region against it.

| code | full name | also called (in papers and datasets) | atlas | datasets |
|---|---|---|---|---|
| `front` | front part of the example tissue | anterior | - | smith2024, jones2025 |
| `back` | back part of the example tissue | posterior | - | smith2024 |

Notes:
- When a dataset's region is broader or narrower than a claim's, map it in `config/contexts.toml` `[region_map]`
  and say so here; the rating counts other-region evidence as indirect.
