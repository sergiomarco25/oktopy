# Region glossary

One row per region code of `config/schema.toml` and `config/harmonise.toml`. Papers name regions in many ways
(Brodmann areas, gyri, abbreviations); claims and datasets must agree on what a code means, or a claim is tested in the
wrong tissue. Fill it before extracting claims, and check every claim's region against it.

| code | full name | also called (in papers and datasets) | Brodmann area / atlas | datasets |
|---|---|---|---|---|
| `example` | Example region | the names papers use | area | datasets with it |

Notes:
- A Brodmann area is not a synonym for a gyrus name: check the area (for example, area 10 is frontopolar cortex,
  not dorsolateral prefrontal cortex).
- When a dataset's region is broader or narrower than a claim's, map it in `config/contexts.toml` `[region_map]`
  and say so here; the rating counts other-region evidence as indirect.
