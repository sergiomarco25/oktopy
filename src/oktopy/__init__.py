"""oktopy: the field-independent engine behind a re-tested claim store.

Published conclusions of a field's single-cell studies become structured claims; each claim is re-tested at donor
level in its own dataset and in every other harmonised dataset, pooled across datasets, rated (verdict and certainty)
and served as a SQLite store, a knowledge graph, an atlas page and an MCP server.

Modules: config (the project file), schema, claims (merge), contexts, signatures, harmonise, pseudobulk, de_cache
(edgeR, R), retest (claim tests), claim_meta (pooling), build_store, rating, atlas, kg (knowledge graph), registry,
calibration, power and multiverse (robustness), release, release_diff, public_findings, linked_data, audit,
mcp_server and mcp_http (serving), newproject and example (starting a project), runner (all steps in order).
Everything about one field (its datasets, papers, cell types, regions, contrasts, wording) comes from the project's
config files and claim files; the engine names none of it (tests/test_engine_boundary.py). A project adds its own
package for what is specific to it: dataset loaders and any analyses of its own.
"""
__version__ = "0.1.0"
